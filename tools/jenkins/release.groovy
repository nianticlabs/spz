#!/usr/bin/env groovy

@Library('algpipeline') _

def setRunOnMain() {
  if (env.BRANCH_NAME in ['main']) {
    return [githubPush()]
  }
  return []
}

properties([
  buildDiscarder(logRotator(daysToKeepStr: '90', numToKeepStr: '50', artifactNumToKeepStr: '3')),
  // only run release packaging for changes on the 'main' branch
  // (i.e., do not trigger packaging when git tags are pushed to the repo)
  pipelineTriggers(
    setRunOnMain()
  )
])

// zap_venv is for cleaning and rebuilding the python virtual environment -- you'll want to
// set this to 'true' for a Jenkins Replay when you've updated a git+ssh python dependency
// (if there is no module version bump, pip will keep the stale one around by default)
def zap_venv = true

// Build for Linux and macOS across Python 3.10, 3.11, 3.12
def profiles = [
  // Linux
  [host:'ubuntu', name: 'Ubuntu-Python3.10', label: 'builder&&(Ubuntu22||RedHat8)', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useGcc(11), zap_py_venv: zap_venv, conda_file: './tools/conda310.yaml'],
  [host:'ubuntu', name: 'Ubuntu-Python3.11', label: 'builder&&(Ubuntu22||RedHat8)', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useGcc(11), zap_py_venv: zap_venv, conda_file: './tools/conda311.yaml'],
  [host:'ubuntu', name: 'Ubuntu-Python3.12', label: 'builder&&(Ubuntu22||RedHat8)', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useGcc(11), zap_py_venv: zap_venv, conda_file: './tools/conda312.yaml'],
  // macOS
  [host:'mac', name: 'MacOS-Python3.10', label: 'Xcode_16_0&&mac', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useXcode('16.0'), zap_py_venv: zap_venv, conda_file: './tools/conda310.yaml'],
  [host:'mac', name: 'MacOS-Python3.11', label: 'Xcode_16_0&&mac', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useXcode('16.0'), zap_py_venv: zap_venv, conda_file: './tools/conda311.yaml'],
  [host:'mac', name: 'MacOS-Python3.12', label: 'Xcode_16_0&&mac', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useXcode('16.0'), zap_py_venv: zap_venv, conda_file: './tools/conda312.yaml'],
]

smartStage("Load utils.groovy") {
    node('builder') {
    utils = scmLoad(scm, 'tools/jenkins/utils.groovy')
    }
}

env.bump_tag = ''
env.tag = ''

def build_spz_wheel(profile) {
    smartStage(profile.name) {
    nodeTimeout(profile.label, [time: profile.timeout, unit: profile.timeout_unit]) {
      withCredentials([
        // These are credentials available on Jenkins. This will make them avaialble as env vars for build.
        usernamePassword(credentialsId: 'aws_access_techtransfer', usernameVariable: 'AWS_ACCESS_KEY_ID', passwordVariable: 'AWS_SECRET_ACCESS_KEY'),
        usernamePassword(credentialsId: 's3dibot_artifactory_adobe', passwordVariable: 'ARTIFACTORY_API_KEY', usernameVariable: 'ARTIFACTORY_USER'),
        usernamePassword(credentialsId: 's3dibot-artifactory-uw2-apikey', passwordVariable: 'ARTIFACTORY_UW2_API_KEY', usernameVariable: 'ARTIFACTORY_UW2_USER'),
      ]) {
        withSshCredentials() {
          try {
              def venv_name = '.venv'

              smartStage("Setup-${profile.name}") {
                smartCleanWs() // ensure there aren't any old .whl builds hanging around
                printEnvVarsAndJobParams()
                checkout scm
                utils.setupCondaEnvironment(venv_name, profile)
              } // setup

              smartStage("Wheel-${profile.name}") {
                def wheel_version = "v${env.bump_tag}"
                def release_mode = "release"
                utils.pythonWheelOps(venv_name, wheel_version, release_mode, profile)
              } // wheel
            currentBuild.result = "SUCCESS"
          } catch(Exception e) {
            println("+++++ build_spz_wheel() for ${profile.name} on ${NODE_NAME} caught Exception:\n    ${e}")
            currentBuild.result = "FAILURE"
            // archive any logs lying around before the exception
            utils.archiveLogs(profile.host)
            throw e   // Need this so that Jenkins job can declare the correct verdict
          } finally {
            println("+++++ currentBuild.result=${currentBuild.result}")

            if (currentBuild.result != "FAILURE") {
              println("+++++ Cleaning up ${WORKSPACE} on ${NODE_NAME}")
              smartCleanWs()
            }

            printEnvVarsAndJobParams()
          } // finally
        } // withSshCredentials
      } // with credentials
    } // nodeTimeout
  } // stage
} // build_spz_wheel

timestamps {
  // Cancel old jobs on the same branch
  abortPreviousBuild(isPRBuild())
  smartStage("Tag release") {
    node('builder') {
      withSshCredentials() {
        checkout scm
        env.tag = sh(returnStdout: true, script: "git ls-remote --tags | awk -F'/' '{print $3}' | sort -V | tail -n 1 | sed 's/^v//' | sed 's/\^{}$//'").trim()
        env.bump_tag = sh(returnStdout: true, script: '''echo "${tag%.*}.$((${tag##*.}+1))"''').trim()
        sh(script: '''git tag v${bump_tag} && git push --tags''')
      }
    }
  }
  parallel profiles.collectEntries { profile ->
    [ "${profile.name}" : {
        build_spz_wheel(profile)
    } ]
  }
}