#!/usr/bin/env groovy

@Library('algpipeline') _

// This file is adapted from https://git.corp.adobe.com/3di/python-scaffold
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

// Please be mindful of the limited resources on Jenkins when setting the timeout and issuing retries for jobs.
// If you need to replay a job on only one platform, please comment the other platforms out for the replay to conserve build resources
def profiles = [
  [host:'mac', name: 'MacOS', label: 'Xcode_16_0', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useXcode('16.0')],
  [host:'ubuntu', name: 'Ubuntu' ,  label: 'builder&&(Ubuntu22||RedHat8)', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useGcc(11)],
]

smartStage("Load utils.groovy") {
    node('builder') {
    utils = scmLoad(scm, 'tools/jenkins/utils.groovy')
    }
}

env.bump_tag = ''
env.tag = ''

smartStage("Compute and Push Tag") {
  node('builder') {
    withCredentials([
      string(credentialsId: 'git-ghec-api-s3dibot-apitoken', variable: 'GIT_TOKEN'),
    ]) {
      sshagent(['gerrit-ssh', 's3dibot_gitcorp_ssh']) {
      checkout scm
      env.tag = cmd(returnStdout: true, script: "git ls-remote --tags https://x-access-token:${GIT_TOKEN}@github.com/Adobe-3di/neural-assets.git | awk -F'/' '{print \$3}' | sort -V | tail -n 1 | sed 's/^v//' | sed 's/\\^{}\$//'").trim()
      env.bump_tag = cmd(returnStdout: true, script: '''echo "v${tag%.*}.$((${tag##*.}+1))"''')
      cmd(script: '''git tag ${bump_tag} && git push https://x-access-token:${GIT_TOKEN}@github.com/Adobe-3di/neural-assets.git --tags''')
      }
    }
  }
}

def build_neural_assets_wheel(profile, bump_tag) {
    smartStage(profile.name) {
    nodeTimeout(profile.label, [time: profile.timeout, unit: profile.timeout_unit]) {
      withCredentials([
        // These are credentials available on Jenkins. This will make them avaialble as env vars for build.
        usernamePassword(credentialsId: 'aws_access_techtransfer', usernameVariable: 'AWS_ACCESS_KEY_ID', passwordVariable: 'AWS_SECRET_ACCESS_KEY'),
        usernamePassword(credentialsId: 's3dibot_artifactory_adobe', passwordVariable: 'ARTIFACTORY_API_KEY', usernameVariable: 'ARTIFACTORY_USER'),
        usernamePassword(credentialsId: 's3dibot-artifactory-uw2-apikey', passwordVariable: 'ARTIFACTORY_UW2_API_KEY', usernameVariable: 'ARTIFACTORY_UW2_USER'),
        string(credentialsId: 'git-ghec-api-s3dibot-apitoken', variable: 'GIT_TOKEN'),
      ]) {
        sshagent(['gerrit-ssh', 's3dibot_gitcorp_ssh']) {
          try {
              def venv_name = '.venv'

              smartStage("Setup-${profile.host}") {
                smartCleanWs() // ensure there aren't any old .whl builds hanging around
                printEnvVarsAndJobParams()
                checkout scm
                utils.setupCondaEnvironment(venv_name, profile)
              } // setup

              smartStage("Wheel-${profile.host}") {
                def wheel_version = bump_tag
                def release_mode = "release"
                utils.pythonWheelOps(venv_name, wheel_version, release_mode, profile)
              } // wheel
            currentBuild.result = "SUCCESS"
          } catch(Exception e) {
            println("+++++ build_neural_assets_wheel() for ${profile.name} on ${NODE_NAME} caught Exception:\n    ${e}")
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
        } // sshagent
      } // with credentials
    } // nodeTimeout
  } // stage
} // build_neural_assets_wheel

timestamps {
  // Cancel old jobs on the same branch
  abortPreviousBuild(isPRBuild())
  parallel profiles.collectEntries { profile ->
    [ "${profile.name}" : {
        build_neural_assets_wheel(profile, env.bump_tag)
    } ]
  }
}