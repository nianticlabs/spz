#!/usr/bin/env groovy

@Library('algpipeline') _

// This file is adapted from https://git.corp.adobe.com/3di/python-scaffold
properties([
  buildDiscarder(logRotator(daysToKeepStr: '15', numToKeepStr: '75'))
])

def zap_venv = false
def clean_workspace = false

// Please be mindful of the limited resources on Jenkins when setting the timeout and issuing retries for jobs.
// If you need to replay a job on only one platform, please comment the other platforms out for the replay to conserve build resources
def profiles = [
  [host:'mac', name: 'MacOS', label: 'Xcode_16_0&&mac', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useXcode('16.0')],
  [host:'ubuntu', name: 'Ubuntu' ,  label: 'builder&&(Ubuntu22||RedHat8)', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useGcc(11)],
]

def wheel_profile = [name:'linux', platform: 'linux', label: 'builder&&linux', zap_py_venv: false, timeout: [time: '60', unit: 'MINUTES']]

smartStage("Load utils.groovy") {
    node('builder') {
        utils = scmLoad(scm, 'tools/jenkins/utils.groovy')
    }
}

def build_py_project_scaffold(profile) {
  smartStage(profile.name) {
    nodeTimeout(profile.label, [time: profile.timeout, unit: profile.timeout_unit]) {
      withCredentials([
          usernamePassword(credentialsId: 'aws_access_techtransfer', usernameVariable: 'AWS_ACCESS_KEY_ID', passwordVariable: 'AWS_SECRET_ACCESS_KEY'),
          usernamePassword(credentialsId: 's3dibot_artifactory_adobe', passwordVariable: 'ARTIFACTORY_API_KEY', usernameVariable: 'ARTIFACTORY_USER'),
          usernamePassword(credentialsId: 's3dibot-artifactory-uw2-apikey', passwordVariable: 'ARTIFACTORY_UW2_API_KEY', usernameVariable: 'ARTIFACTORY_UW2_USER')
      ]) {
          withSshCredentials() {
          try {
            def venv_name = '.venv'

            smartStage("Setup-${profile.host}") {
              if (profile.clean_workspace) {
                smartCleanWs()  // from algpipeline with smart retries
              }
              printEnvVarsAndJobParams()
              checkout scm
              utils.setupCondaEnvironment(venv_name, profile)
              utils.runInConda(name: venv_name, script: "python tools/artifacts.py test_data_manifests/unit_test_neural_assets.yaml clean")
              utils.runInConda(name: venv_name, script: "python tools/artifacts.py test_data_manifests/unit_test_neural_assets.yaml pull")
            } // Setup

            smartStage("Lint-${profile.host}") {
              utils.lintPython(venv_name)
            } // LintPython

            smartStage("Test-${profile.host}") {
              utils.testPython(venv_name, profile)
            } // TestPython

            smartStage("BuildTypeScript") {
              utils.runInConda(name: venv_name, script: "npx tsc --build")
            } // BuildTypeScript
            currentBuild.result="SUCCESS"

          } catch(Exception e) {
            println("${profile.name} on ${NODE_NAME} caught Exception:\n    ${e}")
            currentBuild.result = "FAILURE"
            throw e   // Need this so that Jenkins job can declare the correct verdict
          } finally {
            println("currentBuild.result=${currentBuild.result}")
            if (currentBuild.result != "FAILURE") {
              println("Cleaning up ${WORKSPACE} on ${NODE_NAME}")
              cleanWs()
            }
            // dump env vars again to show newly added ones during the run
            printEnvVarsAndJobParams()
          } // finally
        } // sshagent
      } // with credentials
    } // nodeTimeout()
  } // smartStage(profile.name)
} // Build

timestamps {
  // Cancel old jobs on the same branch
  abortPreviousBuild(isPRBuild())
  parallel profiles.collectEntries { profile ->
    [ "${profile.name}" : {
        build_py_project_scaffold(profile)
    } ]
  }

  smartStage("Build wheel") {
    nodeTimeout(wheel_profile.label, [time: wheel_profile.timeout.time, unit: wheel_profile.timeout.unit]) {
      withCredentials([
        // These are credentials available on Jenkins. This will make them avaialble as env vars for build.
        usernamePassword(credentialsId: 'aws_access_techtransfer', usernameVariable: 'AWS_ACCESS_KEY_ID', passwordVariable: 'AWS_SECRET_ACCESS_KEY'),
        usernamePassword(credentialsId: 's3dibot_artifactory_adobe', passwordVariable: 'ARTIFACTORY_API_KEY', usernameVariable: 'ARTIFACTORY_USER'),
        usernamePassword(credentialsId: 's3dibot-artifactory-uw2-apikey', passwordVariable: 'ARTIFACTORY_UW2_API_KEY', usernameVariable: 'ARTIFACTORY_UW2_USER'),
      ]) {
        sshagent(['gerrit-ssh', 's3dibot_gitcorp_ssh']) {
          try {
            smartCleanWs() // ensure there aren't any old .whl builds hanging around
            checkout scm
            def venv_name = '.venv'
            utils.setupCondaEnvironment(venv_name, wheel_profile)
            def wheel_version = "0.0.0.dev+pr${CHANGE_ID}build${BUILD_NUMBER}"
            def release_mode = "dev"
            utils.pythonWheelOps(venv_name, wheel_version, release_mode, wheel_profile)
          } catch(Exception e) {
            println("+++++ Build wheel on ${NODE_NAME} caught Exception:\n    ${e}")
            currentBuild.result = "FAILURE"
            throw e   // Need this so that Jenkins job can declare the correct verdict
          } finally {
            println("+++++ currentBuild.result=${currentBuild.result}")
            printEnvVarsAndJobParams()
          }
        }  // withSshCredentials()
      } // with credentials
    } // nodeTimeout
  }
}