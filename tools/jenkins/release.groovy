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
  // Python 3.11 builds
  [host:'mac', name: 'MacOS-Python3.11', python_version: '3.11', label: 'Xcode_16_0&&mac', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useXcode('16.0'), zap_py_venv: zap_venv],
  [host:'ubuntu', name: 'Ubuntu-Python3.11', python_version: '3.11', label: 'builder&&(Ubuntu22||RedHat8)', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useGcc(11), zap_py_venv: zap_venv],
  [host:'windows', name: 'Windows-Python3.11', python_version: '3.11', label: 'builder&&win', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.vs(2022, 'x64', '14.38'), zap_py_venv: zap_venv],

  // Python 3.12 builds
  [host:'mac', name: 'MacOS-Python3.12', python_version: '3.12', label: 'Xcode_16_0&&mac', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useXcode('16.0'), zap_py_venv: zap_venv],
  [host:'ubuntu', name: 'Ubuntu-Python3.12', python_version: '3.12', label: 'builder&&(Ubuntu22||RedHat8)', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.useGcc(11), zap_py_venv: zap_venv],
  [host:'windows', name: 'Windows-Python3.12', python_version: '3.12', label: 'builder&&win', timeout: '45', timeout_unit: 'MINUTES', toolchain: cmd.vs(2022, 'x64', '14.38'), zap_py_venv: zap_venv],
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
    sshagent(['gerrit-ssh', 's3dibot_gitcorp_ssh']) {
    checkout scm
    // Get the current branch name and extract the version
    def current_branch = cmd(returnStdout: true, script: "git rev-parse --abbrev-ref HEAD").trim()
    if (current_branch == "HEAD") {
        // We're in detached HEAD state, try to find the branch containing this commit
        def commit_hash = cmd(returnStdout: true, script: "git rev-parse HEAD").trim()
        def branch_cmd
        if (isUnix()) {
            branch_cmd = "git branch -a --contains ${commit_hash} | grep 'remotes/origin/adobe-v' | head -n 1 | sed 's/remotes\\/origin\\///'"
        } else {
            // Windows PowerShell command for branch detection
            branch_cmd = """powershell -Command "\$branches = git branch -a --contains ${commit_hash} | Select-String 'remotes/origin/adobe-v'; if (\$branches) { \$branches[0].ToString().Trim().Replace('remotes/origin/', '') } else { Write-Output '' }" """
        }
        current_branch = cmd(returnStdout: true, script: branch_cmd).trim()
        if (!current_branch) {
            error "Could not determine branch name from commit ${commit_hash}"
        }
    }
    def version_match = current_branch =~ /adobe-v(\d+\.\d+)\.\d+/
    if (!version_match) {
        error "Current branch '${current_branch}' does not match expected format 'adobe-vX.Y.Z'"
    }
    def major_minor = version_match[0][1]

    // Find the latest tag matching the major.minor version
    def tag_cmd
    if (isUnix()) {
        tag_cmd = "git ls-remote --tags | grep 'v${major_minor}' | awk -F'/' '{print \$3}' | sort -V | tail -n 1 | sed 's/^v//' | sed 's/\\^{}\$//'"
    } else {
        // Windows PowerShell command with null handling
        tag_cmd = """powershell -Command "\$tags = git ls-remote --tags | Select-String 'v${major_minor}'; if (\$tags) { \$tags | ForEach-Object { \\$_.ToString().Split('/')[-1] } | Sort-Object { [version]\\$_.TrimStart('v').TrimEnd('^{}') } | Select-Object -Last 1 | ForEach-Object { \\$_.TrimStart('v').TrimEnd('^{}') } } else { Write-Output '' }" """
    }
    env.tag = cmd(returnStdout: true, script: tag_cmd).trim()

    // If no matching tag exists, start with .0
    if (!env.tag) {
    env.tag = "${major_minor}.0"
    }

    // Only increment the patch number
    env.bump_tag = cmd(returnStdout: true, script: '''echo "v${tag%.*}.$((${tag##*.}+1))"''')
    cmd(script: '''git tag ${bump_tag} && git push --tags''')
    }
  }
}

def build_spz_wheel(profile, bump_tag) {
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
                def wheel_version = bump_tag
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
  parallel profiles.collectEntries { profile ->
    [ "${profile.name}" : {
        build_spz_wheel(profile, env.bump_tag)
    } ]
  }
}