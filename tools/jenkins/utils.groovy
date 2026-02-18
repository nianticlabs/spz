def setupUvEnvironment(profile) {
  cmd 'uv venv --python 3.12'
  cmd 'uv pip install -r tools/requirements.txt'
  cmd script: 'uv pip install -e ".[tests]" --upgrade', buildEnv: profile.toolchain
}

def testPython(profile) {
  try {
    uvr(script: "pytest --junit-xml=.tmp_test_out/${profile.name}/test_report.xml --log-file=.tmp_test_out/${profile.name}/test_log.txt", buildEnv: profile.toolchain)
  } finally {
    archiveArtifacts artifacts: ".tmp_test_out/${profile.name}/test_log.txt"
    junit testResults: ".tmp_test_out/${profile.name}/test_report.xml"
  }
}

def pythonWheelOps(wheel_version, release_mode, profile, is_pr=false) {
  withEnv(["PIP_EXTRA_INDEX_URL=https://$ARTIFACTORY_UW2_USER:$ARTIFACTORY_UW2_API_KEY@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-tech-transfer-3di-release/simple https://$ARTIFACTORY_UW2_USER:$ARTIFACTORY_UW2_API_KEY@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-adobeshared-release/simple"]) {
    def repo_url_tech_transfer = "https://artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-tech-transfer-3di-${release_mode}"
    def repo_url_adobeshared = "https://artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-adobeshared-${release_mode}"
    def script = "boa_toolkit package build --wheel_dir wheelhouse --version ${wheel_version}"
    echo "Running: ${script}"
    uvr(script: script, buildEnv: profile.toolchain)
    uvr(script: "boa_toolkit package check --wheel_dir wheelhouse --smoke_test tools/smoke_test.py", buildEnv: profile.toolchain)
    if (!is_pr) {
        echo "Uploading ${repo_url_adobeshared}..."
        uvr(script: "boa_toolkit package upload --wheel_dir wheelhouse --pypi_url ${repo_url_adobeshared}", buildEnv: profile.toolchain)
    }
    echo "Uploading ${repo_url_tech_transfer}..."
    uvr(script: "boa_toolkit package upload --wheel_dir wheelhouse --pypi_url ${repo_url_tech_transfer}", buildEnv: profile.toolchain)
  }
}

def archiveLogs(host) {
  archiveArtifacts artifacts: ".tmp_test_out/${host}/test_log.txt", allowEmptyArchive: true
  junit testResults: ".tmp_test_out/${host}/test_report.xml"
}

def withSlackNotificationsOnFailure(Closure body) {
  // change the channel here if needed
  slackNotify("#tech-transfer-build-status", [skip: ["SUCCESS", "UNSTABLE"], blocksBuilder: { ctx ->
    slackNotify.getDefaultBlocks(ctx).plus([
      [
        "type": "section",
        "text": [
          "type": "mrkdwn",
          "text": "Build failed" // optionally add @ldap here to get the notification
        ]
      ]
    ])
  }]) {
    body()
  }
}

return this