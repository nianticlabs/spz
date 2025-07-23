def setupCondaEnvironment(venv_name, profile) {
    def debug_flag = '' // set to '-vvv' for detailed debugging info

    echo "setupCondaEnvironment() venv_name=${venv_name}"
    if (profile.zap_py_venv) {
        def file_exists = fileExists(venv_name)
        if (file_exists) {
            runInConda(name: venv_name, label: "Delete existing venv", script: "conda deactivate && conda env remove --prefix ${WORKSPACE}/${venv_name} -y")
            echo "INFO: Conda venv ${venv_name} has been removed."
        }
        else {
            echo "INFO: Conda venv ${venv_name} does not exist. No need to delete the venv."
        }
    }

    // Build conda create command with optional Python version override
    def extra_params = debug_flag

    runInConda.createEnv(name: venv_name, file: "./tools/conda.yaml", extra_params: extra_params)

    runInConda(name: venv_name, label: "Check tools version",
    script: [
            "conda --version",
            "python --version"
    ])

    def setup_script = nodeUtils.shyIsUnix() ? '. tools/setup.sh' : 'powershell tools\\setup.ps1'
    if (profile.host == 'mac') {
        setup_script = 'source useXcode 16.0 MacOSX && ' + setup_script
    }
    // set env var for debugging the script: CONDA_DEBUG_FLAG=-vvv if needed
    env_list = ["CONDA_DEBUG_FLAG=${debug_flag}"]
    withEnv(env_list) {
        runInConda(name: venv_name, label: "Setup virtual env", script: setup_script)
    }
}

def lintPython(venv_name) {

    def linter_failures = []
    try {
        runInConda(name: venv_name, script: "python ./tools/lint.py ruff --verbose")
    } catch (Exception e) {
        echo "Got exception: pylint ${e}"
        linter_failures.add("Running ruff")
    }

    try {
        runInConda(name: venv_name, script: "python ./tools/lint.py mypy --verbose")
    } catch (Exception e) {
        echo "Got exception: mypy ${e}"
        linter_failures.add("Running mypy")
    }

    if (!linter_failures.isEmpty()) {
        String msg = "The following linter phases have failed\n"
        for (String m in linter_failures) {
        msg = msg + m + "\n"
        }
        msg = msg + "To see the faulty files, see the corresponding section of each linter phase.\n"
        error(msg)
    }
}

def testPython(venv_name, profile) {
    try {
    runInConda(name: venv_name, script: "python tests/main.py --xml=.tmp_test_out/${profile.name}/test_report.xml --log=.tmp_test_out/${profile.name}/test_log.txt")
    } finally {
    archiveArtifacts artifacts: ".tmp_test_out/${profile.name}/test_log.txt"
    junit testResults: ".tmp_test_out/${profile.name}/test_report.xml"
}}

def pythonWheelOps(venv_name, wheel_version, release_mode, profile, is_pr=false) {

    withEnv(["PIP_EXTRA_INDEX_URL=https://$ARTIFACTORY_UW2_USER:$ARTIFACTORY_UW2_API_KEY@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-tech-transfer-3di-release/simple https://$ARTIFACTORY_UW2_USER:$ARTIFACTORY_UW2_API_KEY@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-adobeshared-release/simple"]) {
        def repo_url_tech_transfer = "https://artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-tech-transfer-3di-${release_mode}"
        def repo_url_adobeshared = "https://artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-adobeshared-${release_mode}"
        try {
            def script = "boa_toolkit package build --wheel_dir wheelhouse --version ${wheel_version}"
            if (profile.host == 'mac') {
                script = 'source useXcode 16.0 MacOSX && ' + script
            }
            echo "Running: ${script}"
            runInConda(name: venv_name, script: script)
        } catch (Exception e) {
            echo "Got exception: ${e}"
            throw e
        }
        try {
            runInConda(name: venv_name, script: "boa_toolkit package check --wheel_dir wheelhouse --smoke_test tools/smoke_test.py")
        } catch (Exception e) {
            echo "Got exception: ${e}"
            throw e
        }
        try {
            if (!is_pr) {
                echo "Uploading ${repo_url_adobeshared}..."
                runInConda(name: venv_name, script: "boa_toolkit package upload --wheel_dir wheelhouse --pypi_url ${repo_url_adobeshared}")
            }
            echo "Uploading ${repo_url_tech_transfer}..."
            runInConda(name: venv_name, script: "boa_toolkit package upload --wheel_dir wheelhouse --pypi_url ${repo_url_tech_transfer}")
        } catch (Exception e) {
            echo "Got exception: ${e}"
            throw e
        }
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