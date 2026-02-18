def setupUvEnvironment(profile) {
  // Verify uv and uvr are available on this node
  cmd 'uv --version'
  uvr(script: 'python --version')

  def python_version = profile.python_version ?: '3.12'
  if (profile.zap_py_venv || !fileExists('.venv')) {
    cmd "uv venv --python ${python_version}"
  }
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
  // Clean and prepare wheelhouse directory
  echo "Preparing wheelhouse directory..."
  if (nodeUtils.shyIsUnix()) {
    cmd("mkdir -p wheelhouse && rm -f wheelhouse/*.whl")
  } else {
    cmd("if not exist wheelhouse mkdir wheelhouse")
    cmd("del /Q wheelhouse\\*.whl 2>nul || echo No existing wheels to clean")
  }

  // Define repository names based on release mode
  def tech_transfer_index = (release_mode == "release") ? "tech-transfer-release" : "tech-transfer-dev"

  // Step 1: Build the wheel
  echo "Building wheel version ${wheel_version}..."

  // For scikit-build-core, we modify pyproject.toml to use a static version
  // CMake project() VERSION must be numeric-only, so we can't use dev versions there
  // By changing 'dynamic = ["version"]' to 'version = "..."', scikit-build-core
  // will use the static version instead of reading from CMakeLists.txt
  // If wheel_version matches the version already in CMakeLists.txt, skip the override
  def cmake_version = cmd(returnStdout: true, script: "grep -oP 'VERSION\\s+\\K[0-9]+\\.[0-9]+\\.[0-9]+' CMakeLists.txt | head -1").trim()
  if (wheel_version != cmake_version) {
    echo "Modifying pyproject.toml to set version ${wheel_version}..."
    if (nodeUtils.shyIsUnix()) {
      cmd("sed -i 's/dynamic = \\[\"version\"\\]/version = \"${wheel_version}\"/' pyproject.toml")
    } else {
      cmd("powershell -Command \"(Get-Content pyproject.toml) -replace 'dynamic = \\[\\\"version\\\"\\]', 'version = \\\"${wheel_version}\\\"' | Set-Content pyproject.toml\"")
    }
  }

  cmd(script: "uv build --wheel --out-dir wheelhouse", buildEnv: profile.toolchain)

  // Step 2: Validate the wheel was built
  def wheelFiles
  if (nodeUtils.shyIsUnix()) {
    wheelFiles = cmd(returnStdout: true, script: 'ls wheelhouse/*.whl 2>/dev/null | wc -l').trim().toInteger()
  } else {
    wheelFiles = cmd(returnStdout: true, script: 'powershell -Command "(Get-ChildItem wheelhouse\\*.whl -ErrorAction SilentlyContinue | Measure-Object).Count"').trim().toInteger()
  }

  if (wheelFiles == 0) {
    error "No wheel files found in wheelhouse/ directory after build"
  }

  // Get wheel filename for logging
  def wheel_file
  if (nodeUtils.shyIsUnix()) {
    wheel_file = cmd(returnStdout: true, script: 'ls wheelhouse/*.whl | head -1 | xargs basename').trim()
  } else {
    wheel_file = cmd(returnStdout: true, script: 'powershell -Command "Get-ChildItem wheelhouse\\*.whl | Select-Object -First 1 -ExpandProperty Name"').trim()
  }
  echo "Built wheel: ${wheel_file}"

  // Step 3: Smoke test - install wheel in isolated environment and test import
  echo "Running smoke test..."
  if (nodeUtils.shyIsUnix()) {
    uvr(script: "uv run --with wheelhouse/*.whl --no-project --index-strategy unsafe-best-match -- python tools/smoke_test.py", buildEnv: profile.toolchain)
  } else {
    uvr(script: "uv run --with wheelhouse/${wheel_file} --no-project --index-strategy unsafe-best-match -- python tools/smoke_test.py", buildEnv: profile.toolchain)
  }

  // Step 4: Publish wheels to repositories
  if (!is_pr) {
    // For release builds, also publish to adobeshared
    echo "Publishing to adobeshared-release..."
    if (nodeUtils.shyIsUnix()) {
      cmd("uv publish --index adobeshared-release wheelhouse/*.whl")
    } else {
      cmd('for %%f in (wheelhouse\\*.whl) do uv publish --index adobeshared-release %%f')
    }
  } else {
    echo "Skipping adobeshared publish (PR build)"
  }

  // Always publish to tech-transfer (dev or release based on release_mode)
  echo "Publishing to ${tech_transfer_index}..."
  if (nodeUtils.shyIsUnix()) {
    cmd("uv publish --index ${tech_transfer_index} wheelhouse/*.whl")
  } else {
    cmd('for %%f in (wheelhouse\\*.whl) do uv publish --index ' + tech_transfer_index + ' %%f')
  }

  echo "Wheel operations completed successfully"
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