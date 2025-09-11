#!/bin/bash

# Setup pypi index urls
export PIP_EXTRA_INDEX_URL="https://$ARTIFACTORY_UW2_USER:$ARTIFACTORY_UW2_API_KEY@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-tech-transfer-3di-release/simple \
https://$ARTIFACTORY_UW2_USER:$ARTIFACTORY_UW2_API_KEY@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-adobeshared-release/simple"

# ask conda to show folder name instead of tag or abs path for envs in non default location.
conda config --set env_prompt '({name}) '

# you can set debug conda setup issues by running this way:
#   CONDA_DEBUG_FLAG=-vvv . tools/setup.sh
# Optional first argument: path to conda.yaml. Defaults to tools/conda.yaml
CONDA_FILE_PATH="${1:-tools/conda.yaml}"

conda env update -q $CONDA_DEBUG_FLAG --prefix .venv/ --file "$CONDA_FILE_PATH" || exit -1
conda activate .venv/ || exit -1
pip install -e ".[tests]" --upgrade
conda activate .venv/ || exit -1
echo "[setup] done!"

