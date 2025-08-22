#!/bin/bash

# ask conda to show folder name instead of tag or abs path for envs in non default location.
conda config --set env_prompt '({name}) '

# you can set debug conda setup issues by running this way:
#   CONDA_DEBUG_FLAG=-vvv . tools/setup.sh
conda env update -q $CONDA_DEBUG_FLAG --prefix .venv/ --file "tools/conda.yaml" || exit -1
conda activate .venv/ || exit -1
pip install -e ".[tests]" --upgrade
conda activate .venv/ || exit -1
echo "[setup] done!"

