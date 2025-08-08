# Setup pypi index urls
$env:PIP_EXTRA_INDEX_URL="https://" + $env:ARTIFACTORY_UW2_USER + ":" + $env:ARTIFACTORY_UW2_API_KEY + "@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-tech-transfer-3di-release/simple " + `
"https://" + $env:ARTIFACTORY_UW2_USER + ":" + $env:ARTIFACTORY_UW2_API_KEY + "@artifactory-uw2.adobeitc.com/artifactory/api/pypi/pypi-adobeshared-release/simple "

#Activate conda env
conda config --set ssl_verify no
# ask conda to show folder name instead of tag or abs path for envs in non default location.
conda config --set env_prompt '({name}) '
if( -not $? ) { exit -1 }

# you can set debug conda setup issues by running this way:
#   $ENV:CONDA_DEBUG_FLAG="-vvv"\n& tools\setup.ps1
conda env update $ENV:CONDA_DEBUG_FLAG --prefix .venv\ --file "tools\conda.yaml"
if( -not $? ) { exit -1 }
conda activate .venv\
if( -not $? ) { exit -1 }

# Install package in dev mode
pip install -U -e ".[tests]"
if( -not $? ) { exit -1 }


conda activate .venv\
if( -not $? ) { exit -1 }
