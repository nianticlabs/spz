
Setup your dev env with

    uv venv --python 3.12
    uv pip install -r tools/requirements.txt
    uv pip install -e ".[tests]" --upgrade

Prefix all your Python commands with `uv run` (or activate the venv in your shell)

    uv run python --version
