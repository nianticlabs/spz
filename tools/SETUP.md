# Developer Setup

## Prerequisites

- **Git** installed and configured
- **Adobe VPN / Adobe Internal Network** connection

## Set up Artifactory Credentials

You need credentials for Cloud (UW2) Artifactory to access internal package indices.

- `ARTIFACTORY_UW2_USER` — your LDAP username
- `ARTIFACTORY_UW2_API_KEY` — fetch from https://artifactory-uw2.adobeitc.com/

To get your API key: click the link, select `SAML SSO` to log in, then click your username (upper right) > Edit Profile > API Key.

### Mac / Linux

Add to your shell profile (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
export ARTIFACTORY_UW2_USER=...
export ARTIFACTORY_UW2_API_KEY=...
```

Then add to `~/.netrc`:

```
machine artifactory-uw2.adobeitc.com
    login YOUR_LDAP_USERNAME
    password YOUR_UW2_API_KEY
```

### Windows

Set environment variables (add to your PowerShell profile for persistence):

```powershell
$env:ARTIFACTORY_UW2_USER="..."
$env:ARTIFACTORY_UW2_API_KEY="..."
```

Then add to `%USERPROFILE%\_netrc`:

```
machine artifactory-uw2.adobeitc.com
    login YOUR_LDAP_USERNAME
    password YOUR_UW2_API_KEY
```

## Install uv

uv is a fast Python package manager and dependency resolver.

**Mac / Linux:**

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

## Set up the project

```bash
uv venv --python 3.12
uv pip install -r tools/requirements.txt
uv pip install -e ".[tests]" --upgrade
```

## Running commands

Prefix all Python commands with `uv run` (or activate the `.venv` in your shell):

```bash
uv run python --version
```

## Running tests

```bash
uv run pytest
```
