# mate_ba

[![Copier Badge][copier-badge]][copier-url]
[![Pixi Badge][pixi-badge]][pixi-url]
![License][license-badge]
[![CI Badge][ci-badge]][ci-url]
[![conda-forge Badge][conda-forge-badge]][conda-forge-url]
[![PyPI Badge][pypi-badge]][pypi-url]
[![Python version Badge][pypi-version-badge]][pypi-version-url]

Decorrelates data, allowing for accurate calculation of molecular brightness achieving better precision with fewer data points. This protocol of data analysis provides insights into the effects of temporal data correlation on molecular brightness calculations, a key parameter in fluorescence imaging and molecular dynamics. Temporal correlations can distort brightness estimates and significantly increase data volume and computational demands, we underscore the necessity of data decorrelation strategies to enhance accuracy and efficiency

## Install

Using [pixi][pixi-url],
install from PyPI with:

```sh
pixi add --pypi mate-ba
```

or install the latest development version from GitHub with:

```sh
pixi add --pypi mate-ba@https://github.com/EstradaLab-BiophotonicsGroup/mate-ba.git
```

Otherwise,
use `pip` or your `pip`-compatible package manager:

```sh
pip install mate-ba  # from PyPI
pip install git+https://github.com/EstradaLab-BiophotonicsGroup/mate-ba.git  # from GitHub
```

## Development

This project is managed by [pixi][pixi-url].
You can install it for development using:

```sh
git clone https://github.com/EstradaLab-BiophotonicsGroup/mate-ba
cd mate-ba
pixi run pre-commit-install
```

Pre-commit hooks are used to lint and format the project.

### Testing

Run tests using:

```sh
pixi run test
```

### Publishing to PyPI

When a tagged commit is pushed to GitHub,
the GitHub Action defined in `.github/workflows/ci.yml`
builds and publishes the package to PyPI.

Tag a commit and push the tags with:

```sh
git tag <my-tag>
git push --tags
```

Trusted publishing must be enabled once in [PyPI Publishing](https://pypi.org/manage/account/publishing/).
Fill the following values in the form:

```
PyPI Project Name: mate-ba
            Owner: EstradaLab-BiophotonicsGroup
  Repository name: mate-ba
    Workflow name: ci.yml
 Environment name: pypi
```

[ci-badge]: https://img.shields.io/github/actions/workflow/status/EstradaLab-BiophotonicsGroup/mate-ba/ci.yml
[ci-url]: https://github.com/EstradaLab-BiophotonicsGroup/mate-ba/actions/workflows/ci.yml
[conda-forge-badge]: https://img.shields.io/conda/vn/conda-forge/mate-ba?logoColor=white&logo=conda-forge
[conda-forge-url]: https://prefix.dev/channels/conda-forge/packages/mate-ba
[copier-badge]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/copier-org/copier/master/img/badge/badge-black.json
[copier-url]: https://github.com/copier-org/copier
[license-badge]: https://img.shields.io/badge/license-GNU-blue
[pixi-badge]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json
[pixi-url]: https://pixi.sh
[pypi-badge]: https://img.shields.io/pypi/v/mate-ba.svg?logo=pypi&logoColor=white
[pypi-url]: https://pypi.org/project/mate-ba
[pypi-version-badge]: https://img.shields.io/pypi/pyversions/mate-ba?logoColor=white&logo=python
[pypi-version-url]: https://pypi.org/project/mate-ba
