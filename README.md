# mate_ba

[![Copier Badge][copier-badge]][copier-url]
[![Pixi Badge][pixi-badge]][pixi-url]
![License][license-badge]

Decorrelates data, allowing for accurate calculation of molecular brightness achieving better precision with fewer data points. This protocol of data analysis provides insights into the effects of temporal data correlation on molecular brightness calculations, a key parameter in fluorescence imaging and molecular dynamics. Temporal correlations can distort brightness estimates and significantly increase data volume and computational demands, we underscore the necessity of data decorrelation strategies to enhance accuracy and efficiency


## Development

This project is managed by [pixi][pixi-url].
You can install it for development using:

```sh
git clone https://github.com/EstradaLab-BiophotonicsGroup/mate-ba
cd mate-ba
pixi run pre-commit-install
```

Pre-commit hooks are used to lint and format the project.

## Usage

```sh
pixi run python '.\Code\Fig 1 and 2 - Molecular brightness using N data.py'
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
