============
Contributing
============

Contributions are welcome.

Submit contributions as a pull request on GitHub.


Standards
---------
- All contributions must be PEP complaint.
- All modules, classes, and functions must have Sphinx style doc strings.
- All parameters and variables should have type annotations.
- A reasonable amount of code block documentation is expected as well.

Tests should be written first before implementing functionality. See https://en.wikipedia.org/wiki/Test-driven_development


Repo Maintainers
----------------

To create a new release for pypi and conda:

1. Clone this repository
2. Make sure all tests pass first by running :code:`./setup.py test`
3. Use git or the GitHub interface to tag the repository with the next version. Precede the version number with a 'v'. Follow symantic versioning rules https://semver.org/:

    Given a version number MAJOR.MINOR.PATCH, increment the:

    MAJOR version when you make incompatible API changes,

    MINOR version when you add functionality in a backwards compatible manner, and

    PATCH version when you make backwards compatible bug fixes.

4. Create a release via the GitHub interface for the new tag.
5. Build the package by running :code:`./setup.py build sdist`
6. Upload the new package to pypi by running :code:`twine upload -u brinkmanlab dist/<file with latest version number>.tar.gz`

Bioconda will detect the upload and automatically update its recipe. A Bioconda maintainer will manually approve the update.

