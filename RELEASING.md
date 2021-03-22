# Release LingPy

## Release Types
  
* a major release is one where the number after the dot increases (2.5 to 2.6), major releases are accompanied by a new release of the documentation at http://lingpy.org
* a minor release is one where the second number increases (2.5 to 2.5.1), minor releases may be accompanied by an update of the documentation at https://lingpy.github.io

## Release Procedure

- Do platform test via tox:
  ```
  tox -r
  ```

- Make sure docs can be built:
  ```
  cd ..
  git clone https://github.com/lingpy/doc
  cd doc/
  sphinx-apidoc -f ../lingpy/src/lingpy -o source/reference
  make html
  cd ../lingpy
  ```

- Update the version number, by removing the trailing `.dev0` in:
  - `setup.py`
  - `src/lingpy/__init__.py`

- Create the release commit:
  ```shell
  git commit -a -m "release <VERSION>"
  ```

- Create a release tag:
  ```
  git tag -a v<VERSION> -m"<VERSION> release"
  ```

- Release to PyPI (see https://github.com/di/markdown-description-example/issues/1#issuecomment-374474296):
  ```shell
  rm dist/*
  python setup.py sdist
  twine upload dist/*
  rm dist/*
  python setup.py bdist_wheel
  twine upload dist/*
  ```

- Push to github:
  ```
  git push origin
  git push --tags
  ```

- Change version for the next release cycle, i.e. incrementing and adding .dev0
  - `setup.py`
  - `src/lingpy/__init__.py`

- Commit/push the version change:
  ```shell
  git commit -a -m "bump version for development"
  git push origin
  ```


### Major Release

* make sure the checks for the last PR have passed and coverage for tests has not decreased
* make sure the documentation is up to date
* release the code via GitHub
* modify the metadata on Zenodo
* upload data to PyPi
* announce the release in different channels (LINGUISTList, etc.)

### Minor Release

* make sure the checks for the last PR have passed and coverage for tests has not decreased
* make sure the documentation is up to date
* release the code via GitHub
* modify the metadata on Zenodo
* upload data to PyPi

## Major Releases Planned

With version 2.6, we plan to release 2.7 in 2021(?) with some new functions,
including new approaches from new future members on cognate detection and some
basic approaches to phylogenetic reconstruction which we did not have yet the
time to include. You can check out what we plan by comparing the [Milestone
2.7](https://github.com/lingpy/lingpy/issues?q=is%3Aopen+is%3Aissue+milestone%3Alinpgy2.7)
on GitHub.

With version 2.7 we will start working on a complete relaunch of LingPy which
will no longer guarantee backwards-compatibility and which will also reduce
certain functions which LingPy offered in the past, as they will be done by
alternative software packages. This new package, which we would ideally
finalize in 2022 will be called LingPy3 and appear as version 3.0.
