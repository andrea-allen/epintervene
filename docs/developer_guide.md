# Developing EpIntervene
Guide for contributing features and updating the Python package.

## First: your local codebase
1. Make sure you have your own up-to-date version of the codebase locally, via `git clone https://github.com/andrea-allen/epintervene.git` to start, or if you already have the codebase locally, pull the latest changes down to your local environment.
2. To pull down the latest changes use `git pull --rebase` to ensure any local changes you have avoid conflicts

## Changing the code
1. Change the code to your liking
2. Use `sandbox.py`, `interventions_sandbox.py`, or another personal local file to test out your new features
3. When you're happy, now it's time to push the code and update the package:

## Pushing your new code and updating the package
#### Note: Updating the package will update the software for everyone who uses it and downloads the new version. Make sure you've thoroughly tested your code, experimented with it, and are happy with the new changes.
1. Update the `setup.py` and `setup.cfg` files with a new version:
   1. Set the line `VERSION = 'X.Y.z'` by appropriately updating the version number: If a change is a new MAJOR version (rare) bump `X` by one. If there is a new MINOR version (a new feature) bump `Y` by one. If a small change or fix is made that does not hinder existing functions or add new ones, bump `z` by one.
2. Commit your code locally:
```
git add your_new_or_changed_file
git add setup.py
git add setup.cfg
git commit -m "Short one-line description of your new changes"
```
2. Push your new commit with the new version
```
git push
```
3. This version will now be pushed to the remote repo and it will go through the pipeline. Check out the `Actions` tab on the Github repository to watch the pipeline progress. The code will be published to the TestPyPi index but NOT the PyPi index (the actual package) yet. HOWEVER, make sure the pipeline run is GREEN before continuing! If it's not green, check out the error in the logs and see what went wrong. Go back to step 1 to fix things.
4. TAGGING to update the package: As the last step, you need to create a local tag for your new version, and push that tag to the remote in order to trigger a new release to the `epintervene` package. To make a new tag, in the command line
```
git tag -a vX.Y.z -m "vX.Y.z"
```
matching the exact version number you JUST updated the code to. Then, you need to PUSH the tag. To do so,
```
git push origin vX.Y.z
```
with the name of your new tag.

Congrats! You've now updated the package. Make sure to keep an eye on https://github.com/andrea-allen/epintervene/actions to make sure the pipeline is triggered and the release is green. Then go to https://pypi.org/project/epintervene/ and see that your latest version has been published.