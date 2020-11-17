#!/bin/bash
command -v git >/dev/null 2>&1 || { echo >&2 "I require git but it's not installed.  Aborting."; exit 1; }


command -v travis >/dev/null 2>&1 || { echo >&2 "I require travis (https://github.com/travis-ci/travis.rb) but it's not installed.  Aborting."; exit 1; }
echo "Next the deploy key will be encrypted for Travis CI to use. Also we will enable the newly created repository in Travis CI."
travis login -u sophiasage
travis encrypt-file .travis_ci_gh_pages_deploy_key --add before_script -r sophiasage/cn_hyperarr
travis enable -r sophiasage/cn_hyperarr
echo "Done!"

echo "Finally we add the encrypted key and commit changes to GitHub. You may be asked for your Personal Access Token one last time."

git add .travis_ci_gh_pages_deploy_key.enc
git add .travis.yml
git commit -m "Deploy built documentation to GitHub"
git push

