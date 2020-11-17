#!/bin/bash
curl -s https://api.github.com/repos/sophiasage/cn_hyperarr | grep "Not Found" > /dev/null || { echo >&2 "The repository cn_hyperarr already exists in sophiasage.  Aborting."; exit 1; }

echo "Trying to create a new repository on github.com."
echo "You will be asked for the GitHub password corresponding to the user sophiasage".
echo "Instead of your password, you MUST use a Personal Access Token !!!"
echo "( See here how to get one:"
echo "  https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line )"
echo "sophiasage/cn_hyperarr"
echo "Computation of congruence normality for hyperplane arrangements"

curl -s -u 'sophiasage' https://api.github.com/user/repos -d '{"name":"cn_hyperarr","description":"Computation of congruence normality for hyperplane arrangements"}' >> install.log

echo "Repository successfully created."
