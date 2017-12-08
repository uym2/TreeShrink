# ! /bin/bash

echo "Installing BMS package in R"

R CMD INSTALL dependencies/BMS_0.3.3.tar.gz


u=""

while [[ $# -gt 1 ]]; do
    key="$1"

case $key in
-U|--user)
u="--user"
shift # past argument
;;
esac
shift
done

echo "Installing treeshrink"


python setup.py develop $u
