#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR/src
find ./ -iname "*.pyc" -exec rm {} \;
cd ..
while true; do
    read -p "Install Rail-RNA in ~/rail-rna? (y/n): " yn
    case $yn in
        [Yy]* ) rm -rf ~/rail-rna; cp -rf ./src ~/rail-rna ; break;;
        [Nn]* ) exit;;
        * ) echo "Please enter y or n.";;
    esac
done
echo 'Installed Rail-RNA.'
echo 'If prompted, give permission to copy "rail-rna" executable to /usr/local/bin.'
sudo cat >/usr/local/bin/rail-rna <<EOF
#!/usr/bin/env bash

/usr/bin/env python ~/rail-rna \$@
EOF
sudo chmod +x /usr/local/bin/rail-rna
echo '"rail-rna" is now an executable in PATH. Enter "rail-rna -h" to get started.'