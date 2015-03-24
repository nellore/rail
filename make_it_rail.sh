#!/usr/bin/env bash
## Makes Rail-RNA package downloaded by EC2 nodes in Elastic MapReduce job
## and Rail-RNA installer executable. The former is placed in packages/ and
## the latter is placed in releases/. The version number in src/version.py
## is respected.
PACKAGES=packages
RELEASES=releases
CWD=$(pwd)
cd $(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd src
VER=$(python -c "import version; print version.version_number,")
# Change version number in rail-rna.txt
cd rna/driver
echo -ne "\xe2\x88\x80 Rail-RNA v${VER}\n" >rail-rna.txt
cd ../..
# Create installer
INTERMEDIATE=rail-rna_installer.zip
rm -rf ../rail-rna*installer*
zip ../rail-rna_installer.zip $(find . -not -name \*.pyc -not -name .DS\_Store | xargs)
cd ..
mkdir -p $RELEASES
TARGET=${RELEASES}/install_rail-rna-${VER}
rm -rf $TARGET
echo '#!/usr/bin/env python' | cat - $INTERMEDIATE >$TARGET
rm -rf $INTERMEDIATE
chmod 755 $TARGET
# Create package
ARNAME=rail-rna-${VER}.tar.gz
rm -rf $ARNAME
tar czvf ${ARNAME} --exclude '*.pyc' --exclude '*.tar.gz' --exclude '.DS_Store' src lib
mkdir -p $PACKAGES
mv $ARNAME ${PACKAGES}/$ARNAME
cd $PACKAGES
FULLPACK=$(pwd)
# Return to original dir
cd $CWD
echo "Run this to upload to S3."
echo "s3cmd put --acl-public ${FULLPACK}/${ARNAME} s3://rail-emr/bin/"