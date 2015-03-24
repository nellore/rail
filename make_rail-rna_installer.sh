#!/usr/bin/env bash
CWD=$(pwd)
cd $(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd src
VER=$(python -c "import version; print version.version_number,")
INTERMEDIATE=rail-rna_installer.zip
rm -rf ../rail-rna*installer*
zip ../rail-rna_installer.zip $(find . -not -name \*.pyc -not -name .DS\_Store | xargs)
cd ..
mkdir -p releases
TARGET=releases/install_rail-rna-${VER}
rm -rf $TARGET
echo '#!/usr/bin/env python' | cat - $INTERMEDIATE >$TARGET
rm -rf $INTERMEDIATE
chmod 755 $TARGET
cd $CWD