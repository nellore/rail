#!/usr/bin/env bash
CWD=$(pwd)
cd $(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd src
VER=$(python -c "import version; print version.version_number,")
INTERMEDIATE=rail-rna_installer.zip
rm -rf ../rail-rna*installer*
zip ../rail-rna_installer.zip $(find . -name \*.py | xargs)
cd ..
mkdir -p installers
TARGET=installers/rail-rna-${VER}_installer
rm -rf $TARGET
echo '#!/usr/bin/env python' | cat - $INTERMEDIATE >$TARGET
rm -rf $INTERMEDIATE
chmod +x $TARGET
cd $CWD