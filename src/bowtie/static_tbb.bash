#!/bin/bash

set -e

if [[ $# -lt 1 ]]; then
    echo "$0 bowtie|bowtie2 [branch]"
    exit 1
fi

bowtie_version=$1
bowtie_branch=$2
bowtie_commit=$3

source /hbb_exe/activate
mkdir /mybin
echo  'res=`echo $@ | sed "s/-L.*$//"`; echo $res; echo $res; /opt/rh/devtoolset-2/root/usr/bin/ar $res;' > /mybin/ar
chmod +x /mybin/ar && export PATH=/mybin:$PATH

set -x

TBB_ROOT=/io
TBB_LIB=/io/lib
if [[ ! -d /io/lib || ! -d /io/include ]]; then
    curl -LO https://github.com/01org/tbb/archive/2017_U5.tar.gz
    TBB_ROOT="/`tar -tzf /2017_U5.tar.gz | head -1`"
    tar xzf 2017_U5.tar.gz && pushd $TBB_ROOT && make extra_inc=big_iron.inc && popd
    TBB_LIB=`ls -d $TBB_ROOT/build/linux*release`
    cp -r $TBB_ROOT/include /io
    mkdir io/lib && cp $TBB_LIB/*.a /io/lib
fi

export CPLUS_INCLUDE_PATH=$TBB_ROOT/include
export LIBRARY_PATH=$TBB_LIB:$LIBRARY_PATH
export LD_LIBRARY_PATH=$TBB_LIB:$LD_LIBRARY_PATH

yum -y install git zip unzip

[[ -e /io/$bowtie_version ]] && echo "Error, clone already there" && exit 1
git clone https://github.com/BenLangmead/$bowtie_version.git -- $bowtie_version
pushd $bowtie_version
[[ -n "${bowtie_commit}" ]] && git reset --hard "${bowtie_commit}"
makefile_var="LIBS"
if [[ $bowtie_version == "bowtie" ]]; then makefile_var=EXTRA_FLAGS; fi
if [[ ! -z $bowtie_branch ]]; then git checkout $bowtie_branch; fi
sed -e "/^${makefile_var}/s/${makefile_var} =\(.*\)/${makefile_var} = $\(LDFLAGS\) \1/"\
    -e s/tbbmalloc_proxy/tbbmalloc/g -i'' Makefile
cp -r ../$bowtie_version /io

make ${bowtie_version}-`[ $bowtie_version == "bowtie2" ] && echo "pkg" || echo "bin.zip"`
libcheck ${bowtie_version}-{align,build,inspect}-*
cp *.zip /io
rm -rf /io/$bowtie_version /io/include /io/lib
