Scripts for building and uploading binary packages for Bowtie and Bowtie 2 for use with Rail-RNA.

The Building process uses the [Holy Build Box](https://hub.docker.com/r/phusion/holy-build-box-64/) image for [Docker](https://www.docker.com).  All binaries built this way ought to be quite portable across Linux distros.

The `master.sh` script drives the build process, ultimately placing two `bowtie*.zip` files in the current directory.

The `upload.sh` script uploads the zips to the Rail-RNA "requester pays" bucket on S3 and makes them publicly readable.
