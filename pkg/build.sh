#!/bin/bash
cd ..
set -o errexit
rm -rf fpop-release
cp -r pkg fpop-release
PKG_TGZ=$(R CMD build fpop-release|grep building|sed "s/.*\(fpop.*.tar.gz\).*/\1/")
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
