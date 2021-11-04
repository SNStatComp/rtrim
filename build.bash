#!/bin/bash

# Bash script to build the rtrim R package
# version 1: Mark van der Loo
# version 2: Patrick Bogaart
# - added generation of install script 

R=R
CHECKARG="--as-cran"
while [ $# -gt 0 ] ; do
  case "$1" in 
    -dev)
       R=Rdev
       shift 1 ;;
    *)
       CHECKARG="$CHECKARG $1"
       shift 1 ;;
  esac
done

echo "######## Removing building information..."
rm -rf output


echo "######## Generate documentation..."
$R -q -f roxygen.R

echo "######## Building package in output..."
mkdir output
cd output
$R CMD build ../pkg
echo "######## Testing package with $CHECKARG ..."
for x in *.tar.gz 
do 
    $R CMD check $CHECKARG $x
done

echo "######## Creating installation script..."
TARGET=install.R
PRE="install.packages(\""
TGZ=$(eval ls -1 *.tar.gz | head -1)
POST="\", repos=NULL, type=\"source\")"
echo "# R script to install rtrim package from source" > $TARGET
echo $PRE$TGZ$POST >> $TARGET

echo "**BUILT USING $R"
$R --version



