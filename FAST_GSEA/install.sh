#!/bin/bash
# Not really added in shell environment
# Retrieve kernel name
KERNEL=$(uname -s)
if [ $KERNEL = "Linux" ]; then
  echo "Checking for CPU architecture"
else
  echo "Unsupported OS"
fi
# Retrieve CPU architecture
ARCH=$(uname -m)
if [ $ARCH = "i386" ];then
  ARCH="32"
elif [ $ARCH = "i486" ];then
  ARCH="32"
elif [ $ARCH = "i586" ];then
  ARCH="32"
elif [ $ARCH = "i686" ];then
  ARCH="32"
elif [ "$ARCH" = "ppc64 ppc" ];then
  ARCH="64"
elif [ "$ARCH" = "x86_64" ];then
  ARCH="64"

else
  echo "Unsoported Architecture, please try to install manually"
  exit 1
fi
# Install MiniConda (i.e lightest Conda) if needed
read -p "Do you already have Conda installed on your system?(y/n)? " answer
case ${answer:0:1} in
    y|Y )
        echo "Ok"
    ;;
    * )
	if [ "$ARCH" = "64" ]
          then
            wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
            bash Miniconda2-latest-Linux-x86_64.sh
            rm Miniconda2-latest-Linux-x86_64.sh*
        else
          wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86.sh
          bash Miniconda2-latest-Linux-x86.sh
          rm Miniconda2-latest-Linux-x86.sh*
	fi
    ;;
esac
# Generate an alias for the main script
echo 'alias fastGSEA="python $(pwd)'/src/fastGSEA.py'"' >> ~/.bashrc
# Create a new environment from packages.yml and activate it
conda env create -f packages.yml
