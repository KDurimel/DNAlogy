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
  echo "32 bits detected"
elif [ $ARCH = "i486" ];then
  ARCH="32"
  echo "32 bits detected"
elif [ $ARCH = "i586" ];then
  ARCH="32"
  echo "32 bits detected"
elif [ $ARCH = "i686" ];then
  ARCH="32"
  echo "32 bits detected"
elif [ "$ARCH" = "ppc64 ppc" ];then
  ARCH="64"
  echo "64 bits detected"
elif [ "$ARCH" = "x86_64" ];then
  ARCH="64"
  echo "64 bits detected"

else
  echo "Unsuported Architecture, please try to install manually"
  exit 1
fi
# Install MiniConda (i.e lightest Conda) if needed
read -p "Do you already have Conda installed on your system?(y/n) " answer
case ${answer:0:1} in
    y|Y )
        echo "Ok"
    ;;
    * )
  echo "(Mini)Conda download will start in 5 seconds. Type CTRL+C for exit."
  sleep 5
  if [ "$ARCH" = "64" ]
          then
            wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
            bash Miniconda2-latest-Linux-x86_64.sh -y
            rm Miniconda2-latest-Linux-x86_64.sh*
        else
          wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86.sh
          bash Miniconda2-latest-Linux-x86.sh -y
          rm Miniconda2-latest-Linux-x86.sh*
  fi
    ;;
esac
# Generate an alias for the main script
read -p "Do you want to add fastGSEA (recommended) in your shell environment?(y/n) " answer
case ${answer:0:1} in
    y|Y )
      echo '# Added by fastGSEA installer' >> ~/.bashrc
      ABSPATH=$(readlink -f $0)
      ABSDIR=$(dirname $ABSPATH)
      chmod +x $ABSDIR/src/*
      echo 'export PATH=$PATH:'$ABSDIR'/src/' >> ~/.bashrc
      source ~/.bashrc
  echo "Ok, fastGSEA added and sourced."
  sleep 1
        # Create a new environment from packages.yml and activate it
  echo 'fastGSEA environment build will start in 5 seconds...'
  sleep 5
  conda env create -f $ABSDIR/packages.yml
  sleep 1
  echo 'fastGSEA sucessfully installed. Type "source activate gsea_env" to start using it.'
  exec bash # reload bashrc not only in this script context
    ;;
    * )
       echo "Ok, so please install packages.yml (look at the readme) and source gsea_env"
    ;;
esac