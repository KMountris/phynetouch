#!/bin/sh

# Instructions for building gcc 8.x from source.

# This gcc build script is free software; you can redistribute it and/or modify
# it under the terms of the MIT license.


#======================================================================
# User configuration
#======================================================================

# Provide the version of gcc being built (e.g. 9.1.0)
#
gcc_version=11

# Additional makefile options.  E.g., "-j 4" for parallel builds.  Parallel
# builds are faster, however it can cause a build to fail if the project
# makefile does not support parallel build.
#make_flags="-j 2"
make_flags=""

# Source and build directories for GCC
#
root_dir=$(pwd)
source_dir=${root_dir}/gcc-${gcc_version}-source
build_dir=${root_dir}/gcc-${gcc_version}-build
install_dir=${root_dir}/gcc-${gcc_version}-install

arch_flags="-march=x86-64"
build_target=x86_64-unknown-linux-gnu


# Set the languages that you want GCC to support
# list of comma seperated entries: e.g., c,c++
# available: ada,c,c++,d,fortran,go,jit,lto,objc,obj-c++
#
lang=c,c++,fortran


#======================================================================
# Modify only if you know what you are doing
#======================================================================

packageversion="$(whoami)-$(hostname -s)"


# Get GCC source files for specified version
git clone git://gcc.gnu.org/git/gcc.git ${source_dir}


# Download the GCC compilation prerequisites
cd ${source_dir}
git checkout remotes/origin/releases/gcc-${gcc_version}
./contrib/download_prerequisites



#======================================================================
# Directory creation
#======================================================================


# Compile GCC
mkdir "$install_dir"
mkdir "$build_dir"
cd "${build_dir}"

CC=gcc
CXX=g++
OPT_FLAGS="-O2 $gflags -Wall  $arch_flags"
CC="$CC" CXX="$CXX" CFLAGS="$OPT_FLAGS" \
    CXXFLAGS="`echo " $OPT_FLAGS " | sed 's/ -Wall / /g'`" \
    $source_dir/configure --prefix=${install_dir} \
    --enable-languages=${lang} \
    --disable-multilib \
    --enable-bootstrap \
    --enable-shared \
    --enable-threads=posix \
    --enable-checking=release \
    --with-system-zlib \
    --enable-__cxa_atexit \
    --disable-libunwind-exceptions \
    --enable-linker-build-id \
    --disable-vtable-verify \
    --with-default-libstdcxx-abi=new \
    --enable-libstdcxx-debug  \
    --without-included-gettext  \
    --enable-plugin \
    --disable-initfini-array \
    --disable-libgcj \
    --enable-plugin  \
    --with-tune=generic \
    --build=${build_target} \
    --target=${build_target} \
    --host=${build_target} \
    --with-pkgversion="$packageversion"


make BOOT_CFLAGS="$OPT_FLAGS" ${make_flags} bootstrap
make check
make install


# Create a shell script that users can source to bring gcc into shell
# environment
cat << EOF > ${install_dir}/activate-gcc-${gcc_version}
# source this script to bring gcc ${gcc_version} into your environment
export PATH=${install_dir}/bin:\$PATH
export LD_LIBRARY_PATH=${install_dir}/lib:${install_dir}/lib64:\$LD_LIBRARY_PATH
export MANPATH=${install_dir}/share/man:\$MANPATH
export INFOPATH=${install_dir}/share/info:\$INFOPATH
EOF

#end
