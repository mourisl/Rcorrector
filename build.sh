#!/bin/sh

tar -xzf jellyfish-2.1.3.tar.gz
cd jellyfish-2.1.3
./configure
make

cd ../
make
