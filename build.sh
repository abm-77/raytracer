#!/bin/bash
clear
echo "Building..."
if [ ! -d ./build ]
then
    mkdir ./build
fi

cd ./build

clang++ -std=c++14 -g -o ../bin/application ../source/main.cpp 


echo "Done!"