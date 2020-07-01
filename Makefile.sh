#/usr/bin/bash
#This file compiles all C++ scripts to executable programs
#Please make sure boost library has been installed

echo "This file compiles all C++ scripts to executable programs"
echo "Please make sure boost library has been installed"
echo "Start compling files..."

for file in *cpp
do
    echo "Compiling "$file
    g++ $file -o ${file%.*}
done
