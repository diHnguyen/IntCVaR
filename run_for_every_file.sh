#!/bin/bash

for i in {1..20}; 
do
    /Applications/Julia-1.0.app/Contents/Resources/julia/bin/julia /Users/din/Desktop/IntCVaR/BasicAlg_4min.jl $i
    #echo "Hello $I"
    # echo myCoolCommand $file $file.result;
    # change "echo myCoolCommand ..." to whatever you want to run
    # like, "mauve $file $file.result;"
done