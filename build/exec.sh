#!/bin/bash

make
./ensitsp $1 $2 $3 -s > fichiersvg.svg
if test $? != 3; then
    eog fichiersvg.svg &
fi