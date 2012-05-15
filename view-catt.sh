#!/bin/sh
DIR=$PWD
which gradsc >& /dev/null
if [ $? -ne 0 ] ; then
        echo " grads is not installed or not in PATH"
        exit 8
fi
gradsc -clb "run figcatt.gs" >& grads.out
if [ $? -ne 0 ] ; then
        exit 8
fi
\rm -rf co pm25 temp topo src1
