#!/bin/bash

wd=$PWD

cmd="rm -rf ${wd}/DBFiles ${wd}/Plots"
echo $cmd
$cmd
