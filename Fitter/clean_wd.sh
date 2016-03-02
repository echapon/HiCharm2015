#!/bin/bash

wd=$PWD

cmd="rm -rf ${wd}/Output ${wd}/DataSet"
echo $cmd
$cmd
