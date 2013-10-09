#!/bin/csh

rm output
mkdir foobar
rm foobar/*

condor_submit condorscript
