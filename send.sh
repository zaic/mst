#!/bin/sh

rm -f src.tar.bz2
tar cjf src.tar.bz2 makefile makefile.h *.cpp *.h
scp src.tar.bz2 xeonphi:./z/mst
