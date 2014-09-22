#!/bin/sh

rm -f src.tar.xz
tar cJf src.tar.xz *.cpp *.h
scp src.tar.xz xeonphi:./z/mst
