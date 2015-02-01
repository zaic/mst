#!/bin/sh

rm -f src.tar.xz
tar cJf src.tar.xz *.cpp *.h third_party_generator/*.h third_party_generator/*.cpp third_party_generator/Makefile
scp src.tar.xz xeonphi:./z/mst
