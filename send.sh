#!/bin/sh

rm -f src.tar.xz
tar cJf src.tar.xz *.cpp *.h third_party_generator
scp src.tar.xz xeonphi:./z/mst
