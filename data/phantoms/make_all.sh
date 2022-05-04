#!/bin/sh

make phantoms

./phantoms $@ --mode=diamond
./phantoms $@ --mode=sphere
./phantoms $@ --mode=cube
./phantoms $@ --mode=pyramid
./phantoms $@ --mode=bell
./phantoms $@ --mode=shepp-logan
./phantoms $@ --mode=shepp-logan_bold
./phantoms $@ --mode=stark
./phantoms $@ --mode=dorn
