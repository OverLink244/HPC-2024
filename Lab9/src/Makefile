Version:=winograd_V1.1.2
FastFlags:=-fopenmp -O3 -mfma -Ofast -march=x86-64 -mtune=generic
OriFlags:= -fopenmp -O3 
File:=small
Mode:=0
Data:=_data

all:
	gcc -std=c11  $(FastFlags) -g driver.c $(Version).c -o winograd -l openblas
	
	# gcc -std=c11 -D__DEBUG -fopenmp -O0 -g driver.c winograd.c -o winograd -l openblas

run:
	make
	export OMP_NUM_THREADS=8
	mkdir -p ./$(Data)
	./winograd $(File).conf $(Mode) > ./$(Data)/$(Version)_$(File).m