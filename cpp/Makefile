CC=g++
#CPPFLAGS  = -O0 -g -Wall -fPIC -ftree-vectorize -std=c++11
CPPFLAGS=-std=c++11 -O3 -ffast-math -ftree-vectorize 
CPPLFLAGS= $(CPPFLAGS) -shared -Wl,-soname

all: STM

STM:
	$(CC) $(CPPFLAGS) -c -fPIC ProbeSTM_spd.cpp
	$(CC) $(CPPLFLAGS),ProbeSTM_spd_lib.so -o ProbeSTM_spd_lib.so ProbeSTM_spd.o
clean:
	rm *.o *.so *.pyc
