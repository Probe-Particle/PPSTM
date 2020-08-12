# C++ compilers: CC - linux; MCC - Mac OS - we recommend to update g++ compiler on your MAC#
CC=g++
MCC=g++-9
# flagss Linux; OMPFLAGS used only for OMP compiltion of #
OMPFLAGS=-fopenmp
CPPFLAGS=-std=c++11 -O3 -ffast-math -ftree-vectorize 
CPPLFLAGS= $(CPPFLAGS) -shared -Wl,-soname
# flags Max Os #
MCPPFLAGS=-std=c++11 -O3 -ffast-math -ftree-vectorize
MCPPLFLAGS= $(MCPPFLAGS) -shared -Wl,-install_name


all: STM IO

STM:
	$(CC) $(CPPFLAGS) -c -fPIC ProbeSTM_spd.cpp
	$(CC) $(CPPLFLAGS),ProbeSTM_spd_lib.so -o ProbeSTM_spd_lib.so ProbeSTM_spd.o
STMPAR:
	$(CC) $(OMPFLAGS) $(CPPFLAGS) -c -fPIC ProbeSTM_spd.cpp
	$(CC) $(OMPFLAGS) $(CPPLFLAGS),ProbeSTM_spd_lib.so -o ProbeSTM_spd_lib.so ProbeSTM_spd.o
MSTM:
	$(MCC) $(MCPPFLAGS) -c -fPIC ProbeSTM_spd.cpp
	$(MCC) $(MCPPLFLAGS),ProbeSTM_spd_lib.so -o ProbeSTM_spd_lib.so ProbeSTM_spd.o
MSTMPAR:
	$(MCC) $(OMPFLAGS) $(MCPPFLAGS) -c -fPIC ProbeSTM_spd.cpp
	$(MCC) $(OMPFLAGS) $(MCPPLFLAGS),ProbeSTM_spd_lib.so -o ProbeSTM_spd_lib.so ProbeSTM_spd.o

IO:
	$(CC) $(CPPFLAGS) -c -fPIC IO.cpp
	$(CC) $(CPPLFLAGS),IO_lib.so -o IO_lib.so IO.o
MIO:
	$(MCC) $(MCPPFLAGS) -c -fPIC IO.cpp
	$(MCC) $(MCPPLFLAGS),IO_lib.so -o IO_lib.so IO.o
clean:
	rm *.o *.so *.pyc
