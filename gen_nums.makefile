# Makefile for gen_nums
# Generated by CMake Version 2.55
#    Build date Dec 11 2017 19:22:26
#    Ryan Koehler, ryan@verdascend.com
# Tue Jun  5 15:31:47 2018
#

TARGET = gen_nums
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    autil.o    dfutil.o   mutil.o    putil.o    sutil.o     \
          sysinfo.o  util.o     gen_nums.o 

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}

autil.o:  autil.c    prim.h     

dfutil.o: dfutil.c   prim.h     

mutil.o:  mutil.c    prim.h     

putil.o:  putil.c    prim.h     

sutil.o:  sutil.c    prim.h     

sysinfo.o: sysinfo.c  prim.h     

util.o:   util.c     prim.h     

gen_nums.o: gen_nums.c prim.h     gen_nums.h 
