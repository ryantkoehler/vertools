# Makefile for cmake
# Generated by CMake Version 2.55
#    Build date May 29 2016 18:40:42
#    Ryan Koehler, ryan@verdascend.com
# Sun Jun  5 11:03:01 2016
#

TARGET = cmake
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    putil.o    sutil.o    sysinfo.o  util.o     cmake.o    

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}

putil.o:  putil.c    prim.h     

sutil.o:  sutil.c    prim.h     

sysinfo.o: sysinfo.c  prim.h     

util.o:   util.c     prim.h     

cmake.o:  cmake.c    prim.h     
