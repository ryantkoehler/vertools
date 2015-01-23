# Makefile for chardump
# Generated by CMake Version 2.55
#    Build date Jan  4 2015 10:16:12
#    Ryan Koehler, ryan@verdascend.com
# Thu Jan 22 17:00:04 2015
#

TARGET = chardump
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    putil.o    sutil.o    sysinfo.o  util.o     chardump.o 

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}

putil.o:  putil.c    prim.h     

sutil.o:  sutil.c    prim.h     

sysinfo.o: sysinfo.c  prim.h     

util.o:   util.c     prim.h     

chardump.o: chardump.c prim.h     
