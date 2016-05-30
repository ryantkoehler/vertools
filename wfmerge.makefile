# Makefile for wfmerge
# Generated by CMake Version 2.55
#    Build date Mar  4 2016 11:08:00
#    Ryan Koehler, ryan@verdascend.com
# Sun May 29 18:40:48 2016
#

TARGET = wfmerge
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    autil.o    dfutil.o   mutil.o    putil.o    sutil.o     \
          sysinfo.o  util.o     wfutil.o   wfmerge.o  

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}

autil.o:  autil.c    prim.h     

dfutil.o: dfutil.c   prim.h     

mutil.o:  mutil.c    prim.h     

putil.o:  putil.c    prim.h     

sutil.o:  sutil.c    prim.h     

sysinfo.o: sysinfo.c  prim.h     

util.o:   util.c     prim.h     

wfutil.o: wfutil.c   prim.h     score.h    wfutil.h   wordfreq.h 

wfmerge.o: wfmerge.c  prim.h     dna.h      wfutil.h   wfmerge.h  
