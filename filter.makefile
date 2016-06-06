# Makefile for filter
# Generated by CMake Version 2.55
#    Build date May 29 2016 18:40:42
#    Ryan Koehler, ryan@verdascend.com
# Sun Jun  5 11:02:51 2016
#

TARGET = filter
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    autil.o    dfutil.o   mutil.o    putil.o    sutil.o     \
          sysinfo.o  util.o     numlist.o  score.o    stat.o      \
          wordlist.o filter.o   

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}

autil.o:  autil.c    prim.h     

dfutil.o: dfutil.c   prim.h     

mutil.o:  mutil.c    prim.h     

putil.o:  putil.c    prim.h     

sutil.o:  sutil.c    prim.h     

sysinfo.o: sysinfo.c  prim.h     

util.o:   util.c     prim.h     

numlist.o: numlist.c  prim.h     numlist.h  

score.o:  score.c    prim.h     score.h    

stat.o:   stat.c     prim.h     score.h    numlist.h  stat.h     

wordlist.o: wordlist.c prim.h     wordlist.h 

filter.o: filter.c   prim.h     wordlist.h filter.h   
