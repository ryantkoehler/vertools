# Makefile for lassoo
# Generated by CMake Version 2.55
#    Build date Mar  4 2016 11:08:00
#    Ryan Koehler, ryan@verdascend.com
# Fri Mar  4 11:15:29 2016
#

TARGET = lassoo
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    autil.o    dfutil.o   mutil.o    putil.o    sutil.o     \
          util.o     sysinfo.o  numlist.o  wordlist.o bitpool.o   \
          score.o    lassoo.o   

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}

autil.o:  autil.c    prim.h     

dfutil.o: dfutil.c   prim.h     

mutil.o:  mutil.c    prim.h     

putil.o:  putil.c    prim.h     

sutil.o:  sutil.c    prim.h     

util.o:   util.c     prim.h     

sysinfo.o: sysinfo.c  prim.h     

numlist.o: numlist.c  prim.h     numlist.h  

wordlist.o: wordlist.c prim.h     wordlist.h 

bitpool.o: bitpool.c  prim.h     bitpool.h  

score.o:  score.c    prim.h     score.h    

lassoo.o: lassoo.c   prim.h     score.h    bitpool.h  numlist.h   \
          lassoo.h   
