# Makefile for bit_util
# Generated by CMake Version 2.55
#    Build date Jan  4 2015 10:16:12
#    Ryan Koehler, ryan@verdascend.com
# Thu Jan 22 17:00:16 2015
#

TARGET = bit_util
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    autil.o    dfutil.o   mutil.o    putil.o    sutil.o     \
          util.o     sysinfo.o  numlist.o  wordlist.o table.o     \
          table_io.o table_str.o table_val.o bitpool.o  bit_util.o 

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

table.o:  table.c    prim.h     table.h    

table_io.o: table_io.c prim.h     table.h    

table_str.o: table_str.c prim.h     table.h    

table_val.o: table_val.c prim.h     table.h    

bitpool.o: bitpool.c  prim.h     bitpool.h  

bit_util.o: bit_util.c prim.h     bitpool.h  table.h    bit_util.h 
