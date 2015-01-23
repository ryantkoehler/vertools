# Makefile for scoretab
# Generated by CMake Version 2.55
#    Build date Jan  4 2015 10:16:12
#    Ryan Koehler, ryan@verdascend.com
# Thu Jan 22 17:00:00 2015
#

TARGET = scoretab
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    autil.o    dfutil.o   mutil.o    putil.o    sutil.o     \
          util.o     sysinfo.o  numlist.o  wordlist.o score.o     \
          sctab_ga.o score_io.o sctab_mix.o sctab_part.o sctab_sc.o  \
          sctab_str.o stat.o     histogram.o mut_info.o table.o     \
          table_io.o table_ops.o table_stat.o table_str.o table_val.o  \
          table_wf.o scoretab.o 

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

score.o:  score.c    prim.h     score.h    

sctab_ga.o: sctab_ga.c prim.h     numlist.h  wordlist.h table.h     \
          scoretab.h 

score_io.o: score_io.c prim.h     score.h    

sctab_mix.o: sctab_mix.c prim.h     numlist.h  wordlist.h score.h     \
          table.h    scoretab.h 

sctab_part.o: sctab_part.c prim.h     score.h    table.h    scoretab.h 

sctab_sc.o: sctab_sc.c prim.h     score.h    table.h    scoretab.h 

sctab_str.o: sctab_str.c prim.h     score.h    table.h    stat.h      \
          scoretab.h 

stat.o:   stat.c     prim.h     score.h    numlist.h  stat.h     

histogram.o: histogram.c prim.h     numlist.h  stat.h     histogram.h 

mut_info.o: mut_info.c prim.h     numlist.h  stat.h     table.h     \
          histogram.h mut_info.h 

table.o:  table.c    prim.h     table.h    

table_io.o: table_io.c prim.h     table.h    

table_ops.o: table_ops.c prim.h     table.h    score.h    

table_stat.o: table_stat.c prim.h     table.h    

table_str.o: table_str.c prim.h     table.h    

table_val.o: table_val.c prim.h     table.h    

table_wf.o: table_wf.c prim.h     table.h    

scoretab.o: scoretab.c prim.h     score.h    stat.h     table.h     \
          mut_info.h scoretab.h 
