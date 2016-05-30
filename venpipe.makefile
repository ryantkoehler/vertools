# Makefile for venpipe
# Generated by CMake Version 2.55
#    Build date Mar  4 2016 11:08:00
#    Ryan Koehler, ryan@verdascend.com
# Sun May 29 18:40:54 2016
#

TARGET = venpipe
CC = gcc
CFLAGS = -Wall -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
SYSLIBS = -lm
OBJS =    autil.o    dfutil.o   mutil.o    putil.o    sutil.o     \
          sysinfo.o  util.o     dna.o      dna_char.o dna_file.o  \
          dna_in.o   dna_out.o  dna_seqs.o dna_snp.o  snp_char.o  \
          ven_engy.o ven_io.o   ven_salt.o ven_str.o  venpipe.o   \
          libRNA.a   

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) ${SYSLIBS}

autil.o:  autil.c    prim.h     

dfutil.o: dfutil.c   prim.h     

mutil.o:  mutil.c    prim.h     

putil.o:  putil.c    prim.h     

sutil.o:  sutil.c    prim.h     

sysinfo.o: sysinfo.c  prim.h     

util.o:   util.c     prim.h     

dna.o:    dna.c      prim.h     dna.h      

dna_char.o: dna_char.c prim.h     dna.h      

dna_file.o: dna_file.c prim.h     dna.h      

dna_in.o: dna_in.c   prim.h     dna.h      

dna_out.o: dna_out.c  prim.h     dna.h      

dna_seqs.o: dna_seqs.c prim.h     dna.h      

dna_snp.o: dna_snp.c  prim.h     dna.h      

snp_char.o: snp_char.c prim.h     dna.h      

ven_engy.o: ven_engy.c prim.h     venpipe.h  utils.h    fold_vars.h  \
          fold.h     part_func.h 

ven_io.o: ven_io.c   prim.h     venpipe.h  

ven_salt.o: ven_salt.c prim.h     venpipe.h  

ven_str.o: ven_str.c  prim.h     venpipe.h  

venpipe.o: venpipe.c  prim.h     venpipe.h  
