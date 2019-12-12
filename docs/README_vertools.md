### README
VerTools 1.5
9/18/18 RTK

Source files for vertools collection; 
This package has only source (no examples, tests, etc)



-----------------------------------------------------------------------------
### Building tools
To build programs, either call *make* for each one (using <prog>.makefile) 
or call the script that calls *make* in a loop:

    scripts/make_all_makes.sh

To copy / move resulting binaries to some directory (e.g. "~/vertools/"):

    find . -maxdepth 1 -type f -perm -100 -exec cp {} ~/vertools/ \;

    find . -maxdepth 1 -type f -perm -100 -exec mv {} ~/vertools/ \;

A script that does the same (maybe easier to call?) is this "mv_binaries.sh".
This finds binaries in the current directory and moves them to a specified
output dir. For example (from top level vertools dir after make_all_makes.sh):

    mv_binaries.sh ~/vertools/
    

**Note** that some programs will FAIL TO BUILD unless depencencies are met; These require extra libraries to link.

***venpipe*** needs the Vienna RNA library; Download this from here: https://www.tbi.univie.ac.at/RNA/

Also see [README_ViennaRNA_RTK.md](/README_ViennaRNA_RTK.md)

A script to download Vienna, build the (whole) package, then copy the library
to link with vertools "get_build_Vienna.sh". This should be run from the top
level vertools dir like this:

    scripts/get_build_Vienna.sh

***plotmat*** needs the gd graphics library; Download this from wherever...  depending on platform?


-----------------------------------------------------------------------------
### Tool descriptions
(Very) brief description of tools (*update 6/24/16 RTK*)
* ***alphcont*** Alphabet content. Used to calculate fraction of bases, base types, max runs, etc.
* ***bit_util*** Bit utility. Used to manipulate bitfield files
* ***blastout*** BLAST output parser. Used to summarize BLAST hits, get stats, etc.
* ***bpaste*** Big paste. *Legacy (may be in scripts?)* ... TODO: retire
* ***ccase*** Change case. *Legacy (may be in scripts?)* ... TODO: retire
* ***chardump*** Character dump. Reports / cleans 'weird' characters in text files (like dos end-of-line)
* ***cmake*** C makefile generator. Given list of C files, create makefile (Build tool only)
* ***comp_seq*** Compare sequences. Calculate, display various alignments; Similarity, complimentarity, 3'end, etc... TODO: Replace alignment code (e.g. Smith-waterman)
* ***dna_util*** DNA utility. Sequence file manipulations and stats. Also calculates sequence complexity.
* ***filepick*** File pick. Select subsection(s) of file based on keywords (may be in scripts)
* ***filter*** Filter file. Select lines based on numbers, random sample.
* ***fname*** File name. *Legacy (may be in scripts?)* ... TODO: retire
* ***gen_clis*** Generate C file list. This generates separate C lists suitable for cmake (Build tool only)
* ***gen_nums*** Generate numbers. Generate numbers with various patterns, including random.
* ***gen_prot*** Generate C function prototypes (Build tool only)
* ***gen_seq*** Generate sequences. Random or systematic sequence generation with various constraints.
* ***getsnps*** Get SNPs. Extract SNP sequences annotated in various ways.
* ***lassoo*** Lassoo. Diversity library tool updated since publication. Uses bitwise item representations.
* ***num*** Number tool. *Legacy (may be in scripts?)* ... TODO: retire
* ***numstat*** Number stats. Summary stats and histograms from numbers.
* ***pick_seq*** Pick sequences. Tool to pick diverse sequences, for example different barcodes.
* ***plotmat*** Plot matrix. Heatmat jpg file from number data file.
* ***scoretab*** Score table. Number table file manipulations, including applying score transforms and row product.
* ***seqtweak*** Sequence tweak. Modify sequences in various ways (e.g. 'mutate' starting sequence)
* ***shuffle*** Shuffle. Shuffle lines in file.
* ***tm_util*** Tm utility. Sequence thermodynamic calculations, including length and thermo value bounded extraction (e.g. candidate primers)
* ***venpipe*** Vienna pipeline. Wrapper for Vienna RNA folding library (can be used with DNA parameters, salt correction too)
* ***wfmerge*** Word frequence merge. Merge word frequency tables.
* ***wordfreq*** Word frequency. Calculate the frequency of N-mer words in set of sequences.

