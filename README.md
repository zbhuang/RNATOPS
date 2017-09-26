RNATOPS version 1.2

GENERAL INFORMATION
-------------------
RNATOPS is a program that searches for RNA secondary structures based on the notion of a structure graph to specify the consensus structure of an RNA family.

This program is implemented in C++. It has been compiled and tested on several systems, including Desktop Linux computers, a Linux cluster, and a SUN workstation running SunOS 5.1.

Authors:
  Zhibin Huang, Yong Wu, Joseph Robertson, Yingfeng Wang
  Department of Computer Science
  University of Georgia
  Athens, GA 30602

Copyright (C) 2008 University of Georgia. All Rights Reserved.

IMPROVEMENTS OVER RNATOPS version V1.1
--------------------------------------
Support hmm-filter-based searches to improve efficiency. 

[-hmmfilter] with this option, search produces hits based on an automatically selected profile hmm (as filter) and 
using the whole structure profile on these hit sequences, the program produces final search results directly.

WHAT IS CONTAINED IN THIS PACKAGE
---------------------------------
The software package of RNATOPS and one example are contained in this package.

INSTALLATION
------------
Compile RNATOPS.

Type "make" to build an executable program called "rnatops".

USAGE
-----
1. How to run this program
Usage: rnatops <-tf training_file> <-gf genome_file> [other options]

2. Command line parameter list:
[-tf]			training_file
[-gf]			genome_file
[-hmmfilter]	search with automatic-hmmfilter
[-pcnt]			pseudocount_value [DEFAULT=0.001]
[-k]			k_value [DEFAULT=10]
[-th]			threshold_value_for_candidate_hit [DEFAULT=0.0]
[-no]			number of overlapping nts between two stems when overlap is allowed [DEFAULT=2]
[-mc]			to merge candidates in preprocessing
[-ms]			when taking the merge strategy, take the candidate with the (m)ax score 
                or the one with the (s)hortest length (m|s) [DEFAULT=s]
[-ns]			number of shift positions allowed in the merge strategy [DEFAULT=0]
[-ni]			number of insertions allowed in a null loop [DEFAULT=3]
[-pcv]			weight for prior base pairing matrix [DEFAULT=0.4]
[-pf]			prior_file_name [DEFAULT='./base_pair_prior.txt']
[-pv]			pcoeff value [DEFAULT=2.0]
[-st]			split threshold [DEFAULT=6]
[-js]			stepsize in jump strategy [DEFAULT=1]
[-jt]			score_filtering threshold in jump strategy [DEFAULT=0.0]
[-r]			to also search reverse complement strand
[-ps]			to print structure alignment info
[-d]			to print out the debug info
[-time]			to print out the time info
[-pscore]		to print score info for stem and loop
[-dfilter]		to print out the debug info for the filter based search
[-ddpinput]		to print out the debug info of the input in dp search
[-ddptree]		to print out the debug info of the tree in dp search
[-ddpwin]		to print out the debug info of the window in dp search
[-ddptd]		to print out the debug info of the td in dp search
[-ddpsearch]	to print out the debug info of the dp search

3. Command line parameter explanation:
[-tf]	training_file
File of RNA training data in the pasta format to produce the RNA structure model.
Two examples (one pk-free and one pk) are included in the directory of "example".

[-gf]	genome_file
File to be searched for the structure modeled with the training data.
Current version requires that genome_file only contain ACGT nucleotide characters.
See Section of "PREPARATION OF THE GENOME FILE".

[-hmmfilter]	search with automatic-hmmfilter
This option is to do the automatically selected hmm filter-based fast search,
If not specified this option, then the program will do the whole-structure search.

[-pcnt] pseudocount_value [DEFAULT=0.001]
The value of pseudocount.
Default value is 0.001.

[-k]	k_value [DEFAULT=10]
The number of candidates for each stem in the structure to be searched in genome_file.
Default value is 10.

[-th]	threshold_value_for_candidate_hit [DEFAULT=0.0]
the threshold value which is used to determine whether the current structure is hit or not.
Default value is 0.0.

[-no]	number of overlapping nts between two stems when overlap is allowed [DEFAULT=2]
The number of overlapping nucleotides allowed between adjacent candidates.
Default value is 2.

[-mc]	to merge candidates in preprocessing
If you use this option, the strategy of merging similar candidates will be taken.
By default, it is used when doing the structure search.

[-ms]	when taking merge strategy, take the candidate with (m)ax score or one with the (s)hortest length (m|s) [DEFAULT=s]
When using the -mc option, merging candidates, then you need to choose which merging strategy you will choose.
Choosing the representive candidate with max score or shortest length.
m: max score;  s: shortest length.
Default value is s.

[-ns]	number of shift positions allowed in the merge strategy [DEFAULT=0]
When merging candidates, the candidates within the number of shift positions will be merged together.
Default value is 0.

[-ni]	number of insertions allowed in null loop [DEFAULT=3]
For the null loop model, the number of inserted nucleotides will be allowed. 
Default value is 3.

[-pcv]	prior weight value [DEFAULT=0.4]
The weight for prior base pairing frequency matrix. 
Default value is 0.4.

[-pf]	prior_file [DEFAULT='./base_pair_prior.txt']
The file containing prior pair frequency matrix (5x5). 
base_pair_prior.txt.
(Note: User can put this file any places as long as he needs to specify the absolute path of this file or he can use his own matrix file)

[-pv]	pcoeff value [DEFAULT=2.0]
The weight of distance penalty. 
Default value is 2.0.

[-st]	split threshold [DEFAULT=6]
The standard deviation value to group the training sequences.  
Default value is 6.

[-js]	stepsize in jump strategy [DEFAULT=1]
This program supports the skip-and-jump strategy.
The scanning window is shifted by the 'stepsize' nucleotides to speed up search.
2 or 3 for the stepsize is suggested. 

[-jt]	score_filtering threshold in jump strategy [DEFAULT=0.0]
Default value is 0.0.

[-r]	to search the complementary strand 
This option can be used to search the complementary genome strand.

[-ps]	to print out info of the folded structure and structure alignment.
This option can give you more info of the candidate hit the program finds.
Definitions of folded structure and structure alignment are given in "EXPLANATION OF THE RESULTS AND OUTPUT" in this file.

[-time]		to print out the time info
This option can give you the time info about your running RNATOPS.

[-pscore]	to print score info for stem and loop
This option can give you the score info for every stem and loop.

[-dfilter]	to print out the debug info for the filter based search
[-ddpinput]	to print out the debug info of the input in dp search
[-ddptree]	to print out the debug info of the tree in dp search
[-ddpwin]	to print out the debug info of the window in dp search
[-ddptd]	to print out the debug info of the td in dp search
[-ddpsearch]	to print out the debug info of the dp search
[-d]		to print out the debug info
If you want to debug the program, you can use -d option.
We suggest you do not debug this program until you have a good understanding of the code 
because it will produce a large volume of data especially for the the tree decomposition based dynamic programming.

4.  How to use RNATOPS to do the search.
Generally user can do two kinds of search: automatic hmm-filter-based search and whole-structure search in a stream line or direct whole-structure search. 

The data directory 'example' contains one sample using the option of "-hmmfilter" and no "-hmmfilter" option:

the script for running the program.

./rnatops -tf ./example/t_02.txt -gf ./example/revised-g.txt -hmmfilter -r >./example/wholesearchresult_t.txt
./rnatops -tf ./example/t_02_twopastaline.txt -gf ./example/revised-g.txt -hmmfilter -r >./example/wholesearchresult_t_twopastaline.txt
./rnatops -tf ./example/tmRNA_43_training_24_02.txt -gf ./example/RC_NC_008787.fasta -hmmfilter -r >./example/wholesearchresult_RC_NC_008787.txt
./rnatops -tf ./example/tmRNA_43_training_24_02_twopastaline.txt -gf ./example/RC_NC_008787.fasta -hmmfilter -r >./example/wholesearchresult_RC_NC_008787_twopastaline.txt

EXPLANATION OF THE RESULTS AND OUTPUT
-------------------------------------
The following information of every found structure candidate (of a score above the threshold) is generated as a result: 
	search direction;
	total alignment score, begin and end positions in the genome; 
	folded structure of the candidate; structural alignment with the model;

The following is an example of running result from RC_NC_008787.fasta.

------------------------------------------------------------
- Whole Structure Result                                   -
-                                                          -
- Profile file : tmRNA_43_training_24_02.txt               -
- Profile length : 582                                     -
-                                                          -
- Filter type : HMM                                        -
- Filter info : positions from 535    to 581               -
-                                                          -
- Genome file : RC_NC_008787.fasta                         -
- Number of sequences : 1                                  -
- Total length of sequences : 1616554                      -
-                                                          -
- Search parameters setting:                               -
- Pseudocount = 0.001                                      -
- Num of stem candidates = 10                              -
- Score threshold for hits = 0                             -
- Num of nt overlap between stems = 2                      -
- Candidate representatives only = Yes                     -
- Shortest candidate representatives = Yes                 -
- IShiftNumMergeCand = No                                  -
- Nts allowed in null loops = 3                            -
- Pcoeff = 2                                               -
- Search with jump = Yes                                   -
- Search step size = 1                                     -
- Search reversed complement sequence = Yes                -
------------------------------------------------------------

Whole structure search hit 1
----------------------------
RC_NC_008787
Plus search result
Hit Positions: 337938-338291
Alignment score = 17.2457

Folded structure

         1 YYYYYY..........WWWWWWWWWW.....VVVVVUUUUTTTTTTAAA**GGGGG.AAA 61
    337939 GGAGCGACTTGGCTTCGACAGGAGTAAGTCTGCTTAGATGGCATGTCGCTTTGGGCAAAG 337999

        61 AA......GGGGGGG.........................................IIII 121
    337999 CGTAAAAAGCCCAAATAAAATTAAACGCAAACAACGTTAAATTCGCTCCTGCTTACGCTA 338059

       121 IIHHHH........HHH*IIII*MMMMMM.LLLL*JJJ.BBBBBBJJJ*LLLLMMMMMMM 181
    338059 AAGCTGCGTAAGTTCAGTTGAGCCTGAAATTTAAGTCATACTATCTAGCTTAATTTTCGG 338119

       181 ........BBBBBB..OOOOOOOO**CCCC...OOOO.OOOOOO...CCCCCC....QQQ 241
    338119 TCATTTTTGATAGTGTAGCCTTGCGTTTGACAAGCGTTGAGGTGAAATAAAGTCTTAGCC 338179

       241 QQQQQQQ.......P**DDDDD...P**QQQQQQQQ.........DDDiDDDD...TTTT 301
    338179 TTGCTTTTGAGTTTTGGAAGATGAGCGAAGTAGGGTGAAGTAGTCATCTTTGCTAAGCAT 338239

       301 TT...UUUU....VVVVV.WWWWWWWWWW.XXXXX.......XXXXXYYYY.YY 355
    338239 GTAGAGGTCTTTGTGGGATTATTTTTGGACAGGGGTTCGATTCCCCTCGCTTCC 338293

Structure alignment

         1 YYYYYYY..........WWWWWWWWWW.....VVVVVUUUUTTTTTTAAA**GGGGG.AA 61
    337939 GG-AGCGACUUGGCUUCGACAGGAGUAAGUCUGCUUAGAUGGCAUGUCGCUUUGGGCAAA 337999

        61 AAA.....-GGGGGGG..............................------------.. 121
    337999 GCGUAAAAAGCCCAAAUAA-AAUUAAACGCAAACAACGUUAAAU--UCGCUCCUGCUUAC 338059

       121 IIIIIIHHHH.........HHH*IIII*MMMMMM-LLLL*JJJ.BBBBBBJJJ*LLLLMM 181
    338059 GCUAAAGCUGCGUAAGU-UCAGUUGAGCCUGAAAUUUAAGUCAUACUAUCUAGCUUAAUU 338119

       181 MMMMM........BBBBBB..--OOOOOOOO**CCCC...OOOOrOOOOOO.......rC 241
    338119 UUCGGUCAUUUUUGAUAGU--GUAGCCUUGCGUUUGACAAGCGUUGAGGUG---AA--AU 338179

       241 CCCCC..--QQQQQQQQQQ------.P**DDDDD...P**QQQQQQQQ....-......D 301
    338179 AAAGUCUUAGCCUUGCUUUUGAGUUUUGGAAGAUGAGCGAAGUAGGGUGAAGU--AGUCA 338239

       301 DD-DDDD...TTTTTT.-.UUUU----VVVVV.....WWWWWWWWWW.XXXXX....... 361
    338239 UCUUUGCUAAGCAUGUAGAGGUCUUUGUGGGA----UUAUUUUUGGACAGGGGUUCGAUU 338299

       361 XXXXXYYYYYYY 373
    338299 CCCCUCGCUUCC 338311

------------------------------------------------------------
- Searched done : with RNATOPS V1.2                        -
- By RNA-Informatics @ UGA                                 -
- Total no of hits: 1                                      -
- Total time used 0.273942        hours                    -
- Time: 23:09:40 EST 2008-11-21                            -
------------------------------------------------------------

Note:
The folded structure is the predicted secondary structure of the found RNA candidate 
in the genome, where a base pair is annotated by a pair of letter. (See referenced pasta paper). 

The structure alignment shows the optimal alignment between the model and the 
RNA sequence found on the genome. It contains additional information of matches,
deletions, and insertions computed by the RNATOPS program.

CONTACT INFORMATION
-------------------
For suggestions, questions, requests and bug report:
please contact Liming Cai at cai@cs.uga.edu.
