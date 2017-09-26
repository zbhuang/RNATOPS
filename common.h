#ifndef _COMMON_H_
#define _COMMON_H_

#include "constdef.h"
#include "timecalc.h"

#include "GenomeData.h"

#include <iomanip>
//using std::setfill;
//using std::setw;

#define SEARCH_HMMFILTER			0
#define	SEARCH_SUBSTRUCTUREFILTER	1
#define	SEARCH_WHOLESTRUCTURE		2
#define	SEARCH_STREAMLINE			3

#define GENOME_TAG_HIT_EXTENSION		"Extension of the hit"
#define GENOME_TAG_HIT_EXTENSION_TAIL	"--------"
#define GENOME_TAG_TOTAL_TIME_USED		"Total time used"
#define GENOME_TAG_TOTAL_NUM_HIT		"Total no of hits"

#include <string.h>	//gcc3.4

int isRNABase(char nucleotide);
int isRNABase(int nucleotide);
int chartonum(char nucleotide);
char numtochar(int base);

bool isCanonicalPair(char basex, char basey);

int isFileName( char *str);

void printPairEmP( double prob[5][5] );
void printPairEmP( int pairprob[5][5] );
void printPairHz( int pairprob[5][5], double prob[5][5]);
void printBaseEmP( double prob[4] );
void printSingleEmProb( double baseprob[4] );
int obtainSequence(  char * filename, char* &bufseq );

//calculate normalization of two matrix
void normalizeMatrix( double orgm[5][5], double priorm[5][5], double h);

char * strrpl(char *s, const char s1, const char s2);
//char * strdel(char *s, const char s1);
char * strappend(char *s, const char s1);
char * strextract(char *s, int start_pos, int end_pos);
string extractfilename(string s);
string printYN(int iFlag);
string printDirection(int iFlag);
void printGenomeSegments(char * c_seq, int num);

//
int isHmmFilterFile(char * filename);
/*
void printValuesForParameters(double	pseudocount,
							  int		top_k,
							  double	threshold,
							  int		num_nts_overlap_between_stem,
							  int		iFlagMergeCandInPreprocess,
							  int		iFlagCandwithShortestLength,
							  int		iShiftNumMergeCand,
							  int		iAllowedNullLoopInsNum,
							  double	pcoeff,
							  int		iJumpStrategy,
							  int		iStepSize,
							  double	dScoreThresholdInJump,
							  int		iFlagPrintStrAlignInfo,
							  int		iFlagSearchReverse,
							  int		iFlagPrintDebugInfo,
							  int		iFlagPrintDebugInputInfo,
							  int		iFlagPrintDebugTreeInfo,
							  int		iFlagPrintDebugDPWinInfo,
							  int		iFlagPrintDebugTDInfo,
							  int		iFlagPrintDebugDPSearchInfo);
*/
void printValuesForParameters(double	pseudocount,
							  int		top_k,
							  double	threshold,
							  int		num_nts_overlap_between_stem,
							  int		iFlagMergeCandInPreprocess,
							  int		iFlagCandwithShortestLength,
							  int		iShiftNumMergeCand,
							  int		iAllowedNullLoopInsNum,
							  double	pcoeff,
							  int		iJumpStrategy,
							  int		iStepSize,
							  int		iFlagSearchReverse);

void printCurrentTime();
void printOneStarLine();
void printAdditionalInfo();
void printFilterSearchResultHeader(string	trainfile,
								   int		profilelength,
								   int		filterauto, 
								   int		filtertype, 
								   int		filterbegin, int filterend, 
								   string	genomefile,	//char * genomefile,
								   int		genomenum,
								   long		total_genome_length);
void printFilterHitIndex(int index);
void printHitPos(int start, int end);
void printHitScore(double score);
void printFilterHitAlignment(char *seq, char *align);
void printFilterHitExtenPos(int start, int end);
void printFilterHitExtenNtsHeader();
//void printFilterHitExtenNtsTail();
void printTotalHitNum(int hit_num);
void printTotalTimeInfo(double dTime);
void printSearchEndingTimeInfo();
void printWholeSearchResultHeader(string	trainfile,
								  int		profilelength,
								  string	genomefile,
								  int		genomenum,
								  long		total_genome_length,
								  double	pseudocount,
								  int		top_k,
								  double	threshold,
								  int		num_nts_overlap_between_stem,
								  int		iFlagMergeCandInPreprocess,
								  int		iFlagCandwithShortestLength,
								  int		iShiftNumMergeCand,
								  int		iAllowedNullLoopInsNum,
								  double	pcoeff,
								  int		iJumpStrategy,
								  int		iStepSize,
								  int		iFlagSearchReverse,
								  int		streamline,
								  int		searchtype,
								  int		filterbegin,
								  int		filterend);
int	getProfileLength(char * trainfile);
void printWholeSearchHitIndex(int index);
void printEValue(double evalue);
void formatStringOutput(string illustration, string line_one, string line_two, string line_three, int start, int end, int flag);
void formatOneOutput(string line_one_unit, string line_two_unit, string line_thr_unit, int index, int start, int flag);
void formatStringOutput2(char * illustration, char * line_one, char * line_two, char * line_three, int start, int end, int flag);
void formatOneOutput2(char * line_one_unit, char * line_two_unit, char * line_thr_unit, int index, int start, int flag);

#endif
