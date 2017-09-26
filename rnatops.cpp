#include "DataLoader.h"
#include "ScfgBuilder.h"
#include "CYK.h"
#include "HMMBuilder.h"
#include "Viterbi.h"
#include "FViterbi.h"
#include "CandScanner.h"
#include "DynamicP.h"

#include "GenomeData.h"

#include "StructureSearch.h"
#include "common.h"

void Usage();

int main(int argc, char *argv[])
{
	char *	c_run_mode;
	int		run_mode	= RUNNINGMODE_STANDALONE;
	int		searchtype	= SEARCH_WHOLESTRUCTURE;
	char *	filter_file;
	int		filterbegin	= -1;
	int		filterend	= -1;
	char *	realtrainfilename;
	char *	realgenomefilename;

	char	cOffsetMethod='L';  //'L': compute as a whole from left-end; 'S': compute as a summation from left-end
	char *	train_file;
	char *	genome_file;
	double	pseudocount = 0.001;
	int		top_k = 10;
	double	threshold = 0.0;
    int		num_nts_overlap_between_stem = 2;
	int		iFlagMergeCandInPreprocess = 1;	//0;		//
	int		iFlagCandwithShortestLength = 1;
	int		iShiftNumMergeCand = 0;
	int		iAllowedNullLoopInsNum = 3;

	int		iFlagPrintDebugInfo = 0;	//print out the debug info
	int		iFlagPrintDebugFilterInfo	= 0;	//
	int		iFlagPrintDebugInputInfo	= 0;	//the debug info of the flag in dp search
	int		iFlagPrintDebugTreeInfo		= 0;	//the debug info of the tree in dp search
	int		iFlagPrintDebugDPWinInfo	= 0;	//the debug info of the window in dp search
	int		iFlagPrintDebugTDInfo		= 0;	//the debug info of the td in dp search
	int		iFlagPrintDebugDPSearchInfo	= 0;	//the debug info of the dp search
	double	priorh = 0.4;
	char *	prior_file;
	double	pcoeff = 2.0;
	int		iSplitThreshold = 6;
	int		iJumpStrategy = 1;	//jump strategy
	int		iStepSize = 1;
	double	dScoreThresholdInJump = 0.0;
	int		iFlagPrintStrAlignInfo = 1;	//0;
	int		iFlagSearchReverse = GENOME_SEARCH_DIRECTION_ONE;
	int		iFlagShowTimeInfo = 1;	//0;
	int		iFlagPrintScoreInfo = 0;

	int		iFlagAutoFilter = 0;		//by default in standalone mode, filter function will be called.
	int		iFlagUserDefineFilter = 0;

	int		iFlagStreamLine = 0;	//1;//

	train_file = '\0';
	genome_file = '\0';
	prior_file = "./base_pair_prior.txt";
	c_run_mode = "standalone";//"webserver";
	realtrainfilename	= NULL;
	realgenomefilename	= NULL;
	filter_file = '\0';

	if(argc == 1)
		Usage();

	int i;
    // process arguments
    for (i=1; i<argc; i++) {
        if (argv[i][0] == '-') {
            if (!strcmp(argv[i], "-rmode")) {
                if ((i+1) >= argc) {	Usage();	}
                c_run_mode = argv[++i];
				if(strcmp(c_run_mode, "webserver") && strcmp(c_run_mode, "standalone")) {	
					Usage();	
				}
            } else if (!strcmp(argv[i], "-tf")) {
                if ((i+1) >= argc) {	Usage();	}
                train_file = argv[++i];
            } else if (!strcmp(argv[i], "-gf")) { 
                if ((i+1) >= argc) {	Usage();	}
                genome_file = argv[++i];
//			} else if (!strcmp(argv[i], "-rtfn")) {
//				if ((i+1) >= argc) {	Usage();	}
//				realtrainfilename = argv[++i];
//			} else if (!strcmp(argv[i], "-rgfn")) { 
//				if ((i+1) >= argc) {	Usage();	}
//				realgenomefilename = argv[++i];
			} else if (!strcmp(argv[i], "-pcnt")) { 
                if ((i+1) >= argc) {	Usage();	}
                pseudocount = atof(argv[++i]);
            } else if (!strcmp(argv[i], "-k")) { 
                if ((i+1) >= argc) {	Usage();	}
                top_k = atoi(argv[++i]);
            } else if (!strcmp(argv[i], "-th")) { 
                if ((i+1) >= argc) {	Usage();	}
                threshold = atof(argv[++i]);
            } else if (!strcmp(argv[i], "-no")) { 
                if ((i+1) >= argc) {	Usage();	}
                num_nts_overlap_between_stem = atoi(argv[++i]);
            } else if (!strcmp(argv[i], "-mc")) {	//merge candidate in preprocessing
                if (i >= argc) {	Usage();	}
				iFlagMergeCandInPreprocess = 1;
            } else if (!strcmp(argv[i], "-ms")) {	//merge strategy, the candidate with (m)ax score or one with the (s)hortest length
                if ((i+1) >= argc) {	Usage();	}
				if(argv[++i][0] == 'm')
	                iFlagCandwithShortestLength = 0;
            } else if (!strcmp(argv[i], "-ns")) { 
                if ((i+1) >= argc) {	Usage();	}
                iShiftNumMergeCand = atoi(argv[++i]);
            } else if (!strcmp(argv[i], "-ni")) { 
                if ((i+1) >= argc) {	Usage();	}
                iAllowedNullLoopInsNum = atoi(argv[++i]);
            } else if (!strcmp(argv[i], "-pcv")) { 
                if ((i+1) >= argc) {	Usage();	}
                priorh = atof(argv[++i]);
            } else if (!strcmp(argv[i], "-pf")) { 
                if ((i+1) >= argc) {	Usage();	}
                prior_file = argv[++i];
            } else if (!strcmp(argv[i], "-pv")) { 
                if ((i+1) >= argc) {	Usage();	}
                pcoeff = atof(argv[++i]);
            } else if (!strcmp(argv[i], "-st")) { 
                if ((i+1) >= argc) {	Usage();	}
                iSplitThreshold = atoi(argv[++i]);
			} else if (!strcmp(argv[i], "-js")) { 
                if ((i+1) >= argc) {	Usage();	}
                iStepSize = atoi(argv[++i]);
            } else if (!strcmp(argv[i], "-jt")) { 
                if ((i+1) >= argc) {	Usage();	}
                dScoreThresholdInJump = atof(argv[++i]);
//			} else if (!strcmp(argv[i], "-ps")) {	//print out info of folded structure and structure alignment
//				if (i >= argc) {	Usage();	}
//				iFlagPrintStrAlignInfo = 1;
			} else if (!strcmp(argv[i], "-pscore")) {	//print out info of folded structure and structure alignment
				if (i >= argc) {	Usage();	}
				iFlagPrintScoreInfo = 1;
			} else if (!strcmp(argv[i], "-hmmfilter")) {
				if (i >= argc) {	Usage();	}
				searchtype = SEARCH_HMMFILTER;
				iFlagAutoFilter = 1;
				iFlagStreamLine	= 1;
//			} else if (!strcmp(argv[i], "-streamline")) {
//				if (i >= argc) {	Usage();	}
//				iFlagStreamLine = 1;
//			} else if (!strcmp(argv[i], "-substructurefilter")) {	//print out info of folded structure and structure alignment
//				if (i >= argc) {	Usage();	}
//				searchtype = SEARCH_SUBSTRUCTUREFILTER;
//			} else if (!strcmp(argv[i], "-autofilter")) {
//				if (i >= argc) {	Usage();	}
//				iFlagAutoFilter = 1;
//			} else if (!strcmp(argv[i], "-filterfile")) {
//				if ((i+1) >= argc) {	Usage();	}
//				filter_file = argv[++i];
//			} else if (!strcmp(argv[i], "-userdefinedfilterfile")) {
//				if ((i+1) >= argc) {	Usage();	}
//				iFlagUserDefineFilter = 1;
//				filter_file = argv[++i];
//			} else if (!strcmp(argv[i], "-filterbegin")) { 
//				if ((i+1) >= argc) {	Usage();	}
//				filterbegin = atoi(argv[++i]);
//			} else if (!strcmp(argv[i], "-filterend")) { 
//				if ((i+1) >= argc) {	Usage();	}
//				filterend = atoi(argv[++i]);
//			} else if (!strcmp(argv[i], "-time")) {	//print out the time info
//				if (i >= argc) {	Usage();	}
//				iFlagShowTimeInfo = 1;
            } else if (!strcmp(argv[i], "-d")) {	//print out the debug info
				if (i >= argc) {	Usage();	}
				iFlagPrintDebugInfo = 1;
            } else if (!strcmp(argv[i], "-dfilter")) {	//print out the debug info
				if (i >= argc) {	Usage();	}
				iFlagPrintDebugFilterInfo = 1;
            } else if (!strcmp(argv[i], "-ddpinput")) {	//print out the debug info
				if (i >= argc) {	Usage();	}
				iFlagPrintDebugInputInfo = 1;
            } else if (!strcmp(argv[i], "-ddptree")) {	//print out the debug info
				if (i >= argc) {	Usage();	}
				iFlagPrintDebugTreeInfo = 1;
            } else if (!strcmp(argv[i], "-ddpwin")) {	//print out the debug info
				if (i >= argc) {	Usage();	}
				iFlagPrintDebugDPWinInfo = 1;
            } else if (!strcmp(argv[i], "-ddptd")) {	//print out the debug info
				if (i >= argc) {	Usage();	}
				iFlagPrintDebugTDInfo = 1;
            } else if (!strcmp(argv[i], "-ddpsearch")) {	//print out the debug info
				if (i >= argc) {	Usage();	}
				iFlagPrintDebugDPSearchInfo = 1;
            } else if (!strcmp(argv[i], "-r")) {	//reverse complement strand search
				if (i >= argc) {	Usage();	}
				iFlagSearchReverse = GENOME_SEARCH_DIRECTION_BOTH;
            } else { 
                Usage();
            }
        } else {
            Usage();
        } 
    }

	if(!strcmp(c_run_mode, "webserver"))
		run_mode = RUNNINGMODE_WEBSERVER;
	else
		run_mode = RUNNINGMODE_STANDALONE;

	if(iFlagUserDefineFilter == 1) {
		int iSearchType = isHmmFilterFile(filter_file);
		if(iSearchType == 0)
			searchtype = SEARCH_SUBSTRUCTUREFILTER;
		else if(iSearchType == 1)
			searchtype = SEARCH_HMMFILTER;
	}

	int		genome_num		= 0;
	int *	genome_length	= NULL;
	int *	genome_type		= NULL;
	int *	genome_direction= NULL;
	int *	genome_extleft	= NULL;
	int *	genome_extright	= NULL;
	char **	genome_name		= NULL;
	char **	genome_sequence		= NULL;
	char **	genome_reversebuf	= NULL;
	double ** genome_baseFreq	= NULL;

	float fFilterTime = 0.0;

	GenomeData genomeObj;
	StructureSearch search;

	string str_train_file;
	string str_genome_file;
	if(realtrainfilename == NULL)
	{
		str_train_file = extractfilename(train_file);
	} else {
		str_train_file = extractfilename(realtrainfilename);
	}
				
	if(realgenomefilename == NULL) {
		str_genome_file = extractfilename(genome_file);
	} else {
		str_genome_file = extractfilename(realgenomefilename);
	}

	int iFlagFilterAuto = 1;
	if(searchtype == SEARCH_HMMFILTER)	//hmmfilter search
	{
		//load the original genome files
//		genomeObj.loadMultipleGenomes(genome_file);
		genomeObj.loadMultipleGenomes2(genome_file);
		genome_num		= genomeObj.getMultipleGenomeNum();
		genome_name		= genomeObj.getMultipleGenomeName();
		genome_sequence	= genomeObj.getMultipleGenomeDataC();
		genome_length	= genomeObj.getMultipleGenomeLength();
		genome_type		= genomeObj.getMultipleGenomeType();
		genome_direction= genomeObj.getMultipleSearchDirect();	//20080905
		genome_baseFreq = genomeObj.getMultipleGenomeBaseFreq();

		//support rc hmm filter search 20080905
		if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_BOTH) {
			genomeObj.generateMultipleGenomeReverseData();
			genome_reversebuf = genomeObj.getMultipleGenomeReverseDataC();
		}

		//just calculate filterbegin and filterend	20080927
		int j;
		if(iFlagAutoFilter == 1)
		{
			//calling hmm filter selection
			vector <string> cc;
			string thestandard="TACGU";//"_ACGU";

			FilterSelection filterselection(train_file);
			cc=filterselection.getFilter();

			if (filterselection.getfilterflag()== 0) {
				cout<<"No filter is selected, pls checking the filter selection. Thanks."<<endl;
				exit(1);
			}
			else
			{
				filterbegin = filterselection.getbegposi();
				filterend	= filterselection.getendposi();

				if(iFlagPrintDebugFilterInfo)
				{
//					for (i=0; i<(int)cc.size(); i++)
//					{
//						for (j=0; j<(int)cc[0].length(); j++)
//							cout<<cc[i][j];
//						cout<<endl;
//					}
					cout<<"filterbegin="<<filterbegin<<endl;
					cout<<"filterend  ="<<filterend<<endl;
				}

			}
		}
		else
		{
			//get one of the training sequences
			DataLoader trainLoader;
			trainLoader.inputData(train_file);
			int ** trainingSeqs = trainLoader.getPointerOfSeqs();
			int length = trainLoader.getSeqLength();
			string one_training_seq;
			for(i=0; i<length; i++)
				one_training_seq += numtochar(trainingSeqs[0][i]);

			//
			DataLoader filterLoader;
			filterLoader.inputData(filter_file);
			int ** filteringSeqs = filterLoader.getPointerOfSeqs();
			length = filterLoader.getSeqLength();
			string one_filtering_seq;
			for(i=0; i<length; i++)
				one_filtering_seq += numtochar(filteringSeqs[0][i]);

			filterbegin = one_training_seq.find(one_filtering_seq);
			filterend	= filterbegin + one_filtering_seq.length();

			if(iFlagPrintDebugFilterInfo)
			{
				cout<<"filterbegin="<<filterbegin<<endl;
				cout<<"filterend  ="<<filterend<<endl;
			}
		}

		if(iFlagStreamLine)
		{
			//print whole structure search result header. 20081005
			//whole structure search
			printWholeSearchResultHeader(	str_train_file,
											getProfileLength(train_file),
											str_genome_file,	//genomefile,
											genome_num,
											genomeObj.getTotalGenomeLength(),
											pseudocount,
											top_k,
											threshold,
											num_nts_overlap_between_stem,
											iFlagMergeCandInPreprocess,
											iFlagCandwithShortestLength,
											iShiftNumMergeCand,
											iAllowedNullLoopInsNum,
											pcoeff,
											iJumpStrategy,
											iStepSize,
											iFlagSearchReverse,
											iFlagStreamLine,
											SEARCH_HMMFILTER,
											filterbegin,
											filterend);
		} else {
			//print filter search result header. 20080909
			printFilterSearchResultHeader(	str_train_file,
											getProfileLength(train_file),
											iFlagFilterAuto, 
											SEARCH_HMMFILTER, 
											filterbegin,
											filterend, 
											str_genome_file,	//genome_file,
											genome_num,
											genomeObj.getTotalGenomeLength());
		}

		//hmm filter search
		search.search(	top_k,
						threshold,
						num_nts_overlap_between_stem,
						iFlagMergeCandInPreprocess,
						iFlagCandwithShortestLength,
						iShiftNumMergeCand,
						iAllowedNullLoopInsNum,
						pcoeff,
						iJumpStrategy,
						iStepSize,
						dScoreThresholdInJump,
						iFlagSearchReverse,
						iFlagShowTimeInfo,
						iFlagPrintStrAlignInfo,
						iFlagPrintScoreInfo,
						iFlagPrintDebugInfo,
						iFlagPrintDebugFilterInfo,
						iFlagPrintDebugInputInfo,
						iFlagPrintDebugTreeInfo,
						iFlagPrintDebugDPWinInfo,
						iFlagPrintDebugTDInfo,
						iFlagPrintDebugDPSearchInfo,
						cOffsetMethod,
						iSplitThreshold,
						train_file,
						iFlagAutoFilter,
						filter_file,
						filterbegin,
						filterend,
						pseudocount,
						priorh,
						prior_file,
						genome_num,
						genome_type,
						genome_direction,
						genome_length,
						genome_name,
						genome_sequence,
						genome_reversebuf,
						genome_baseFreq,
						genome_extleft,
						genome_extright,
						fFilterTime,		//
						run_mode,
						SEARCH_HMMFILTER,
						iFlagStreamLine);	//
	} 
	else 
	{
		//common part of sub-structurefilter and whole-structure search
//		genomeObj.loadMultipleGenomes(genome_file);
		genomeObj.loadMultipleGenomes2(genome_file);

		genome_num		= genomeObj.getMultipleGenomeNum();
		genome_name		= genomeObj.getMultipleGenomeName();
		genome_sequence	= genomeObj.getMultipleGenomeDataC();
		genome_length	= genomeObj.getMultipleGenomeLength();
		genome_type		= genomeObj.getMultipleGenomeType();
		genome_direction= genomeObj.getMultipleSearchDirect();	//20080905
		genome_extleft	= genomeObj.getMultipleGenomeExtLeft();
		genome_extright	= genomeObj.getMultipleGenomeExtRight();
		genome_baseFreq	= genomeObj.getMultipleGenomeBaseFreq();

		fFilterTime = genomeObj.getTotalFilterSearchTime();	//

		if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_BOTH) {
			genomeObj.generateMultipleGenomeReverseData();
			genome_reversebuf = genomeObj.getMultipleGenomeReverseDataC();
		}

		if(searchtype == SEARCH_WHOLESTRUCTURE) 
		{
			if(genome_type[0] == GENOME_SEGMENT) {
				iFlagSearchReverse = GENOME_SEARCH_DIRECTION_ONE;
			}
			//whole structure search
			printWholeSearchResultHeader(	str_train_file,
											getProfileLength(train_file),
											str_genome_file,	//genomefile,
											genome_num,
											genomeObj.getTotalGenomeLength(),
											pseudocount,
											top_k,
											threshold,
											num_nts_overlap_between_stem,
											iFlagMergeCandInPreprocess,
											iFlagCandwithShortestLength,
											iShiftNumMergeCand,
											iAllowedNullLoopInsNum,
											pcoeff,
											iJumpStrategy,
											iStepSize,
											iFlagSearchReverse,
											0,
											0,
											0,
											0);
		} else {
			//first calculate the filterbegin and filterend	20080927
			//
			if(iFlagAutoFilter == 1) {
				//to be developped	20080927
			} else {
				//get one of the training sequences
				DataLoader trainLoader;
				trainLoader.inputData(train_file);
				int ** trainingSeqs = trainLoader.getPointerOfSeqs();
				int length = trainLoader.getSeqLength();
				string one_training_seq;
				for(i=0; i<length; i++)
					one_training_seq += numtochar(trainingSeqs[0][i]);

				//
				DataLoader filterLoader;
				filterLoader.inputData(filter_file);
				int ** filteringSeqs = filterLoader.getPointerOfSeqs();
				length = filterLoader.getSeqLength();
				string one_filtering_seq;
				for(i=0; i<length; i++)
					one_filtering_seq += numtochar(filteringSeqs[0][i]);

				filterbegin = one_training_seq.find(one_filtering_seq);
				filterend	= filterbegin + one_filtering_seq.length();

				if(iFlagPrintDebugFilterInfo)
				{
					cout<<"filterbegin="<<filterbegin<<endl;
					cout<<"filterend  ="<<filterend<<endl;
				}
			}
			//print filter search result header. 20080909
			printFilterSearchResultHeader(	str_train_file,
											getProfileLength(train_file),
											iFlagFilterAuto, 
											SEARCH_SUBSTRUCTUREFILTER, 
											filterbegin,
											filterend, 
											str_genome_file,	//genome_file,
											genome_num,
											genomeObj.getTotalGenomeLength());
		}
		
		//structure search
		search.search(top_k,
					threshold,
					num_nts_overlap_between_stem,
					iFlagMergeCandInPreprocess,
					iFlagCandwithShortestLength,
					iShiftNumMergeCand,
					iAllowedNullLoopInsNum,
					pcoeff,
					iJumpStrategy,
					iStepSize,
					dScoreThresholdInJump,
					iFlagSearchReverse,
					iFlagShowTimeInfo,
					iFlagPrintStrAlignInfo,
					iFlagPrintScoreInfo,
					iFlagPrintDebugInfo,
					iFlagPrintDebugFilterInfo,
					iFlagPrintDebugInputInfo,
					iFlagPrintDebugTreeInfo,
					iFlagPrintDebugDPWinInfo,
					iFlagPrintDebugTDInfo,
					iFlagPrintDebugDPSearchInfo,
					cOffsetMethod,
					iSplitThreshold,
					train_file,
					iFlagAutoFilter,
					filter_file,	//
					filterbegin,	//
					filterend,		//
					pseudocount,
					priorh,
					prior_file,
					genome_num,
					genome_type,
					genome_direction,
					genome_length,
					genome_name,
					genome_sequence,
					genome_reversebuf,
					genome_baseFreq,
					genome_extleft,
					genome_extright,
					fFilterTime,	//
					run_mode,
					searchtype,		
					iFlagStreamLine);	//
	}

    return 0;
}
void Usage() {
	cout<<endl<<"Usage: rnatops	[-tf] training_file [-gf] genome_file [other options]"<<endl
		<<"[-tf]	training_file"<<endl
		<<"[-gf]	genome_file"<<endl
//		<<"[-rmode]	running mode (webserver|standalone)[DEFAULT=standalone]"<<endl
		<<"[-hmmfilter]				to do automatic-hmmfilter search"<<endl
//		<<"[-substructurefilter]	to do substructurefilter search"<<endl
//		<<"[-autofilter]	to do filter search, automatically construct the filter-file"<<endl
//		<<"[-filterfile]	filter-file for searching"<<endl
//		<<"[-userdefinedfilterfile]	user-defined-filter-file for searching"<<endl
//		<<"[-filterbegin]	begin position for filter in training file"<<endl
//		<<"[-filterend]		end position for filter in training file"<<endl
//		<<"[-streamline]	to do hmmfilter search and whole structure search in a stream line"<<endl
//		<<"[-rtfn]	real training file name"<<endl
//		<<"[-rgfn]	real genome file name"<<endl
		<<"[-pcnt]	pseudocount_value [DEFAULT=0.001]"<<endl
		<<"[-k]		k_value [DEFAULT=10]"<<endl
		<<"[-th]	threshold_value_for_candidate_hit [DEFAULT=0.0]"<<endl
		<<"[-no]	number of nts in stem when overlap is allowed [DEFAULT=2]"<<endl
		<<"[-mc]	to merge candidate in preprocessing"<<endl
		<<"[-ms]	when taking the merge strategy, take the candidate with the (m)ax score or the one with the (s)hortest length (m|s) [DEFAULT=s]"<<endl
		<<"[-ns]	number of shift pos allowed in merge strategy [DEFAULT=0]"<<endl
		<<"[-ni]	number of insertion allowed in null loop [DEFAULT=3]"<<endl
		<<"[-pcv]	prior coef value [DEFAULT=0.4]"<<endl
		<<"[-pf]	prior_file (DEFAULT='./base_pair_prior.txt')"<<endl
		<<"[-pv]	pcoeff value [DEFAULT=2.0]"<<endl
		<<"[-st]	split threshold [DEFAULT=6]"<<endl
		<<"[-js]	stepsize in jump strategy [DEFAULT=1]"<<endl
		<<"[-jt]	score_filtering threshold in jump strategy [DEFAULT=0.0]"<<endl
		<<"[-r]		to also search reverse complement strand"<<endl
//		<<"[-ps]	to print structure alignment info"<<endl
		<<"[-d]		to print out the debug info"<<endl
//		<<"[-time]	to print out the time info"<<endl
		<<"[-pscore]	to print score info for stem and loop"<<endl	//
		<<"[-dfilter]	to print out the debug info for the filter based search result that generates the whole search hit"<<endl
		<<"[-ddpinput]	to print out the debug info of the input in dp search"<<endl
		<<"[-ddptree]	to print out the debug info of the tree in dp search"<<endl
		<<"[-ddpwin]	to print out the debug info of the window in dp search"<<endl
		<<"[-ddptd]		to print out the debug info of the td in dp search"<<endl
		<<"[-ddpsearch]	to print out the debug info of the dp search"<<endl;
	exit(1);
}
