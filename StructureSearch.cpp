// StructureSearch.cpp: implementation of the StructureSearch class.
//
//////////////////////////////////////////////////////////////////////

#include "StructureSearch.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

StructureSearch::StructureSearch()
{

}

StructureSearch::~StructureSearch()
{

}

/*
Input parameters:
				train data
				genome data
				hmm model
*/
void StructureSearch::hmm_search(int	top_k,
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
								int		iFlagSearchReverse,
								int		iFlagShowTimeInfo,
								int		iFlagPrintStrAlignInfo,
								int		iFlagPrintScoreInfo,		//
								int		iFlagPrintDebugInfo,
								int		iFlagPrintDebugInputInfo,
								int		iFlagPrintDebugTreeInfo,
								int		iFlagPrintDebugDPWinInfo,
								int		iFlagPrintDebugTDInfo,
								int		iFlagPrintDebugDPSearchInfo,
								char	cOffsetMethod,
								int		iSplitThreshold,
								char *	train_file,
								int		iFlagAutoFilter,
								char *	filter_file,	//
								int		filterbegin,	//
								int		filterend,		//
								double	pdcnt,
								double	priorh,
								char *	prior_file,
								int		genome_num,
								int	*	genome_length,
								char **	genome_name,
								char **	genome_sequence,
								char **	rc_genome_sequence,
								double **genome_baseFreq,
								int		runmode,
								int		iFlagPrintFilterInfo,	//
								int		iFlagStreamLine)
{
	int		iterate = 0;

	float	total_cpu_time = 0.0;
	float	hmm_search_time = 0.0;
	int		search_count = 0;	
	//num_whole_structure_hit 20080906
	int		num_whole_structure_hit = 0;	

	for(iterate=0; iterate<genome_num; iterate++)
	{
//		search_count = 0;

		//nc search
		hmm_search_time = 0.0;
		this->hmm_search_oneway(
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
							dScoreThresholdInJump,
							GENOME_SEARCH_DIRECTION_PLUS,
							iFlagShowTimeInfo,
							iFlagPrintStrAlignInfo,
							iFlagPrintScoreInfo,		//
							iFlagPrintDebugInfo,
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
							pdcnt,
							priorh,
							prior_file,
							genome_num,
							genome_length[iterate],
							genome_name[iterate],
							genome_sequence[iterate],
//							rc_genome_sequence[iterate],
							genome_baseFreq[iterate][0],	//A
							genome_baseFreq[iterate][1],	//C
							genome_baseFreq[iterate][2],	//G
							genome_baseFreq[iterate][3],	//U
							runmode,
							iFlagPrintFilterInfo,	//
							iFlagStreamLine,
							hmm_search_time,
							search_count);
		total_cpu_time += hmm_search_time;

		if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_BOTH)
		{
			hmm_search_time = 0.0;
			this->hmm_search_oneway(
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
								dScoreThresholdInJump,
								GENOME_SEARCH_DIRECTION_MINUS,	//iFlagSearchReverse,
								iFlagShowTimeInfo,
								iFlagPrintStrAlignInfo,
								iFlagPrintScoreInfo,		//
								iFlagPrintDebugInfo,
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
								pdcnt,
								priorh,
								prior_file,
								genome_num,
								genome_length[iterate],
								genome_name[iterate],
//								genome_sequence[iterate],
								rc_genome_sequence[iterate],
								genome_baseFreq[iterate][3],	//U
								genome_baseFreq[iterate][2],	//G
								genome_baseFreq[iterate][1],	//C
								genome_baseFreq[iterate][0],	//A
								runmode,
								iFlagPrintFilterInfo,	//
								iFlagStreamLine,
								hmm_search_time,
								search_count);
			total_cpu_time += hmm_search_time;
		}
	}

	if(iFlagShowTimeInfo){
		cout<<endl;
		printAdditionalInfo();
		printTotalHitNum(search_count);	//20080915
		printTotalTimeInfo(total_cpu_time);
		printSearchEndingTimeInfo();
		printOneStarLine();
	}
}

void StructureSearch::hmm_search_oneway(	
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
					int		iFlagSearchReverse,
					int		iFlagShowTimeInfo,
					int		iFlagPrintStrAlignInfo,
					int		iFlagPrintScoreInfo,		//
					int		iFlagPrintDebugInfo,
					int		iFlagPrintDebugInputInfo,
					int		iFlagPrintDebugTreeInfo,
					int		iFlagPrintDebugDPWinInfo,
					int		iFlagPrintDebugTDInfo,
					int		iFlagPrintDebugDPSearchInfo,
					char	cOffsetMethod,
					int		iSplitThreshold,
					char *	train_file,
					int		iFlagAutoFilter,
					char *	filter_file,	//
					int		filterbegin,	//
					int		filterend,		//
					double	pdcnt,
					double	priorh,
					char *	prior_file,
					int		genome_num,
					int	 	genome_length,
					char *	genome_name,
					char *	genome_sequence,
//					char *	rc_genome_sequence,
					double	genome_baseFreq_A,
					double	genome_baseFreq_C,
					double	genome_baseFreq_G,
					double	genome_baseFreq_U,
					int		runmode,
					int		iFlagPrintFilterInfo,	//
					int		iFlagStreamLine,
					float &	hmm_search_time,
					int	 &	count)
{
	int		i, j;
	clock_t start = 0, end = 0;
	float total_cpu_time = 0.0;
	float structure_search_time = 0.0;

	if(iFlagShowTimeInfo)
		start = clock();

	double	pdcntgap= 0.001;
	int		numlp	= 0;
	double	probv = 0.0;

//	int		pre_hit_num = 0;	//for streamline search

	DataLoader trainLoader;
	trainLoader.setOffsetMethod(cOffsetMethod);
	trainLoader.setSplitThreshold(iSplitThreshold);
	trainLoader.inputData(train_file);

	DataLoader filterLoader;
	if(runmode == RUNNINGMODE_WEBSERVER)
		filterLoader.inputData(filter_file);
	else	//RUNNINGMODE_STANDALONE
	{
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

//				if(iFlagPrintFilterInfo)
//				{
//					for (i=0; i<(int)cc.size(); i++)
//					{
//						for (j=0; j<(int)cc[0].length(); j++)
//							cout<<cc[i][j];
//						cout<<endl;
//					}
//					cout<<"filterbegin="<<filterbegin<<endl;
//					cout<<"filterend  ="<<filterend<<endl;
//				}

				filterLoader.inputData(cc);
			}
		}
		else
		{
			//get one of the training sequences
			int ** trainingSeqs = trainLoader.getPointerOfSeqs();
			int length = trainLoader.getSeqLength();
			string one_training_seq;
			for(i=0; i<length; i++)
				one_training_seq += numtochar(trainingSeqs[0][i]);

			//
			filterLoader.inputData(filter_file);
			int ** filteringSeqs = filterLoader.getPointerOfSeqs();
			length = filterLoader.getSeqLength();
			string one_filtering_seq;
			for(i=0; i<length; i++)
				one_filtering_seq += numtochar(filteringSeqs[0][i]);

			filterbegin = one_training_seq.find(one_filtering_seq);
			filterend	= filterbegin + one_filtering_seq.length();

//			if(iFlagPrintFilterInfo)
//			{
//				cout<<"filterbegin="<<filterbegin<<endl;
//				cout<<"filterend  ="<<filterend<<endl;
//			}
		}
	}

	int algpos[2], loopNo=0;
	int winsize = 0;
	int original_winsize = filterLoader.getScanWinLength(3, '+');

	//should call other functions to get these two values
	int extension_left	= 0;
	int extension_right = 0;
	int extension_pos_start = 0;
	int extension_pos_end	= 0;

	int	hitlength = 0;

	filterLoader.setBaseFreq(genome_baseFreq_A, 
							genome_baseFreq_C, 
							genome_baseFreq_G, 
							genome_baseFreq_U);

	HMMBuilder loopbd;
	loopbd.setPseudocnt(pdcnt, pdcntgap);
	loopbd.buildHMM(filterLoader);
	Loop *plp=loopbd.getAllLoops(numlp);
    
	winsize = original_winsize*2;     //note: make the real windowsize is long enough
	loopNo	= 0;
	probv	= INVLDPROB;
	for(i=0; i<2; i++) {
		algpos[i]=-1;
	}
   
	char *pStr = new char[winsize+1];
    
	double	optimal_score;
	int		pre_end_pos = -1;
	int		optimal_pos_start = -1, optimal_pos_end = -1;	//for hmm-filter
	optimal_score = SMALLEST;

	FVTBSearch vtber(1);
	vtber.setHitThreshold(0.0);
	char *retAlg=NULL, *retSeq=NULL;
	char *pre_retAlg=NULL, *pre_retSeq=NULL;	//

	//store the subgenome based on hmmfilter result
	vector <string> subgenome_based_on_hmmfilter;

	char * c_subgenome_seq = NULL;
	int tmp_count = 0;	//just for the dfilter

	for(i=0; i<=genome_length-winsize; i++)
	{
		//1. set the search sequence
		for(j=0; j<winsize; j++)
			pStr[j] = genome_sequence[i+j];
		pStr[j] = '\0';
		//2. search for all candidates
		FVTBCandidate * fcd = vtber.SearchOneNoMerge( &plp[loopNo], pStr );
        
		if(fcd == NULL)  //no hit
			continue;
		else 
		{             //find hit
			fcd->getPosition(algpos);
			fcd->getBothShortSeqs(&retSeq, &retAlg);
			algpos[0]=algpos[0]+i;
			algpos[1]=algpos[1]+i;
			probv = fcd->getScore();

			//merging the result
			if(optimal_score == SMALLEST)
			{
				//the first time we meet the candidate hit.
				if(optimal_score < probv) {
					optimal_score		= probv;
					optimal_pos_start	= algpos[0];
					optimal_pos_end		= algpos[1];
					pre_end_pos			= algpos[1];	//

					pre_retAlg = new char[strlen(retAlg)+1];
					strcpy(pre_retAlg, retAlg);
					pre_retSeq = new char[strlen(retSeq)+1];
					strcpy(pre_retSeq, retSeq);
				}
			}
			else
			{
				//not the first time of meeting the candidate hit
				if(algpos[0] < pre_end_pos)	//pointing to the same hit
				{
					//belongs to the same candidate hit
					if(optimal_score < probv) {
						optimal_score		= probv;
						optimal_pos_start	= algpos[0];
						optimal_pos_end		= algpos[1];
						pre_end_pos			= algpos[1];

						if( pre_retSeq != NULL ) {   delete [ ] pre_retSeq;  pre_retSeq = NULL;  }
						if( pre_retAlg != NULL ) {   delete [ ] pre_retAlg;  pre_retAlg = NULL;  }
						pre_retAlg = new char[strlen(retAlg)+1];
						strcpy(pre_retAlg, retAlg);
						pre_retSeq = new char[strlen(retSeq)+1];
						strcpy(pre_retSeq, retSeq);
					}
				}
				else	//other candidate hit is coming
				{
					hitlength = optimal_pos_end - optimal_pos_start + 1;
					trainLoader.getBothExtLength(filterbegin, filterend, hitlength, extension_left, extension_right);
					count++;
/*
					//20080909
					if(iFlagPrintFilterInfo) {
						printFilterHitIndex(count);
					}
*/					
					if(optimal_pos_start > extension_left)
					{
						extension_pos_start = optimal_pos_start-extension_left;
						extension_pos_end	= optimal_pos_end+extension_right;
						c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
/*
						if(iFlagPrintFilterInfo)
						{
							cout<<">"<<genome_name<<"("<<optimal_pos_start<<"-"<<optimal_pos_end<<")"<<"["<<extension_pos_start<<"-"<<extension_pos_end<<"]"<<endl;
							if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_PLUS)
								cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
							else
								cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
							printHitPos(optimal_pos_start, optimal_pos_end);
						}
*/
					}
					else
					{
						extension_pos_start = 0;
						extension_pos_end	= optimal_pos_end+extension_right;
						c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
/*							
						if(iFlagPrintFilterInfo)
						{
							cout<<">"<<genome_name<<"("<<optimal_pos_start<<"-"<<optimal_pos_end<<")"<<"["<<extension_pos_start<<"-"<<extension_pos_end<<"]"<<endl;
							if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_PLUS)
								cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
							else
								cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
							printHitPos(optimal_pos_start, optimal_pos_end);
						}
*/
					}
/*
					if(iFlagPrintFilterInfo)
					{
						subgenome_based_on_hmmfilter.push_back(string(c_subgenome_seq));

						printHitScore(optimal_score);
						printFilterHitAlignment(pre_retSeq, pre_retAlg);

						printFilterHitExtenPos(extension_pos_start, extension_pos_end);
						printFilterHitExtenNtsHeader();
						printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE);
					}
*/
					//do the whole structure search
					if(iFlagStreamLine) {
						count -= 1;
						if(iFlagPrintFilterInfo) {
							tmp_count = count;
						}
						end = clock();
						float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
						cpu_time = cpu_time/3600;
						total_cpu_time += cpu_time;

						structure_search_time = 0.0;
						this->structure_search_based_on_one_hmmfilterresult(top_k,
																			threshold,
																			num_nts_overlap_between_stem,
																			iFlagMergeCandInPreprocess,	//1,//
																			iFlagCandwithShortestLength,
																			iShiftNumMergeCand,
																			iAllowedNullLoopInsNum,
																			pcoeff,
																			iJumpStrategy,
																			iStepSize,
																			dScoreThresholdInJump,
																			iFlagSearchReverse,	//PLUS | MINUS
																			iFlagShowTimeInfo,
																			iFlagPrintStrAlignInfo,
																			iFlagPrintScoreInfo,		//
																			iFlagPrintDebugInfo,
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
																			-1,//filterbegin,	//
																			-1,//filterend,		//
																			pdcnt,
																			priorh,
																			prior_file,
//																			hitlength,
																			genome_name,
																			c_subgenome_seq,
																			optimal_pos_start,
																			extension_left,
																			optimal_pos_end,
																			extension_right,
																			runmode,
																			count,	//pre_hit_num,			//
																			structure_search_time);	//
						total_cpu_time += structure_search_time;
						start = clock();
//						count += pre_hit_num;	//20080906
						if(iFlagPrintFilterInfo) {
							if(tmp_count == count) {
								//no whole search hit found
								count = tmp_count + 1;
							} else if(count == tmp_count + 1) {
								//if whole structure hit is found, output its mapping hmm filter hit

								printHitScore(optimal_score);
								printFilterHitAlignment(pre_retSeq, pre_retAlg);

								printFilterHitExtenPos(extension_pos_start, extension_pos_end);
								printFilterHitExtenNtsHeader();
								printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE);

								if(iFlagPrintDebugInfo) {
									cout<<endl;
								}
							}
						}
					}

					if(c_subgenome_seq != NULL) {
						delete [] c_subgenome_seq;
						c_subgenome_seq = NULL;
					}

					if( retSeq != NULL ) {   delete [ ] retSeq;  retSeq = NULL;  }
					if( retAlg != NULL ) {   delete [ ] retAlg;  retAlg = NULL;  }

					optimal_score		= probv;
					optimal_pos_start	= algpos[0];
					optimal_pos_end		= algpos[1];
					pre_end_pos			= algpos[1];
				}
			}
			delete [] fcd;
		}
	}
	delete [] pStr;		pStr = NULL;

	//20080709 segmentation fault
	if(optimal_pos_start != -1 && optimal_pos_end != -1)
	{
		hitlength = optimal_pos_end - optimal_pos_start + 1;
		trainLoader.getBothExtLength(filterbegin, filterend, hitlength, extension_left, extension_right);
		count++;
/*
		//20080909
		if(iFlagPrintFilterInfo) {
			printFilterHitIndex(count);
		}
*/
		if(optimal_pos_start > extension_left)
		{
			extension_pos_start = optimal_pos_start-extension_left;
			extension_pos_end	= optimal_pos_end+extension_right;
			c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
/*				
			if(iFlagPrintFilterInfo)
			{
				cout<<">"<<genome_name<<"("<<optimal_pos_start<<"-"<<optimal_pos_end<<")"<<"["<<extension_pos_start<<"-"<<extension_pos_end<<"]"<<endl;
				if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_PLUS)
					cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
				else
					cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
				printHitPos(optimal_pos_start, optimal_pos_end);
			}
*/
		}
		else
		{
			extension_pos_start = 0;
			extension_pos_end	= optimal_pos_end+extension_right;
			c_subgenome_seq = strextract(genome_sequence, 0, optimal_pos_end+extension_right);
/*				
			if(iFlagPrintFilterInfo)
			{
				cout<<">"<<genome_name<<"("<<optimal_pos_start<<"-"<<optimal_pos_end<<")"<<"["<<extension_pos_start<<"-"<<extension_pos_end<<"]"<<endl;
				if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_PLUS)
					cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
				else
					cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
				printHitPos(optimal_pos_start, optimal_pos_end);
			}
*/
		}
/*
		if(iFlagPrintFilterInfo)
		{
			subgenome_based_on_hmmfilter.push_back(string(c_subgenome_seq));

			printHitScore(optimal_score);
			printFilterHitAlignment(pre_retSeq, pre_retAlg);

			printFilterHitExtenPos(extension_pos_start, extension_pos_end);
			printFilterHitExtenNtsHeader();
			printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE);
		}
*/
		//do the whole structure search
		if(iFlagStreamLine) {
			count -= 1;
			if(iFlagPrintFilterInfo) {
				tmp_count = count;
			}
			end = clock();
			float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
			cpu_time = cpu_time/3600;
			total_cpu_time += cpu_time;

			structure_search_time = 0.0;
			this->structure_search_based_on_one_hmmfilterresult(top_k,
																threshold,
																num_nts_overlap_between_stem,
																iFlagMergeCandInPreprocess,	//1,//
																iFlagCandwithShortestLength,
																iShiftNumMergeCand,
																iAllowedNullLoopInsNum,
																pcoeff,
																iJumpStrategy,
																iStepSize,
																dScoreThresholdInJump,
																iFlagSearchReverse,	//PLUS | MINUS
																iFlagShowTimeInfo,
																iFlagPrintStrAlignInfo,
																iFlagPrintScoreInfo,		//
																iFlagPrintDebugInfo,
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
																-1,//filterbegin,	//
																-1,//filterend,		//
																pdcnt,
																priorh,
																prior_file,
//																hitlength,
																genome_name,
																c_subgenome_seq,
																optimal_pos_start,
																extension_left,
																optimal_pos_end,
																extension_right,
																runmode,
																count,	//pre_hit_num,			//
																structure_search_time);	//
			total_cpu_time += structure_search_time;
			start = clock();
//			count += pre_hit_num;	//20080906
			if(iFlagPrintFilterInfo) {
				if(tmp_count == count) {
					//no whole search hit found
					count = tmp_count + 1;
				} else if(count == tmp_count + 1)
				{
					//if whole structure hit is found, output its mapping hmm filter hit

					printHitScore(optimal_score);
					printFilterHitAlignment(pre_retSeq, pre_retAlg);

					printFilterHitExtenPos(extension_pos_start, extension_pos_end);
					printFilterHitExtenNtsHeader();
					printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE);

					if(iFlagPrintDebugInfo) {
						cout<<endl;
					}
				}
			}
		}

		if(c_subgenome_seq != NULL) {
			delete [] c_subgenome_seq;
			c_subgenome_seq = NULL;
		}
	}

	if( retSeq != NULL ) {   delete [ ] retSeq;  retSeq = NULL;  }
	if( retAlg != NULL ) {   delete [ ] retAlg;  retAlg = NULL;  }

	if( pre_retSeq != NULL ) {   delete [ ] pre_retSeq;  pre_retSeq = NULL;  }
	if( pre_retAlg != NULL ) {   delete [ ] pre_retAlg;  pre_retAlg = NULL;  }

	loopbd.freeAllLoops();

	if(iFlagShowTimeInfo)
	{
		end = clock();
		float cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
		cpu_time = cpu_time/3600;
		total_cpu_time += cpu_time;
		hmm_search_time += total_cpu_time;	//20080906
	}
}

void StructureSearch::structure_search_based_on_one_hmmfilterresult(int		top_k,
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
																	int		iFlagSearchReverse,
																	int		iFlagShowTimeInfo,
																	int		iFlagPrintStrAlignInfo,
																	int		iFlagPrintScoreInfo,		//
																	int		iFlagPrintDebugInfo,
																	int		iFlagPrintDebugInputInfo,
																	int		iFlagPrintDebugTreeInfo,
																	int		iFlagPrintDebugDPWinInfo,
																	int		iFlagPrintDebugTDInfo,
																	int		iFlagPrintDebugDPSearchInfo,
																	char	cOffsetMethod,
																	int		iSplitThreshold,
																	char *	train_file,
																	int		iFlagAutoFilter,
																	char *	filter_file,	//
																	int		filterbegin,	//
																	int		filterend,		//
																	double	pdcnt,
																	double	priorh,
																	char *	prior_file,
//																	int		genome_length,
																	char *	genome_name,
																	char *	genome_segment,
																	int		optimal_pos_start,
																	int		extension_left,
																	int		optimal_pos_end,
																	int		extension_right,
																	int		runmode,
																	int	 &	pre_hit_num,
																	float &	structure_hmm_search_time)
{
	clock_t start_all, end_all;
	float total_time = 0.0;
	start_all = clock();

	int j, genome_length;
	//prepare for the whole structure search
	int		inside_genome_num = 1;
				
	int	*	inside_genome_type;
	inside_genome_type = new int[1];
	inside_genome_type[0] = GENOME_SEGMENT;

	int *	inside_genome_extleft;
	inside_genome_extleft = new int[1];

	int *	inside_genome_extright;
	inside_genome_extright = new int[1];

	if(optimal_pos_start > extension_left) {
		inside_genome_extleft[0] = optimal_pos_start-extension_left;
		inside_genome_extright[0] = optimal_pos_end+extension_right;
	} else {
		inside_genome_extleft[0] = 0;
		inside_genome_extright[0] = optimal_pos_end+extension_right;
	}

//	genome_length = inside_genome_extright[0] - inside_genome_extleft[0];// + 1;
	genome_length = inside_genome_extright[0] - inside_genome_extleft[0] + 1;

	int	*	inside_genome_length;
	inside_genome_length = new int[1];
	inside_genome_length[0] = genome_length;//hitlength;

//	cout<<genome_segment<<endl;

	char **	inside_genome_name;
	inside_genome_name = new char*[1];
	inside_genome_name[0] = new char[strlen(genome_name)+1];//genome_name[iterate])+1];
	strcpy(inside_genome_name[0], genome_name);//genome_name[iterate]);

	char **	inside_genome_sequence;
	inside_genome_sequence = new char*[1];
	inside_genome_sequence[0] = new char[strlen(genome_segment)+1];//c_subgenome_seq)+1];
	strcpy(inside_genome_sequence[0], genome_segment);//c_subgenome_seq);

	char ch;

	double ** inside_genome_baseFreq;
	inside_genome_baseFreq = new double*[1];
	inside_genome_baseFreq[0] = new double[4];
	int numbase[5]={0,0,0,0,0};
	for(j=0; j<genome_length; j++)
	{
		ch = genome_segment[j];//c_subgenome_seq[j];
		if(	ch=='a' || ch=='A' ||
			ch=='c' || ch=='C' ||
			ch=='g' || ch=='G' ||
			ch=='u' || ch=='U' ||
			ch=='t' || ch=='T') 
		{
			numbase[chartonum(ch)]++;
		}
	}
	int allbasenum = numbase[1]+numbase[2]+numbase[3]+numbase[4];
	for(j=1; j<5; j++) {
		inside_genome_baseFreq[0][j-1]=(double)numbase[j]/allbasenum;
	}  

	end_all = clock();
	total_time = ((float) (end_all - start_all)) / CLOCKS_PER_SEC;
	total_time = total_time/3600;

	float structure_search_time = 0.0;
	//here the whole structure search should be executed. 20080803
	//this->structure_search();
	this->structure_search(	top_k,
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
							iFlagPrintScoreInfo,		//
							iFlagPrintDebugInfo,
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
							pdcnt,
							priorh,
							prior_file,
							inside_genome_num,
							inside_genome_type,
							NULL,
							inside_genome_length,
							inside_genome_name,
							inside_genome_sequence,
							NULL,	//inside_rc_genome_sequence,
							inside_genome_baseFreq,
							inside_genome_extleft,
							inside_genome_extright,
							0.0,	//filtertime
							runmode,	//runmode,
							SEARCH_WHOLESTRUCTURE,	//SEARCH_STREAMLINE);//
							1,				//streamline_flag
							pre_hit_num,
							structure_search_time);	//
	
	total_time += structure_search_time;
	structure_hmm_search_time += total_time;

	if(inside_genome_type != NULL) {
		delete [] inside_genome_type;
		inside_genome_type = NULL;
	}

	if(inside_genome_length != NULL) {
		delete [] inside_genome_length;
		inside_genome_length = NULL;
	}

	if(inside_genome_name != NULL) {
		delete [] inside_genome_name[0];
		inside_genome_name[0] = NULL;
		delete [] inside_genome_name;
		inside_genome_name = NULL;
	}

	if(inside_genome_sequence != NULL) {
		delete [] inside_genome_sequence[0];
		inside_genome_sequence[0] = NULL;
		delete [] inside_genome_sequence;
		inside_genome_sequence = NULL;
	}

	if(inside_genome_baseFreq != NULL) {
		delete [] inside_genome_baseFreq[0];
		inside_genome_baseFreq[0] = NULL;
		delete [] inside_genome_baseFreq;
		inside_genome_baseFreq = NULL;
	}

	if(inside_genome_extleft != NULL) {
		delete [] inside_genome_extleft;
		inside_genome_extleft = NULL;
	}

	if(inside_genome_extright != NULL) {
		delete [] inside_genome_extright;
		inside_genome_extright = NULL;
	}
}
/*
Compared to hmm_search, here trainLoader is totally different from the previous trainLoader
Here filterLoader is the same as the trainLoader in hmm_search.
*/
void StructureSearch::structure_search(	int		top_k,
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
										int		iFlagSearchReverse,
										int		iFlagShowTimeInfo,
										int		iFlagPrintStrAlignInfo,
										int		iFlagPrintScoreInfo,		//
										int		iFlagPrintDebugInfo,
										int		iFlagPrintDebugInputInfo,
										int		iFlagPrintDebugTreeInfo,
										int		iFlagPrintDebugDPWinInfo,
										int		iFlagPrintDebugTDInfo,
										int		iFlagPrintDebugDPSearchInfo,
										char	cOffsetMethod,
										int		iSplitThreshold,
										char *	train_file,
										int		iFlagAutoFilter,
										char *	filter_file,	//
										int		filterbegin,	//
										int		filterend,		//
										double	pdcnt,
										double	priorh,	//
										char *	prior_file,	//
										int		genome_num,
										int	*	genome_type,
										int	*	genome_direction,
										int	*	genome_length, 
										char **	genome_name,
										char **	genome_sequence, 
										char **	rc_genome_sequence, 
										double ** genome_baseFreq,
										int *	genome_extleft,
										int *	genome_extright,
										float	filtertime,
										int		runmode,
										int		searchtype,
										int		streamline_flag,
										int	&	pre_hit_num,		//
										float & struct_search_time)	//
{
/*
	//for the debug 20081010
	///////////////////////////////////////////
	cout<<"top_k = "<<top_k
		<<"threshold = "<<threshold<<endl
		<<"num_nts_overlap_between_stem = "<<num_nts_overlap_between_stem<<endl
		<<"iFlagMergeCandInPreprocess = "<<iFlagMergeCandInPreprocess<<endl
		<<"iFlagCandwithShortestLength = "<<iFlagCandwithShortestLength<<endl
		<<"iShiftNumMergeCand = "<<iShiftNumMergeCand<<endl
		<<"iAllowedNullLoopInsNum = "<<iAllowedNullLoopInsNum<<endl
		<<"pcoeff = "<<pcoeff<<endl
		<<"iJumpStrategy = "<<iJumpStrategy<<endl
		<<"iStepSize = "<<iStepSize<<endl
		<<"dScoreThresholdInJump = "<<dScoreThresholdInJump<<endl
		<<"iFlagSearchReverse = "<<iFlagSearchReverse<<endl
		<<"iFlagShowTimeInfo = "<<iFlagShowTimeInfo<<endl
		<<"iFlagPrintStrAlignInfo = "<<iFlagPrintStrAlignInfo<<endl
		<<"iFlagPrintScoreInfo = "<<iFlagPrintScoreInfo<<endl
		<<"iFlagPrintDebugInfo = "<<iFlagPrintDebugInfo<<endl
		<<"iFlagPrintDebugInputInfo = "<<iFlagPrintDebugInputInfo<<endl
		<<"iFlagPrintDebugTreeInfo = "<<iFlagPrintDebugTreeInfo<<endl
		<<"iFlagPrintDebugDPWinInfo = "<<iFlagPrintDebugDPWinInfo<<endl
		<<"iFlagPrintDebugTDInfo = "<<iFlagPrintDebugTDInfo<<endl
		<<"iFlagPrintDebugDPSearchInfo = "<<iFlagPrintDebugDPSearchInfo<<endl
		<<"cOffsetMethod = "<<cOffsetMethod<<endl
		<<"iSplitThreshold = "<<iSplitThreshold<<endl
		<<"train_file = "<<train_file<<endl
		<<"iFlagAutoFilter = "<<iFlagAutoFilter<<endl;

	if(filter_file != NULL)
		cout<<"filter_file = "<<filter_file<<endl;
	
	cout<<"filterbegin = "<<filterbegin<<endl
		<<"filterend = "<<filterend<<endl
		<<"pdcnt = "<<pdcnt<<endl
		<<"priorh = "<<priorh<<endl
		<<"prior_file = "<<prior_file<<endl;

	cout<<"genome_num = "<<genome_num<<endl;
	for(int i=0; i<genome_num; i++)
	{
		cout<<"genome_type["<<i<<"] = "<<genome_type[i]<<endl;
//		cout<<"genome_direction["<<i<<"] = "<<genome_direction[i]<<endl;
		cout<<"genome_length["<<i<<"] = "<<genome_length[i]<<endl;
		cout<<"genome_name["<<i<<"] = "<<genome_name[i]<<endl;
		cout<<"genome_sequence["<<i<<"] = "<<genome_sequence[i]<<endl;
//		cout<<"rc_genome_sequence["<<i<<"] = "<<rc_genome_sequence[i]<<endl;
		for(int j=0; j<4; j++)
		{
			cout<<"genome_baseFreq["<<i<<"]["<<j<<"] = "<<genome_baseFreq[i][j]<<endl;
		}
		cout<<"genome_extleft["<<i<<"] = "<<genome_extleft[i]<<endl
			<<"genome_extright["<<i<<"] = "<<genome_extright[i]<<endl;
	}

	cout<<"filtertime = "<<filtertime<<endl
		<<"runmode = "<<runmode<<endl
		<<"searchtype = "<<searchtype<<endl
		<<"streamline_flag = "<<streamline_flag<<endl
		<<"pre_hit_num = "<<pre_hit_num<<endl
		<<"struct_search_time = "<<struct_search_time<<endl;
	////////////////////////////////////////////////////////
*/
	int iterate = 0;

	//Tree decomposition Part
	Tree_bag *	overall_root;
	DynamicP dp;

	dp.setSearchParams(	top_k,				//[k value]
						threshold,			//[threshold]
						num_nts_overlap_between_stem,	//[number of nts in stem when overlap is allowed]
						iFlagMergeCandInPreprocess,		//whether taking the merge-candidate strategy in preprocessing
						iFlagCandwithShortestLength,
						iShiftNumMergeCand,			//[num of shift when merge candidate in preprocessing]
						iAllowedNullLoopInsNum,		//[number of insertion allowed in null loop]
						pcoeff,				//pcoeff
						iJumpStrategy,		//[skip-and-jump strategy]
						iStepSize,			//[stepsize in skip-and-jump strategy]
						dScoreThresholdInJump,	//[score_filtering threshold in jump strategy]
						iFlagPrintStrAlignInfo,	//[structure alignment info -s|-n]
						iFlagPrintScoreInfo,	//
						iFlagPrintDebugInfo,	//[debug info -n|-d]
						iFlagPrintDebugInputInfo,
						iFlagPrintDebugTreeInfo,
						iFlagPrintDebugDPWinInfo,
						iFlagPrintDebugTDInfo,
						iFlagPrintDebugDPSearchInfo);

	dp.setPreHitNum(pre_hit_num);	//

	DataLoader trainLoader;
	trainLoader.setOffsetMethod(cOffsetMethod);
	trainLoader.setSplitThreshold(iSplitThreshold);

	if(searchtype == SEARCH_SUBSTRUCTUREFILTER)	//sub-structure filter search
	{
		//filter file is the training file for sub-structure filter search
		trainLoader.inputData(filter_file);

		if(runmode == RUNNINGMODE_WEBSERVER)
		{
			//do nothing
		}
		else	//RUNNINGMODE_STANDALONE
		{
			if(iFlagAutoFilter == 1)
			{
				//calling substructure filter selection
//				filterbegin = filterselection.getbegposi();
//				filterend	= filterselection.getendposi();
			}
			else
			{
				//get one of the training sequences
				DataLoader fulltrainLoader;
				fulltrainLoader.setOffsetMethod(cOffsetMethod);
				fulltrainLoader.setSplitThreshold(iSplitThreshold);
				fulltrainLoader.inputData(train_file);

				int ** trainingSeqs = fulltrainLoader.getPointerOfSeqs();
				int length = fulltrainLoader.getSeqLength();
				string one_training_seq;
				int i;
				for(i=0; i<length; i++)
					one_training_seq += numtochar(trainingSeqs[0][i]);

				//
//				trainLoader.inputData(filter_file);	//
				int ** filteringSeqs = trainLoader.getPointerOfSeqs();
				length = trainLoader.getSeqLength();
				string one_filtering_seq;
				for(i=0; i<length; i++)
					one_filtering_seq += numtochar(filteringSeqs[0][i]);

				filterbegin = one_training_seq.find(one_filtering_seq);
				filterend	= filterbegin + one_filtering_seq.length();
			}
		}
	}
	else	//whole-structure search
	{
		trainLoader.inputData(train_file);
	}

	//num of pastaline	20081015
	dp.setNumPastaLine(trainLoader.getPastaLines());

	int numOfTrainSeq = trainLoader.getNumOfTrainSeq();
	int length = trainLoader.getSeqLength();
	int ** trainingSeqs = trainLoader.getPointerOfSeqs();

	dp.setTrainingSeqs(numOfTrainSeq, length, trainingSeqs);
	if(dp.getFlagPrintDebugInputInfo())
		dp.printTrainingSeqs();
					
	dp.setStems(trainLoader.getNumstem(), trainLoader.getStemarray());
	if(dp.getFlagPrintDebugInputInfo())
		dp.printStems();

	dp.setLoops(trainLoader.getNumloop(), trainLoader.getLooparray());
	if(dp.getFlagPrintDebugInputInfo())
		dp.printLoops();

	//Tree Decomposition Generator
	dp.build_tree();

	dp.identify_loopneighbor();
	dp.build_node_mapping();

	//add the pParent in every tree node.
	overall_root = dp.getMyTDRoot();
	overall_root->pParent = overall_root;
	dp.setTreeBagParent(overall_root);

	//re-arrange the node-list
	overall_root->shared_node_num = overall_root->nodenum;
	dp.reArrangeChildNodelist(overall_root);

	if(dp.getFlagPrintDebugTreeInfo()) {
		printf("<tree decomposition>\n");
		dp.printTreebags(overall_root);
		printf("</tree decomposition>\n");
	}

	dp.buildStemIdxArray();
	if(dp.getFlagPrintDebugTDInfo())
		dp.printStemIdxArray();
	dp.allocateStemCandIdxArray();

	dp.preprocess_tree(overall_root);

	dp.buildTreeNodePath(overall_root);

	int winsize = 0, minsize = 0;
	float total_searching_time = 0.0;
	for(iterate=0; iterate<genome_num; iterate++)
	{
		if(genome_type[iterate] != GENOME_NONE)
		{
			trainLoader.setBaseFreq(genome_baseFreq[iterate][0], 
									genome_baseFreq[iterate][1], 
									genome_baseFreq[iterate][2], 
									genome_baseFreq[iterate][3]);

			//
			if(dp.getFlagPrintDebugDPWinInfo()) {
				cout<<genome_baseFreq[iterate][0]<<" | "
					<<genome_baseFreq[iterate][1]<<" | "
					<<genome_baseFreq[iterate][2]<<" | "
					<<genome_baseFreq[iterate][0]<<endl;
			}

			winsize = trainLoader.getScanWinLength(3, '-');
			minsize = trainLoader.getMiniWinLength( );

			int numsm=0,  numlp=0;

			ScfgBuilder	scfgBuilder;
			scfgBuilder.setPseudocnt(pdcnt, 0.001);
			//coefficient of prior matrix
			scfgBuilder.setPriorCoef(priorh);
			scfgBuilder.setPriorFilePath(prior_file);
			scfgBuilder.buildSCFG(trainLoader);
			Stem *	pStemModel = scfgBuilder.getAllStems(numsm);

			HMMBuilder	loopBuilder;
			loopBuilder.setPseudocnt(pdcnt, 0.001);
			loopBuilder.buildHMM(trainLoader);
			Loop *	pLoopModel = loopBuilder.getAllLoops(numlp);

			//load the random/genome sequence file
			dp.setGenomeSequence(genome_length[iterate], genome_sequence[iterate]);

			//
			//dp.setFlagSearchReverseGenome(iFlagSearchReverse);

			//load the reverse genome sequence.
			if(genome_type[iterate] == ORIGINAL_GENOME)
			{
				if(iFlagSearchReverse == GENOME_SEARCH_DIRECTION_BOTH) {
					dp.setFlagSearchReverseGenome(iFlagSearchReverse);
					dp.setRGenomeSequence(genome_length[iterate], rc_genome_sequence[iterate]);
				}
			} else {
				//GENOME_SEGMENT. filter search result
//				dp.setFlagSearchReverseGenome(genome_direction[iterate]);
				dp.setFlagSearchReverseGenome(iFlagSearchReverse);
			}

			dp.searchPK(winsize, 
						minsize, 
						pStemModel, 
						pLoopModel, 
						overall_root, 
						genome_extleft[iterate], 
						genome_name[iterate], 
						searchtype, 
//						fulltrainLoader, 
						train_file,			//
						cOffsetMethod,		//
						iSplitThreshold,	//
						filterbegin, 
						filterend, 
						genome_baseFreq[iterate]);

			loopBuilder.freeAllLoops();
			scfgBuilder.freeAllStems();

			if(iFlagShowTimeInfo == 1) {
//				printf("Time for Preprocessing:%f hours\n", dp.getCpuTimePreprocessing());
//				printf("Time for DP:%f hours\n", dp.getCpuTimeDP());

				if(searchtype == SEARCH_SUBSTRUCTUREFILTER) {
//					cout<<endl<<"Time for total filter search "<<dp.getCpuTimeAll()<<" hours"<<endl;
//					printCurrentTime();

					cout<<endl;
					printAdditionalInfo();
					printTotalHitNum(dp.getNumHit());	//20080915
					printTotalTimeInfo(dp.getCpuTimeAll());

					printSearchEndingTimeInfo();
//					cout<<"*: time: ";
//					printCurrentTime();
//					cout<<right<<setw(NUM_NTS_PER_LINE-49)<<"*"<<endl;

					printOneStarLine();

				}
				else{
//					if(dp.getNumHit() == 0) {
//						cout<<endl<<"---------------"<<endl;
//						cout<<"No Hit"<<endl;
//						cout<<genome_name[iterate]<<endl;
//					}
					total_searching_time += dp.getCpuTimeAll();
//					printf("Time for the current genome segment %f hours\n", dp.getCpuTimeAll());
				}
			}
		}
//		else
//		{
//			//No results have been found
//			cout<<endl
//				<<"---------------"<<endl
//				<<genome_name[iterate]<<endl
//				<<"No results have been found"<<endl;
//		}
	}

	total_searching_time += filtertime;
	struct_search_time = total_searching_time;	//
	
	if(dp.getNumHit() == 0)
	{
		if(!streamline_flag)
		{
//			cout<<endl<<"No Hit has been found"<<endl;	//20080915
		}
	}

	pre_hit_num += dp.getNumHit();

	if(!streamline_flag)
	{
		if(searchtype == SEARCH_WHOLESTRUCTURE) {	//20080813
//			cout<<endl<<"--------------------------------------------------"<<endl;
//			printf("Time for the whole search process %f hours\n", total_searching_time);
//			cout<<"Finished at ";printCurrentTime();

			cout<<endl;
			printAdditionalInfo();
			printTotalHitNum(pre_hit_num);	//20080915
			printTotalTimeInfo(total_searching_time);

			printSearchEndingTimeInfo();
//			cout<<"*: time: ";
//			printCurrentTime();
//			cout<<right<<setw(NUM_NTS_PER_LINE-49)<<"*"<<endl;
			printOneStarLine();
		}
	}

	dp.freeTreeNodePath();
	dp.freeStemIdxArray();
	dp.freeStemCandIdxArray();

	dp.free_node_mapping();	//node_mapping[] is needed in the scan_genome function
							//opposite to the previous dp.allocate_node_mapping();
	dp.removeTreeBagParent(overall_root);
	dp.postprocess_tree(overall_root);
	dp.free_tree(overall_root);
}

void StructureSearch::search(int	top_k,
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
							int		iFlagSearchReverse,
							int		iFlagShowTimeInfo,
							int		iFlagPrintStrAlignInfo,
							int		iFlagPrintScoreInfo,
							int		iFlagPrintDebugInfo,
							int		iFlagPrintDebugFilterInfo,
							int		iFlagPrintDebugInputInfo,
							int		iFlagPrintDebugTreeInfo,
							int		iFlagPrintDebugDPWinInfo,
							int		iFlagPrintDebugTDInfo,
							int		iFlagPrintDebugDPSearchInfo,
							char	cOffsetMethod,
							int		iSplitThreshold,
							char *	train_file,
							int		iFlagAutoFilter,
							char *	filter_file,	//
							int		filterbegin,	//
							int		filterend,		//
							double	pdcnt,
							double	priorh,
							char *	prior_file,
							int		genome_num,
							int	*	genome_type,
							int	*	genome_direction,
							int	*	genome_length, 
							char **	genome_name,
							char **	genome_sequence, 
							char **	rc_genome_sequence, 
							double ** genome_baseFreq,
							int *	genome_extleft,
							int *	genome_extright,
							float	filtertime,
							int		runmode,
							int		searchtype,
							int		iFlagStreamLine)
{
	if(searchtype == SEARCH_HMMFILTER) {
		this->hmm_search(top_k,
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
						iFlagPrintScoreInfo,		//
						iFlagPrintDebugInfo,
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
						pdcnt,
						priorh,
						prior_file,
						genome_num,
//						genome_type,
						genome_length,
						genome_name,
						genome_sequence,
						rc_genome_sequence,
						genome_baseFreq,
//						genome_extleft,
//						genome_extright,
//						filtertime,
						runmode,
						iFlagPrintDebugFilterInfo,	//
						iFlagStreamLine);			//
	} else {
		int pre_hit_num = 0;
		float time = 0.0;
		this->structure_search(top_k,
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
								pdcnt,
								priorh,
								prior_file,
								genome_num,
								genome_type,
								genome_direction,
								genome_length,
								genome_name,
								genome_sequence,
								rc_genome_sequence,
								genome_baseFreq,
								genome_extleft,
								genome_extright,
								filtertime,
								runmode,
								searchtype,
								0,	//streamline flag
								pre_hit_num,	//pre_hit_num
								time);			//search time
	}
}


