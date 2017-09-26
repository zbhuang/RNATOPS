// StructureSearch.h: interface for the StructureSearch class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_STRUCTURESEARCH_H__94292247_F0A1_479D_83C9_93D1EADD5AFA__INCLUDED_)
#define AFX_STRUCTURESEARCH_H__94292247_F0A1_479D_83C9_93D1EADD5AFA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "ScfgBuilder.h"
#include "HMMBuilder.h"
#include "common.h"
#include "FViterbi.h"
#include "DynamicP.h"
#include "GenomeData.h"
#include "FilterSelection.h"

#define RUNNINGMODE_WEBSERVER	0
#define RUNNINGMODE_STANDALONE	1

#define SEARCH_HMMFILTER			0
#define	SEARCH_SUBSTRUCTUREFILTER	1
#define	SEARCH_WHOLESTRUCTURE		2
#define	SEARCH_STREAMLINE			3

class StructureSearch  
{
private:
	void hmm_search(int		top_k,
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
					int	 *	genome_length,
					char **	genome_name,
					char **	genome_sequence,
					char **	rc_genome_sequence,
					double **genome_baseFreq,
					int		runmode,
					int		iFlagPrintDebugFilterInfo,	//
					int		iFlagStreamLine);

	void hmm_search_oneway(	
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
						double	genome_baseFreq_A,
						double	genome_baseFreq_C,
						double	genome_baseFreq_G,
						double	genome_baseFreq_U,
						int		runmode,
						int		iFlagPrintDebugFilterInfo,	//
						int		iFlagStreamLine,
						float & hmm_search_time,
						int	 &	count);

	void structure_search_based_on_one_hmmfilterresult(	int		top_k,
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
//														int		genome_length,
														char *	genome_name,
														char *	genome_segment,
														int		optimal_pos_start,
														int		extension_left,
														int		optimal_pos_end,
														int		extension_right,
														int		runmode,
														int	&	total_hit,
														float & time);

	void structure_search(	int		top_k,
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
							int		streamline_flag,
							int	&	pre_hit_num,			//
							float & struct_search_time);	//

public:
	StructureSearch();
	virtual ~StructureSearch();

	void search(	int		top_k,
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
					int		searchmode,
					int		searchtype,
					int		iFlagStreamLine);

};

#endif // !defined(AFX_STRUCTURESEARCH_H__94292247_F0A1_479D_83C9_93D1EADD5AFA__INCLUDED_)
