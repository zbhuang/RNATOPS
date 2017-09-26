#include "common.h"

//-------------------------------
int isRNABase(char nucleotide)
{
    int result;
    nucleotide = toupper( nucleotide);
    switch( nucleotide )
    {
        //case GAP:
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'U':
               result=1;
               break;
        default:
            	 //printf("There is invalid residue: %c\n", nucleotide);
            	 result=0;
    }
    return result;
}

int isRNABase(int nucleotide)
{
    int result;
    nucleotide = toupper( nucleotide);
    switch( nucleotide )
    {
        //case GAP:
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
               result=1;
               break;
        default:
            	 //printf("There is invalid residue: %c\n", nucleotide);
            	 result=0;
    }
    return result;
}

int chartonum(char nucleotide)
{
    int result;
    nucleotide = toupper( nucleotide);
    switch( nucleotide )
    {
        case GAP:  result=0;
                   break;
        case 'A':  result=1;
                   break;
        case 'C':  result=2;
                   break;
        case 'G':  result=3;
                   break;
        case 'T':
        case 'U':  result=4;
                   break;
        default:   
            	 printf("There is invalid residue: %c, which is treated as gap.\n", nucleotide);
            	 result=0;
    }
    return result;
}

char numtochar(int base)
{
    char result;
    switch( base )
    {
        case 0:  result=GAP;
                   break;
        case 1:  result='A';
                   break;
        case 2:  result='C';
                   break;
        case 3:  result='G';
                   break;
        case 4:  result='U';
                   break;
        default: printf("There is invalid residue: %d, which is treated as gap\n", base);
                 result=' ';
    }
    return result;
}

//-------------------------------
bool isCanonicalPair(char basex, char basey)
{
    char bex=' ';
    char bey=' ';
    bool result=false;
    
    bex = toupper( basex);
    bey = toupper( basey); 
    
    if( bex=='U' && bey=='A' )
        result=true;
    else if( bex=='A' && bey=='U' )
        result=true;
    else if( bex=='G' && bey=='C' )
        result=true;
    else if( bex=='C' && bey=='G' )
        result=true;
    else if( bex=='U' && bey=='G' )
        result=true;
    else if( bex=='G' && bey=='U' )
        result=true;
    else
        result=false;
    return result;
}

int isFileName( char *str)
{
    int i=0, len=0;
    int flag=0;
    char base=' ';
    
    len=strlen( str );
    for( i=0; i<len; i++ )
    {
        base=toupper( str[i] );
        if( base=='\r' || base=='\n')
            break;
        if( base!='A' && base!='C' && base!='G' && base!='U' && base!='T' && base!='-' )
        {
            flag=1;
            break;
        }
    }
    return flag;
}

void printPairEmP( double pairprob[5][5] )
{
    int i=0, j=0;
    char base[5]={'-', 'A','C', 'G', 'U'};
    printf("  ");
    for(i=0; i<5; i++ )
    {
        printf("%9c", base[i] );
    }
    printf("\n");
    for(i=0; i<5; i++ )
    {
        printf("%5c", base[i]);
         for(j=0; j<5; j++)
        {
            printf("%9.5f", pairprob[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printPairEmP( int pairprob[5][5] )
{
    int i=0, j=0;
    char base[5]={'-', 'A','C', 'G', 'U'};
    printf("     ");
    for(i=0; i<5; i++ )
    {
        printf("%9c", base[i] );
    }
    printf("\n");
    for(i=0; i<5; i++ )
    {
        printf("%5c", base[i]);
         for(j=0; j<5; j++)
        {
            printf("%9d", pairprob[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printPairHz( int pairprob[5][5], double prob[5][5])
{
    int i=0, j=0;
    char base[5]={'-', 'A','C', 'G', 'U'};
    printf("# of Pairs:                   ");
    printf("        Percentage: \n");
    printf("     ");
    for(i=0; i<5; i++ )
        printf("%5c", base[i] );
    printf("   ");
    for(i=0; i<5; i++ )
        printf("%8c", base[i] );
    printf("\n");
    
    for(i=0; i<5; i++ )
    {
        printf("%5c", base[i]);
        for(j=0; j<5; j++)
            printf("%5d", pairprob[i][j]);
        printf("      ");
        for(j=0; j<5; j++)
            printf("%8.4f", prob[i][j]);
        printf("\n");
    }
    printf("\n");    
}

void printBaseEmP( double baseprob[4] )
{
    int i=0;
    char base[4]={'A','C', 'G', 'U'};
    printf("  ");
    for(i=0; i<4; i++ )
    {
        printf("%9c", base[i] );
    }
    printf("\n");
    printf("%5c", ' ');
    for(i=0; i<4; i++ )
    {
        printf("%9.5f", baseprob[i]);
    }
    printf("\n");
    printf("\n");
}

void printSingleEmProb( double baseprob[4] )
{
    int i=0;
    printf("                    <");
    for(i=0; i<4; i++ )
        printf("%9.5f", baseprob[i]);
    printf(" >\n");
}

//--------------------------------------------------

int obtainSequence(  char * filename, char* &bufseq )
{
    #ifdef DEBUG_DEV
        printf("obtainSequence\n");
    #endif
    int i=0, templen=0;
    char ch;
    FILE * fp=NULL;

    //----Open the search file
    fp = fopen(filename, "r");
    if( fp == NULL )
    {
        printf(" Fail to open %s\n", filename);
        return 0;
    }

    //get the size of str: len
    fseek( fp, 0, SEEK_SET);
    while ( !feof(fp) )
    {
        ch=fgetc (fp);
        if( isRNABase(ch) )
            templen++;
    }
    //allocate the memory
    bufseq= new char[templen+1];
    bufseq[templen] = '\0';
    //get the str
    fseek( fp, 0, SEEK_SET);
    while ( !feof(fp) )
    {
        ch=fgetc (fp);
        if( isRNABase(ch) )
            bufseq[i++]=ch;
    }
    if( i!=templen )
    {
        printf("There's something wrong with target file.\n");
        return 0;
    }
    fclose(fp);
    return templen;
}

//(1-h)*orgm + h*priorm
void normalizeMatrix( double orgm[5][5], double priorm[5][5], double h)
{
    int i=0, j=0;
    double sumvalue=0.0;
    double ht=0.0, hp=0.0;
    if( h>1 )
        hp=1;
    else if( h<0)
        hp=0;
    else
        hp=h;
    
    ht=1-hp;    
    //print the original matrix
    //printf("Orignal matrix:\n");
    //printPairEmP( orgm );

    //printf("Prior matrix:\n");
    //printPairEmP( priorm );
    
    for(i=0; i<5; i++ )
    {
        for(j=0; j<5; j++ )
        {
            orgm[i][j]=ht*orgm[i][j] + hp*priorm[i][j];
            sumvalue += orgm[i][j];
        }
    }
    
    for(i=0; i<5; i++ )
    {
        for(j=0; j<5; j++ )
        {
            orgm[i][j]=orgm[i][j]/sumvalue;
        }
    }
    //print the new matrix
    //printf("New matrix:\n");
    //printPairEmP( orgm );
    //printf("\n\n");
    
}
char * strrpl(char *s, const char s1, const char s2)
{
	int m;
	int iLength = strlen(s);
	char *newChar = new char[iLength+1];
	for(m=0; m<iLength; m++)
	{
		if(s[m] == s1)
			newChar[m] = s2;
		else
			newChar[m] = s[m];
	}
	newChar[m] = '\0';
	return newChar;
}
/*
char * strdel(char *s, const char s1)
{
	int m;
	int iLength = strlen(s);
	int count = 0;
	char *newChar = new char[iLength+1];
	for(m=0; m<iLength; m++)
	{
		if(s[m] != s1)
			newChar[count++] = s[m];
	}
	newChar[count] = '\0';
	return newChar;
}
*/
char * strappend(char *s, const char s1)
{
    int m;
	int length = strlen(s);
	char * new_sequence = new char[length + 2];
	for(m=0; m<length; m++)
	{
		new_sequence[m] = s[m];
	}
	new_sequence[m++] = s1;
	new_sequence[m] = '\0';
	return new_sequence;
}
char * strextract(char *s, int start_pos, int end_pos)
{
    int m;
//	int length = strlen(s);
	char * new_sequence = new char[end_pos - start_pos + 2];
//	for(m=start_pos; m<end_pos; m++) {
	for(m=start_pos; m<=end_pos; m++) {
		new_sequence[m-start_pos] = s[m];
	}
	new_sequence[m-start_pos] = '\0';
	return new_sequence;
}
string extractfilename(string s)
{
	int b_onlyfilename = 1;
	int ipos = s.find_first_of("/");
	if(ipos != -1)
		b_onlyfilename = 0;
	if(b_onlyfilename == 1)
		return s;
	else {
		ipos = s.find_last_of("/");
		return s.substr(ipos+1, s.length());
	}
}
string printYN(int iFlag)
{
	if(iFlag)
		return string("Yes");
	else
		return string("No");
}
string printDirection(int iFlag)
{
	if(iFlag == GENOME_SEARCH_DIRECTION_BOTH)
		return string("Yes");
	else
		return string("No");
}
/*
void printGenomeSegments(char * c_seq, int num)
{
	string s_seq(c_seq);
	int length = s_seq.length();
	for(int i=0; i*num<length; i++)
	{
		cout<<s_seq.substr(i*num, num)<<endl;
	}
}
*/
void printGenomeSegments(char * uSeq, int num)
{
	char * tSeq = strrpl(uSeq, 'U', 'T');
	string s_seq(tSeq);
	int length = s_seq.length();
	for(int i=0; i*num<length; i++) {
		cout<<s_seq.substr(i*num, num)<<endl;
	}
	if(tSeq != NULL) {
		delete [] tSeq;
		tSeq = NULL;
	}
}
//------------------------------------------------------------------
int isHmmFilterFile(char * filename)
{
	int iResult = 1;
	FILE *  pFile = fopen(filename, "r");
    if( pFile == NULL ) {
        printf(" Fail to open %s\n", filename);
        exit( 0 );
    }
	fseek(pFile, 0, SEEK_SET);

	const int line_max_num = 200;
	char c_line[line_max_num-1];
	fgets(c_line, line_max_num, pFile);
//	string s_line(c_line);
//	for(int i=0; i<s_line.length() && (iResult==0); i++)
	for(unsigned int i=0; i<strlen(c_line) && (iResult==1); i++)
	{
		if(c_line[i] != '.' && c_line[i] != '\n')
			iResult = 0;
	}
	fclose(pFile);

	return iResult;
}
//deprecated
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
							  int		iFlagPrintDebugDPSearchInfo)
{
	cout<<"-------------Parameter settings-------------"<<endl;
	cout<<"pseudocount="<<pseudocount<<endl;
	cout<<"k="<<top_k<<endl;
	cout<<"Threshold="<<threshold<<endl;
	cout<<"Number of overlap between stems="<<num_nts_overlap_between_stem<<endl;
	cout<<"Take the merge-candidate strategy in preprocessing="<<printYN(iFlagMergeCandInPreprocess)<<endl;
	cout<<"About the merge strategy, take the candidate with the (s)hortest length="<<printYN(iFlagCandwithShortestLength)<<endl;
	cout<<"iShiftNumMergeCand="<<printYN(iShiftNumMergeCand)<<endl;
	cout<<"iAllowedNullLoopInsNum="<<iAllowedNullLoopInsNum<<endl;
	cout<<"pcoeff="<<pcoeff<<endl;
	cout<<"iJumpStrategy="<<printYN(iJumpStrategy)<<endl;
	cout<<"iStepSize="<<iStepSize<<endl;
	cout<<"iScoreThresholdInJump="<<dScoreThresholdInJump<<endl;
	cout<<"Print structure alignment option="<<printYN(iFlagPrintStrAlignInfo)<<endl;
	cout<<"reverse complement search="<<printYN(iFlagSearchReverse)<<endl;
//	cout<<"Debug option="<<printYN(iFlagPrintDebugInfo)<<endl;
//	cout<<"Print the debug info of input-checking in dp search="<<printYN(iFlagPrintDebugInputInfo)<<endl;
//	cout<<"Print the debug info of tree in dp search="<<printYN(iFlagPrintDebugTreeInfo)<<endl;
//	cout<<"Print the debug info of window size in dp search="<<printYN(iFlagPrintDebugDPWinInfo)<<endl;
//	cout<<"Print the debug info of td in dp search="<<printYN(iFlagPrintDebugTDInfo)<<endl;
//	cout<<"Print the debug info of the dp search="<<printYN(iFlagPrintDebugDPSearchInfo)<<endl;
	cout<<"--------------------------------------------"<<endl;
}
*/
void printCurrentTime()
{
	time_t     now;
	struct tm  *ts;
	char       buf[80];
 
    // Get the current time
	now = time(NULL);
 
    // Format and print the time, "ddd yyyy-mm-dd hh:mm:ss zzz"
	ts = localtime(&now);
//	strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", ts);
	strftime(buf, sizeof(buf), "%H:%M:%S %Z %Y-%m-%d", ts);
    cout<<left<<setw(40)<<buf;	//<<endl;
}

void printSearchEndingTimeInfo()
{
	cout<<"* Time: ";
	printCurrentTime();
	cout<<right<<setw(NUM_NTS_PER_LINE-48)<<"*"<<endl;
}

void printOneStarLine()
{
	cout<<"************************************************************"<<endl;
}
void printAdditionalInfo()
{
	printOneStarLine();
	cout<<"* Searched done : with RNATOPS V1.2                        *"<<endl
		<<"* By RNA-Informatics @ UGA                                 *"<<endl;
}

void printFilterSearchResultHeader(string	trainfile,
								   int		profilelength,
								   int		filterauto, 
								   int		filtertype, 
								   int		filterbegin, int filterend, 
								   string	genomefile,	//char * genomefile,
								   int		genomenum,
								   long		total_genome_length)
{
	printOneStarLine();
	cout<<"* Filtering Result                                         *"<<endl;
	cout<<"*                                                          *"<<endl;
	cout<<"* Profile file : "<<left<<setw(NUM_NTS_PER_LINE - 18)<<trainfile<<"*"<<endl;
	cout<<"* Profile length : "<<left<<setw(NUM_NTS_PER_LINE - 20)<<profilelength<<"*"<<endl;
	cout<<"*                                                          *"<<endl;
	if(filterauto)
		cout<<"* Filter generation: automatic seleted                     *"<<endl;
	else
		cout<<"* Filter generation: manual seleted                        *"<<endl;
	if(filtertype == SEARCH_HMMFILTER)
		cout<<"* Filter type : HMM                                        *"<<endl;
	else
		cout<<"* Filter type : substructure                               *"<<endl;
	cout<<"* Filter info : positions from "<<left<<setw(6)<<filterbegin<<" to "<<left<<setw(NUM_NTS_PER_LINE-42)<<filterend<<"*"<<endl;
//	cout<<"* : regions from "<<left<<setw(6)<<filterbegin<<" to "<<left<<setw(NUM_NTS_PER_LINE-28)<<filterend<<"*"<<endl;
	cout<<"* Genome file : "<<left<<setw(NUM_NTS_PER_LINE-17)<<genomefile<<"*"<<endl;
	cout<<"* Number of sequences : "<<left<<setw(NUM_NTS_PER_LINE-25)<<genomenum<<"*"<<endl;
	cout<<"* Total length of sequences : "<<left<<setw(NUM_NTS_PER_LINE-31)<<total_genome_length<<"*"<<endl;
	printOneStarLine();
}

void printFilterHitIndex(int index)
{
	cout<<endl<<"Filtering hit "<<index<<endl
		<<"---------------"<<endl;
}
void printHitPos(int start, int end)
{
	cout<<"Hit Positions: "<<start<<"-"<<end<<endl;
}
void printFilterHitAlignment(char *uSeq, char *align)
{
	char * tSeq = strrpl(uSeq, 'U', 'T');
	cout<<"Alignment to the filter"<<endl
		<<tSeq<<endl
		<<align<<endl;
	if(tSeq != NULL) {
		delete [] tSeq;
		tSeq = NULL;
	}
}
void printHitScore(double score)
{
	cout<<"Alignment score = "<<score<<endl;
}
void printFilterHitExtenPos(int start, int end)
{
	cout<<"Extension positions: "<<start<<"-"<<end<<endl;
}
void printFilterHitExtenNtsHeader()
{
	cout<<GENOME_TAG_HIT_EXTENSION<<endl;
}

//void printFilterHitExtenNtsTail()
//{
//	cout<<GENOME_TAG_HIT_EXTENSION_TAIL<<endl;
//}

void printTotalHitNum(int hit_num)
{
	cout<<"* "<<GENOME_TAG_TOTAL_NUM_HIT<<": "<<left<<setw(15)<<hit_num<<right<<setw(NUM_NTS_PER_LINE-35)<<"*"<<endl;
}

void printTotalTimeInfo(double dTime)
{
	cout<<"* "<<GENOME_TAG_TOTAL_TIME_USED<<" "<<left<<setw(15)<<dTime<<" hours "<<right<<setw(NUM_NTS_PER_LINE-40)<<"*"<<endl;
}

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
								  int		filtertype, 
								  int		filterbegin, 
								  int		filterend)
{
	printOneStarLine();
	cout<<"* Whole Structure Result                                   *"<<endl;
	cout<<"*                                                          *"<<endl;
	cout<<"* Profile file : "<<left<<setw(NUM_NTS_PER_LINE - 18)<<trainfile<<"*"<<endl;
	cout<<"* Profile length : "<<left<<setw(NUM_NTS_PER_LINE - 20)<<profilelength<<"*"<<endl;
	cout<<"*                                                          *"<<endl;
	if(streamline)	//adding some filter info
	{
		//
		if(filtertype == SEARCH_HMMFILTER)
			cout<<"* Filter type : HMM                                        *"<<endl;
		else
			cout<<"* Filter type : substructure                               *"<<endl;
		cout<<"* Filter info : positions from "<<left<<setw(6)<<filterbegin<<" to "<<left<<setw(NUM_NTS_PER_LINE-42)<<filterend<<"*"<<endl;
		cout<<"*                                                          *"<<endl;
	}
	cout<<"* Genome file : "<<left<<setw(NUM_NTS_PER_LINE - 17)<<genomefile<<"*"<<endl;
	cout<<"* Number of sequences : "<<left<<setw(NUM_NTS_PER_LINE - 25)<<genomenum<<"*"<<endl;
	cout<<"* Total length of sequences : "<<left<<setw(NUM_NTS_PER_LINE - 31)<<total_genome_length<<"*"<<endl;
	cout<<"*                                                          *"<<endl;
	printValuesForParameters(pseudocount,
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
							iFlagSearchReverse);
	printOneStarLine();
}
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
							  int		iFlagSearchReverse)
{
	cout<<"* Search parameters setting:                               *"<<endl;
	cout<<"* Pseudocount = "<<left<<setw(NUM_NTS_PER_LINE - 17)<<pseudocount<<"*"<<endl;
	cout<<"* Num of stem candidates = "<<left<<setw(NUM_NTS_PER_LINE - 28)<<top_k<<"*"<<endl;
	cout<<"* Score threshold for hits = "<<left<<setw(NUM_NTS_PER_LINE - 30)<<threshold<<"*"<<endl;
	cout<<"* Num of nt overlap between stems = "<<left<<setw(NUM_NTS_PER_LINE - 37)<<num_nts_overlap_between_stem<<"*"<<endl;
	cout<<"* Candidate representatives only = "<<left<<setw(NUM_NTS_PER_LINE - 36)<<printYN(iFlagMergeCandInPreprocess)<<"*"<<endl;
	cout<<"* Shortest candidate representatives = "<<left<<setw(NUM_NTS_PER_LINE - 40)<<printYN(iFlagCandwithShortestLength)<<"*"<<endl;
	cout<<"* IShiftNumMergeCand = "<<left<<setw(NUM_NTS_PER_LINE - 24)<<printYN(iShiftNumMergeCand)<<"*"<<endl;
	cout<<"* Nts allowed in null loops = "<<left<<setw(NUM_NTS_PER_LINE - 31)<<iAllowedNullLoopInsNum<<"*"<<endl;
	cout<<"* Pcoeff = "<<left<<setw(NUM_NTS_PER_LINE - 12)<<pcoeff<<"*"<<endl;
	cout<<"* Search with jump = "<<left<<setw(NUM_NTS_PER_LINE - 22)<<printYN(iJumpStrategy)<<"*"<<endl;
	cout<<"* Search step size = "<<left<<setw(NUM_NTS_PER_LINE - 22)<<iStepSize<<"*"<<endl;
	cout<<"* Search reversed complement sequence = "<<left<<setw(NUM_NTS_PER_LINE - 41)<<printDirection(iFlagSearchReverse)<<"*"<<endl;
}

int getProfileLength(char * trainfile)
{
	int length = 0;
	FILE *  pFile = fopen(trainfile, "r");
    if( pFile == NULL ) {
        printf(" Fail to open %s\n", trainfile);
        exit( 0 );
    }
	fseek(pFile, 0, SEEK_SET);

	const int line_max_num = 1000;
	char c_line[line_max_num-1];
	fgets(c_line, line_max_num, pFile);
	fclose(pFile);

	//Yf Wang report this bug. 20081029
//	length = strlen(c_line);
	string str_line(c_line);
	length = str_line.find_first_of("\n");
	return length;
}

void printWholeSearchHitIndex(int index)
{
	cout<<endl<<"Whole structure search hit "<<index<<endl
		<<"----------------------------"<<endl;
}

void printEValue(double evalue)
{
	cout<<"E-value = "<<evalue<<endl;
}

void formatStringOutput(string illustration, string line_one, string line_two, string line_thr, int start, int end, int flag)
{
	//NUM_NTS_PER_LINE
	cout<<endl<<illustration<<endl;

	string line_one_unit;
	string line_two_unit;
	string line_thr_unit;
	int length = line_one.length();
	for(int i=0; i*NUM_NTS_PER_LINE<length; i++) {
		line_one_unit.assign(line_one.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		if(flag)
			line_two_unit.assign(line_two.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		line_thr_unit.assign(line_thr.substr(i*NUM_NTS_PER_LINE, NUM_NTS_PER_LINE));
		formatOneOutput(line_one_unit, line_two_unit, line_thr_unit, i, start, flag);
	}
}

void formatOneOutput(string line_one_unit, string line_two_unit, string line_thr_unit, int index, int start, int flag)
{
	int length_one = line_one_unit.length();
	int length_two = 0;
	if(flag)
		length_two = line_two_unit.length();
	int length_thr = line_thr_unit.length();
	if(flag) {
		cout<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_two_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_two)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	} else {
		cout<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	}
}

/*
 * To solve the error of "terminate called after throwing an instance of 'std::out_of_range'"
 * 20081031
 */
void formatStringOutput2(char * illustration, char * line_one, char * line_two, char * line_thr, int start, int end, int flag)
{
	int i, j;
	//NUM_NTS_PER_LINE
	cout<<endl<<illustration<<endl;

	char * line_one_unit = NULL;
	char * line_two_unit = NULL;
	char * line_thr_unit = NULL;
	int length = strlen(line_one);
	for(i=0; i*NUM_NTS_PER_LINE<length; i++) {
		if((i+1)*NUM_NTS_PER_LINE > length) {
			line_one_unit = new char[length - i*NUM_NTS_PER_LINE + 1];
			for(j=0; j<(length - i*NUM_NTS_PER_LINE); j++)
				line_one_unit[j] = line_one[i*NUM_NTS_PER_LINE+j];
			line_one_unit[j] = '\0';

			line_thr_unit = new char[length - i*NUM_NTS_PER_LINE + 1];
			for(j=0; j<(length - i*NUM_NTS_PER_LINE); j++)
				line_thr_unit[j] = line_thr[i*NUM_NTS_PER_LINE+j];
			line_thr_unit[j] = '\0';
		} else {
			line_one_unit = new char[NUM_NTS_PER_LINE + 1];
			for(j=0; j<NUM_NTS_PER_LINE; j++)
				line_one_unit[j] = line_one[i*NUM_NTS_PER_LINE+j];
			line_one_unit[j] = '\0';

			line_thr_unit = new char[NUM_NTS_PER_LINE + 1];
			for(j=0; j<NUM_NTS_PER_LINE; j++)
				line_thr_unit[j] = line_thr[i*NUM_NTS_PER_LINE+j];
			line_thr_unit[j] = '\0';
		}
		if(flag) {
			if((i+1)*NUM_NTS_PER_LINE > length) {
				line_two_unit = new char[length - i*NUM_NTS_PER_LINE + 1];
				for(j=0; j<(length - i*NUM_NTS_PER_LINE); j++)
					line_two_unit[j] = line_two[i*NUM_NTS_PER_LINE+j];
				line_two_unit[j] = '\0';
			} else {
				line_two_unit = new char[NUM_NTS_PER_LINE + 1];
				for(j=0; j<NUM_NTS_PER_LINE; j++)
					line_two_unit[j] = line_two[i*NUM_NTS_PER_LINE+j];
				line_two_unit[j] = '\0';
			}
		}
		formatOneOutput2(line_one_unit, line_two_unit, line_thr_unit, i, start, flag);
		if(line_one_unit != NULL) {
			delete [] line_one_unit;
			line_one_unit = NULL;
		}
		if(line_two_unit != NULL) {
			delete [] line_two_unit;
			line_two_unit = NULL;
		}
		if(line_thr_unit != NULL) {
			delete [] line_thr_unit;
			line_thr_unit = NULL;
		}
	}
}

void formatOneOutput2(char * line_one_unit, char * line_two_unit, char * line_thr_unit, int index, int start, int flag)
{
	int length_one = strlen(line_one_unit);
	int length_two = 0;
	if(flag)
		length_two = strlen(line_two_unit);
	int length_thr = strlen(line_thr_unit);
	if(flag) {
		cout<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_two_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_two)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	} else {
		cout<<endl
			<<right<<setw(10)<<(index*NUM_NTS_PER_LINE+1)<<" "<<line_one_unit<<" "<<(index*NUM_NTS_PER_LINE+1+length_one)<<endl
			<<right<<setw(10)<<(start+1+index*NUM_NTS_PER_LINE)<<" "<<line_thr_unit<<" "<<(start+1+index*NUM_NTS_PER_LINE+length_thr)<<endl;
	}
}
