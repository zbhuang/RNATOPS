#include "Viterbi.h"

//********************************************************************************
//*******************    implementation of VTBCandidate class         ***************
//********************************************************************************
VTBCandidate::VTBCandidate( )
{
    #ifdef DEBUG_DEV
        printf("VTBCandidate::VTBCandidate\n");
    #endif
    initAllMembers( );
}

VTBCandidate::~VTBCandidate( )
{
    #ifdef DEBUG_DEV
        printf("VTBCandidate::~VTBCandidate\n");
    #endif
    freeAllocMemory( );
}

void VTBCandidate::initAllMembers( )
{
    #ifdef DEBUG_DEV
        printf("VTBCandidate::initAllMembers\n");
    #endif
    int i=0;
    
    pos[0]=-1;
    pos[1]=-1;
    pos[2]=-1;
    for(i=0; i<2; i++ )
        loc[i]=INVALID;
    
     prob=INVLDPROB;     //probability
     cadno=-1;      //-1, only one;
     finseq=NULL;   //final sequence with gap
     finstr=NULL;   //final structure
     finlen=0;
}

void VTBCandidate::initVTBCandidates( )
{
    #ifdef DEBUG_DEV
        printf("VTBCandidate::initVTBCandidates\n");
    #endif
    freeAllocMemory( );
    initAllMembers( );
}

void VTBCandidate::freeAllocMemory( )
{
    #ifdef DEBUG_DEV
        printf("VTBCandidate::freeAllocMemory\n");
    #endif
    if( finseq!=NULL )
    {
        delete [] finseq;
        finseq=NULL;
    }
    if( finstr!=NULL )
    {
        delete [] finstr;
        finstr=NULL;
    }
}

int VTBCandidate::hasData( )
{
    #ifdef DEBUG_DEV
        printf("VTBCandidate::hasData\n");
    #endif
    if( loc[0]==INVALID )
        return 0;
    else
        return 1;
}

double VTBCandidate::getLog(  )
{
    return prob;
}

void VTBCandidate::printAlgnedSeq( )
{
    #ifdef DEBUG_DEV
        printf("VTBCandidate::printAlgnedSeq\n");
    #endif
    if( cadno != -1 )
    {
        printf("No.%d:  ", cadno);
    }
    printf("Max Log-Prob: %10.6f at [%d, %d]\n", prob, loc[0], loc[1]);
    //printf("Result:\n");
    printf("   ");
    printf("%s\n", finseq );
    printf("   ");
    printf("%s\n\n", finstr );
}


//********************************************************************************
//*******************    implementation of VTBSearch class         ***************
//********************************************************************************
VTBSearch::VTBSearch(  )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBSearch\n");
    #endif
    searchSeq=NULL;
    seqlen=0;
    numstates=0;
  	Vt=NULL;
  	Tb=NULL;
    numtops=1;
    allowNumIns=0;
    tops=new VTBCandidate[ numtops ];
}

VTBSearch::VTBSearch( int nummax )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBSearch( parameter )\n");
    #endif
    searchSeq=NULL;
    seqlen=0;
    numstates=0;
  	Vt=NULL;
  	Tb=NULL;
    numtops=nummax;
    allowNumIns=0;
    tops=new VTBCandidate[ numtops ];
}

VTBSearch::~VTBSearch(  )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::~VTBSearch\n");
    #endif
    freeSearchSeq( );
    free3DTables( );
    if(tops!=NULL)
    {
       delete [ ] tops;
       tops=NULL;
    }
}

void VTBSearch::initTopCands( )
{
    int i=0;
    for(i=0; i<numtops; i++ )
    {
        tops[i].initVTBCandidates( );
    }
}

int VTBSearch::malloc3DTables(  )
{
    #ifdef DEBUG_DEV
       printf("VTBSearch::malloc3DTables\n");
    #endif
    
    int i=0, j=0;
    
    Vt = new double** [ seqlen+2  ];
    if( Vt==NULL )
        return 0;
    
    Tb = new char** [ seqlen+2  ];
    if( Tb==NULL )
        return 0;
    
    for( i=0; i<seqlen+2; i++ )
    {
        Vt[i]=new double* [ numstates ];
        if( Vt[i]==NULL )
            return 0;
        Tb[i]=new char* [ numstates ];
        if( Tb[i]==NULL )
            return 0;
    }
    
    for( i=0; i<seqlen+2; i++ )
    {
        for( j=0; j<numstates; j++ )
        {
            Vt[i][j]=new double [ 3 ];
            if( Vt[i][j]==NULL )
                return 0;
            Tb[i][j]=new char [ 3 ];
            if( Tb[i][j]==NULL )
                return 0;
        }
    }
    return 1;    
}

void VTBSearch::freeSearchSeq(  )
{
    #ifdef DEBUG_DEV
       printf("VTBSearch::freeSearchSeq\n");
    #endif
    if( searchSeq!=NULL )
    {
        delete [ ] searchSeq;
        searchSeq=NULL;
    }
}

void VTBSearch::free3DTables(  )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::free3DTables\n");
    #endif
    int i=0, j=0;
    //for Vt
    if( Vt!=NULL )
    {
        for( i=0; i<seqlen+2; i++ )
        {
            if( Vt[i]!=NULL )
            {
                for( j=0; j<numstates; j++ )
                {
                    if( Vt[i][j]!=NULL )
                    {
                        delete [ ] Vt[i][j];
                        Vt[i][j]=NULL;
                    }
                }
                delete [ ] Vt[i];
                Vt[i]=NULL;
            }
        }
        delete [ ] Vt;
        Vt=NULL;
    }

    //for Tb
    if( Tb!=NULL )
    {
        for( i=0; i<seqlen+2; i++ )
        {
            if( Tb[i]!=NULL )
            {
                for( j=0; j<numstates; j++ )
                {
                    if( Tb[i][j]!=NULL )
                    {
                        delete [ ] Tb[i][j];
                        Tb[i][j]=NULL;
                    }
                }
                delete [ ] Tb[i];
                Tb[i]=NULL;
            }
        }
        delete [ ] Tb;
        Tb=NULL;
    }
}
    
void VTBSearch::init3DTables(  )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::init3DTables\n");
    #endif
    
    int i=0, j=0, k=0;
    for(i=0; i<seqlen+2; i++)
    {
        for(j=0; j<numstates; j++)
        {
            for(k=0; k<3; k++)
            {
                Vt[i][j][k]=INVLDPROB;
                Tb[i][j][k]=' ';
            }
        }
    }
}

int VTBSearch::obtainSequence(  char * filename )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::obtainSequence\n");
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
    searchSeq= new int[templen];
    //get the str
    fseek( fp, 0, SEEK_SET);
    while ( !feof(fp) )
    {
        ch=fgetc (fp);
        if( isRNABase(ch) )
            searchSeq[i++]=chartonum(ch);
    }
    if( i!=templen )
    {
        printf("Wrong with length in target file.\n");
        return 0;
    }
    fclose(fp);
    return templen;
}

int VTBSearch::VTBSearchFile( char *filename, Loop *curloop )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBSearchFile\n");
    #endif
    
    if( curloop->getLoopLength( )<=0 )
        return 0;
    if( !curloop->isLogProb( ) )
        curloop->logProbforLoop( );
    
    int status=0;

    if( searchSeq!=NULL )
    {
        delete [ ] searchSeq;
        searchSeq=NULL;
    }
    //freeVTBMemory( );
    
    numstates = curloop->numstates;
    seqlen = obtainSequence( filename );
    if( seqlen == 0)
        return 0;    
    printSearchSeq( );

    //allocate P and G memory
    status = malloc3DTables( );
    if( status==0)
    {
        printf("Fail to allocate memory in malloc3DTables\n");
        return 0;
    }
    VTBimplement( curloop );
    free3DTables( );
    //printResultSeq( );
    return 1;
}

int	VTBSearch::update3DandSeq( int newlen,  int newnumstate )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::update3DandSeq\n");
    #endif
    int change=0;
    
    if( newlen!= seqlen || newnumstate!= numstates )
    {
        free3DTables( );
        change = 1;
    }

    if( newlen!= seqlen )
    {
        seqlen = newlen;
        if( searchSeq!=NULL )
        {
            delete [ ] searchSeq;
            searchSeq=NULL;
        }
        searchSeq=new int[seqlen];
    }
    
    if( newnumstate!= numstates )
    {
        numstates = newnumstate;
    }
    
    if( change )
    {
        if( malloc3DTables( )==0 )
        {
            printf("Fail to allocate memory in malloc3DTables\n");
            return 0;
        }
    }
    return 1;
}

int VTBSearch::VTBSearchString( char *instr, int newlen, Loop *curloop )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBSearchString\n");
    #endif
    
    int i=0, retval=0;
    int *tempstr = NULL;
    
    if( newlen < 0 )
        return 0;
    
    //copy searchSeq
    tempstr = new int[newlen];
    for( i=0; i<newlen; i++ )
         tempstr[i]=chartonum( instr[i] );
  
    retval = VTBSearchString( tempstr, newlen, curloop );
    
    delete [ ] tempstr;
    return retval;
}

int VTBSearch::VTBSearchString( int *instr, int newlen, Loop *curloop )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBSearchString\n");
    #endif
    
    if( !curloop->isLogProb( ) )
        curloop->logProbforLoop( );
        
    if( curloop->getLoopLength( )<=0 )
        return 0;    
    
    int i=0;
    int newnumstate=0;
    
    newnumstate = curloop->numstates;
    
    //1.check the dimension of 3D
    update3DandSeq( newlen, newnumstate );
    
    //2.copy searchSeq
    for( i=0; i<seqlen; i++ )
    {
         searchSeq[i]= instr[i];
    }
    
    //3.do VTB profile
    VTBimplement( curloop );
    
    //printResultSeq( );
    return 1;
}

double VTBSearch::VTBSearchMax( Loop *curloop, char *instr, int start, int &end, char **retpath, char **retSeq )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBSearchMax\n");
    #endif
    int i=0, curlen=0;
    int *subseq=NULL;
    double theprob=INVLDPROB;
    
    //process the passing parameters
    curlen=end-start+1;
    if( curlen < 0)
        curlen=0;
    subseq= new int[ curlen ];
    for( i=0; i<curlen; i++ )
        subseq[i]= chartonum( instr[start+i] );
    
    //Loop model is empty
    if( curloop->getLoopLength( )<=0 )
    {
        if( curlen <= allowNumIns )
            theprob=0;
        else
            theprob=INVLDPROB;
    }
    else  //Loop model is not empty
    {
        VTBSearchString( subseq, curlen, curloop );
        theprob = tops->getLog( );  // update the variable prob
        end=tops->loc[1];      //end position
        if( retpath != NULL )  //traceback path
        {
            *retpath=new char[ tops->finlen+1 ];
            strcpy( *retpath, tops->finstr );
        }
        if( retSeq != NULL )
        {
            *retSeq=new char[ tops->finlen+1 ];
            strcpy( *retSeq, tops->finseq );
        }
    }
    delete [ ] subseq;
    return theprob;
}

VTBCandidate * VTBSearch::getCandHead( )
{
    return  tops;
}

int VTBSearch::setAllowedInsNum( int num )
{
    allowNumIns = num;
    return num;
}

void VTBSearch::printSearchSeq( )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::printSearchSeq\n");
    #endif
    int i=0;
    printf("\n======================= Search Sequence =========================\n");
    printf("Target Sequences of length %d:\n", seqlen );
    printf("   ");
    for(i=0; i<seqlen; i++)
    {
        printf("%c", numtochar(searchSeq[i]));
    }
    printf("\n\n");
}

void VTBSearch::printResultSeq( )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::printResultSeq\n");
    #endif
    int i=0;
    for( i=0; i<numtops; i++ )
    {
        if( tops[i].hasData( ) )
        {
            tops[i].printAlgnedSeq( );
        }
    }
}

void VTBSearch::output3Dtables(  )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::output3Dtables\n");
    #endif
    int i=0, j=0;
    printf("Match table( x:seq=%d+1; y:states=%d )\n", seqlen, numstates );
    for( i=0; i<numstates; i++ )
    {
        if( i==0 )
            printf("B0:  ");
        else if(i==numstates-1 )
            printf("E%d:  ", i);
        else
            printf("M%d:  ", i);
        for( j=0; j<seqlen+2; j++ )
        {
            if( Vt[j][i][ HMM_M ] > INVLDPROB )
                printf(" %10.6f(%c) ", Vt[j][i][ HMM_M ], Tb[j][i][ HMM_M ]);
            else
                printf(" %10s(%c) ", "-INF", Tb[j][i][ HMM_M ]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("Insertion table( x:seq=%d+1; y:states=%d )\n", seqlen, numstates );
    for( i=0; i<numstates; i++ )
    {
        if( i==0 )
            printf("B0:  ");
        else if(i==numstates-1 )
            printf("E%d:  ", i);
        else
            printf("M%d:  ", i);
        for( j=0; j<seqlen+2; j++ )
        {
            if( Vt[j][i][ HMM_I ] > INVLDPROB )
                printf(" %10.6f(%c) ", Vt[j][i][ HMM_I ], Tb[j][i][ HMM_I ]);
            else
                printf(" %10s(%c) ", "-INF", Tb[j][i][ HMM_I ]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("Deletion table( x:seq=%d+1; y:states=%d )\n", seqlen, numstates );
    for( i=0; i<numstates; i++ )
    {
        if( i==0 )
            printf("B0:  ");
        else if(i==numstates-1 )
            printf("E%d:  ", i);
        else
            printf("M%d:  ", i);
            
        for( j=0; j<seqlen+2; j++ )
        {
            if( Vt[j][i][ HMM_D ] > INVLDPROB )
                printf(" %10.6f(%c) ", Vt[j][i][ HMM_D ], Tb[j][i][ HMM_D ]);
            else
                printf(" %10s(%c) ", "-INF", Tb[j][i][ HMM_D ]);
        }
        printf("\n");
    }
    printf("\n");
}


// Pick the maximum value from a and b
double VTBSearch::max(double fromI, double fromM )
{
    double temp;
	
    if (fromI > fromM )
        temp = fromI;
    else
        temp = fromM;
    return temp;
}

// Pick the maximum value from a, b and c
double VTBSearch::max(double fromD, double fromI, double fromM )
{
    double temp;
	
    if ( fromD > fromI )
        temp=fromD;
    else
        temp=fromI;
    
    if ( temp > fromM )
        return temp;
    else
        return fromM;
}

// Give the type of the max one
char VTBSearch::maxstate (double fromI, double fromM )
{
    char type;
    if (fromI > fromM )
        type = 'I';
    else
        type = 'M';
  	return type;
}

// Give the type of the max one
char VTBSearch::maxstate (double fromD, double fromI, double fromM )
{
    char type;
    double  temp=0.0;
    if ( fromD > fromI )
    {
        temp = fromD;
        type = 'D';
    }
    else
    {
        temp = fromI;
        type = 'I';
    }
    if ( temp > fromM )
        return type;
    else
        return 'M';
}

//insert Top N number of positions
double VTBSearch::insertCurPos( int idi, int idj, int idk,  double cprob )
{
    int i=0, j=0;
    for( i=0; i<numtops; i++ )
    {
        if( cprob > tops[i].prob )
        {
            for( j=numtops-1; j>i; j--)
            {
                tops[j].pos[0]=tops[j-1].pos[0];
                tops[j].pos[1]=tops[j-1].pos[1];
                tops[j].pos[2]=tops[j-1].pos[2];
                tops[j].prob=tops[j-1].prob;
            }
            tops[i].pos[0]=idi;
            tops[i].pos[1]=idj;
            tops[i].pos[2]=idk;
            tops[i].prob=cprob;
            break;
        }
    }
    return tops[ numtops-1 ].prob;
}

int VTBSearch::getGlobalProbs(  )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::getGlobalProbs\n");
    #endif
    
    int i=0, k=0, num=0;
    int floor=numstates-1;
    double thres=INVLDPROB, curprob=INVLDPROB;
    
    //[seq][state][ D=0/I=1/M=2 ]
    i=seqlen+1;  
    for(k=0; k<3; k++ )
    {
        curprob = Vt[i][floor][k];
        if( curprob >= thres && curprob!=INVLDPROB)
        {
            thres=insertCurPos( i, floor,  k,  curprob );
      	}
    }
    //count the final number of positions
    for(i=0; i<numtops; i++ )
    {
        if( tops[i].pos[0]!=-1)
            num++;
    }
    return num;
}

// get the n candidates with max probability
int VTBSearch::getTopnMaxProbs(  )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::getTopnMaxProbs\n");
    #endif
    
    int i=0, k=0, num=0;
    int floor=numstates-1;
    double thres=INVLDPROB, curprob=INVLDPROB;
    
    for(i=0; i<seqlen+2; i++)    //[seq][state][ D=0/I=1/M=2 ]
    {
        for(k=0; k<3; k++ )
        {
            curprob = Vt[i][floor][k];
            if( curprob >= thres && curprob!=INVLDPROB)
            {
                thres=insertCurPos( i, floor,  k,  curprob );
            }
      	}
    }
    //count the final number of positions
    for(i=0; i<numtops; i++ )
    {
        if( tops[i].pos[0]!=-1)
            num++;
    }
    return num;
}

int VTBSearch::isInRemSeq( int *seq, int size,  int num)
{
    int i=0, found=0;
    for( i=0; i<size; i++ )
    {
        if( num == seq[i] )
        {
            seq[i]=INVALID;
            found++;
        }
    }
    return found;
}

int VTBSearch::VTBTraceBack( VTBCandidate &cad )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBTraceBack\n");
    #endif

    int i=0, j=0, k=0;
    int alld=0, numdel=0;
    int newlen=0; 
    int templen=5*seqlen+numstates;
    int *rem=new int[ templen ];
    int size=0;
    char *align=new char[ templen+1 ];
    
    char curstate=' ';
    
    for( i=0; i<templen; i++)
    {
        align[i]='.';
        rem[i]=INVALID;
    }
    align[i]='\0';
    
    //***** [seq][state][ D=0/I=1/M=2 ]
    i=cad.pos[0];
    j=cad.pos[1];
    k=cad.pos[2];
    
    curstate = Tb[i][j][k];
    i--;
    j--;
    cad.loc[1]=i-1;
    cad.loc[0]=i-1;
    while( i!=0 || j!=0 )
    { 
        cad.loc[0]=i-1;
        //printf(" j=%d, i=%d, state=%c\n", j, i, curstate);
        switch( curstate )
        {
            case 'D':  //0
                      curstate = Tb[i][j][HMM_D];
                      rem[ size++ ]=i-1;
                      j--;
                      break;
            case 'I':  //1
                      curstate = Tb[i][j][HMM_I];
                      align[i-1] = AlgHmmIns;
                      i--;
                      break;
            case 'M':  //2
                      curstate = Tb[i][j][HMM_M];
                      align[i-1] = AlgHmmMatch;
                      i--;
                      j--;
                      break;
            case ' ': 
                      break;
            default:
                      printf(" Unknown nonterminal:\"%c\" id: %d!\n", curstate, k);
                      return 0;
        }
    }
    //***********process the output string
    newlen=seqlen + size;
    cad.finlen=newlen;
    cad.finstr= new char[newlen+1];
    cad.finseq= new char[newlen+1];
    //for( i=0; i<size; i++ )
    //    printf("%d  ", rem[i] );
    //printf("\n" );
    for( i=0, k=0; i<newlen; i++)
    {
        numdel=isInRemSeq( rem, size,  i-1-alld);
        if( numdel )
        {
            //printf("k=%d, numdel=%d\n", k, numdel);
            for(j=0; j<numdel; j++ )
            {
                cad.finseq[i+j]='-';
                cad.finstr[i+j]=AlgHmmDel;
            }
            i+=numdel-1;
            alld+=numdel;
        }
        else
        {
            cad.finseq[i]= numtochar(searchSeq[k]);
            cad.finstr[i]= align[k];
            k++;
        }
    }
    cad.finseq[i]='\0';
    cad.finstr[i]='\0';
    delete [ ] rem;
    delete [ ] align;
    return 1;
}

void VTBSearch::VTBimplement( Loop *curLoop )
{
    #ifdef DEBUG_DEV
        printf("VTBSearch::VTBimplement\n");
    #endif
    
    int i=0;
    int nt=0, st=0;
    int residue=-1;
    double	fromM=0, fromI=0, fromD=0;
    int	nowtops=0;
    
    HMMGram *preM = NULL,	*preI = NULL,	*preD = NULL;
    HMMGram *curM = NULL,	*curI = NULL,	*curD = NULL;    
    HMMGram *beg = curLoop->pLoopGram;
    
    //****** initialization stage ****** [ seq ][ state ][ D=0/I=1/M=2 ]
    initTopCands( );
    init3DTables( );
    Vt[0][0][ HMM_M ]=0;
    
    // initialize the first column in D
    Vt[0][1][ HMM_D ] = Vt[0][0][ HMM_M ] + beg->trprob[ HMM_D ];
    Tb[0][1][ HMM_D ] = 'M';
    curD=beg->nextD;
    for (i=2; i<numstates; i++)//  2-----numstates-1
    {
       	Vt[0][i][ HMM_D ] = Vt[0][i-1][ HMM_D ] + curD->trprob[ HMM_D ];
       	Tb[0][i][ HMM_D ] = 'D';
        preD = curD;
       	curD=curD->nextD;
    }
    
    if( seqlen == 0)
    {
        //M  seqlen+1,  numstates
        if( preD!=NULL )
        {
            fromD = Vt[seqlen][numstates-2][ HMM_D ] + preD->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = fromD;
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = 'D';
        }
        else
        {
            fromM = Vt[seqlen][numstates-2][ HMM_M ] + beg->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = fromM;
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = 'M';
        }
    }
    else   //seqlen>0
    {
        // initialize the first row in I
        curI=beg->nextI;
        Vt[1][0][ HMM_I ] = Vt[0][0][HMM_M] + beg->trprob[HMM_I] + curI->emission[ searchSeq[0]-1 ];
        Tb[1][0][ HMM_I ] = 'M';
        for (i=2; i<seqlen+1; i++)
        {
           	Vt[i][0][ HMM_I ] = Vt[i-1][0][HMM_I] + curI->trprob[HMM_I] + curI->emission[ searchSeq[i-1]-1 ];
           	Tb[i][0][ HMM_I ] = 'I';
        }

        //****** main loop stage ****** [seq][state][ D=0/I=1/M=2 ]
        if( beg->nextM->state != 'E' )
        {
             for ( nt=1; nt<seqlen+1; nt++ )
            {
                preM = beg;
                preI = beg->nextI;
                curM = beg->nextM;
                curD = beg->nextD;
                curI = curM->nextI;
            
                residue = searchSeq[nt-1]-1;
                st=1;
                // M--match
                fromM = Vt[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
    	           fromI = Vt[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];             
    	           Vt[nt][st][ HMM_M ] = max( fromI, fromM ) + curM->emission[ residue ];
    	           Tb[nt][st][ HMM_M ] = maxstate( fromI, fromM );

                // I--insertion
                fromM = Vt[nt-1][st][ HMM_M ] + curM->trprob[ HMM_I ];
                fromI = Vt[nt-1][st][ HMM_I ] + curI->trprob[ HMM_I ];
                fromD = Vt[nt-1][st][ HMM_D ] + curD->trprob[ HMM_I ];
                Vt[nt][st][ HMM_I ] = max( fromD, fromI, fromM ) + curI->emission[ residue ];
                Tb[nt][st][ HMM_I ] = maxstate( fromD, fromI, fromM );
        
                // D--deletion
                fromM = Vt[nt][st-1][ HMM_M ] + preM->trprob[ HMM_D ];
                fromI = Vt[nt][st-1][ HMM_I ] + preI->trprob[ HMM_D ];
                Vt[nt][st][ HMM_D ] = max( fromD, fromI, fromM );
                Tb[nt][st][ HMM_D ] = maxstate( fromI, fromM );

    	           preM = curM;
    	           preD = curD;
    	           preI = curI;
    	           curM = preM->nextM;
    	           curD = preM->nextD;
                curI = curM->nextI;
      	
                for ( st=2; st<numstates-1;  st++ )
                {
                     // M--match
                    fromM = Vt[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
                     fromI = Vt[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];
                     fromD = Vt[nt-1][st-1][ HMM_D ] + preD->trprob[ HMM_M ];
                     Vt[nt][st][ HMM_M ] = max(fromD, fromI, fromM ) + curM->emission[ residue ];
                     Tb[nt][st][ HMM_M ] = maxstate(fromD, fromI, fromM);

                    // I--insertion
                    fromM = Vt[nt-1][st][ HMM_M ] + curM->trprob[ HMM_I ];
                    fromI = Vt[nt-1][st][ HMM_I ] + curI->trprob[ HMM_I ];
                    fromD = Vt[nt-1][st][ HMM_D ] + curD->trprob[ HMM_I ];
                    Vt[nt][st][ HMM_I ] = max(fromD, fromI, fromM ) + curI->emission[ residue ];
                    Tb[nt][st][ HMM_I ] = maxstate(fromD, fromI, fromM);
        
                    // D--deletion
                    fromM = Vt[nt][st-1][ HMM_M ] + preM->trprob[ HMM_D ];
                    fromI = Vt[nt][st-1][ HMM_I ] + preI->trprob[ HMM_D ];
                    fromD = Vt[nt][st-1][ HMM_D ] + preD->trprob[ HMM_D ];
                    Vt[nt][st][ HMM_D ] = max(fromD, fromI, fromM );
                    Tb[nt][st][ HMM_D ] = maxstate(fromD, fromI, fromM);

                     preM = curM;
                     preD = curD;
                     preI = curI;
             
                     curM = curM->nextM;
                     curD = preM->nextD;
                    curI = curM->nextI;
                }
        
                st=numstates-1;
                // E--end state
                fromM = Vt[nt-1][st-1][ HMM_M ] + preM->trprob[ HMM_M ];
    	           fromI = Vt[nt-1][st-1][ HMM_I ] + preI->trprob[ HMM_M ];
    	           fromD = Vt[nt-1][st-1][ HMM_D ] + preD->trprob[ HMM_M ];
    	           Vt[nt][st][ HMM_M ] = max(fromD, fromI, fromM );
                 Tb[nt][st][ HMM_M ] = maxstate(fromD, fromI, fromM);
            }
            //last column & last row
            //M  seqlen+1,  numstates,
            fromM = Vt[seqlen][numstates-2][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = Vt[seqlen][numstates-2][ HMM_I ] + preI->trprob[ HMM_M ];
            fromD = Vt[seqlen][numstates-2][ HMM_D ] + preD->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = max(fromD, fromI, fromM );
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = maxstate(fromD, fromI, fromM);
        }
        else
        {
            preM = beg;
            preI = beg->nextI;
            // curM = beg->nextM;
            // E--end state
            fromM = Vt[seqlen-1][0][ HMM_M ] + preM->trprob[ HMM_M ];
             fromI = Vt[seqlen-1][0][ HMM_I ] + preI->trprob[ HMM_M ];
             Vt[seqlen][1][ HMM_M ] = max( fromI, fromM );
            Tb[seqlen][1][ HMM_M ] = maxstate( fromI, fromM);
        
            //last column
            //M  seqlen+1,  numstates,
            fromM = Vt[seqlen][numstates-2][ HMM_M ] + preM->trprob[ HMM_M ];
            fromI = Vt[seqlen][numstates-2][ HMM_I ] + preI->trprob[ HMM_M ];
            Vt[ seqlen+1 ][numstates-1][ HMM_M ] = max(fromI, fromM );
            Tb[ seqlen+1 ][numstates-1][ HMM_M ] = maxstate(fromI, fromM);
        }
    }
    //output3Dtables(  );

    //****** find the maximal prob ****** [seq][state][ D=0/I=1/M=2 ]
    //nowtops=getTopnMaxProbs(  );  //get maximum prob
    nowtops=getGlobalProbs(  );     //global alignment
    for(i=0; i<nowtops; i++)
    {
         tops[i].cadno=i+1;
        VTBTraceBack( tops[ i ] );
    }
}
