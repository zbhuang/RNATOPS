// DynamicP.cpp: implementation of the DynamicP class.
//
//////////////////////////////////////////////////////////////////////
#pragma warning( disable: 4786 )

#include "DynamicP.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DynamicP::DynamicP()
{
	pStem			= NULL;
	pLoop			= NULL;
	pCongraph		= NULL;
	pTrainingSeqs	= NULL;
	node_mapping	= NULL;
	//
	stemstate		= NULL;

	seqlength		= 0;
	numtrainingseqs = 0;

	numnode			= 0;
	nodestart		= 0;

	genomeSequenceLength = 0;
	cGenomeSequence = NULL;
	cRGenomeSequence = NULL;

	scanLength = 0;
	scanSequence = NULL;

	candidateStems = NULL;

	r_scores = NULL;

	this->numhit = 0;
	this->pre_hit_num = 0;	//

	iTheWholeTreeNodeIndex			= 0;
	iOptimalTheWholeTreeNodeIndex	= 0;
	iFlagSearchReverseGenome		= 0;

	flagPrintDebugInfo			= 0;
	flagPrintDebugInputInfo		= 0;
	flagPrintDebugTreeInfo		= 0;
	flagPrintDebugDPWinInfo		= 0;
	flagPrintDebugTDInfo		= 0;
	flagPrintDebugDPSearchInfo	= 0;
}

DynamicP::~DynamicP()
{
	this->free_pTrainingSeqs();

	this->free_pStem();

	this->free_pLoop();

	this->free_pCongraph();

	this->free_node_mapping();

	//free memeory for genome sequence
	this->freeGenomeSequence();
	this->freeRGenomeSequence();

}

void DynamicP::setLeadingNum(int num)
{
	this->iLeading = num;
}

int DynamicP::getLeadingNum()
{
	return this->iLeading;
}

void DynamicP::free_pTrainingSeqs()
{
	int i;
	if(pTrainingSeqs != NULL) {
		for(i=0; i<this->getNumtrainingseqs(); i++) {
			delete [] pTrainingSeqs[i];
			pTrainingSeqs[i] = NULL;
		}
		delete [] pTrainingSeqs;
		pTrainingSeqs = NULL;
	}
}

void DynamicP::free_pStem()
{
	if(pStem != NULL)
		delete [] pStem;
	pStem = NULL;
}
void DynamicP::free_pLoop()
{
	if(pLoop != NULL)
		delete [] pLoop;
	pLoop = NULL;
}
void DynamicP::free_pCongraph()
{
	int i;
	if(pCongraph != NULL) {
		for(i=0; i<this->getNumnode(); i++) {
			delete [] pCongraph[i];
			pCongraph[i] = NULL;
		}
		delete [] pCongraph;
		pCongraph = NULL;
	}
}
void DynamicP::free_node_mapping()
{
	if(node_mapping != NULL)
		delete [] node_mapping;
	node_mapping = NULL;
}

void DynamicP::freeGenomeSequence()
{
	if(cGenomeSequence != NULL)
		delete [] cGenomeSequence;
	cGenomeSequence = NULL;
}

void DynamicP::freeRGenomeSequence()
{
	if(cRGenomeSequence != NULL)
		delete [] cRGenomeSequence;
	cRGenomeSequence = NULL;
}

int DynamicP::getSeqlength()
{
	return this->seqlength;
}

int DynamicP::getNumtrainingseqs()
{
	return this->numtrainingseqs;
}

void DynamicP::setTrainingSeqs(int numtrainseqs, int trainingseqlength, int ** trainingseqs)
{
	int i, j;

	this->seqlength = trainingseqlength;
	this->numtrainingseqs = numtrainseqs;

	this->pTrainingSeqs = new int* [numtrainseqs];
	if(this->pTrainingSeqs == 0) {
		printf("Fail to allocate memory!\n");
		exit(0);
    }

	for(i=0; i<numtrainseqs; i++)
	{
		this->pTrainingSeqs[i] = new int[trainingseqlength];
		for(j=0; j<trainingseqlength; j++)
		{
			this->pTrainingSeqs[i][j] = trainingseqs[i][j];
		}
	}
}

void DynamicP::printTrainingSeqs()
{
	int i, j;
	int num = this->getNumtrainingseqs();
	int length = this->getSeqlength();

	for(i=0; i<num; i++)
	{
		for(j=0; j<length; j++)
			cout<<numtochar(this->pTrainingSeqs[i][j]);
		cout<<endl;
	}
}

double DynamicP::getThreshold()
{
	return this->threshold;
}
void DynamicP::setThreshold(double dThreshold)
{
	this->threshold = dThreshold;
}
//genome sequence
void DynamicP::setGenomeSequence(int ilength, char * cGenomeSequence)
{
	this->genomeSequenceLength = ilength;
	this->cGenomeSequence = new char[ilength+1];
	strcpy(this->cGenomeSequence, cGenomeSequence);
}
int	DynamicP::getGenomeSequenceLength()
{
	return this->genomeSequenceLength;
}
char * DynamicP::getGenomeSequence()
{
	return this->cGenomeSequence;
}
void DynamicP::printGenomeSequence()
{
	cout<<"printGenomeSequence()"<<endl;
	cout<<this->cGenomeSequence<<endl;
}

//reverse genome sequence
void DynamicP::setRGenomeSequence(int ilength, char * cRGenomeSequence)
{
	this->cRGenomeSequence = new char[ilength+1];
	strcpy(this->cRGenomeSequence, cRGenomeSequence);
}
char * DynamicP::getRGenomeSequence()
{
	return this->cRGenomeSequence;
}
void DynamicP::printRGenomeSequence()
{
	cout<<"printRGenomeSequence()"<<endl;
	cout<<this->cRGenomeSequence<<endl;
}

void DynamicP::setNumstems(int numstems)
{
	this->numstems = numstems;
}

void DynamicP::setNumloops(int numloops)
{
	this->numloops = numloops;
}
int DynamicP::getNumstems()
{
	return this->numstems;
}

int DynamicP::getNumloops()
{
	return this->numloops;
}

StemInfo * DynamicP::getPStem()
{
	return this->pStem;
}

LoopInfo * DynamicP::getPLoop()
{
	return this->pLoop;
}

void DynamicP::setStems(int numstems, StemLocation *stemarray)
{
	int i;

	this->setNumstems(numstems);

	this->pStem = new StemInfo[numstems];

    for(i=0; i<numstems; i++) {
		this->pStem[i].charid	= stemarray[i].charid;	//
		if(this->getNumPastaLine() == TWO_PASTALINE)
		{
			//support stemid index in two pastalines i.e. A1, A2, a1, a2. 20081015
			string stemstrid	= stemarray[i].getStemStrID();
			if(this->pStem[i].charid == stemstrid[0])
				this->pStem[i].charid_idx	= stemstrid[1];
			else
				cout<<"Something wrong in parsing stemid and stemstrid"<<endl;
		} else {
			this->pStem[i].charid_idx	= NULL;
		}
		this->pStem[i].start1	= stemarray[i].Pos[0];
		this->pStem[i].end1		= stemarray[i].Pos[1];
		this->pStem[i].start2	= stemarray[i].Pos[2];
		this->pStem[i].end2		= stemarray[i].Pos[3];
		if(i == 0) {
			this->iStemStartPosInStructureLine = stemarray[i].Pos[0];
		}
		if(i == (numstems-1)) {
			this->iStemEndPosInStructureLine = stemarray[i].Pos[3];
		}
    }
}

void DynamicP::printStems()
{
	cout<<"----Stem info----"<<endl;
	int i;
	int num = this->getNumstems();
	cout<<"- num of stem: "<<num<<" -"<<endl;
	StemInfo *	_pStem = this->getPStem();
    for(i=0; i<num; i++) {
		cout<<_pStem[i].start1<<"\t"
			<<_pStem[i].end1<<"\t"
			<<_pStem[i].start2<<"\t"
			<<_pStem[i].end2<<endl;
    }
}

void DynamicP::setLoops(int numloops, LoopLocation *looparray)
{
	int i;
	this->setNumloops(numloops);
	this->pLoop = new LoopInfo[numloops];

	for(i=0; i<numloops; i++) {
		this->pLoop[i].start= looparray[i].Pos[0];
		this->pLoop[i].end	= looparray[i].Pos[1];
	}
}

void DynamicP::printLoops()
{
	cout<<"----Loop info----"<<endl;
	int i;
	int num = this->getNumloops();
	cout<<"- num of loop: "<<num<<" -"<<endl;
	LoopInfo *	_pLoop = this->getPLoop();
    for(i=0; i<num; i++) {
		cout<<_pLoop[i].start<<"\t"
			<<_pLoop[i].end<<endl;
    }
}

int	DynamicP::getNumnode()
{
	return this->numnode;
}

int	DynamicP::getNodestart()
{
	return this->nodestart;
}

//construct the graph to represent the consensus structure.
void DynamicP::build_tree()
{
	std::vector< TRnaGraphNode<int> > v;
	std::map< int, std::pair<int, TRnaGraphNode<int> > > m;
	int i, j;
	int numstem = this->getNumstems();
	int numloop = this->getNumloops();
	int numnode = 2*numstem + 2;
	this->numnode = numnode;

	//allocate memory for pCongraph
	this->pCongraph = new int* [numnode];
    for(i=0; i<numnode; i++) {
		this->pCongraph[i] = new int[numnode];
    }

    for(i=0; i<numnode; i++){
		for(j=0; j<numnode; j++){
			this->pCongraph[i][j] = NONCONNECTED;
		}
    }

    this->nodestart = 1;
    for(i=0; i<numstem; i++)
	{
		//Assign the nodes i+1 and i+2 to the stem i.     
		this->pStem[i].lg_id = nodestart;  
		m.insert(std::pair<int,std::pair<int, TRnaGraphNode<int> > >(pStem[i].start1, std::pair<int, TRnaGraphNode<int> >(pStem[i].start2,TRnaGraphNode<int>(-1,-1,-1,pStem[i].lg_id))));
		this->pStem[i].rg_id = nodestart+1;
		m.insert(std::pair<int,std::pair<int, TRnaGraphNode<int> > >(pStem[i].start2, std::pair<int, TRnaGraphNode<int> >(pStem[i].start1,TRnaGraphNode<int>(-1,-1,-1,pStem[i].rg_id))));
		this->pCongraph[nodestart][nodestart+1] = NONDIRECTED;
		this->pCongraph[nodestart+1][nodestart] = NONDIRECTED;
		this->nodestart += 2;
    }

	i = 0;
	std::map<int,std::pair<int, TRnaGraphNode<int> > >::iterator it = m.begin();
	while( it != m.end() )
	{
		it->second.second.ID(i);
		it->second.second.Next(i+1);
		m.find(it->second.first)->second.second.Pair(i);
		it++;
		i++;
	}

	it = m.begin();
	while( it != m.end() )
	{
		v.push_back( it->second.second );
		it++;
	}
/*
	cout << "<--------------" << endl;
	for( i = 0; i < v.size(); i++ )
	{
		cout << v[i].ID() << ' ' << v[i].Next() << ' ' << v[i].Pair() << ' ' << v[i].Payload() << endl;
	}
	cout << "-------------->" << endl;
*/
	TRnaGraphTd< int, Tree_bag > rgtd( v );
	rgtd.Construct();
//	cout << "###############" << endl;
//	rgtd.Print();
	//get the root of the tree-bag
	this->pTDRoot = rgtd.Decomposition();

//	this->IterateTD(this->pTDRoot);
//	this->IterateTD_Payload_topdown(this->pTDRoot);
//	this->IterateTD_Payload_bottomup(this->pTDRoot);
	mypTDRoot = this->copy_treedecomp(this->pTDRoot);

/*
	//20080630
	if(this->pTDRoot != NULL) {
		delete this->pTDRoot;
		this->pTDRoot = NULL;
	}
*/
	//20080702
	rgtd.Release();

//	this->printTreebags(mypTDRoot);

    //Construct the edges in the graph, go through all the
    //loops in the consensus structure and make connections between nodes.
	int source = -1, sink = -1;
	int loop_start, loop_end, start_found, end_found;

    for(i=0; i<numloop; i++)
	{
		loop_start	= this->pLoop[i].start;
		loop_end	= this->pLoop[i].end;
		start_found	= 0;
		end_found	= 0;

		for(j=0; j<numstem; j++)
		{
			//Go through the list of stems to check for the corresponding stem; 
			if(loop_start == this->pStem[j].end1+1) {
				source		= this->pStem[j].lg_id;
				start_found = 1;
			} else if(loop_start == this->pStem[j].end2+1) {
				source		= this->pStem[j].rg_id;
				start_found = 1;
			}
			if(loop_end == this->pStem[j].start1-1) {
				sink		= pStem[j].lg_id;
				end_found	= 1;
			} else if(loop_end == this->pStem[j].start2-1) {
				sink		= pStem[j].rg_id;
				end_found	= 1;
			}
		}//for j.
		if(start_found == 0) {
			source = 0;
		}

		if(end_found == 0) {
			sink = this->getNodestart();
		}

		this->pCongraph[source][sink] = DIRECTED;
	}//for i.
/*
	//Now start checking if there are contiguous stems without a loop in between.
	int s1_start1, s1_end1, s1_start2, s1_end2;
	int s2_start1, s2_end1, s2_start2, s2_end2;

	for(i=0; i<numstem; i++)
	{
		s1_start1 = this->pStem[i].start1;
		s1_end1 = this->pStem[i].end1;
		s1_start2 = this->pStem[i].start2;
		s1_end2 = this->pStem[i].end2;

		for(j=0; j<numstem; j++)
		{
			s2_start1	= this->pStem[j].start1;
			s2_end1		= this->pStem[j].end1;
			s2_start2	= this->pStem[j].start2;
			s2_end2		= this->pStem[j].end2;

			if(s1_end1+1 == s2_start1)
			{
				//If the end of the stem1 is the starting of the start2.      
				source	= this->pStem[i].lg_id;
				sink	= this->pStem[j].lg_id;
				this->pCongraph[source][sink] = WEIGHTZERO;  //Set the weight of the edge to be zero.
			}
			if(s2_end2+1 == s1_start2)
			{
				//If the end of the right stem2 is the starting of the right half of stem1.
				source	= this->pStem[j].rg_id;
				sink	= this->pStem[i].rg_id;
				this->pCongraph[source][sink] = WEIGHTZERO;  //Set the weight of the edge to be zero.
			}
			if(s2_end1+1 == s1_start2)
			{
				//If the end of the left stem2 is the starting of the right half of stem1. 
				source	= this->pStem[j].lg_id;
				sink	= this->pStem[i].rg_id;
				this->pCongraph[source][sink] = WEIGHTZERO; //Set the weight of the edge to be zero.
			}
			if(s1_end2+1 == s2_start1)
			{
				//If the end of the right stem1 is the starting of the left half of stem2. 
				source	= this->pStem[i].rg_id;
				sink	= this->pStem[j].lg_id;
				this->pCongraph[source][sink] = WEIGHTZERO; //Set the weight of the edge to be zero.
			}//if
		}//for j.
	}//for i.
*/
	//Finally, the source and sink are checked to see if they should be connected to the stems.     
	for(i=0; i<numstem; i++)
	{
//		if(this->pStem[i].start1 == 0)
		if(this->pStem[i].start1 == this->iStemStartPosInStructureLine)
		{
			//If there is a stem starting with the source.   
			this->pCongraph[0][this->pStem[i].lg_id] = WEIGHTZERO;
		}

//		if(this->pStem[i].end2 == this->getSeqlength()-1) {
		if(this->pStem[i].end2 == this->iStemEndPosInStructureLine) {
			//If there is a stem ending with the sink.
			this->pCongraph[this->pStem[i].rg_id][nodestart] = WEIGHTZERO;
		}
	}
}

void DynamicP::print_tree()
{
	cout<<"----Graph Info----"<<endl;
	int i, j;
	int numnode = this->getNumnode();
	for(i=0; i<numnode; i++)
	{
		for(j=0; j<numnode; j++)
			cout<<this->pCongraph[i][j]<<" ";
		cout<<endl;
	}
	cout<<"numnode="<<numnode<<endl;
}

//Find the left stem id and the right stem id (allloops.left_id/right_id/lg_id/rg_id)
//int left_id;	- The identity of the stem where the left boundary point belongs -1 if it is the source.      
//int right_id; - The identity of the stem where the right boundary point belongs -1 if it is the sink. 
//int lg_id;	- The left graph identity for the left end of the loop.     
//int rg_id;	- The right graph identity for the right end of the loop.
void DynamicP::identify_loopneighbor()
{
	int numloop		= this->getNumloops();
	int numstem		= this->getNumstems();
	int nodestart	= this->getNodestart();

	int i, j;
	int l_start, l_end;
	int l_sid, r_sid;
	int l_gid = -1, r_gid = -1;

	for(i=0; i<numloop; i++)
	{
		//Go through every loop in the structure.
		l_start = this->pLoop[i].start;
		l_end	= this->pLoop[i].end;

		l_sid = -1;
		r_sid = -1;

		for(j=0; j<numstem; j++)
		{
			//Go through every stem in the secondary structure to check for the ending points of loops. 
            if(this->pStem[j].end1+1 == l_start || this->pStem[j].end2+1 == l_start)
			{
				l_sid = j;
				if(this->pStem[j].end1+1 == l_start)
					l_gid = this->pStem[j].lg_id;
				else
					l_gid = this->pStem[j].rg_id;
			}
            if(l_end+1 == this->pStem[j].start1 || l_end+1 == this->pStem[j].start2)
			{
				r_sid = j;
				if(this->pStem[j].start1 == l_end+1)
					r_gid = this->pStem[j].lg_id;
				else
					r_gid = this->pStem[j].rg_id;
            }
		}//for j.
		this->pLoop[i].left_stemid	 = l_sid;
		this->pLoop[i].right_stemid = r_sid;
		this->pLoop[i].lg_id	 = l_gid;
		this->pLoop[i].rg_id	 = r_gid;

		if(l_sid == -1) {
			//If the leftmost end of the loop is the source.
			this->pLoop[i].lg_id = 0;
		}
		if(r_sid == -1) {
			//If the rightmost end of the loop is the sink. 
			this->pLoop[i].rg_id = nodestart;
		}
	}//for i.
}

void DynamicP::printLoopNeighbor()
{
	cout<<"----Loop Neighbor Info----"<<endl;
	int numloop = this->getNumloops();

	int i;
	for(i=0; i<numloop; i++) {
        cout<<"this->pLoop["<<i<<"].left_stemid is: "<<this->pLoop[i].left_stemid<<endl;
        cout<<"this->pLoop["<<i<<"].right_stemid is: "<<this->pLoop[i].right_stemid<<endl;
	}

	for(i=0; i<numloop; i++) {
        cout<<"this->pLoop["<<i<<"].lg_id is: "<<this->pLoop[i].lg_id<<endl;
        cout<<"this->pLoop["<<i<<"].rg_id is: "<<this->pLoop[i].rg_id<<endl;
	}
}

//identify the loop id given the left stem id and right stem id
//it is used in the printing out the whole structure alignment info
//including stem and loop
int DynamicP::identifyLoopId(int left_gid, int right_gid)
{
	int		iLoopId = -1;
	int		i;
	int		numloop = this->getNumloops();
	bool	bFound = false;

	for(i=0; i<numloop && !bFound; i++) {
		if((this->pLoop[i].lg_id == left_gid) 
		&& (this->pLoop[i].rg_id == right_gid))
		{
			iLoopId = i;
			bFound	= true;
		}
	}
	return iLoopId;
}
void DynamicP::build_node_mapping()
{
	int nodestart = this->getNodestart();
	node_mapping = new C_graphnode[nodestart+1];

    //Initialize the type and mapping_id of the source and the sink.
    node_mapping[0].node_type = 0;
    node_mapping[0].mapping_id = 0;

	//do the initialization
	if(nodestart>1) {
		for(int i=0; i<nodestart; i++) {
			node_mapping[i].node_type = 1;
		}
	}

    node_mapping[nodestart].node_type = 0;
    node_mapping[nodestart].mapping_id = 0;
}

void DynamicP::allocateStemState(int num)
{
	int i;
	stemstate = new int[num];
	for(i=0; i<num; i++) {
		stemstate[i] = 0;
	}
}

void DynamicP::deallocateStemState()
{
	delete [] stemstate;
	stemstate = NULL;
}

//Compute the number of stems.    
void DynamicP::get_stemcomponents(int num, int& stop) 
{
	int numstem = this->getNumstems();

    int i, j;
    int s1_start1, s1_end1, s1_start2, s1_end2;
    int s2_start1, s2_end1, s2_start2, s2_end2;

    int maxnum, maxstem = 0;

	//Store the number of stems that cross with each given stem in the consensus structure.
    int * numcross = new int[numstem];  
    for(i=0; i<numstem; i++) {
		numcross[i]=0;
    }
 
    for(i=0; i<numstem; i++)
	{
		//Go through all possible stems. 
		if(stemstate[i] == num)
		{
			s1_start1	= this->pStem[i].start1;
			s1_end1		= this->pStem[i].end1;
			s1_start2	= this->pStem[i].start2;
			s1_end2		= this->pStem[i].end2;

			for(j=0; j<numstem; j++)
			{
				//Go through the stems with an integer index greater than i. 
				if(stemstate[j] == num)
				{
					s2_start1	= this->pStem[j].start1; 
					s2_end1		= this->pStem[j].end1;
					s2_start2	= this->pStem[j].start2;
					s2_end2		= this->pStem[j].end2;
					if(s1_start1 < s2_start1 && s1_end2 < s2_end2 && s2_end1 < s1_start2 ) {
						//Increase the number of crossing stems of stem i by 1.
						numcross[i]++;
					} else if(s1_start1 > s2_start1 && s1_end2 > s2_end2 && s1_end1<s2_start2) {
						//Increase the number of crossing stems of stem i by 1.
						numcross[i]++;
					}
				}//if stemstate[j]
			}//for j.
		}//if stemstate[i]
	}//for i.

    //Find the stem with the maximum number of crossing stems. 
    maxnum = 0;

    for(i=0; i<numstem; i++) {
		if(maxnum < numcross[i]) {
			maxnum	= numcross[i];
			maxstem	= i;
        }
	}

	if(maxnum == 0) {
		//Stop the procedure if there is no crossing stems in the array.   
		stop = 1;
	} else {
		stemstate[maxstem] = num+1;
	}

	delete [] numcross;
	numcross = NULL;
}
/*
//Separate the crossing stems into different components those will be considered 
//when the tree decomposition of the graph needs to be constructed.        
int DynamicP::split_stemcomponents()
{
	int numstem = this->getNumstems();

	int		i;
	int		s_id, mappedstem;
	int		stop = 0;

	//this->allocateStemState(numstem);

	s_id = 0;
	mappedstem = 0;

	while(mappedstem < numstem)
	{
		while(stop == 0) {
			this->get_stemcomponents(s_id, stop);
        }

        for(i=0; i<numstem; i++)
		{
			if(stemstate[i] == s_id)
			{
				mappedstem++;
				node_mapping[this->pStem[i].lg_id].node_type	= 1;	//set the node type to be 1.
				node_mapping[this->pStem[i].lg_id].mapping_id	= s_id;
				node_mapping[this->pStem[i].rg_id].node_type	= 1;	//set the node type to be 1.
				node_mapping[this->pStem[i].rg_id].mapping_id	= s_id;
			}//if
		}//for
		s_id++;
		stop=0; ///get stop back to 0
	}//while

	//this->deallocateStemState();
	//s_id is the number of trees we have obtained for the mapping.
	return s_id;
}
*/
//Print out the nodes in one tree bag. 
void DynamicP::printTreenode(Tree_bag *node)
{
    //int i;
    Node_list *nodehead;

	if(node == NULL)
		return;

    nodehead = node->nhead;
    printf("{");
    while(nodehead) {
		//Go through the list of nodes, print out the head information.
		printf("%d ", nodehead->g_node);
		nodehead = nodehead->next;
    }
    printf("}");
}

//Print out the tree_bag, given the root of the tree node. 
void DynamicP::printTreebags(Tree_bag *node)
{
	int i;
	//node_list *nodehead;
	if(node == NULL)
		return;
	printTreenode(node);

	printf(":=");
	for(i=0; i<node->childnum; i++)
		printTreenode(node->pChild[i]);

	printf("\n");

	for(i=0; i<node->childnum; i++) {
		//printf("node->ptnum is: %d\n", node->ptnum);
		printTreebags(node->pChild[i]);
	}
}

//Search for a given node in a linked list. 
//It returns 1 if it is found and returns 0 if it is not.
int DynamicP::searchlist(Node_list *n_head, int nodeid) 
{
	Node_list *nhead;
	nhead = n_head;
	while(nhead) {
		if(nhead->g_node == nodeid) {
			return 1;
		}
		nhead = nhead->next;
	}
	return 0; 
}

//according to the current treebag, rearrange the child's nodelist
void DynamicP::reArrangeChildNodelist(Tree_bag *root)
{
	int i = 0;

	if(root == NULL)
		return;

	Node_list *root_nhead = root->nhead;
//	print_treenodelist(root_nhead);

	Node_list *child_nhead;

//	int num = 0;
	int i_node;
	for(i=0; i<root->childnum; i++) {
		child_nhead = (root->pChild[i])->nhead;

		root->pChild[i]->shared_node_num = 0;
//		cout<<"\t";
//		print_treenodelist(child_nhead);

		Node_list *new_child_nhead = NULL;
		Node_list *new_child_nhead_tail = NULL;
		Node_list *new_child_nhead_head = NULL;

//		cout<<"\t";
		while(child_nhead) {
			i_node = child_nhead->g_node;
//			cout<<" "<<i_node<<" ";
			if(this->searchlist(root_nhead, i_node) == 1) {
				root->pChild[i]->shared_node_num++;
				new_child_nhead = new Node_list[1];
				new_child_nhead->g_node		= i_node;
				new_child_nhead->next		= NULL;
				if(new_child_nhead_tail != NULL) {
					new_child_nhead_tail->next	= new_child_nhead;
					new_child_nhead_tail		= new_child_nhead;
				} else {
					new_child_nhead_tail		= new_child_nhead;
				}
				if(new_child_nhead_head == NULL)
					new_child_nhead_head	= new_child_nhead;
			}
			child_nhead = child_nhead->next;
		}

//		cout<<"\t";
		child_nhead = (root->pChild[i])->nhead;
		while(child_nhead) {
			i_node = child_nhead->g_node;
//			cout<<" "<<i_node<<" ";
			if(this->searchlist(root_nhead, i_node) == 0) {
				new_child_nhead = new Node_list[1];
				new_child_nhead->g_node		= i_node;
				new_child_nhead->next		= NULL;
				if(new_child_nhead_tail != NULL) {
					new_child_nhead_tail->next	= new_child_nhead;
					new_child_nhead_tail		= new_child_nhead;
				} else {
					new_child_nhead_tail		= new_child_nhead;
				}
				if(new_child_nhead_head == NULL)
					new_child_nhead_head	= new_child_nhead;
			}
			child_nhead = child_nhead->next;
		}

		(root->pChild[i])->nhead = new_child_nhead_head;

//		cout<<endl;
	}
//	this->printTreebags(root);

	for(i=0; i<root->childnum; i++) {
		reArrangeChildNodelist(root->pChild[i]);
	}
}

void DynamicP::setTreeBagParent(Tree_bag *root) {
	int i;
	if(root == NULL)
		return;

	for(i=0; i<root->childnum; i++) {
		(root->pChild[i])->pParent = root;
	}

	for(i=0; i<root->childnum; i++) {
		setTreeBagParent(root->pChild[i]);
	}
}
void DynamicP::printTreeBagParent(Tree_bag *root) {
	int i;
	if(root == NULL)
		return;

	for(i=0; i<root->childnum; i++) {
		printTreenode(root->pChild[i]);
		cout<<"'s parent=";
		printTreenode((root->pChild[i])->pParent);
		cout<<endl;
	}

	for(i=0; i<root->childnum; i++) {
		printTreeBagParent(root->pChild[i]);
	}
}

void DynamicP::removeTreeBagParent(Tree_bag *root) {
	int i;
	if(root == NULL)
		return;

	for(i=0; i<root->childnum; i++) {
		removeTreeBagParent(root->pChild[i]);
	}

	root->pParent = NULL;

//	for(i=0; i<root->childnum; i++) {
//		(root->pChild[i])->pParent = NULL;
//	}

}
//Printing out the conformational graph (only the connectivity so far).  
void DynamicP::print_outstructuregraph()
{
	int		numnode = this->getNumnode();

    int		i, j;

    printf("<graph>\n"); 
    for(i=0; i<numnode; i++)
	{
        for(j=0; j<numnode; j++) {
			printf("%d ", this->pCongraph[i][j]);
		}
		printf("\n");
	}
	printf("</graph>\n"); 
}

//Convert between all the vertices in the graph and the stems in the secondary structure.
void DynamicP::get_conversion_gtos()
{
	int numstem = this->getNumstems();
	int i;
	int lg_id, rg_id;

	//Allocate enough enough memory to store the information .
	gidtosid = new Sid_info[2*numstem + 2];

	gidtosid[0].s_id			= -1;	//The starting node (source).
	gidtosid[0].left_or_right	= -1;	//

	gidtosid[nodestart].s_id			= -2;	//The ending node (sink).
	gidtosid[nodestart].left_or_right	= -1;	//

	for(i=0; i<numstem; i++)
	{
		//Now, go through all the stems in the secondary structure and store the conversion relationship in the table.
        lg_id = this->pStem[i].lg_id;
        rg_id = this->pStem[i].rg_id;

        gidtosid[lg_id].left_or_right = 0;
        gidtosid[lg_id].s_id = i;

        gidtosid[rg_id].left_or_right = 1;
        gidtosid[rg_id].s_id = i;
    }
}

void DynamicP::print_gidtosid()
{
	printf("--- print out gidtosid array ---\n");
	int numstem = this->getNumstems();
	int i;

	for(i=0; i<2*numstem+2; i++) {
		cout<<"gidtosid["<<i<<"].s_id="<<gidtosid[i].s_id<<", left_or_right="<<gidtosid[i].left_or_right<<endl;
	}
}

void DynamicP::allocate_dptable(Tree_bag *root)
{
	if(root == NULL)
		return;

//	this->printTreenode(root);	cout<<endl;
	int		i, nodenum, stnum, iTotal;
	nodenum = root->nodenum;
	stnum	= root->stnum;

	//Now start allocating memory needed for storing intermediate computational results.
	root->enum_arrays = new int[stnum];
	root->enum_arrayo = new int[nodenum];

	if(root->pParent == root)	//the decision of being the root
	{
		iTotal = 1;
	} else {
		//Applying the blocking technique proposed by Dr.Cai
		iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(root->shared_node_num)));
	}
	
	root->t_head = new New_node_combination[iTotal];

	for(i=0; i<root->childnum; i++)
		this->allocate_dptable(root->pChild[i]); 
}

void DynamicP::initDPTable(Tree_bag *root)
{
	int		i, nodenum, iTotal;
	nodenum = root->nodenum;

	if(root->pParent == root)	//the decision of being the root
	{
		iTotal = 1;
	} else {
		//Applying the blocking technique proposed by Dr.Cai
//		iTotal = floor(pow(double(this->iLeading), double(root->shared_node_num)));
		iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(root->shared_node_num)));
	}

	for(i=0; i<iTotal; i++) {
		root->t_head[i].score			= SMALLEST;
		root->t_head[i].optimal			= -1;
		root->t_head[i].optimalIndex	= -1;
	}

	for(i=0; i<root->childnum; i++)
		this->initDPTable(root->pChild[i]); 
}

//Construct the list of stems or half stems for all the nodes in the tree.
void DynamicP::construct_stemlist(Tree_bag *root) 
{
	int		node1, node2 = -1;
	int		node_cnt1, node_cnt2; 
	int		found_pair;
	int		stem_num;
	Strunit *	head = NULL, *tail = NULL;
	Node_list *	n_head1 = NULL, *n_head2 = NULL;
     
	n_head1 = root->nhead;
	head	= NULL;

	while(n_head1) {
		n_head1->visited = 0;	//Initialize the visited bit to be zero.
		n_head1 = n_head1->next;
	}

	node_cnt1	= 0;
	stem_num	= 0;
	n_head1		= root->nhead;

	while(n_head1)
	{
		if(n_head1->visited == 0)
		{
			node1 = n_head1->g_node;
			//just for trace. 
//			printf("node1=%d\n", node1);
			found_pair	= 0;
			n_head2		= n_head1->next; 
			node_cnt2	= node_cnt1+1;
             
			while(n_head2 != NULL && n_head2->visited != 0) {
				n_head2 = n_head2->next;
				node_cnt2++;
			}

			///if(n_head2!=NULL && n_head2->visited==0){
			if((n_head2 != NULL && n_head2->visited == 0) 
				|| (n_head2 == NULL && n_head1->visited == 0))
			{
                //If node2 has not been visited yet.
                while(n_head2)
				{
					node2 = n_head2->g_node;
					//just for trace. 
//					printf("\tnode2=%d\n", node2);
					if(gidtosid[node1].s_id == gidtosid[node2].s_id)
					{
						//If node1 and node2 belongs to the same stem.
						found_pair = 1;
						n_head2->visited = 1;
						break;
					}
					node_cnt2++;
					n_head2 = n_head2->next;
				}

                if(head == NULL && gidtosid[node1].s_id >= 0)
				{
					head = new Strunit[1];	//allocate memory
					if(found_pair == 1)
					{
						head->type = 1; 
						if(gidtosid[node1].left_or_right == 0)
						{ 
							head->g_id1 = node1;
							head->g_id2 = node2;
							head->left_nid = node_cnt1;
							head->right_nid = node_cnt2;
							head->stem_id = gidtosid[node1].s_id;
						} else {
							head->g_id1 = node2;
							head->g_id2 = node1;
							head->left_nid = node_cnt2;
							head->right_nid = node_cnt1;
							head->stem_id = gidtosid[node1].s_id;
						}
						stem_num++;
					} else {
						head->type = 2;
						head->g_id1 = node1;
						head->g_id2 = NEGATIVE_ONE;		//20081019
						head->left_nid = node_cnt1;
						head->right_nid= NEGATIVE_ONE;	//20081019
						head->stem_id = gidtosid[node1].s_id;
						stem_num++;
					}
					head->next = NULL;
					tail = head;
				}//if
                else if(gidtosid[node1].s_id>=0)
				{
					tail->next = new Strunit[1];	//allocate memory
					tail = tail->next;

					if(found_pair == 1)
					{
						tail->type = 1;
						if(gidtosid[node1].left_or_right == 0)
						{
							tail->g_id1 = node1;
							tail->g_id2 = node2;
							tail->left_nid = node_cnt1;
							tail->right_nid = node_cnt2;
							tail->stem_id = gidtosid[node1].s_id;
						} else {
							tail->g_id1 = node2;
							tail->g_id2 = node1;
							tail->left_nid = node_cnt2;
							tail->right_nid = node_cnt1;
							tail->stem_id = gidtosid[node1].s_id;
						}
						stem_num++;
					} else {
						tail->type = 2;
						tail->g_id1 = node1;
						tail->g_id2 = NEGATIVE_ONE;		//20081019
						tail->left_nid = node_cnt1;
						tail->right_nid= NEGATIVE_ONE;	//20081019
						tail->stem_id = gidtosid[node1].s_id;
						stem_num++;
					}
					tail->next = NULL;
				}//else if
			}//if(n_head2!=NULL && n_head2->visited==0)
		}//if
		node_cnt1++;
		n_head1 = n_head1->next;
	}//while
	root->stnum = stem_num; 
	root->su_stem = head;
}

void DynamicP::tree_stemlist(Tree_bag *root)
{
	int		i;
	this->construct_stemlist(root);

	for(i=0; i<root->childnum; i++)
		this->tree_stemlist(root->pChild[i]);
}

//Construct the list of loops for a given node in the tree.
//At the end, also incorporate stem information into the data structure.
void DynamicP::construct_looplist(Tree_bag *root)
{
	int numloop = this->getNumloops();
	int		i;
	int		node1, node2;
	int		loop_num;
	Node_list * n_head1 = NULL, * n_head2 = NULL;
	Strunit * head = NULL, * tail = NULL;
	Strunit * s_head = NULL;

	n_head1 = root->nhead;
	loop_num = 0;
	head = NULL;

	while(n_head1)
	{
		node1 = n_head1->g_node;
		//just for trace. 
//		printf("node1=%d\n", node1);
		n_head2 = n_head1->next;
        while(n_head2)
		{
			node2 = n_head2->g_node;
			//just for trace. 
//			printf("\tnode2=%d\n", node2);
			if(this->pCongraph[node1][node2] == DIRECTED || this->pCongraph[node2][node1] == DIRECTED)
			{
				for(i=0; i<numloop; i++)
				{
					if((node1 == this->pLoop[i].lg_id && node2 == this->pLoop[i].rg_id)
						|| (node1 == this->pLoop[i].rg_id && node2 == this->pLoop[i].lg_id))
					{
						if(head == NULL)
						{
							head = new Strunit[1];	//allocate memory
							head->type = 0;
							if(this->pCongraph[node1][node2] == DIRECTED) {
								head->g_id1 = node1;
								head->g_id2 = node2;
							} else {
								head->g_id1 = node2;
								head->g_id2 = node1;
							}
							head->loop_id = i;
							head->next = NULL;
							tail = head;
						} else {
							tail->next = new Strunit[1];	//allocate memory
							tail = tail->next;
							tail->type = 0;

							if(this->pCongraph[node1][node2] == DIRECTED){
								tail->g_id1 = node1;
								tail->g_id2 = node2;
							} else {
								tail->g_id1 = node2;
								tail->g_id2 = node1;
							}
							tail->loop_id = i;
							tail->next = NULL;
						}
						loop_num++;
					}//if
				}//for i.
			}//if
			n_head2 = n_head2->next;
		}//while
		n_head1 = n_head1->next;
	}//while
	root->lpnum = loop_num;
	root->su_loop = head;

	//Now, go through the linked list to incorporate more information into the data structure.
	while(head) {
		s_head = root->su_stem;
		//Initialize the values of lstem and rstem to be NULL.   
		head->lstem = NULL;
		head->rstem = NULL;

		while(s_head)
		{
			if(s_head->g_id1 == head->g_id1 || s_head->g_id2 == head->g_id1) {
				head->lstem = s_head;   
			}

			if(s_head->g_id2 == head->g_id2 || s_head->g_id1 == head->g_id2) {
				head->rstem = s_head;
			}

			s_head = s_head->next;
		}
		head = head->next;
	}
}

//Generate the loop list for all the nodes in the tree.  
void DynamicP::tree_looplist(Tree_bag *root)
{
    int		i;
/*
	cout<<"----tree_looplist----"<<endl;
	this->printTreenode(root);
	cout<<endl;
*/
	this->construct_looplist(root);

    for(i=0; i<root->childnum; i++)
			this->tree_looplist(root->pChild[i]);
}

//Generate the list of stems that is shared by both a parent node and a child node. 
//Also the array stem_image in the root node indicates whether a
//a particular stem is also contained in a child node or not.  
void DynamicP::generate_sharedstemlist(Tree_bag * root, Tree_bag * child, int c_id)
{
	int		i;
	int		prog_id;
	int		node1, node2;
	int		n_found1, n_found2;
	Strunit		*s_head;
	Node_list	*n_head;
    
	//Check the list of stems.
	s_head = root->su_stem;

	//Now allocate sufficient memory for the array. 
	root->stem_image[c_id] = new int[root->nodenum];	//allocate memory

	for(i=0; i<root->nodenum; i++) {
		root->stem_image[c_id][i] = -1; 
	}

	//Now, go through the list of stems in the stem list.   
	prog_id = 0;

	while(s_head)
	{
		n_found1 = 0;   ///
		n_found2 = 0;
		node1 = -1;
		node2 = -1;

		if(s_head->type == 1) {
           //If the s_head is a whole stem.
           node1 = s_head->g_id1;
           node2 = s_head->g_id2;
           n_found1 = 0;
           n_found2 = 0;
        } else {
           //Otherwise, s_head contains a half stem
			//and only the graph node for that particular half stem is searched in the following procedure.
           node1 = s_head->g_id1;
           n_found1 = 0;
        }

        n_head = child->nhead;

		while(n_head) {
			if(n_head->g_node == node1) {
				n_found1 = 1;
			}
			if(n_head->g_node == node2) {
				n_found2 = 1;
			}
			n_head = n_head->next;
		}

		if(s_head->type == 1 && n_found1 == 1 && n_found2 == 1) {
			//Set the corresponding flag value to be one.
			root->stem_image[c_id][prog_id] = 1;
		}

		if(s_head->type == 2 && n_found1 == 1) {
			root->stem_image[c_id][prog_id] = 1;
		}
		prog_id++;
		s_head = s_head->next;
    }//while
}

//Similarly, Generate the list of loops shared by a child and a root. 
//The loop_image array in the node root will takes the value of either 1 or -1, 
//when the value is 1, it means that the child also contains the loop.   
void DynamicP::generate_sharedlooplist(Tree_bag * root, Tree_bag * child, int c_id)
{
     int	i;
     int	prog_id;
     int	node1, node2;
     int	n_found1, n_found2;
     Strunit *		l_head;
     Node_list *	n_head;

     //Now, allocate memory to the array.
     l_head = root->su_loop;

     root->loop_image[c_id] = new int[root->nodenum];

	//Initialize all the array elements to be -1.
	for(i=0; i<root->nodenum; i++) {
		root->loop_image[c_id][i] = -1;
	}

	//Now go through the list of loops.
	prog_id = 0;

	while(l_head)
	{
        node1 = l_head->g_id1;
        node2 = l_head->g_id2;

        n_found1 = 0;
        n_found2 = 0;

        n_head = child->nhead;

		while(n_head)
		{
			if(n_head->g_node == node1) {
				n_found1 = 1;
			}
			if(n_head->g_node == node2) {
				n_found2 = 1;
			}
			n_head = n_head->next;
		}

		if(n_found1 == 1 && n_found2 == 1) {
			//Set up the loop id to be 1.
			root->loop_image[c_id][prog_id] = 1;
        }

        prog_id++;
        l_head = l_head->next;
	}//while
}

//Generate the shared stem and loop list for all the tree nodes in the tree. 
void DynamicP::tree_sharedstemloop(Tree_bag * root)
{
	int		i;

	//Initialize both arrays to be NULL.
	for(i=0; i<MAXNUMPT; i++) {
		root->stem_image[i] = NULL;
		root->loop_image[i] = NULL;
	}

	for(i=0; i<root->childnum; i++) {
		this->generate_sharedstemlist(root, root->pChild[i], i);
		this->generate_sharedlooplist(root, root->pChild[i], i);
	}

	//Call the function recursively to generate the list for all the tree nodes. 
	for(i=0; i<root->childnum; i++) {
		this->tree_sharedstemloop(root->pChild[i]);
	}
}

//Generate the mapping array between a node and the parent for all the stem nodes in the tree bag node.      
void DynamicP::construct_mappingarray(Tree_bag * parent, Tree_bag * child)
{
	int		i;
	int		c_node, p_node;
	int		c_loc, p_loc;
	int		f_flag, s_count;
	Node_list	*c_head;
	Node_list	*p_head;
	Strunit		*s_head;

	//Allocate sufficient memory space for the array state_array.
	//make state_array as a local variable
	int * state_array = new int[child->nodenum];
	child->set_up = new int[child->stnum];

	//Initialize the state_array to be all -1.
	for(i=0; i<child->nodenum; i++) {
		state_array[i] = -1;
	}

	//Now start looking at every node in the child and
	//check if it is also included in the tree bag node of parent. 
	c_head = child->nhead;
	c_loc = 0;

	while(c_head)
	{
		c_node = c_head->g_node;
        p_loc = 0;
        p_head = parent->nhead;
        f_flag = 0;

        while(p_head)
		{
			p_node = p_head->g_node;
			if(p_node == c_node) {
				f_flag = 1;
				break;
			}
			p_loc++;
			p_head = p_head->next;
		}

		if(f_flag == 1) {
			//If the found flag is found to be one.
			state_array[c_loc] = p_loc;
		}

		c_loc++;
		c_head = c_head->next;
	}//while
 
	//Now, go through the list of stems in the child and generate the array of set_up in the child node 
    //such that a mapping between the stems of the child node and the nodes in the parent node is set up.  
    s_head = child->su_stem;
    s_count = 0;
    while(s_head)
	{
		if(s_head->type == 1)
		{
			if(state_array[s_head->left_nid] != -1) {
				child->set_up[s_count] = state_array[s_head->left_nid];
			} else if(state_array[s_head->right_nid] != -1) {
				child->set_up[s_count] = state_array[s_head->right_nid];
			} else {
				child->set_up[s_count] = -1;
			}
		} else if(s_head->type == 2) {
			if(state_array[s_head->left_nid] != -1){
				child->set_up[s_count] = state_array[s_head->left_nid];
			} else {
				child->set_up[s_count] = -1;
			}
		}
		s_count++;
		s_head = s_head->next;
	}//while

	if(state_array != NULL) {
		delete [] state_array;	
		state_array = NULL;
	}
}

//Constructing the mapping array for all the nodes in the tree recursively. 
void DynamicP::tree_mappingarray(Tree_bag * root) 
{
    int		i;
    int		childnum;

    //If root is found to be NULL.
    if(root == NULL)
		return;

    childnum = root->childnum;
    
	for(i=0; i<childnum; i++)
		this->construct_mappingarray(root, root->pChild[i]);

    for(i=0; i<childnum; i++) {
        this->tree_mappingarray(root->pChild[i]);
    }
}

void DynamicP::print_treebag_set_up(Tree_bag *root)
{
	int i, j;
	if(root == NULL)
		return;

	this->printTreenode(root);
	cout<<endl;
	//print out its child's set_up
    for(i=0; i<root->childnum; i++) {
			cout<<"\t";
			this->printTreenode(root->pChild[i]);
//			cout<<endl;
		for(j=0; j<root->pChild[i]->stnum; j++) {
			cout<<"\t"<<root->pChild[i]->set_up[j]<<" ";
		}
		cout<<endl;
	}

	//recursively call itself
    for(i=0; i<root->childnum; i++) {
        this->print_treebag_set_up(root->pChild[i]);
    }
}

void DynamicP::preprocess_tree(Tree_bag *root)
{
	//1.Perform a conversation between the graph node id and the stem id represented in the secondary structure.  
	this->get_conversion_gtos();
//	this->print_gidtosid();

//	this->free_pStem();

	//3.Construct the list of stems and loops for all nodes in the tree.
//	cout<<"---- this->tree_stemlist(root) ----"<<endl;
	this->tree_stemlist(root);

//	cout<<"---- this->tree_looplist(root) ----"<<endl;
	this->tree_looplist(root);

//	this->printSuStemLoopInfo(root);
//	exit(0);

	this->free_pCongraph();
//	this->free_pLoop();

	//2.Allocate enough memory to perform dynamic programming and store intermediate results.
	this->allocate_dptable(root);

	//4.Construct the list of shared stems and loops for all nodes in the tree.
	this->tree_sharedstemloop(root);

	//5.Construct the arrays used for mapping purposes. 
	this->tree_mappingarray(root);

//	cout<<"print_treebag_set_up"<<endl;
//	this->print_treebag_set_up(root);

}

void DynamicP::free_gtos()
{
	if(gidtosid != NULL)
		delete [] gidtosid;
	gidtosid = NULL;
}
void DynamicP::free_dptable(Tree_bag *root)
{
	int i;
	//20080630
    for(i=0; i<root->childnum; i++)
		this->free_dptable(root->pChild[i]);

	//just for debug. 20080701
//	this->printTreenode(root);
//	cout<<endl;

	if(root->t_head != NULL) {
		delete [] root->t_head;			
		root->t_head = NULL;
	}
	if(root->enum_arrays != NULL) {
		delete [] root->enum_arrays;	
		root->enum_arrays = NULL;
	}
	if(root->enum_arrayo != NULL) {
		delete [] root->enum_arrayo;	
		root->enum_arrayo = NULL;
	}
}

Strunit * DynamicP::find_stemtail(Strunit *cur, int index)
{
	int i;
	Strunit * tail = cur;
	if(tail != NULL) {
		i = 0;
		while(i < index) {
			tail = tail->next;
			i++;
		}
	}
	return tail;
}

void DynamicP::free_stem(Tree_bag *root) 
{
	int i;
	Strunit * head, * tail;
	head = root->su_stem;

	for(i=root->stnum; i>=0; i--)
	{
		tail = this->find_stemtail(head, i);
		if(tail != NULL)
		{
//			cout<<"delete "<<tail->g_id1<<","<<tail->g_id2<<endl;
			if(tail->lstem != NULL) {
				tail->lstem = NULL;
			}
			delete tail;		tail = NULL;
		}
	}
	if(root->su_stem != NULL) {
		root->su_stem = NULL;
	}
}

void DynamicP::free_stemlist(Tree_bag *root)
{
	int		i;

	//starts from the child node
	for(i=0; i<root->childnum; i++)
		this->free_stemlist(root->pChild[i]);

//	this->printTreenode(root);
//	cout<<endl;

	this->free_stem(root);
}

Strunit * DynamicP::find_looptail(Strunit *cur, int index)
{
	int i;
	Strunit * tail = cur;
	if(tail != NULL) {
		i = 0;
		while(i < index) {
			tail = tail->next;
			i++;
		}
	}
	return tail;
}

void DynamicP::free_loop(Tree_bag *root) 
{
	int i;
	Strunit * head, * tail;
	head = root->su_loop;

	for(i=root->lpnum; i>=0; i--)
	{
		tail = this->find_looptail(head, i);
		if(tail != NULL)
		{
//			cout<<"delete "<<tail->g_id1<<","<<tail->g_id2<<endl;
			delete tail;
			tail = NULL;
		}
	}

	if(root->su_loop != NULL) {
		root->su_loop = NULL;
	}
}

void DynamicP::free_looplist(Tree_bag * root)
{
    int		i;

    for(i=0; i<root->childnum; i++)
		this->free_looplist(root->pChild[i]);

//	this->printTreenode(root);
//	cout<<endl;

	this->free_loop(root);
}

void DynamicP::free_sharedstemlist(Tree_bag * root, Tree_bag * child, int c_id)
{
	if(root->stem_image[c_id] != NULL) {
		delete [] root->stem_image[c_id];
		root->stem_image[c_id] = NULL;
	}
}

//Similarly, Generate the list of loops shared by a child and a root. 
//The loop_image array in the node root will takes the value of either 1 or -1, 
//when the value is 1, it means that the child also contains the loop.   
void DynamicP::free_sharedlooplist(Tree_bag * root, Tree_bag * child, int c_id)
{
	if(root->loop_image[c_id] != NULL) {
		delete [] root->loop_image[c_id];  
		root->loop_image[c_id] = NULL;
	}
}

void DynamicP::free_sharedstemloop(Tree_bag * root)
{
	int		i;

	//Call the function recursively to generate the list for all the tree nodes. 
	for(i=0; i<root->childnum; i++)
		this->free_sharedstemloop(root->pChild[i]);

//	this->printTreenode(root);
//	cout<<endl;

	for(i=0; i<root->childnum; i++) {
			this->free_sharedstemlist(root, root->pChild[i], i);
			this->free_sharedlooplist(root, root->pChild[i], i);
	}
}

//Generate the mapping array between a node and the parent for all the stem nodes in the tree bag node.      
void DynamicP::deconstruct_mappingarray(Tree_bag * parent, Tree_bag * child)
{
	if(child->set_up != NULL) {
		delete [] child->set_up;
		child->set_up = NULL;
	}
}

//Constructing the mapping array for all the nodes in the tree recursively. 
void DynamicP::free_mappingarray(Tree_bag * root) 
{
//	cout<<"---- free_mappingarray ----"<<endl;
    int		i;

    //If root is found to be NULL.
    if(root == NULL)
		return;

    for(i=0; i<root->childnum; i++) {
        this->free_mappingarray(root->pChild[i]);
    }

//	this->printTreenode(root);
//	cout<<endl;

	for(i=0; i<root->childnum; i++)
		this->deconstruct_mappingarray(root, root->pChild[i]);
}

//Do not forget to free some pointers in Tree_bag, e.g. * nhead, * t_head, * t_pointers[MAXNUMPT].
void DynamicP::postprocess_tree(Tree_bag *root)
{
	this->free_gtos();

	this->free_dptable(root);

	this->free_stemlist(root);	//Do not forget to delete su_stem/su_loop->lstem/rstem. 

	this->free_looplist(root);

	this->free_sharedstemloop(root);

	this->free_mappingarray(root);
}
//------------------------------------------------------
void DynamicP::print_treenode(Tree_bag * node)
{
    Node_list *nodehead;

    nodehead = node->nhead;
    cout<<"{";
    while(nodehead) {
		//Go through the list of nodes, print out the head information.
		cout<<nodehead->g_node<<" ";
		nodehead = nodehead->next;
    }
    cout<<"}"<<endl;
}

void DynamicP::print_treenodelist(Node_list * node_list)
{
    Node_list *nodehead;

    nodehead = node_list;
    cout<<"{";
    while(nodehead) {
		//Go through the list of nodes, print out the head information.
		cout<<nodehead->g_node<<" ";
		nodehead = nodehead->next;
    }
    cout<<"}"<<endl;
}

void DynamicP::print_tree(Tree_bag * root)
{
	cout<<"---- print_tree ----"<<endl;
    int		i;

    //If root is found to be NULL.
    if(root == NULL)
		return;

    for(i=0; i<root->childnum; i++) {
        this->print_tree(root->pChild[i]);
    }

	this->print_treenode(root);
}

void DynamicP::printSuStemLoopInfo(Tree_bag * root)
{
	int		i;
	Strunit *	s_head;
	Strunit *	l_head;

	for(i=0; i<root->childnum; i++) {
		this->printSuStemLoopInfo(root->pChild[i]);
	}

	cout<<endl<<"compnode:\t";
	this->print_treenode(root);

	//1. printing out the info of su_stem
	cout<<"Info of su_stem"<<endl;
	s_head	= root->su_stem;
	while(s_head)
	{
		cout<<"\tg_id1:"<<s_head->g_id1<<"\t\tg_id2:"<<s_head->g_id2<<endl;
		cout<<"\tleft_nid:"<<s_head->left_nid<<"\tright_nid:"<<s_head->right_nid<<endl;
		s_head = s_head->next;
	}

	//2. printing out the info of su_loop
	cout<<"Info of su_loop"<<endl;
	l_head = root->su_loop;
	while(l_head)
	{
		cout<<"\tg_id1:"<<l_head->g_id1<<"\t\tg_id2:"<<l_head->g_id2<<endl;
//		cout<<"\tleft_nid:"<<l_head->left_nid<<"\tright_nid:"<<l_head->right_nid<<endl;
		l_head = l_head->next;
	}
}

Node_list * DynamicP::find_treebag_nhead_tail(Node_list *cur, int index)
{
	int i;
	Node_list * tail = cur;
	if(tail != NULL) {
		i = 0;
		while(i < index) {
			tail = tail->next;
			i++;
		}
	}
	return tail;
}

void DynamicP::free_treebag_nhead(Node_list *nodehead, int nodenum) 
{
	int i;
	Node_list * head, * tail;
	head = nodehead;

//	cout<<"delete:\t{";
	for(i=nodenum; i>=0; i--)
	{
		tail = this->find_treebag_nhead_tail(head, i);
		if(tail != NULL)
		{
//			cout<<tail->g_node<<" ";
			delete tail;
			tail = NULL;
		}
	}
//	cout<<"}"<<endl;
	//20080701
	nodehead = NULL;
}

void DynamicP::free_treenode(Tree_bag * node)
{
	this->free_treebag_nhead(node->nhead, node->nodenum);

	//the treenode
	if(node != NULL) {
		delete [] node;
		node = NULL;
	}
}

void DynamicP::free_tree(Tree_bag * root)
{
    int		i;
    
	for(i=0; i<root->childnum; i++) {
        this->free_tree(root->pChild[i]);
    }

//	this->printTreenode(root);
//	cout<<endl;

	this->free_treenode(root);
}

void DynamicP::buildTreeNodePath(Tree_bag * root)
{
	//get the whole path of all tree node
	iNumOftheWholeTreeNode = this->countNumOftheWholeTreeNode(root);
	if(this->getFlagPrintDebugTreeInfo())
		cout<<"iNumOftheWholeTreeNode="<<iNumOftheWholeTreeNode<<endl;

	pTheWholeTreeNode = new int[iNumOftheWholeTreeNode];	//
	memset(pTheWholeTreeNode, 0, iNumOftheWholeTreeNode);
	this->getTheWholeTreeNode(root);

	if(this->getFlagPrintDebugTreeInfo()) {
		for(int j=0; j<iTheWholeTreeNodeIndex; j++) {
			cout<<pTheWholeTreeNode[j]<<"-";
		}
	}

	pOptimalTheWholeTreeNode = new int[iNumOftheWholeTreeNode];	//
}
void DynamicP::freeTreeNodePath()
{
	this->iOptimalTheWholeTreeNodeIndex = 0;
	if(pOptimalTheWholeTreeNode != NULL) {
		delete [] pOptimalTheWholeTreeNode;
		pOptimalTheWholeTreeNode = NULL;
	}

	this->iTheWholeTreeNodeIndex = 0;
	if(pTheWholeTreeNode != NULL) {
		delete [] pTheWholeTreeNode;
		pTheWholeTreeNode = NULL;
	}
}

void DynamicP::searchPK(int iWinSize, 
						int iMinSize, 
						Stem * pStemModel, 
						Loop * pLoopModel, 
						Tree_bag * root, 
						int shift_left, 
						char * genome_name, 
						int searchtype, 
//						DataLoader trainLoader, 
						char * train_file,		//
						char cOffsetMethod,		//
						int	 iSplitThreshold,	//
						int filterbegin, 
						int filterend, 
						double * genome_baseFreq)
{
	this->r_scores = new double[ACTUAL_CANDIDATE_NUM];
	this->cWindowSequence = new char[iWinSize+1];

	this->scan_genome(	iWinSize, 
						iMinSize, 
						pStemModel, 
						pLoopModel, 
						1, 
						root, 
						this->getGenomeSequenceLength(), 
						this->getGenomeSequence(), 
						shift_left, 
						genome_name, 
						searchtype, 
//						trainLoader, 
						train_file,			//
						cOffsetMethod,		//
						iSplitThreshold,	//
						filterbegin, 
						filterend, 
						genome_baseFreq,
						GENOME_SEARCH_DIRECTION_PLUS);
	this->freeGenomeSequence();

	if(cWindowSequence != NULL) {
		delete [] this->cWindowSequence;
		this->cWindowSequence = NULL;
	}
	if(r_scores != NULL) {
		delete [] this->r_scores;
		this->r_scores = NULL;
	}

	if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_BOTH)
	{
		this->r_scores = new double[ACTUAL_CANDIDATE_NUM];
		this->cWindowSequence = new char[iWinSize+1];

		this->scan_genome(	iWinSize, 
							iMinSize, 
							pStemModel, 
							pLoopModel, 
							1, 
							root, 
							this->getGenomeSequenceLength(), 
							this->getRGenomeSequence(), 
							shift_left, 
							genome_name, 
							searchtype, 
//							trainLoader, 
							train_file,			//
							cOffsetMethod,		//
							iSplitThreshold,	//
							filterbegin, 
							filterend, 
							genome_baseFreq,
							GENOME_SEARCH_DIRECTION_MINUS);
		this->freeRGenomeSequence();

		if(cWindowSequence != NULL) {
			delete [] this->cWindowSequence;
			this->cWindowSequence = NULL;
		}
		if(r_scores != NULL) {
			delete [] this->r_scores;
			this->r_scores = NULL;
		}
	}

}

double DynamicP::searchCurrentPos(	int iWinSize, Stem *pStemModel, 
									Loop *pLoopModel, int type, 
									Tree_bag * root, int curPos, int genomelength)
{
	int j, k;
	int	iSize = iWinSize;

	this->initDPTable(root);

	if((curPos+iSize) < genomelength)
	{
		for(j=curPos; j<curPos+iSize; j++)
			cWindowSequence[j-curPos] = scanSequence[j];
		cWindowSequence[iSize] = '\0';
	}
	else
	{
		//
		for(j=curPos; j<genomelength; j++)
			cWindowSequence[j-curPos] = scanSequence[j];
		cWindowSequence[genomelength-curPos] = '\0';
	}
//	cout<<cWindowSequence<<endl;

	int		numstem = this->getNumstems();
	int		numloop = this->getNumloops();
	double	profilescore;

	//--------------------------------
	//start computing the time for preprocessing
	start_preprocessing = clock();

	CandCollect	col(numstem, numloop);

	bool bTwoReg = true;
	col.setTwoRegionScan(bTwoReg);

	bool offsetvld = true;
	col.enableOffset( offsetvld );  //enable or disable offset

	double droprate = -1;
	col.setCYKMergeParamter(this->getFlagMergeCandInPreprocess(),
							this->getFlagCandwithShortestLength(),
							this->getShiftNumMergeCand(),
							droprate);

	double timesSD = 3.0;
	col.setCoeffStandardDev(timesSD);	//set the coefficient of standarad deviation

    double plowerbound=1.0;				//start to add penalty out of this bound
	col.setCYKLengthPenalty(plowerbound, this->getPcoeff());	//length penalty parameters

	//---------------------------------------------------------------------------
	if(this->getFlagPrintDebugDPWinInfo())
		cout<<"searchCandids function begins ..."<<endl;
	col.searchCandids(pStemModel, numstem, this->iLeading, INVLDPROB, cWindowSequence);
	if(this->getFlagPrintDebugDPWinInfo())
		cout<<"searchCandids function ends ..."<<endl;

	int iCandidateStemNum;
	this->candidateStems = new CandidateStem[numstem];

	//make sure every stem has candidates, if one of them has none, then move to the next pos
	bNoCandidateInCurrentWindow = 1;
	for(j=0; j<numstem; j++)
	{
		iCandidateStemNum = col.stemgrp[j].num;
		if(iCandidateStemNum > 0)
			bNoCandidateInCurrentWindow *= 1;
		else
			bNoCandidateInCurrentWindow *= 0;
	}

	if(bNoCandidateInCurrentWindow == 1)
	{
		for(j=0; j<numstem; j++)
		{
			iCandidateStemNum = col.stemgrp[j].num;
			candidateStems[j].num = iCandidateStemNum;
			candidateStems[j].s_locinfo = new Stemloc[iCandidateStemNum];	//

			//
			if(this->getFlagPrintStrAlignInfo()) {
				candidateStems[j].leftArm	= new StructureAlign[iCandidateStemNum];
				candidateStems[j].rightArm	= new StructureAlign[iCandidateStemNum];
			}

			for(k=0; k<iCandidateStemNum; k++)
			{
				candidateStems[j].s_locinfo[k].a = col.stemgrp[j].pMember[k].stempos[0][0];
				candidateStems[j].s_locinfo[k].b = col.stemgrp[j].pMember[k].stempos[0][1];
				candidateStems[j].s_locinfo[k].c = col.stemgrp[j].pMember[k].stempos[1][0];
				candidateStems[j].s_locinfo[k].d = col.stemgrp[j].pMember[k].stempos[1][1];
				candidateStems[j].s_locinfo[k].p_score = col.stemgrp[j].pMember[k].prob;

				if(this->getFlagPrintDebugInfo()) {
					cout<<"("<<candidateStems[j].s_locinfo[k].a<<"~"<<candidateStems[j].s_locinfo[k].b<<")"
					<<"-("<<candidateStems[j].s_locinfo[k].c<<"~"<<candidateStems[j].s_locinfo[k].d<<")"
					<<" \t"<<candidateStems[j].s_locinfo[k].p_score<<endl;
				}

				//
				if(this->getFlagPrintStrAlignInfo()) {
					//
					char	*cLeftSequence = NULL, *cRightSequence = NULL, 
							*cLeftStructureAlign = NULL, *cRightStructureAlign = NULL;
					char	*leftAlign = NULL, *rightAlign = NULL;
					//support stemid index in two pastalines 20081015
					char	*leftAlign_charid_idx = NULL, *rightAlign_charid_idx = NULL;
					int		iLeftLength=0, iRightLength=0;

					col.getAlignmentStr(j, k, 
										cLeftSequence, cLeftStructureAlign, iLeftLength, 
										cRightSequence, cRightStructureAlign, iRightLength);

					candidateStems[j].leftArm[k].sequence	= new char[iLeftLength + 1];
					candidateStems[j].leftArm[k].alignment	= new char[iLeftLength + 1];

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
						candidateStems[j].leftArm[k].pastaline_idx	= new char[iLeftLength + 1];

					candidateStems[j].rightArm[k].sequence	= new char[iRightLength + 1];
					candidateStems[j].rightArm[k].alignment	= new char[iRightLength + 1];

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
						candidateStems[j].rightArm[k].pastaline_idx	= new char[iRightLength + 1];

					//replace cLeftStructureAlign/cRightStructureAlign "[" and "]" with this->pStem[j].charid
					leftAlign	= strrpl(cLeftStructureAlign, '[', this->pStem[j].charid);
					rightAlign	= strrpl(cRightStructureAlign, ']', this->pStem[j].charid);
					strcpy(candidateStems[j].leftArm[k].sequence,	cLeftSequence);
					strcpy(candidateStems[j].leftArm[k].alignment,	leftAlign);
					strcpy(candidateStems[j].rightArm[k].sequence,	cRightSequence);
					strcpy(candidateStems[j].rightArm[k].alignment, rightAlign);

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
					{
						leftAlign_charid_idx	= strrpl(cLeftStructureAlign, '[', this->pStem[j].charid_idx);
						rightAlign_charid_idx	= strrpl(cRightStructureAlign, ']', this->pStem[j].charid_idx);
						strcpy(candidateStems[j].leftArm[k].pastaline_idx,	leftAlign_charid_idx);
						strcpy(candidateStems[j].rightArm[k].pastaline_idx, rightAlign_charid_idx);
					}

					if(cLeftSequence != NULL) {
						delete [] cLeftSequence;
						cLeftSequence = NULL;
					}
					if(cRightSequence != NULL) {
						delete [] cRightSequence;
						cRightSequence = NULL;
					}
					if(cLeftStructureAlign != NULL) {
						delete [] cLeftStructureAlign;
						cLeftStructureAlign = NULL;
					}
					if(cRightStructureAlign != NULL) {
						delete [] cRightStructureAlign;
						cRightStructureAlign = NULL;
					}
					if(leftAlign != NULL) {
						delete [] leftAlign;
						leftAlign = NULL;
					}
					if(rightAlign != NULL) {
						delete [] rightAlign;
						rightAlign = NULL;
					}
					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE)
					{
						if(leftAlign_charid_idx != NULL) {
							delete [] leftAlign_charid_idx;
							leftAlign_charid_idx = NULL;
						}
						if(rightAlign_charid_idx != NULL) {
							delete [] rightAlign_charid_idx;
							rightAlign_charid_idx = NULL;
						}
					}
				}//if(this->getFlagPrintStrAlignInfo())
			}//for k
		}//for j
	}
	//-------------------
	col.freeAllCandids();

	if(bNoCandidateInCurrentWindow == 1)
	{
		end_preprocessing = clock();
		double	cpu_time_preprocessing = ((double) (end_preprocessing - start_preprocessing)) / CLOCKS_PER_SEC;
		cpu_time_hours_preprocessing += cpu_time_preprocessing/3600;

		start_dp = clock();
		this->compute_profilescore(curPos, root, cWindowSequence, pLoopModel);
		profilescore = this->getMaxProfilescore(root);
		end_dp = clock();
		double	cpu_time_dp = ((double) (end_dp - start_dp)) / CLOCKS_PER_SEC;
		cpu_time_hours_dp += cpu_time_dp/3600;
	} else {
		profilescore = SMALLEST;
	}
	return profilescore;
}
void DynamicP::freeMemForCandidateStems()
{
	int		j, k;
	int		numstem = this->getNumstems();
	//start: free the memory for candidateStems
	for(j=0; j<numstem; j++)
	{
		//
		if(this->getFlagPrintStrAlignInfo()) 
		{
			if((this->candidateStems[j].leftArm	!= NULL)
				&& (this->candidateStems[j].rightArm != NULL))
			{
				for(k=0; k<candidateStems[j].num; k++) {
					delete [] candidateStems[j].leftArm[k].sequence;
					delete [] candidateStems[j].leftArm[k].alignment;
					delete [] candidateStems[j].rightArm[k].sequence;
					delete [] candidateStems[j].rightArm[k].alignment;
					candidateStems[j].leftArm[k].sequence	= NULL;
					candidateStems[j].leftArm[k].alignment	= NULL;
					candidateStems[j].rightArm[k].sequence	= NULL;
					candidateStems[j].rightArm[k].alignment	= NULL;
				}
				delete [] this->candidateStems[j].leftArm;
				delete [] this->candidateStems[j].rightArm;
				this->candidateStems[j].leftArm	 = NULL;
				this->candidateStems[j].rightArm = NULL;
			}
		}
		if(this->candidateStems[j].s_locinfo != NULL)
		{
			delete [] this->candidateStems[j].s_locinfo;
			this->candidateStems[j].s_locinfo = NULL;
		}
	}
	delete [] this->candidateStems;
	this->candidateStems = NULL;
	//end:  free the memory for candidateStems
}
void DynamicP::scan_genome(	int		iWinSize, 
							int		iMinSize, 
							Stem *	pStemModel, 
							Loop *	pLoopModel, 
							int		type, 
							Tree_bag * root, 
							int		genomelength, 
							char *	genome_sequence, 
							int		shift_left, 
							char *	genome_name, 
							int		search_type, 
//							DataLoader trainLoader,
							char *	train_file,		//
							char	cOffsetMethod,		//
							int		iSplitThreshold,	//
							int		filterbegin, 
							int		filterend,
							double * genome_baseFreq,
							int		search_direction)
{
	int		iSize = iWinSize;
	int		numstem = this->getNumstems();
//	int		numloop = this->getNumloops();
	double	profilescore;
	int		segstartpos, segendpos;	//the absolute value of the ending location of a sequence segment

	//1. First of all, decide to take the random sequence or the genome sequence.
	if(type == 0)
	{
		//Firstly, deal with random sequence then computing the threshold
	} else {
		//scan the genome sequence.
		this->scanLength = genomelength;
		this->scanSequence = new char[scanLength+1];
		strcpy(scanSequence, genome_sequence);
	}

	if(this->getFlagPrintDebugDPWinInfo()) {
		cout<<"scanLength="<<scanLength<<endl;
		cout<<"window size="<<iSize<<endl;
	}

	memset(cWindowSequence,0,iSize+1);

	//skip-and-jump strategy. 
//	int		scan_start = 0, scan_end = 0;	//start/end position of scanning
	int		stepsize = this->getStepSize();	//incremented stepsize

	int		start_pos, cur_start_pos, next_start_pos;
	int		pre_end_pos, cur_end_pos;
	double	optimal_score;
	double	next_optimal_score;
	int		optimal_pos = -1;//, optimal_pos_start, optimal_pos_end;
	int		cur_optimal_pos;

	int		curPos = 0;	//start position
	int		endPos = 0;

	if(search_type == SEARCH_SUBSTRUCTUREFILTER)
		endPos = scanLength - iSize + 1;
	else	//if(search_type == SEARCH_WHOLESTRUCTURE)
		endPos = scanLength - iMinSize + 1;
	
	if(this->getFlagPrintDebugDPWinInfo())
		cout<<"endPos="<<endPos<<endl;

	//compute the starting time
	double	cpu_time_all = 0.0;

	start_all = clock();
	cpu_time_hours_all = 0.0;
	cpu_time_hours_preprocessing = 0.0;
	cpu_time_hours_dp = 0.0;

	optimal_score = SMALLEST;
	pre_end_pos = cur_end_pos = -1;
	start_pos = cur_start_pos = next_start_pos = -1;
	int lastStemId = (stemIdxArray[2*numstem-1]-1)/2;

	int count = 0;

	int extension_left	= 0;
	int extension_right = 0;
	int	hitlength		= 0;
	char * c_subgenome_seq = NULL;

	int hit_pos_start	= 0;
	int hit_pos_end		= 0;
	int extension_pos_start = 0;
	int extension_pos_end	= 0;

//	//just for the test
//	curPos = 79;
//	endPos = 80;

	DataLoader trainLoader;
	trainLoader.setOffsetMethod(cOffsetMethod);
	trainLoader.setSplitThreshold(iSplitThreshold);
	trainLoader.inputData(train_file);

	do
	{
		profilescore = searchCurrentPos(iWinSize, pStemModel, pLoopModel, type, root, curPos, genomelength);
		
		if(this->getFlagPrintDebugDPWinInfo()) {
			cout<<curPos<<"\t";
			cout<<profilescore<<endl;
		}

		if(type == 0) {
			//do nothing
		} else {
			if(profilescore > this->getThreshold()) 
			{
				memset(pOptimalTheWholeTreeNode, 0, iNumOftheWholeTreeNode);
				iOptimalTheWholeTreeNodeIndex = 0;

				this->getOptimalPath(root);

				this->buildStemCandIdxArray();
				this->locateStemStartEndPosIdx();
				segstartpos = candidateStems[0].s_locinfo[iStemStartPosIdx].a;
				segendpos	= candidateStems[lastStemId].s_locinfo[iStemEndPosIdx].d;

				if(this->getJumpStrategy() != 0)
				{
					//do some skip-and-jump strategy
					//after hitting the candidate then use stepsize = 2 
					stepsize = this->getStepSize();

					cur_end_pos =	curPos + segendpos;
					cur_start_pos = curPos + segstartpos;//
					if(optimal_score == SMALLEST)
					{
						//the first time we meet the candidate hit.
						if(optimal_score < profilescore) {
							optimal_score = profilescore;
							optimal_pos = curPos;
						}
					}
					else
					{
						if(abs(cur_end_pos - pre_end_pos) <= CANDIDATE_POSITION_THRESHOLD)
						{
							//belongs to the same candidate hit
							if(optimal_score < profilescore) {
								optimal_score = profilescore;
								optimal_pos = curPos;
							}
							//always keep start_pos as the smallest pos for the current candidate hit
							if(cur_start_pos < start_pos)
								start_pos = cur_start_pos;
						}
						else	//other candidate hit is coming
						{
							next_optimal_score = profilescore;
							next_start_pos = cur_start_pos;

							//output the pos with the highest score, that is the optimal_pos
							if(optimal_score > this->getScoreThresInJump())//CANDIDATE_SCORE_THRESHOLD)
							{
								profilescore = searchCurrentPos(iWinSize, pStemModel, pLoopModel, type, root, optimal_pos, genomelength);
								if(profilescore > threshold)	//we know it should be bigger than zero
								{
									memset(pOptimalTheWholeTreeNode, 0, iNumOftheWholeTreeNode);
									iOptimalTheWholeTreeNodeIndex = 0;

									this->getOptimalPath(root);

									this->buildStemCandIdxArray();

									this->locateStemStartEndPosIdx();
									segstartpos = candidateStems[0].s_locinfo[iStemStartPosIdx].a;
									segendpos	= candidateStems[lastStemId].s_locinfo[iStemEndPosIdx].d;

									r_scores[numhit++] = profilescore;

									if(start_pos == -1)
										start_pos = optimal_pos;

									if(search_type == SEARCH_WHOLESTRUCTURE)
									{
//										cout<<endl<<"---------------"<<endl;
//										cout<<"Hit "<<numhit+this->getPreHitNum()<<endl;
										printWholeSearchHitIndex(numhit+this->getPreHitNum());
										cout<<genome_name<<endl;

										if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_PLUS)
											cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
										else if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_MINUS)
											cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
										else if(search_direction == GENOME_SEARCH_DIRECTION_PLUS)
											cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
										else if(search_direction == GENOME_SEARCH_DIRECTION_MINUS)
											cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;

//										cout<<"hit positions = ("<<(optimal_pos+segstartpos+shift_left)<<"-"<<(optimal_pos+segendpos+shift_left)<<")"<<endl;
//										cout<<"hit score = "<<profilescore<<endl;

										printHitPos(optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);
										printHitScore(profilescore);

//										if(this->getFlagPrintStrAlignInfo()) {
										this->buildFoldedStructureInfo();
										this->printFoldedStructureInfo(segstartpos, segendpos, optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);

										this->buildStructureAlignInfo(pLoopModel);
										this->printStructureAlignInfo(optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);

//										this->printStemCandIdxArray();	//just for the debug 20081021

										if(this->getFlagPrintScoreInfo()) {
											this->printScoreForEveryComponent(pLoopModel, false);
											cout<<"-------------------"<<endl;
											this->printCandidateStemInfo();
										}
									}
									else	//sub-structure search
									{
										//20080911
										printFilterHitIndex(numhit+this->getPreHitNum());

										hitlength = segendpos - segstartpos + 1;
										trainLoader.getBothExtLength(filterbegin, filterend, hitlength, extension_left, extension_right);

										hit_pos_start	= optimal_pos + segstartpos;
										hit_pos_end		= optimal_pos + segendpos;
										if(hit_pos_start > extension_left) {
											extension_pos_start = hit_pos_start - extension_left;
										} else {
											extension_pos_start = 0;
										}
										extension_pos_end	= hit_pos_end + extension_right;
										cout<<">"<<genome_name<<"("<<hit_pos_start<<"-"<<hit_pos_end<<")"<<"["<<extension_pos_start<<"-"<<extension_pos_end<<"]"<<endl;

										if(search_direction == GENOME_SEARCH_DIRECTION_PLUS)
											cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
										else
											cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;

										printHitPos(hit_pos_start, hit_pos_end);
										printHitScore(profilescore);

										//alignment info printed in sub-structure search
//										if(this->getFlagPrintStrAlignInfo()) {
//											this->buildFoldedStructureInfo();
//											this->printFoldedStructureInfo(true);
										this->buildStructureAlignInfo(pLoopModel);
										printFilterHitAlignment(structureSequence, strucutreAlignment);
//										if(this->getFlagPrintScoreInfo())
//											this->printScoreForEveryComponent(pLoopModel, true);
//										}

										printFilterHitExtenPos(extension_pos_start, extension_pos_end);
										printFilterHitExtenNtsHeader();
										c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
										printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE);
										if(c_subgenome_seq != NULL) {
											delete [] c_subgenome_seq;
											c_subgenome_seq = NULL;
										}
//										printFilterHitExtenNtsTail();
									}
								}
								freeMemForCandidateStems();
							}//if(optimal_score > score_threshold)
							optimal_score	= next_optimal_score;
							start_pos		= next_start_pos;
							optimal_pos		= curPos;
						}//else	//other candidate hit is coming
					}
				} else {
					//no jump strategy
					cur_optimal_pos = curPos;
					r_scores[numhit++] = profilescore;

					if(search_type == SEARCH_WHOLESTRUCTURE)
					{
//						cout<<endl<<"---------------"<<endl;
//						cout<<"Hit "<<numhit+this->getPreHitNum()<<endl;
						printWholeSearchHitIndex(numhit+this->getPreHitNum());
						cout<<genome_name<<endl;

						if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_PLUS)
							cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
						else if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_MINUS)
							cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
						else if(search_direction == GENOME_SEARCH_DIRECTION_PLUS)
							cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
						else if(search_direction == GENOME_SEARCH_DIRECTION_MINUS)
							cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;

//						cout<<"hit positions = ("<<(optimal_pos+segstartpos+shift_left)<<"-"<<(optimal_pos+segendpos+shift_left)<<")"<<endl;
//						cout<<"hit score = "<<profilescore<<endl;
						printHitPos(optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);
						printHitScore(profilescore);

//						if(this->getFlagPrintStrAlignInfo()) {
						this->buildFoldedStructureInfo();
						this->printFoldedStructureInfo(segstartpos, segendpos, optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);

						this->buildStructureAlignInfo(pLoopModel);
						this->printStructureAlignInfo(optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);

//						this->printStemCandIdxArray();	//just for the debug 20081021

						if(this->getFlagPrintScoreInfo()) {
							this->printScoreForEveryComponent(pLoopModel, false);
							cout<<"-------------------"<<endl;
							this->printCandidateStemInfo();
						}
					}
					else	//sub-structure search
					{
						//20080911
						printFilterHitIndex(numhit+this->getPreHitNum());

						hitlength = segendpos - segstartpos + 1;
						trainLoader.getBothExtLength(filterbegin, filterend, hitlength, extension_left, extension_right);

						hit_pos_start	= optimal_pos + segstartpos;
						hit_pos_end		= optimal_pos + segendpos;
						if(hit_pos_start > extension_left) {
							extension_pos_start = hit_pos_start - extension_left;
						} else {
							extension_pos_start = 0;
						}
						extension_pos_end	= hit_pos_end + extension_right;
						cout<<">"<<genome_name<<"("<<hit_pos_start<<"-"<<hit_pos_end<<")"<<"["<<extension_pos_start<<"-"<<extension_pos_end<<"]"<<endl;

						if(search_direction == GENOME_SEARCH_DIRECTION_PLUS)
							cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
						else
							cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;

						printHitPos(hit_pos_start, hit_pos_end);
						printHitScore(profilescore);

//						if(this->getFlagPrintStrAlignInfo()) {
//							this->printHeader(segstartpos, segendpos, true);
//							this->buildFoldedStructureInfo();
//							this->printFoldedStructureInfo(true);
						this->buildStructureAlignInfo(pLoopModel);
						printFilterHitAlignment(structureSequence, strucutreAlignment);
//						if(this->getFlagPrintScoreInfo())
//							this->printScoreForEveryComponent(pLoopModel, true);
//						}

						printFilterHitExtenPos(extension_pos_start, extension_pos_end);
						printFilterHitExtenNtsHeader();
						c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
						printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE);
						if(c_subgenome_seq != NULL) {
							delete [] c_subgenome_seq;
							c_subgenome_seq = NULL;
						}
//						printFilterHitExtenNtsTail();
					}
					
					stepsize = 1;
					freeMemForCandidateStems();
				}

				curPos += stepsize;
			}
			else	//if not (profilescore > threshold)
			{
				//in current window candidates are found, then we need to release the memory
				if(bNoCandidateInCurrentWindow == 1)
					freeMemForCandidateStems();	//
				if(this->getJumpStrategy() != 0)
				{
					stepsize = this->getStepSize();
				} else {
					stepsize = 1;
				}
				curPos += stepsize;
			}
		}
		//
		if(this->getJumpStrategy() != 0)
			pre_end_pos = cur_end_pos;

		count++;
		//compute the ending time
		if((count % 100) == 0)
		{
			end_all = clock();
			cpu_time_all = ((double) (end_all - start_all)) / CLOCKS_PER_SEC;
			cpu_time_hours_all += cpu_time_all/3600;
			start_all = clock();
			count = 0;
		}
//	} while(i < scanLength-iSize+1);
	} while(curPos < endPos);

	//to the last position
	if(this->getJumpStrategy() != 0)
	{
		//save the last one 
//		if(optimal_pos>=0 && optimal_pos<scanLength-iSize+1)
		if(optimal_pos >= 0 && optimal_pos < endPos)
		{
			profilescore = searchCurrentPos(iWinSize, pStemModel, pLoopModel, type, root, optimal_pos, genomelength);
			if(profilescore > threshold)
			{
				memset(pOptimalTheWholeTreeNode, 0, iNumOftheWholeTreeNode);
				iOptimalTheWholeTreeNodeIndex = 0;

				this->getOptimalPath(root);

				this->buildStemCandIdxArray();
				this->locateStemStartEndPosIdx();
				segstartpos = candidateStems[0].s_locinfo[iStemStartPosIdx].a;
				segendpos	= candidateStems[lastStemId].s_locinfo[iStemEndPosIdx].d;

				r_scores[numhit++] = profilescore;

				if(start_pos == -1)
					start_pos = optimal_pos;

				if(search_type == SEARCH_WHOLESTRUCTURE)
				{
//					cout<<endl<<"---------------"<<endl;
//					cout<<"Hit "<<numhit+this->getPreHitNum()<<endl;
					printWholeSearchHitIndex(numhit+this->getPreHitNum());
					cout<<genome_name<<endl;

					if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_PLUS)
						cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
					else if(this->getFlagSearchReverseGenome() == GENOME_SEARCH_DIRECTION_MINUS)
						cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;
					else if(search_direction == GENOME_SEARCH_DIRECTION_PLUS)
						cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
					else if(search_direction == GENOME_SEARCH_DIRECTION_MINUS)
						cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;

//					cout<<"hit positions = ("<<(optimal_pos+segstartpos+shift_left)<<"-"<<(optimal_pos+segendpos+shift_left)<<")"<<endl;
//					cout<<"hit score = "<<profilescore<<endl;
					printHitPos(optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);
					printHitScore(profilescore);

//					if(this->getFlagPrintStrAlignInfo()) {
					this->buildFoldedStructureInfo();
					this->printFoldedStructureInfo(segstartpos, segendpos, optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);

					this->buildStructureAlignInfo(pLoopModel);
					this->printStructureAlignInfo(optimal_pos+segstartpos+shift_left, optimal_pos+segendpos+shift_left);

//					this->printStemCandIdxArray();	//just for the debug 20081021

					if(this->getFlagPrintScoreInfo()) {
						this->printScoreForEveryComponent(pLoopModel, false);
						cout<<"-------------------"<<endl;
						this->printCandidateStemInfo();
					}
				}
				else	//sub-structure search
				{
					//20080911
					printFilterHitIndex(numhit+this->getPreHitNum());

					hitlength = segendpos - segstartpos + 1;
					trainLoader.getBothExtLength(filterbegin, filterend, hitlength, extension_left, extension_right);

					hit_pos_start	= optimal_pos + segstartpos;
					hit_pos_end		= optimal_pos + segendpos;
					if(hit_pos_start > extension_left) {
						extension_pos_start = hit_pos_start - extension_left;
					} else {
						extension_pos_start = 0;
					}
					extension_pos_end	= hit_pos_end + extension_right;
					cout<<">"<<genome_name<<"("<<hit_pos_start<<"-"<<hit_pos_end<<")"<<"["<<extension_pos_start<<"-"<<extension_pos_end<<"]"<<endl;

					if(search_direction == GENOME_SEARCH_DIRECTION_PLUS)
						cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_PLUS<<endl;
					else
						cout<<GENOME_TAG_SEARCH_RESULT_DIRECTION_MINUS<<endl;

					printHitPos(hit_pos_start, hit_pos_end);
					printHitScore(profilescore);

//					if(this->getFlagPrintStrAlignInfo()) {
//						this->buildFoldedStructureInfo();
//						this->printFoldedStructureInfo(true);
					this->buildStructureAlignInfo(pLoopModel);
					printFilterHitAlignment(structureSequence, strucutreAlignment);
//					if(this->getFlagPrintScoreInfo())
//						this->printScoreForEveryComponent(pLoopModel, true);
//					}

					printFilterHitExtenPos(extension_pos_start, extension_pos_end);
					printFilterHitExtenNtsHeader();
					c_subgenome_seq = strextract(genome_sequence, extension_pos_start, extension_pos_end);
					printGenomeSegments(c_subgenome_seq, NUM_NTS_PER_LINE);
					if(c_subgenome_seq != NULL) {
						delete [] c_subgenome_seq;
						c_subgenome_seq = NULL;
					}
//					printFilterHitExtenNtsTail();
				}
			}
			freeMemForCandidateStems();
		}//if(cur_optimal_pos>=0)
	}

	end_all = clock();
	cpu_time_all = ((double) (end_all - start_all)) / CLOCKS_PER_SEC;
	cpu_time_hours_all += cpu_time_all/3600;

	if(this->scanSequence != NULL) {
		delete [] this->scanSequence;
		this->scanSequence = NULL;
	}
}

void DynamicP::compute_profilescore(int j, Tree_bag * tr_node, char * cSequence, Loop *pLoopModel)
{
	int		i;
	int		sum;
	int		iValue, iIndex;
	double	score, child_score;
	double	l_score;
//	int		max_end;	//A variable that records the maximum end of a tree_bag.
//	int		min_start;	//

	int		l_id;
	int		l_end;//, l_start;
	int		r_start;//, r_end;
	int		loop_end;
	int		loop_len;

	int		iTmp;
	iTmp =	strlen(cSequence);

	Strunit *	s_head = NULL;
	Strunit *	l_head = NULL;
	Strunit *	tail_head = NULL;

	for(i=0; i<tr_node->childnum; i++) {
		this->compute_profilescore(j, tr_node->pChild[i], cSequence, pLoopModel);
	}

	if(this->getFlagPrintDebugDPSearchInfo()) {
		cout<<"compnode:\t";
		this->print_treenode(tr_node);
		cout<<endl;
	}

	//enum_arrays:	all possible combinations of structural units(stem).
	for(i=0; i<tr_node->stnum; i++) {
		tr_node->enum_arrays[i] = 0;
	}

	//enum_arrayo:	all possible combinations of nodes.
	for(i=0; i<tr_node->nodenum; i++) {
		tr_node->enum_arrayo[i] = 0;
	}

	VTBSearch vs(1);
	vs.setAllowedInsNum(this->getAllowedNullLoopInsNum());

	//Compute the maximum value of max_end that can be reached for this particular location on the genome.
	while(1)
	{
		for(i=0; i<tr_node->childnum; i++) {
			tr_node->prob_score[i] = 0.0;
		}

//		max_end = 0;
//		min_start = iTmp;

		//1. dealing with stems.
		//The list of structural units which include both stems and half stems.
		s_head	= tr_node->su_stem;	
		iIndex	= 0;
		score	= INITSCORE;//0.0;

		while(s_head)
		{
			if(s_head->type == 1)	//type=1, a stem.
			{
				//If the stem unit contains a complete stem.
				iValue			= tr_node->enum_arrays[iIndex];
				score			+= this->candidateStems[s_head->stem_id].s_locinfo[iValue].p_score;

				//Compute the real locations of the stems on the genome sequence.
				s_head->start1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].a;
				s_head->end1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].b;
				s_head->start2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].c;
				s_head->end2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].d;

//				if(min_start > s_head->start1)
//					min_start = s_head->start1;
//				if(max_end < s_head->end2)
//					max_end = s_head->end2;

				//consider the score of parts that is shared with children.
				//save them and substract them when possible.
				for(i=0; i<tr_node->childnum; i++)
				{
					if(tr_node->pChild[i] != NULL && tr_node->stem_image[i][iIndex] == 1)
					{
						//If it is included in the child node.  
						tr_node->prob_score[i] += this->candidateStems[s_head->stem_id].s_locinfo[iValue].p_score;               
					}
				}

				tr_node->enum_arrayo[s_head->left_nid]	= iValue;
				tr_node->enum_arrayo[s_head->right_nid] = iValue;
			}
			else
			{
				//Otherwise, only part of the stem is included.
				iValue	= tr_node->enum_arrays[iIndex];

				if(gidtosid[s_head->g_id1].left_or_right == 0)
				{
					//if it is the left part of the stem.
					s_head->start1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].a;
					s_head->end1	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].b;

//					if(min_start > s_head->start1)
//						min_start = s_head->start1;
//					if(max_end < s_head->end1)
//						max_end = s_head->end1;

				} else {
					//Otherwise it is the right part of the stem.
					s_head->start2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].c;
					s_head->end2	= this->candidateStems[s_head->stem_id].s_locinfo[iValue].d;

//					if(min_start > s_head->start2)
//						min_start = s_head->start2;
//					if(max_end < s_head->end2)
//						max_end = s_head->end2;
				}

				tr_node->enum_arrayo[s_head->left_nid] = iValue;
			}
			iIndex++; 
			s_head = s_head->next;
		}
//		cout<<"Score after considering stems = "<<score<<endl;
//		this->print_enum_arrayo(tr_node->enum_arrayo, tr_node->nodenum);
//		cout<<endl;

		//2. dealing with loops.
		l_head	= tr_node->su_loop;
		l_id	= 0;
		int		l_source, l_sink;

		while(l_head)
		{
			//
			l_source = l_sink = 0;

			if(node_mapping[l_head->g_id1].node_type != 0)
			{
				//If the loop doesn't start with the source. 
				if(gidtosid[l_head->g_id1].left_or_right == 0) {
//					l_start = (l_head->lstem)->start1 + 1;
					l_end	= (l_head->lstem)->end1 + 1;
				} else {
//					l_start = (l_head->lstem)->start2 + 1;
					l_end	= (l_head->lstem)->end2 + 1;
				}

				if(node_mapping[l_head->g_id2].node_type != 0) {
					//If the loop doesn't end with the sink. 
					if(gidtosid[l_head->g_id2].left_or_right == 0) {
						r_start = (l_head->rstem)->start1 - 1;
//						r_end	= (l_head->rstem)->end1 - 1;
					} else {
						r_start = (l_head->rstem)->start2 - 1;
//						r_end	= (l_head->rstem)->end2 - 1;
					}
				} 
				else 
				{
					l_sink = 1;
				}
				if(l_sink == 0) {
                    //modify this part because overlap in parts of stems should be allowed
					if(l_end <= r_start + 1 + this->getNumOfNtsOverlapBetweenStem())
					{
						//If there is no overlap occurs, simply take the profile score of the loop.
						loop_end = r_start;
						loop_len = (r_start - l_end) + 1;

						l_score = vs.VTBSearchMax(&pLoopModel[l_head->loop_id], cSequence, l_end, loop_end, NULL, NULL);
						loop_end += l_end;	//make it the absolute coordinate
						score += l_score;

						//consider the score of parts that is shared with children.
						//save them and substract them when possible.
						for(i=0; i<tr_node->childnum; i++) {
							if(tr_node->pChild[i] != NULL && tr_node->loop_image[i][l_id] == 1) {
								tr_node->prob_score[i] += l_score;
							}
						}
					}//end if(l_end <= r_start + 1)
					else
					{
						score = SMALLEST;
					}
				}//end if(l_sink == 0)
				else
				{
					//reaching the sink part. --->End
					//do nothing, just skip it
				}
			}
			else
			{
				//reaching the source part. Start(0)--->
				l_source = 1;
			}
			l_id++;
			l_head = l_head->next;
		}
//		cout<<"Score after considering loops = "<<score<<endl;

		//3. Now, start compute the probability score of the node;
		if(score > MYTHRESHOLD) {
//			int j = 0;

			for(i=0; i<tr_node->childnum; i++) {
				//Extract the information needed to go bottom up from the child node to a parent.
//				this->print_enum_arrayo(tr_node->enum_arrayo, tr_node->nodenum);
//				cout<<endl;
				child_score = this->extract_profilescore(tr_node->pChild[i], tr_node->enum_arrayo, tr_node->nodenum);
				score		+= child_score - tr_node->prob_score[i];
			}//for i

//			this->print_enum_arrayo(tr_node->enum_arrayo, tr_node->nodenum);
//			cout<<"\t\t\t\tscore="<<score<<endl;
			this->tb_indexingstore(tr_node, tr_node->enum_arrayo, tr_node->nodenum, score);
		}//if

		int real_candidate_num = 0;
		int k;
		for(i=tr_node->stnum-1; i>=0; i--) {
			tr_node->enum_arrays[i]++;
			tail_head = tr_node->su_stem;
			k = i;
			while(k>0) {
				tail_head = tail_head->next;
				k--;
			}
			if(tail_head != NULL)
				real_candidate_num = this->candidateStems[tail_head->stem_id].num;
			if(tr_node->enum_arrays[i] < real_candidate_num) {
				break;
			} else {
				tr_node->enum_arrays[i] = 0;
			}
		}

		if(this->getFlagPrintDebugDPSearchInfo()) {
			printf("Check tr_node->enum_arrays:\t");
			for(i=0; i<tr_node->stnum; i++) {
				printf("%d ", tr_node->enum_arrays[i]);
			}
			printf("\n");
		}

		sum = 0;

		for(i=0; i<tr_node->stnum; i++) {
			sum += tr_node->enum_arrays[i];
		}

		if(sum == 0) {
//			this->printTreeNodeInfo(tr_node);
			//do this after enumerating the whole combination of this tree bag.
//			this->setOptimalDPtable(tr_node);
//			this->printOptimalDPtable(tr_node);
			break;
		}
	}
}

//Extract the information needed to go bottom up from the child node to a parent. 
//double TreeDecomposition::extract_profilescore(Tree_bag *child, int *enumerate_array, int num, int * s_pos, int * e_pos)
double DynamicP::extract_profilescore(Tree_bag *child, int *enumerate_array, int num)
{
	int		i;
	int		sum, s_count;
	double	maxscore, score;
	Strunit		*s_head = NULL;
	Node_list	*n_list = NULL;
	
	Strunit *	tail_head = NULL;

	maxscore = SMALLEST;

	if(this->getFlagPrintDebugDPSearchInfo())
		cout<<"\t\t-------- extract_profilescore -------"<<endl;

	//Now we go through the set up list to enumerate all possible combinations.
	for(i=0; i<child->stnum; i++)
    {
		if(child->set_up[i] == -1) {
			//If it is not preset by the parent, set it to be zero. 
			child->enum_arrays[i] = 0;
		} else {
			child->enum_arrays[i] = enumerate_array[child->set_up[i]];
		}
	}
//	//for the purpose of debugging
//	cout<<"\t\t";
//	this->print_enum_arrays(child);

	//Initialize the aray of enum_arrayo to be zero. 
	for(i=0; i<child->nodenum; i++) {
		child->enum_arrayo[i] = 0;
	}

	//Now start doing the enumeration.
	while(1)
	{
		//Now based on the generated array, start doing the look up and conversion; 
        s_head = child->su_stem;
        s_count = 0;

        while(s_head) {
			if(s_head->type == 1) {
				//If it is a complete loop.
				child->enum_arrayo[s_head->left_nid]	= child->enum_arrays[s_count];
				child->enum_arrayo[s_head->right_nid]	= child->enum_arrays[s_count];
			} else {
				child->enum_arrayo[s_head->left_nid]	= child->enum_arrays[s_count];
			}
			s_count++;
			s_head = s_head->next;
		}

		n_list = child->nhead;

		if(this->getFlagPrintDebugDPSearchInfo()) {
			printf("\t\tChild node includes:\t");
			this->printTreenode(child);
			cout<<"\t\t";
			this->print_enum_arrayo(child->enum_arrayo, child->nodenum);
			cout<<endl;
		}

		this->tb_indexingfetch(child, child->enum_arrayo, child->nodenum, &score);

		if(this->getFlagPrintDebugDPSearchInfo()) {
			cout<<"\t\tScore of this child node is "<<score<<endl;
		}

        if(maxscore < score) {
			maxscore = score;
		}
		int real_candidate_num = 0;
		int k=0;
		for(i=child->stnum-1; i>=0; i--) {
			if(child->set_up[i] == -1) {
				child->enum_arrays[i]++;
				tail_head = child->su_stem;
				k = i;
				while(k>0) {
					tail_head = tail_head->next;
					k--;
				}
				if(tail_head != NULL)
					real_candidate_num = this->candidateStems[tail_head->stem_id].num;
				if(child->enum_arrays[i] < real_candidate_num) {
					break;
				} else {
					child->enum_arrays[i] = 0;
				}
			}
		}

		sum = 0;
		for(i=0; i<child->stnum; i++) {
			if(child->set_up[i] == -1) {
				sum += child->enum_arrays[i];
			}
		}
		if(sum == 0) break;
	}//while

	return maxscore;
}

//after enumerating the whole combination of the tree bag,
//we need to check which combination is the optimal.
void DynamicP::setOptimalDPtable(Tree_bag *root)
{
	int		i, j, index = 0;
	int		nodenum, iShared, iTotalShared, iNotShared, iTotalNotShared;
	nodenum		= root->nodenum;
	iShared		= root->shared_node_num;
	iNotShared	= nodenum - iShared;

	double	optimal_value;
	double	dTmp;
	int		optimalIndex	= -1;

//	iTotalShared = floor(pow(double(this->iLeading), double(iShared)));
//	iTotalNotShared = floor(pow(double(this->iLeading), double(iNotShared)));
	iTotalShared = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(iShared)));
	iTotalNotShared = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(iNotShared)));
	for(i=0; i<iTotalShared; i++) {
		//for every shared node
		optimal_value	= SMALLEST;
		optimalIndex	= -1;
		for(j=0; j<iTotalNotShared; j++) {
			index	= i*iTotalNotShared + j;
//			cout<<index<<endl;
			dTmp	= root->t_head[index].score;
			if(optimal_value < dTmp) {
				optimal_value	= dTmp;
				optimalIndex	= index;
			}
		}
		if(optimalIndex > -1)
			root->t_head[optimalIndex].optimal = 1;
	}
}
void DynamicP::printOptimalDPtable(Tree_bag *root)
{
	this->printTreenode(root);
	cout<<endl;

	int		i, j, index = 0;
	int		nodenum, iShared, iTotalShared, iNotShared, iTotalNotShared;
	nodenum		= root->nodenum;
	iShared		= root->shared_node_num;
	iNotShared	= nodenum - iShared;

//	double	optimal_value	= SMALLEST;
//	int		optimalIndex	= -1;

//	char *	cKNary	= new char[nodenum];
//	memset(cKNary, 0, nodenum);

//	iTotalShared = floor(pow(double(this->iLeading), double(iShared)));
	iTotalShared = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(iShared)));
	for(i=0; i<iTotalShared; i++) {
		//for every shared node
//		iTotalNotShared = floor(pow(double(this->iLeading), double(iNotShared)));
		iTotalNotShared = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(iNotShared)));
		for(j=0; j<iTotalNotShared; j++) {
			index	= i*iTotalNotShared + j;
			if(root->t_head[index].optimal == 1) {
//				itoa(index, cKNary, this->iLeading);
//				cout<<"\t"<<cKNary;
				cout<<"\t"<<index;
				cout<<"("<<index<<")"<<"\t"<<root->t_head[index].score<<endl;
			}
		}
	}
//	if(cKNary != NULL) {
//		delete [] cKNary;
//		cKNary = NULL;
//	}
}

double DynamicP::getMaxProfilescore(Tree_bag *root)
{
	double	dOptimalScore = SMALLEST;
	dOptimalScore = root->t_head[0].score;
	return dOptimalScore;
}

void DynamicP::tb_indexingfetch(Tree_bag * tr_node, int * index_array, int num_index, double * score)
{
	int		i, index;
	int shared_num = tr_node->shared_node_num;

	//calculate the index
	index = 0;
	for(i=0; i<shared_num; i++) {
//		index += tr_node->enum_arrayo[i]*floor(pow(double(this->iLeading), double(shared_num-1-i)));
		index += tr_node->enum_arrayo[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(shared_num-1-i)));
	}

//	//calculate the index
//	for(i=0; i<num_index; i++) {
//		index += index_array[i]*floor(pow(double(this->iLeading), double(num_index-1-i)));
//	}

	if(tr_node->t_head[index].optimal == 1) {
		*score	= tr_node->t_head[index].score;
	} else {
//		cout<<"something wrong when setting optimal..."<<endl;
		*score	= SMALLEST;
	}
}

void DynamicP::tb_indexingstore(Tree_bag *tr_node, int *index_array, int num_index, double score)
{
	int		i, index;
	int		optimal_index;

	int		nodenum = tr_node->nodenum;
	int		shared_num = tr_node->shared_node_num;

	if(tr_node->pParent == tr_node) {
		//root
		index = 0;
		if(score > tr_node->t_head[index].score) {
			optimal_index = 0;
			for(i=0; i<nodenum; i++) {
//				optimal_index += tr_node->enum_arrayo[i]*floor(pow(double(this->iLeading), double(nodenum-1-i)));
				optimal_index += tr_node->enum_arrayo[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(nodenum-1-i)));
			}
			tr_node->t_head[index].score		= score;
			tr_node->t_head[index].optimal		= 1;
			tr_node->t_head[index].optimalIndex	= optimal_index;
		}
	} else {
		//calculate the index
		index = 0;
		for(i=0; i<shared_num; i++) {
//			index += tr_node->enum_arrayo[i]*floor(pow(double(this->iLeading), double(shared_num-1-i)));
			index += tr_node->enum_arrayo[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(shared_num-1-i)));
		}

		if(score > tr_node->t_head[index].score) {
			optimal_index = 0;
			for(i=0; i<nodenum; i++) {
//				optimal_index += tr_node->enum_arrayo[i]*floor(pow(double(this->iLeading), double(nodenum-1-i)));
				optimal_index += tr_node->enum_arrayo[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(nodenum-1-i)));
			}
			tr_node->t_head[index].score		= score;
			tr_node->t_head[index].optimal		= 1;
			tr_node->t_head[index].optimalIndex	= optimal_index;

//			this->print_enum_arrayo(index_array, num_index);
//			cout<<"\tscore="<<score<<"\tindex="<<index<<"\toptimalIndex="<<optimal_index;
//			cout<<endl;
		}
	}
}
int	DynamicP::getFlagPrintDebugInfo()
{
	return this->flagPrintDebugInfo;
}

void DynamicP::setFlagPrintDebugInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugInfo = 1;
	else
		this->flagPrintDebugInfo = 0;
}

int	DynamicP::getFlagPrintDebugInputInfo()
{
	return this->flagPrintDebugInputInfo;
}

void DynamicP::setFlagPrintDebugInputInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugInputInfo = 1;
	else
		this->flagPrintDebugInputInfo = 0;
}

int	DynamicP::getFlagPrintDebugTreeInfo()
{
	return this->flagPrintDebugTreeInfo;
}

void DynamicP::setFlagPrintDebugTreeInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugTreeInfo = 1;
	else
		this->flagPrintDebugTreeInfo = 0;
}

int	DynamicP::getFlagPrintDebugDPWinInfo()
{
	return this->flagPrintDebugDPWinInfo;
}

void DynamicP::setFlagPrintDebugDPWinInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugDPWinInfo = 1;
	else
		this->flagPrintDebugDPWinInfo = 0;
}

int	DynamicP::getFlagPrintDebugTDInfo()
{
	return this->flagPrintDebugTDInfo;
}

void DynamicP::setFlagPrintDebugTDInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugTDInfo = 1;
	else
		this->flagPrintDebugTDInfo = 0;
}

int	DynamicP::getFlagPrintDebugDPSearchInfo()
{
	return this->flagPrintDebugDPSearchInfo;
}

void DynamicP::setFlagPrintDebugDPSearchInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintDebugDPSearchInfo = 1;
	else
		this->flagPrintDebugDPSearchInfo = 0;
}

void DynamicP::print_enum_arrays(Tree_bag * tr_node)
{
	cout<<"enum_arrays_info\t";
	int i;
	for(i=0; i<tr_node->stnum; i++) {
		cout<<tr_node->enum_arrays[i]<<" ";
	}
	cout<<endl;

}
void DynamicP::print_enum_arrayo(int *index_array, int num_index)
{
	cout<<"enum_arrayo_info\t";
	int i;
	//enum_arrayo:	all possible combinations of nodes.
	for(i=0; i<num_index; i++) {
		cout<<index_array[i]<<" ";
	}
//	cout<<endl;
}
void DynamicP::print_treebag_num_shared_node(Tree_bag *root)
{
	int i, num;
	if(root == NULL)
		return;

	printTreenode(root);

	num = root->shared_node_num;
	cout<<"\tshared_node_num_with_its_parent = "<<num<<endl;

	for(i=0; i<root->childnum; i++) {
		print_treebag_num_shared_node(root->pChild[i]);
	}
}

void DynamicP::print_stem_image(Tree_bag * root)
{
	int i, j;
	if(root == NULL)
		return;

	printTreenode(root);
	cout<<endl;
	for(i=0; i<root->childnum; i++) {
		for(j=0; j<root->nodenum; j++) {
			cout<<root->stem_image[i][j]<<" ";
		}
		cout<<endl;
	}

	for(i=0; i<root->childnum; i++) {
		print_stem_image(root->pChild[i]);
	}
}
void DynamicP::print_loop_image(Tree_bag * tr_node)
{
}
void DynamicP::printTreeBagAllInfo(Tree_bag *root)
{
	int i, iTotal, nodenum;
	nodenum = root->nodenum;
//	iTotal = floor(pow(double(this->iLeading), double(nodenum)));
	iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(nodenum)));

	if(root == NULL)
		return;
	
	this->printTreenode(root);
	cout<<endl;
//	char * cKNary = new char[nodenum+1];
//	memset(cKNary, 0, nodenum+1);
	cout<<"index\toptimal\tscore\toptimal_c_id"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) {
//		itoa(i, cKNary, this->iLeading);
//		cout<<cKNary;
		cout<<i;
		cout<<"\t"<<root->t_head[i].optimal<<"\t"<<root->t_head[i].score;
//		for(j=0; j<root->ptnum; j++)
//			cout<<"\t"<<root->t_head[i].optimal_c_id[j];
		cout<<endl;
	}
//	if(cKNary != NULL) {
//		delete [] cKNary;
//		cKNary = NULL;
//	}

	//then recursively calls itself
	for(i=0; i<root->childnum; i++) {
		printTreeBagAllInfo(root->pChild[i]);
	}
}
void DynamicP::printChildNum(Tree_bag *root)
{
	int i;
	if(root == NULL)
		return;
	this->printTreenode(root);
	cout<<"\t child num="<<root->childnum<<endl;
	for(i=0; i<root->childnum; i++) {
		printChildNum(root->pChild[i]);
	}
}
void DynamicP::printTreeNodeInfo(Tree_bag *root)
{
	int i, iTotal, nodenum;
	nodenum = root->nodenum;
	if(root->pParent == root)
		iTotal = 1;
	else
//		iTotal = floor(pow(double(this->iLeading), double(root->shared_node_num)));
		iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(root->shared_node_num)));
//	char * cKNary = new char[nodenum+1];
//	memset(cKNary, 0, nodenum+1);
	cout<<"index\toptimal\tscore\toptimalIndex"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) {
		if(root->t_head[i].optimal == 1)
		{
	//		itoa(i, cKNary, this->iLeading);
	//		cout<<cKNary;
			cout<<i;
			cout<<"\t"<<root->t_head[i].optimal<<"\t"<<root->t_head[i].score<<"\t"<<root->t_head[i].optimalIndex;
			cout<<endl;
		}
	}
//	if(cKNary != NULL) {
//		delete [] cKNary;
//		cKNary = NULL;
//	}
}
void DynamicP::printTreeBagChildIdInfo(Tree_bag *root)
{
	int i, iTotal, nodenum;
	nodenum = root->nodenum;
//	iTotal = floor(pow(double(this->iLeading), double(nodenum)));
	iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(nodenum)));

	if(root == NULL)
		return;

//	char * cKNary = new char[nodenum+1];
//	memset(cKNary, 0, nodenum+1);

	this->printTreenode(root);
	cout<<endl;
	cout<<"index\toptimal\tscore\toptimalIndex"<<endl;
	//print info of the current treenode
	for(i=0; i<iTotal; i++) {
		if(root->t_head[i].optimal == 1) {
//			itoa(i, cKNary, this->iLeading);
//			cout<<cKNary;
			cout<<i;
			cout<<"\t"<<root->t_head[i].optimal<<"\t"<<root->t_head[i].score<<"\t"<<root->t_head[i].optimalIndex;
//			for(j=0; j<root->ptnum; j++)
//				cout<<"\t"<<root->t_head[i].optimal_c_id[j];
			cout<<endl;
		}
	}
//	if(cKNary != NULL) {
//		delete [] cKNary;
//		cKNary = NULL;
//	}

	//then recursively calls itself
	for(i=0; i<root->childnum; i++) {
		printTreeBagChildIdInfo(root->pChild[i]);
	}
}

int DynamicP::countNumOftheWholeTreeNode(Tree_bag *root)
{
	int i, nodenum;
	nodenum = root->nodenum;
//	this->printTreenode(root);
//	cout<<"\t"<<nodenum<<endl;

	for(i=0; i<root->childnum; i++) {
		nodenum += countNumOftheWholeTreeNode(root->pChild[i]);
	}
	return nodenum;
}

void DynamicP::getTheWholeTreeNode(Tree_bag *root)
{
	//get the current node's tree node info
	int i, nodenum;
	nodenum = root->nodenum;
	int	* pCurrentNode = new int[nodenum];
	memset(pCurrentNode, 0, nodenum);
	Node_list * nodeList = root->nhead;
	i = 0;
	while(nodeList != NULL) {
		pCurrentNode[i++] = nodeList->g_node;
		nodeList = nodeList->next;
	}

	for(i=0; i<nodenum; i++)
		this->pTheWholeTreeNode[this->iTheWholeTreeNodeIndex++] = pCurrentNode[i];

	if(pCurrentNode != NULL) {
		delete [] pCurrentNode;
		pCurrentNode = NULL;
	}

	for(i=0; i<root->childnum; i++) {
		getTheWholeTreeNode(root->pChild[i]);
	}
}

int DynamicP::findOptimalinRoot(Tree_bag *root)
{
	int i, iTotal, nodenum;
	int index = -1;
	nodenum = root->nodenum;
//	iTotal = floor(pow(double(this->iLeading), double(nodenum)));
	iTotal = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(nodenum)));
	for(i=0; i<iTotal; i++) {
		if(root->t_head[i].optimal == 1) {
			index = i;
		}
	}
	return index;
}

void DynamicP::getOptimalPath(Tree_bag *root)
{
//	cout<<"-- getOptimalPath() --"<<endl;

	int i, count, nodenum;
	nodenum = root->nodenum;

//	this->printTreeNodeInfo(root);

	//starts from the root
	int iRoot = root->t_head[0].optimalIndex;//this->findOptimalinRoot(root);
	if(iRoot != -1) {
		//allocate pCurrentOptimalNode[]
		int	* pCurrentOptimalNode = new int[nodenum];
		memset(pCurrentOptimalNode, 0, nodenum);

//		this->printTreenode(root);
//		cout<<endl;
		//save current optimal index to the pCurrentOptimalNode
		count = 0;
		int optimal_index = iRoot;
		int result = 0;
		while (optimal_index) 
		{
			result = optimal_index % this->iLeading;
			pCurrentOptimalNode[count++] = result;
			optimal_index = optimal_index / this->iLeading;
		}
//		for(i=0; i<count; i++)
//			cout<<pCurrentOptimalNode[i]<<"-";
//		cout<<endl;

		for(i=0; i<nodenum-count; i++)
			this->pOptimalTheWholeTreeNode[this->iOptimalTheWholeTreeNodeIndex++] = 0;

		for(i=0; i<count; i++)
			this->pOptimalTheWholeTreeNode[this->iOptimalTheWholeTreeNodeIndex++] = pCurrentOptimalNode[count-1 - i];

		for(i=0; i<nodenum; i++) {
			pCurrentOptimalNode[i] = this->pOptimalTheWholeTreeNode[i];
		}
		//
//		for(i=0; i<nodenum; i++)
//			cout<<pCurrentOptimalNode[i]<<"-";
//		cout<<endl;

		for(i=0; i<root->childnum; i++) {
			getOptimalTheWholeTreeNode(root->pChild[i], pCurrentOptimalNode, nodenum);
		}

		//de-allocate pCurrentOptimalNode[]
		if(pCurrentOptimalNode != NULL) {
			delete [] pCurrentOptimalNode;
			pCurrentOptimalNode = NULL;
		}
	}
}
int DynamicP::getIdxInNodeList(Node_list *n_head, int nodeid) 
{
	int		index = 0;
	bool	bFound = false;

	Node_list *nhead;
	nhead = n_head;
	while(nhead) {
		if(nhead->g_node == nodeid) {
			bFound = true;
			break;
		}
		nhead = nhead->next;
		index++;
	}
	if(bFound == false)
		index = -1;

	return index;
}

//Calculate the current node's optimal node based on its parent's optimal node
void DynamicP::getOptimalTheWholeTreeNode(Tree_bag *pCurrentNode, int * pParentOptimalNode, int optimal_node_num)
{
	if(pCurrentNode == NULL)
		return;

//	cout<<"-- getOptimalTheWholeTreeNode() --"<<endl;

//	this->printTreenode(pCurrentNode);
//	cout<<endl;

	int i, j, k, count, nodenum;
	nodenum = pCurrentNode->nodenum;

//	for(i=0; i<optimal_node_num; i++)
//		cout<<pParentOptimalNode[i]<<"-";
//	cout<<endl;

	int	* pCurrentOptimalNode = new int[nodenum];

	//build the current node's optimal node from parent's optimal node
	Node_list * nhead;
	nhead = pCurrentNode->nhead;
	int nodeid, index;
	count = 0;
	while(nhead) {
		nodeid = nhead->g_node;
		index = this->getIdxInNodeList(pCurrentNode->pParent->nhead, nodeid);
		if(index != -1) {
			//save them in the pCurrentOptimalNode
			pCurrentOptimalNode[count++] = pParentOptimalNode[index];
		}
		nhead = nhead->next;
	}
	int iShared		= pCurrentNode->shared_node_num;
	int iNotShared	= nodenum - iShared;
	int preindex = 0;
	index = 0;
	if(iShared == count) {
		//calculate the preindex
		//bug
		for(i=0; i<iShared; i++) {
//			cout<<"pCurrentOptimalNode["<<i<<"]="<<pCurrentOptimalNode[i]<<endl;
//			preindex += pCurrentOptimalNode[i]*floor(pow(double(this->iLeading), double(iShared-1-i)));
			preindex += pCurrentOptimalNode[i]*static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(iShared-1-i)));
		}

//		int iTotalNotShared = floor(pow(double(this->iLeading), double(iNotShared)));
		int iTotalNotShared = static_cast<int>(pow(static_cast<double>(this->iLeading), static_cast<double>(iNotShared)));
		int * pLocalOptimalNode = new int[iNotShared];

		//Calculate the following iNotShared bit based on pCurrentOptimalNode and .optimalIndex
		if(pCurrentNode->t_head[preindex].optimal == 1) {
			//it should be one optimal result
			//save i in the pCurrentOptimalNode
			j = 0;
			int ii = pCurrentNode->t_head[preindex].optimalIndex - preindex*iTotalNotShared;
			int result = 0;
			while(ii) {
				result = ii % this->iLeading;
				pLocalOptimalNode[j++] = result;
				ii = ii / this->iLeading;
			}

			for(k=0; k<iNotShared-j; k++) {
				pCurrentOptimalNode[count++] = 0;
			}
			for(k=0; k<j; k++) {
				pCurrentOptimalNode[count++] = pLocalOptimalNode[j-1-k];
			}
		}

		if(pLocalOptimalNode != NULL) {
			delete [] pLocalOptimalNode;
			pLocalOptimalNode = NULL;
		}

//		for(i=0; i<nodenum; i++)
//			cout<<pCurrentOptimalNode[i]<<"-";
//		cout<<endl;

		//copy pCurrentOptimalNode to this->pOptimalTheWholeTreeNode
		for(i=0; i<nodenum; i++) {
			this->pOptimalTheWholeTreeNode[this->iOptimalTheWholeTreeNodeIndex++] = pCurrentOptimalNode[i];
		}
	}

//	for(i=0; i<nodenum; i++)
//		cout<<pCurrentOptimalNode[i]<<"-";
//	cout<<endl;

	for(i=0; i<pCurrentNode->childnum; i++) {
		getOptimalTheWholeTreeNode(pCurrentNode->pChild[i], pCurrentOptimalNode, nodenum);
	}

	if(pCurrentOptimalNode != NULL) {
		delete [] pCurrentOptimalNode;
		pCurrentOptimalNode = NULL;
	}
}

void DynamicP::buildStemIdxArray()
{
	int i, j;
	int stemnum = this->getNumstems();

	stemPosArray = new int[2*stemnum];
	stemIdxArray = new int[2*stemnum];
	for(i=0; i<2*stemnum; i++) {
		stemPosArray[i] = stemIdxArray[i] = 0;
	}

	for(i=0; i<stemnum; i++) {
		j = 2*i;
		stemIdxArray[j]		= j+1;
		stemIdxArray[j+1]	= j+2;
		stemPosArray[j]		= this->pStem[i].start1;
		stemPosArray[j+1]	= this->pStem[i].start2;
	}

	int iTmp;

	for(i=1; i<2*stemnum; i++)
	{
		for(j=2*stemnum-1; j>=i; j--)
		{
			if(stemPosArray[j-1] > stemPosArray[j])
			{
				//swap the element
				iTmp = stemPosArray[j];
				stemPosArray[j] = stemPosArray[j-1];
				stemPosArray[j-1] = iTmp;

				iTmp = stemIdxArray[j];
				stemIdxArray[j] = stemIdxArray[j-1];
				stemIdxArray[j-1] = iTmp;
			}
		}
	}

	if(stemPosArray != NULL) {
		delete [] stemPosArray;
		stemPosArray = NULL;
	}
}

void DynamicP::printStemIdxArray()
{
	int i, j;
	int stemnum = this->getNumstems();

	for(i=0; i<stemnum; i++) {
		j = 2*i;
		cout<<stemIdxArray[j]<<"-"<<stemIdxArray[j+1]<<"-";
	}
	cout<<endl;
}

void DynamicP::freeStemIdxArray()
{
	if(stemIdxArray != NULL)
		delete [] stemIdxArray;
	stemIdxArray = NULL;
}

void DynamicP::allocateStemCandIdxArray()
{
	int stemnum = this->getNumstems();
	stemCandIdxArray = new int[2*stemnum];
}
void DynamicP::buildStemCandIdxArray()
{
	int i, index;
	int iElement;
	int stemnum = this->getNumstems();

	//for every element in stemIdxArray, find its corresponding 
	//optimal index in pOptimalTheWholeTreeNode of pTheWholeTreeNode
	for(i=0; i<2*stemnum; i++) {
		iElement = stemIdxArray[i];
		index = 0;
		while(index<iTheWholeTreeNodeIndex) {
			if(iElement == pTheWholeTreeNode[index])
				break;
			index++;
		}
		stemCandIdxArray[i] = pOptimalTheWholeTreeNode[index];
	}
}

void DynamicP::printStemCandIdxArray()
{
	int i;
	int stemnum = this->getNumstems();

	for(i=0; i<2*stemnum; i++) {
		cout<<stemCandIdxArray[i]<<" ";
	}
	cout<<endl;
}

void DynamicP::freeStemCandIdxArray()
{
	if(stemCandIdxArray != NULL)
		delete [] stemCandIdxArray;
	stemCandIdxArray = NULL;
}

void DynamicP::locateStemStartEndPosIdx()
{
	int stemnum = this->getNumstems();
	iStemStartPosIdx	= stemCandIdxArray[0];
	iStemEndPosIdx		= stemCandIdxArray[2*stemnum-1];
}

int	DynamicP::getFlagMergeCandInPreprocess()
{
	return this->flagMergeCandInPreprocess;
}

void DynamicP::setFlagMergeCandInPreprocess(int pflag)
{
	if(pflag >= 1)
		this->flagMergeCandInPreprocess = 1;
	else
		this->flagMergeCandInPreprocess = 0;
}

int	DynamicP::getFlagCandwithShortestLength()
{
	return this->flagCandwithShortestLength;
}

void DynamicP::setFlagCandwithShortestLength(int pflag)
{
	if(pflag >= 1)
		this->flagCandwithShortestLength = 1;
	else
		this->flagCandwithShortestLength = 0;
}

int	DynamicP::getShiftNumMergeCand()
{
	return this->iShiftNumMergeCand;
}

void DynamicP::setShiftNumMergeCand(int num)
{
	this->iShiftNumMergeCand = num;
}

int	DynamicP::getAllowedNullLoopInsNum()
{
	return this->iAllowedNullLoopInsNum;
}

void DynamicP::setAllowedNullLoopInsNum(int num)
{
	this->iAllowedNullLoopInsNum = num;
}

void DynamicP::setNumOfNtsOverlapBetweenStem(int num)
{
	this->iNumOfNtsOverlapBetweenStem = num;
}

int DynamicP::getNumOfNtsOverlapBetweenStem()
{
	return this->iNumOfNtsOverlapBetweenStem;
}

void DynamicP::IterateTD_topdown(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	unsigned int i = 0;
	cout << pRoot->m_nID<< ':';// << pRoot->m_pPayload << " : ";
	while( i < pRoot->m_vecpNodes.size() )
	{
		cout << pRoot->m_vecpNodes[i]->ID() << ' ';// << pRoot->m_vecpNodes[i]->Payload() << ' ';
		i++;
	}
	cout << endl;
	
	i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		IterateTD_topdown( pRoot->m_vecChildren[i] );
		i++;
	}
}

void DynamicP::IterateTD_Payload_topdown(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	unsigned int i = 0;
	cout << pRoot->m_nID<< ':';// << pRoot->m_pPayload << " : ";
	while( i < pRoot->m_vecpNodes.size() )
	{
		cout << pRoot->m_vecpNodes[i]->Payload() << ' ';
		i++;
	}
	cout << endl;
	
	i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		IterateTD_Payload_topdown( pRoot->m_vecChildren[i] );
		i++;
	}
}

void DynamicP::IterateTD_Payload_bottomup(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	unsigned int i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		IterateTD_Payload_bottomup( pRoot->m_vecChildren[i] );
		i++;
	}

	cout << pRoot->m_nID<< ':';// << pRoot->m_pPayload << " : ";
	i = 0;
	while( i < pRoot->m_vecpNodes.size() )
	{
		cout << pRoot->m_vecpNodes[i]->Payload() << ' ';
		i++;
	}
	cout << endl;
	
}

//create my own tree-bag tree by copying the tree-bag
Tree_bag * DynamicP::copy_treedecomp(TRnaGraphTreeBag<int, Tree_bag>* pRoot)
{
	Tree_bag * cur_node = NULL;
	Node_list * nodelist = NULL;
	unsigned int i = 0;
	//create a tree-bag, and copy the tree-node to its node-list
	cur_node = new Tree_bag[1];
	cur_node->nodenum = pRoot->m_vecpNodes.size();
	cur_node->nhead = new Node_list[1];
	nodelist = cur_node->nhead;

	while(i < (pRoot->m_vecpNodes.size()-1) )
	{
		nodelist->g_node = pRoot->m_vecpNodes[i]->Payload();
		nodelist->next = new Node_list[1];
		nodelist = nodelist->next;
		i++;
	}
	nodelist->g_node = pRoot->m_vecpNodes[i]->Payload();
	nodelist->next = NULL;

	//then deal with its child tree-bag.
	cur_node->childnum = pRoot->m_vecChildren.size();
	i = 0;
	while( i < pRoot->m_vecChildren.size() )
	{
		cur_node->pChild[i] = copy_treedecomp(pRoot->m_vecChildren[i]);
		i++;
	}
	return cur_node;
}

void DynamicP::freeMyTDRoot()
{
	this->free_treenode(this->mypTDRoot);
}
Tree_bag * DynamicP::getMyTDRoot()
{
	return this->mypTDRoot;
}

double DynamicP::getPcoeff()
{
	return this->pcoeff;
}

void DynamicP::setPcoeff(double _pcoeff)
{
	this->pcoeff = _pcoeff;
}

int	DynamicP::getFlagPrintStrAlignInfo()
{
	return this->flagPrintStrAlignInfo;
}

void DynamicP::setFlagPrintStrAlignInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintStrAlignInfo = 1;
	else
		this->flagPrintStrAlignInfo = 0;
}

int	DynamicP::getFlagPrintScoreInfo()
{
	return this->flagPrintScoreInfo;
}

void DynamicP::setFlagPrintScoreInfo(int pflag)
{
	if(pflag >= 1)
		this->flagPrintScoreInfo = 1;
	else
		this->flagPrintScoreInfo = 0;
}

/*
 * The idea is to extract all the information from
 *		pTheWholeTreeNode
 *		stemIdxArray
 *		stemCandidIdxArray
 */
void DynamicP::buildStructureAlignInfo(Loop *pLoopModel)
{
	int m, n;
	int stemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = -1, preRightStemPos = -1;
	int curLeftStemPos = -1, curRightStemPos = -1;
	int pre_gid = 0, cur_gid = 0;

	VTBSearch vtber(1);
	vtber.setAllowedInsNum(this->getAllowedNullLoopInsNum());
	int loopid, start, end;
//	double dLoopScore = 0.0;

	//structure alignment
	memset(structureSequence, 0, MAX_STRUCTURE_ALIGN_LENGTH);
	memset(strucutreAlignment, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	//support stemid index in two pastalines 20081015
	memset(strucutreAlignment_charid_idx, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	char *	tmp_sub_charid_idx = NULL;

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;

			if(stemIdx == 0) {
				strcat(structureSequence, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
				strcat(strucutreAlignment, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment);

				//support stemid index in two pastalines 20081015
				if(this->getNumPastaLine() == TWO_PASTALINE) {
					strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].pastaline_idx);
				}

			} else {
				if(curLeftStemPos > preRightStemPos) 
				{
					if(curLeftStemPos > (preRightStemPos+1))
					{
						//add loop structure alignment info here.

						char * loopAlign	= NULL;
						char * loopSequence	= NULL;

						loopid	= this->identifyLoopId(pre_gid, cur_gid);
						start	= preRightStemPos + 1;
						end		= curLeftStemPos - 1;

						double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, &loopAlign, &loopSequence);
						//add loop sequence/alignment info
						if((this->getAllowedNullLoopInsNum() >= (end-start+1))
							&& (dLoopScore == 0.0))
						{
							for(int n=start; n<=end; n++)
							{
								char *	cTmp = strappend(structureSequence, cWindowSequence[n]);
								strcpy(structureSequence, cTmp);
								delete [] cTmp; cTmp = NULL;
								cTmp = strappend(strucutreAlignment, '.');
								strcpy(strucutreAlignment, cTmp);

								//support stemid index in two pastalines 20081015
								if(this->getNumPastaLine() == TWO_PASTALINE) {
									strcpy(strucutreAlignment_charid_idx, cTmp);
								}
								delete [] cTmp; cTmp = NULL;
							}
						} else {
//							cout<<loopSequence<<" <-> "<<loopAlign<<endl;
							//postprocessing the loopSequence/loopAlign
							//if '-' exists in loopSequence and 'd' in loopalign, then delete 'd' from loopalign
							strcat(structureSequence, loopSequence);
							strcat(strucutreAlignment, loopAlign);

							//support stemid index in two pastalines 20081015
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								strcat(strucutreAlignment_charid_idx, loopAlign);
							}
						}

						if(loopAlign != NULL) {
							delete [] loopAlign;
							loopAlign = NULL;
						}
						if(loopSequence != NULL) {
							delete [] loopSequence;
							loopSequence = NULL;
						}
					}

					//add the next stem sequence/alignment info
					strcat(structureSequence, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
					strcat(strucutreAlignment, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment);

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].pastaline_idx);	//20081015
					}
				} else {
					//overlap
					int	iOverlap = preRightStemPos - curLeftStemPos + 1;
					int	iTmpLength = strlen(strucutreAlignment);
					for(n=0; n<iOverlap; n++) {
						strucutreAlignment[iTmpLength-1-n] = '*';

						//support stemid index in two pastalines 20081015
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							strucutreAlignment_charid_idx[iTmpLength-1-n] = '*';
						}
					}


					iTmpLength = strlen(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
					char * cSeqOther = new char[iTmpLength - iOverlap + 1];
					char * cAliOther = new char[iTmpLength - iOverlap + 1];

					for(n=0; n<(iTmpLength - iOverlap); n++) {
						cSeqOther[n] = candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence[n+iOverlap];
						cAliOther[n] = candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment[n+iOverlap];
					}
					cSeqOther[n] = '\0';
					cAliOther[n] = '\0';
					
					strcat(structureSequence, cSeqOther);
					strcat(strucutreAlignment, cAliOther);

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						tmp_sub_charid_idx = strrpl(cAliOther, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
						strcat(strucutreAlignment_charid_idx, tmp_sub_charid_idx);
						if(tmp_sub_charid_idx != NULL) {
							delete [] tmp_sub_charid_idx;
							tmp_sub_charid_idx = NULL;
						}
					}

					if(cSeqOther != NULL) {
						delete cSeqOther;
						cSeqOther = NULL;
					}
					if(cAliOther != NULL) {
						delete cAliOther;
						cAliOther = NULL;
					}
				}
			}
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;
			if(curLeftStemPos > preRightStemPos) 
			{
				if(curLeftStemPos > (preRightStemPos+1))
				{
					//add loop structure alignment info here.

					char * loopAlign	= NULL;
					char * loopSequence	= NULL;

					loopid	= this->identifyLoopId(pre_gid, cur_gid);
					start	= preRightStemPos + 1;
					end		= curLeftStemPos - 1;
					double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, &loopAlign, &loopSequence);
					//add loop sequence/alignment info
					if((this->getAllowedNullLoopInsNum() >= (end-start+1))
						&& (dLoopScore == 0.0))
					{
						for(int n=start; n<=end; n++)
						{
							char *	cTmp = strappend(structureSequence, cWindowSequence[n]);
							strcpy(structureSequence, cTmp);
							if(cTmp != NULL) {
								delete [] cTmp; cTmp = NULL;
							}
							cTmp = strappend(strucutreAlignment, '.');
							strcpy(strucutreAlignment, cTmp);

							//support stemid index in two pastalines 20081015
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								strcpy(strucutreAlignment_charid_idx, cTmp);
							}
							if(cTmp != NULL) {
								delete [] cTmp; cTmp = NULL;
							}
						}
					} else {
//						cout<<loopSequence<<" <-> "<<loopAlign<<endl;
						strcat(structureSequence, loopSequence);
						strcat(strucutreAlignment, loopAlign);

						//support stemid index in two pastalines 20081015
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							strcat(strucutreAlignment_charid_idx, loopAlign);
						}
					}

					if(loopAlign != NULL) {
						delete [] loopAlign;
						loopAlign = NULL;
					}
					if(loopSequence != NULL) {
						delete [] loopSequence;
						loopSequence = NULL;
					}
				}

				//add the next stem sequence/alignment info
				strcat(structureSequence, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
				strcat(strucutreAlignment, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment);

				//support stemid index in two pastalines 20081015
				if(this->getNumPastaLine() == TWO_PASTALINE) {
					strcat(strucutreAlignment_charid_idx, candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].pastaline_idx);	//20081015
				}
			} else {
				//overlap
				int	iOverlap = preRightStemPos - curLeftStemPos + 1;

				int	iTmpLength = strlen(strucutreAlignment);
				for(n=0; n<iOverlap; n++) {
					strucutreAlignment[iTmpLength-1-n] = '*';
					//support stemid index in two pastalines 20081015
					strucutreAlignment_charid_idx[iTmpLength-1-n] = '*';
				}

				iTmpLength = strlen(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
				char * cSeqOther = new char[iTmpLength - iOverlap + 1];
				char * cAliOther = new char[iTmpLength - iOverlap + 1];

				for(n=0; n<(iTmpLength - iOverlap); n++) {
					cSeqOther[n] = candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence[n+iOverlap];
					cAliOther[n] = candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment[n+iOverlap];
				}
				cSeqOther[n] = '\0';
				cAliOther[n] = '\0';
				strcat(structureSequence, cSeqOther);
				strcat(strucutreAlignment, cAliOther);

				//support stemid index in two pastalines 20081015
				if(this->getNumPastaLine() == TWO_PASTALINE) {
					tmp_sub_charid_idx = strrpl(cAliOther, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
					strcat(strucutreAlignment_charid_idx, tmp_sub_charid_idx);
					if(tmp_sub_charid_idx != NULL) {
						delete [] tmp_sub_charid_idx;
						tmp_sub_charid_idx = NULL;
					}
				}

				if(cSeqOther != NULL) {
					delete cSeqOther;
					cSeqOther = NULL;
				}
				if(cAliOther != NULL) {
					delete cAliOther;
					cAliOther = NULL;
				}
			}

		}
		preLeftStemPos	= curLeftStemPos;
		preRightStemPos	= curRightStemPos;
		
		pre_gid	= cur_gid;
	}//for
}

void DynamicP::printStructureAlignInfo(int real_start, int real_end)
{
	//preprocess the structureAlignment
	char * strAlign_01 = strrpl(strucutreAlignment, 'm', '.');
	char * strAlign_02 = strrpl(strAlign_01, 'd', '.');
	char * strAlign_03 = strrpl(strAlign_02, 'i', '-');
	
	char * strAlign_04 = NULL;
	char * strAlign_05 = NULL;
	char * strAlign_06 = NULL;
	if(this->getNumPastaLine() == TWO_PASTALINE) {
		//preprocess the strucutreAlignment_charid_idx
		strAlign_04 = strrpl(strucutreAlignment_charid_idx, 'm', '.');
		strAlign_05 = strrpl(strAlign_04, 'd', '.');
//		strAlign_06 = strrpl(strAlign_05, 'i', '-');
		strAlign_06 = strrpl(strAlign_05, 'i', '.');

/*
		formatStringOutput(	string("Structure alignment"), 
							string(strAlign_03),
							string(strAlign_06),
							string(structureSequence), 
							real_start,
							real_end,
							1);
*/
		formatStringOutput2("Structure alignment", 
							strAlign_03,
							strAlign_06,
							structureSequence, 
							real_start,
							real_end,
							1);
	} else {
/*
		formatStringOutput(	string("Structure alignment"), 
							string(strAlign_03),
							string(""),
							string(structureSequence), 
							real_start,
							real_end,
							0);
*/
		formatStringOutput2("Structure alignment", 
							strAlign_03,
							"",
							structureSequence, 
							real_start,
							real_end,
							0);
	}


	if(strAlign_01 != NULL) {
		delete [] strAlign_01;
		strAlign_01 = NULL;
	}
	if(strAlign_02 != NULL) {
		delete [] strAlign_02;
		strAlign_02 = NULL;
	}
	if(strAlign_03 != NULL) {
		delete [] strAlign_03;
		strAlign_03 = NULL;
	}
	if(strAlign_04 != NULL) {
		delete [] strAlign_04;
		strAlign_04 = NULL;
	}
	if(strAlign_05 != NULL) {
		delete [] strAlign_05;
		strAlign_05 = NULL;
	}
	if(strAlign_06 != NULL) {
		delete [] strAlign_06;
		strAlign_06 = NULL;
	}
}

/*
 * print the candidate stem info
 * 20081024
 */
void DynamicP::printCandidateStemInfo()
{
	int m, n;
	int stemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = 0, preRightStemPos = 0;
	int curLeftStemPos = 0, curRightStemPos = 0;
	int pre_gid, cur_gid;
	pre_gid = cur_gid = 0;

	for(n=0; n<this->getLeadingNum(); n++)
	{
		if(n == 0)
		{
			cout<<"Result from Preprocessing"<<endl;
			//print out the stem charid
			for(m=0; m<2*stemnum; m++)
			{
				cur_gid = stemIdxArray[m];
				stemIdx = (cur_gid - 1)/2;

				if(cur_gid % 2 != 0) 
				{	
					//odd - left_arm
					cout<<"       "<<this->pStem[stemIdx].charid<<left<<setw(WIDTH-8)<<"";
				} 
				else 
				{	
					//even - right_arm
					cout<<"       "<<this->pStem[stemIdx].charid<<left<<setw(WIDTH-8)<<"";
				}
			}
			cout<<endl;
		}
		//print out the position of one stem
		cout<<left<<setw(2)<<n<<":";
		for(m=0; m<2*stemnum; m++)
		{
			cur_gid = stemIdxArray[m];
			stemIdx = (cur_gid - 1)/2;

			if(cur_gid % 2 != 0) 
			{	
				//odd - left_arm
				curLeftStemPos	= candidateStems[stemIdx].s_locinfo[n].a;
				curRightStemPos = candidateStems[stemIdx].s_locinfo[n].b;
				cout<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
			} 
			else 
			{	
				//even - right_arm
				curLeftStemPos	= candidateStems[stemIdx].s_locinfo[n].c;
				curRightStemPos = candidateStems[stemIdx].s_locinfo[n].d;
				cout<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
			}
		}
		cout<<endl;
		//print out the score for that stem
		cout<<"   ";
		for(m=0; m<2*stemnum; m++)
		{
			cur_gid = stemIdxArray[m];
			stemIdx = (cur_gid - 1)/2;

			if(cur_gid % 2 != 0) 
			{	
				//odd - left_arm
				if(candidateStems[stemIdx].s_locinfo[n].p_score > 1)
					cout<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.1)
					cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.01)
					cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0)
					cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.01)
					cout<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.1)
					cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > -1.0)
					cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else
					cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
			} 
			else 
			{	
				//even - right_arm
				if(candidateStems[stemIdx].s_locinfo[n].p_score > 1)
					cout<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.1)
					cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.01)
					cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0)
					cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.01)
					cout<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.1)
					cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else if(candidateStems[stemIdx].s_locinfo[n].p_score > -1.0)
					cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
				else
					cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[n].p_score<<" ) ";
			}
		}
		cout<<endl;
		if(n == (this->getLeadingNum()-1))
		{
			cout<<endl<<"Result from Dynamic Programming"<<endl;
			//print out the candidate stem index in valid combination
			for(m=0; m<2*stemnum; m++)
			{
				cout<<"       "<<stemCandIdxArray[m]<<left<<setw(WIDTH-8)<<"";
			}
			cout<<endl;
			cout<<"   ";
			for(m=0; m<2*stemnum; m++)
			{
				cur_gid = stemIdxArray[m];
				stemIdx = (cur_gid - 1)/2;

				if(cur_gid % 2 != 0) 
				{	
					//odd - left_arm
					curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
					curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;
					cout<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
				} 
				else 
				{	
					//even - right_arm
					curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
					curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;
					cout<<"["<<right<<setw(3)<<curLeftStemPos<<","<<left<<setw(3)<<curRightStemPos<<"] ";
				}
			}
			cout<<endl;
			//print out the score for that stem
			cout<<"   ";
			for(m=0; m<2*stemnum; m++)
			{
				cur_gid = stemIdxArray[m];
				stemIdx = (cur_gid - 1)/2;

				if(cur_gid % 2 != 0) 
				{	
					//odd - left_arm
					if(candidateStems[stemIdx].s_locinfo[n].p_score > 1)
						cout<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.1)
						cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.01)
						cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0)
						cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.01)
						cout<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.1)
						cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > -1.0)
						cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else
						cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
				} 
				else 
				{	
					//even - right_arm
					if(candidateStems[stemIdx].s_locinfo[n].p_score > 1)
						cout<<"( "<<setw(5)<<setprecision(4)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.1)
						cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0.01)
						cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > 0)
						cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.01)
						cout<<"( "<<setw(5)<<setprecision(0)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > -0.1)
						cout<<"( "<<setw(5)<<setprecision(1)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else if(candidateStems[stemIdx].s_locinfo[n].p_score > -1.0)
						cout<<"( "<<setw(5)<<setprecision(2)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
					else
						cout<<"( "<<setw(5)<<setprecision(3)<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<" ) ";
				}
			}
			cout<<endl;
		}
	}
}

void DynamicP::buildFoldedStructureInfo()
{
	int m, n;
	int start, end;
	int stemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = 0, preRightStemPos = 0;
	int curLeftStemPos = 0, curRightStemPos = 0;
	int pre_gid, cur_gid;
	pre_gid = cur_gid = 0;

	//folded structure
//	memset(strucutreAlignment, 0, MAX_STRUCTURE_ALIGN_LENGTH);
	char	tmpFoldedstrucutre[MAX_STRUCTURE_ALIGN_LENGTH];
	memset(tmpFoldedstrucutre, 0, MAX_STRUCTURE_ALIGN_LENGTH);
    memset(foldedstrucutre, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	//support stemid index in two pastalines 20081015
	char *	tmp_sub_charid_idx = NULL;
	char	tmp_charid_idx[MAX_STRUCTURE_ALIGN_LENGTH];
	memset(tmp_charid_idx, 0, MAX_STRUCTURE_ALIGN_LENGTH);
	memset(foldedstrucutre_charid_idx, 0, MAX_STRUCTURE_ALIGN_LENGTH);

	const int left=0, right=1;
	char * tmpSubFold = NULL;

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;

			if(stemIdx == 0) {
				tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], left);
				strcat(tmpFoldedstrucutre, tmpSubFold);

				//support stemid index in two pastalines 20081015
				if(this->getNumPastaLine() == TWO_PASTALINE)
				{
					tmp_sub_charid_idx = strrpl(tmpSubFold, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
					strcat(tmp_charid_idx, tmp_sub_charid_idx);
					if(tmp_sub_charid_idx != NULL) {
						delete [] tmp_sub_charid_idx;
						tmp_sub_charid_idx = NULL;
					}
				}

				if(tmpSubFold != NULL) {
					delete [] tmpSubFold;
					tmpSubFold = NULL;
				}
			} else {
				if(curLeftStemPos > preRightStemPos) 
				{
					if(curLeftStemPos > (preRightStemPos+1))
					{
						//to loop, just append '.'
						start	= preRightStemPos + 1;
						end		= curLeftStemPos - 1;
						for(n=start; n<=end; n++) {
							strcat(tmpFoldedstrucutre, ".");
							if(this->getNumPastaLine() == TWO_PASTALINE) {
								//support stemid index in two pastalines 20081015
								strcat(tmp_charid_idx, ".");
							}
						}
					}

					//add the next stem sequence/alignment info
					tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], left);
					strcat(tmpFoldedstrucutre, tmpSubFold);

					//support stemid index in two pastalines 20081015
					if(this->getNumPastaLine() == TWO_PASTALINE) 
					{
						tmp_sub_charid_idx = strrpl(tmpSubFold, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
						strcat(tmp_charid_idx, tmp_sub_charid_idx);
						if(tmp_sub_charid_idx != NULL) {
							delete [] tmp_sub_charid_idx;
							tmp_sub_charid_idx = NULL;
						}
					}
					
					if(tmpSubFold != NULL) {
						delete [] tmpSubFold;
						tmpSubFold = NULL;
					}
				} else {
					//overlap
					int	iOverlap = preRightStemPos - curLeftStemPos + 1;
					int	iTmpLength = strlen(tmpFoldedstrucutre);
					for(n=0; n<iOverlap; n++) {
						tmpFoldedstrucutre[iTmpLength-1-n] = '*';
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							//support stemid index in two pastalines 20081015
							tmp_charid_idx[iTmpLength-1-n] = '*';
						}
					}

					iTmpLength = strlen(candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].sequence);
					char * cAliOther = new char[iTmpLength - iOverlap + 1];
					for(n=0; n<(iTmpLength - iOverlap); n++) {
						cAliOther[n] = candidateStems[stemIdx].leftArm[stemCandIdxArray[m]].alignment[n+iOverlap];
					}
					cAliOther[n] = '\0';
					strcat(tmpFoldedstrucutre, cAliOther);
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						//support stemid index in two pastalines 20081015
						tmp_sub_charid_idx = strrpl(cAliOther, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
						strcat(tmp_charid_idx, tmp_sub_charid_idx);
						if(tmp_sub_charid_idx != NULL) {
							delete [] tmp_sub_charid_idx;
							tmp_sub_charid_idx = NULL;
						}
					}
					if(cAliOther != NULL) {
						delete [] cAliOther;
						cAliOther = NULL;
					}
				}
			}
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;
			if(curLeftStemPos > preRightStemPos) 
			{
				if(curLeftStemPos > (preRightStemPos+1))
				{
					//add loop structure alignment info here.
					start	= preRightStemPos + 1;
					end		= curLeftStemPos - 1;
					for(n=start; n<=end; n++) {
						strcat(tmpFoldedstrucutre, ".");
						if(this->getNumPastaLine() == TWO_PASTALINE) {
							//support stemid index in two pastalines 20081015
							strcat(tmp_charid_idx, ".");
						}
					}
				}

				//add the next stem sequence/alignment info
				tmpSubFold = this->buildstemfoldedstructure(stemIdx, stemCandIdxArray[m], right);
				strcat(tmpFoldedstrucutre, tmpSubFold);

				//support stemid index in two pastalines 20081015
				if(this->getNumPastaLine() == TWO_PASTALINE) 
				{
					tmp_sub_charid_idx = strrpl(tmpSubFold, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
					strcat(tmp_charid_idx, tmp_sub_charid_idx);
					if(tmp_sub_charid_idx != NULL) {
						delete [] tmp_sub_charid_idx;
						tmp_sub_charid_idx = NULL;
					}
				}
				
				if(tmpSubFold != NULL) {
					delete [] tmpSubFold;
					tmpSubFold = NULL;
				}
			} else {
				//overlap
				int	iOverlap = preRightStemPos - curLeftStemPos + 1;
				int	iTmpLength = strlen(tmpFoldedstrucutre);
				for(n=0; n<iOverlap; n++) {
					tmpFoldedstrucutre[iTmpLength-1-n] = '*';
					if(this->getNumPastaLine() == TWO_PASTALINE) {
						//support stemid index in two pastalines 20081015
						tmp_charid_idx[iTmpLength-1-n] = '*';
					}
				}

				iTmpLength = strlen(candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].sequence);
				char * cAliOther = new char[iTmpLength - iOverlap + 1];
				for(n=0; n<(iTmpLength - iOverlap); n++) {
					cAliOther[n] = candidateStems[stemIdx].rightArm[stemCandIdxArray[m]].alignment[n+iOverlap];
				}
				cAliOther[n] = '\0';
				strcat(tmpFoldedstrucutre, cAliOther);
				if(this->getNumPastaLine() == TWO_PASTALINE) {
					//support stemid index in two pastalines 20081015
					tmp_sub_charid_idx = strrpl(cAliOther, this->pStem[stemIdx].charid, this->pStem[stemIdx].charid_idx);
					strcat(tmp_charid_idx, tmp_sub_charid_idx);
					if(tmp_sub_charid_idx != NULL) {
						delete [] tmp_sub_charid_idx;
						tmp_sub_charid_idx = NULL;
					}
				}
				if(cAliOther != NULL) {
					delete [] cAliOther;
					cAliOther = NULL;
				}
			}
		}
		preLeftStemPos	= curLeftStemPos;
		preRightStemPos	= curRightStemPos;
		
		pre_gid	= cur_gid;
	}//for
	char * cResultOne = strrpl(tmpFoldedstrucutre, 'l', '.');
	char * cResultTwo = strrpl(cResultOne, 'r', '.');
	strcpy(foldedstrucutre, cResultTwo);
	if(cResultOne != NULL) {
		delete [] cResultOne;
		cResultOne = NULL;
	}
	if(cResultTwo != NULL) {
		delete [] cResultTwo;
		cResultTwo = NULL;
	}

	//support stemid index in two pastalines 20081015
	if(this->getNumPastaLine() == TWO_PASTALINE)
	{
		char * cCharidOne = strrpl(tmp_charid_idx, 'l', '.');
		char * cCharidTwo = strrpl(cCharidOne, 'r', '.');
		strcpy(foldedstrucutre_charid_idx, cCharidTwo);
		if(cCharidOne != NULL) {
			delete [] cCharidOne;
			cCharidOne = NULL;
		}
		if(cCharidTwo != NULL) {
			delete [] cCharidTwo;
			cCharidTwo = NULL;
		}
	}
}
//void DynamicP::printFoldedStructureInfo(bool bDebug)
//{
//	if(bDebug)
//		cout<<"#";
//	cout<<foldedstrucutre<<endl;
//}
void DynamicP::printFoldedStructureInfo(int start, int end, int real_start, int real_end)
{
	//support stemid index in two pastalines 20081015
	if(this->getNumPastaLine() == TWO_PASTALINE)
	{
/*
		formatStringOutput(	string("Folded structure"), 
							string(foldedstrucutre),
							string(foldedstrucutre_charid_idx),		
							string(strextract(cWindowSequence, start, end)), 
							real_start,
							real_end,
							1);
*/
		formatStringOutput2("Folded structure", 
							foldedstrucutre,
							foldedstrucutre_charid_idx,		
							strextract(cWindowSequence, start, end), 
							real_start,
							real_end,
							1);
	} else {
/*
		formatStringOutput(	string("Folded structure"), 
							string(foldedstrucutre),
							string(""),
							string(strextract(cWindowSequence, start, end)), 
							real_start,
							real_end,
							0);
*/
		formatStringOutput2("Folded structure", 
							foldedstrucutre,
							"",
							strextract(cWindowSequence, start, end), 
							real_start,
							real_end,
							0);
	}
}
/*
 * print out some format like the following pasta structure info.
 *
 * 1........10........20........30....
 * ....|....|....|....|....|....|.....
 */
void DynamicP::printHeader(int startPos, int endPos, bool bDebug)
{
	int m;
	int iLength = endPos - startPos + 1;

	if(bDebug)
		cout<<"#";
	for(m=1; m<iLength+1; m++) 
	{
		if((m / 100) > 0)
		{
			if(m%10 == 0)
				cout<<m;
//			else if(m%10 != 1 && m%10 != 2)
//				cout<<".";
			else
			{
				if(m < 1000)
				{
					if(m%10 != 1 && m%10 != 2)
						cout<<".";
				} else if(m < 10000) {
					if(m%10 != 1 && m%10 != 2 && m%10 != 3)
						cout<<".";
				}
			}
		} else {
			if(m==1)
				cout<<m;
			else if(m%10 == 0)
				cout<<m;
			else if(m%10 != 1)
				cout<<".";
		}
	}
	cout<<endl;

	if(bDebug)
		cout<<"#";
	for(m=1; m<iLength+1; m++) {
		if(m % 5 ==0)
			cout<<"|";
		else
			cout<<".";
	}
	cout<<endl;

	if(bDebug)
		cout<<"#";
	cout<<"Folded structure"<<endl;	//

	if(bDebug)
		cout<<"#";
	for(m=startPos; m<endPos+1; m++)
		cout<<cWindowSequence[m];
	cout<<endl;
}
//build the folded structure for every stem based on the sequence/alignment of left/right arm of the stem
char * DynamicP::buildstemfoldedstructure(int stemIdx, int rank, int arm)
{
	int m, n, count, m_num;
	int left_length	= strlen(candidateStems[stemIdx].leftArm[rank].alignment);
	int right_length= strlen(candidateStems[stemIdx].rightArm[rank].alignment);
	char * left_seq					= new char[left_length + 1];;
	char * left_align				= new char[left_length + 1];
	char * right_seq				= new char[right_length + 1];
	char * right_align				= new char[right_length + 1];

	strcpy(left_seq, candidateStems[stemIdx].leftArm[rank].sequence);
	strcpy(left_align, candidateStems[stemIdx].leftArm[rank].alignment);

	strcpy(right_seq, candidateStems[stemIdx].rightArm[rank].sequence);
	strcpy(right_align, candidateStems[stemIdx].rightArm[rank].alignment);

	char * tmp_folded_structure = NULL;
	char * folded_structure		= NULL;
	m_num = 0;
	if(arm == 0)	//left arm
	{
		tmp_folded_structure= new char[left_length + 1];
		folded_structure	= new char[left_length + 1];

		strcpy(tmp_folded_structure, left_align);//candidateStems[stemIdx].leftArm[rank].alignment);
		
		//we need to build its left folded structure 
		//based on its left alignment and right sequence
		for(n=0; n<right_length; n++)
		{
			//from right arm
			//find all the '-'. count the 'M' for every '-'
			if((right_seq[right_length-1-n] != '-') && (right_align[right_length-1-n] == this->pStem[stemIdx].charid))//'M'))	
			{
				//count the num of matching
				m_num++;
			} 
			else
			{
				if(right_seq[right_length-1-n] == '-')
				{
					m_num++;	//including '-'
					//go to the left align and put '.' there
					m = 0;
					count = 0;
					while(count<m_num) {
						if(left_align[m] == this->pStem[stemIdx].charid)//'M')
							count++;
						m++;
					}
					tmp_folded_structure[m-1] = '.';
				}
			}
		}//for
		m = 0;
		for(n=0; n<left_length; n++)
		{
			if(left_seq[n] != '-') {
				folded_structure[m++] = tmp_folded_structure[n];
			}
		}
		folded_structure[m] = '\0';

		if(tmp_folded_structure != NULL)
		{
			delete [] tmp_folded_structure;
			tmp_folded_structure = NULL;
		}
	}
	else	//right arm
	{
		tmp_folded_structure= new char[right_length + 1];
		folded_structure	= new char[right_length + 1];
		strcpy(tmp_folded_structure, right_align);//candidateStems[stemIdx].rightArm[rank].alignment);

		//we need to build its right folded structure 
		//based on its right alignment and left sequence
		for(n=0; n<left_length; n++)
		{
			//from right arm
			//find all the '-'. count the 'M' for every '-'
			if((left_seq[n] != '-') && (left_align[n] == this->pStem[stemIdx].charid))//'M'))	
			{
				//count the num of matching
				m_num++;
			} 
			else
			{
				if(left_seq[n] == '-')
				{
					m_num++;	//including '-'
					//go to the left align and put '.' there
					m = right_length - 1;
					count = 0;
					while(count<m_num) {
						if(right_align[m] == this->pStem[stemIdx].charid)//'M')
							count++;
						m--;
					}
					tmp_folded_structure[m+1] = '.';
				}
			}
		}//for
		m = 0;
		for(n=0; n<right_length; n++)
		{
			if(right_seq[n] != '-') {
				folded_structure[m++] = tmp_folded_structure[n];
			}
		}
		folded_structure[m] = '\0';

		if(tmp_folded_structure != NULL)
		{
			delete [] tmp_folded_structure;
			tmp_folded_structure = NULL;
		}
	}
	//free the memory
	if(left_seq != NULL) {
		delete [] left_seq;
		left_seq = NULL;
	}
	if(left_align != NULL) {
		delete [] left_align;
		left_align = NULL;
	}
	if(right_seq != NULL) {
		delete [] right_seq;
		right_seq = NULL;
	}
	if(right_align != NULL) {
		delete [] right_align;
		right_align = NULL;
	}

	return folded_structure;
}

void DynamicP::printScoreForEveryComponent(Loop *pLoopModel, bool bDebug)
{
	int m;
	int stemIdx, preStemIdx;
	int	stemnum = this->getNumstems();

	int preLeftStemPos = 0, preRightStemPos = 0;
	int curLeftStemPos = 0, curRightStemPos = 0;
	int pre_gid, cur_gid;
	pre_gid = cur_gid = 0;

	VTBSearch vtber(1);
	vtber.setAllowedInsNum(this->getAllowedNullLoopInsNum());
	int loopid, start, end;

	for(m=0; m<2*stemnum; m++)
	{
		cur_gid = stemIdxArray[m];
		stemIdx = (cur_gid - 1)/2;

		if(cur_gid % 2 != 0) 
		{	
			//odd - left_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b;
			//20081021 check the score consistency
//			cout<<"left:   stem["<<stemIdx<<"]'s "<<stemCandIdxArray[m]<<"-th candidate pos info ["<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].a
//				<<", "<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].b
//				<<"]"<<endl;

			if(stemIdx == 0) {
				cout<<"Score for Stem "<<this->pStem[stemIdx].charid<<" : "<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<endl;
			} else {
                loopid	= this->identifyLoopId(pre_gid, cur_gid);
				preStemIdx = (pre_gid - 1)/2;
				if((preRightStemPos+1) <= (curLeftStemPos-1)+1+this->getNumOfNtsOverlapBetweenStem())
				{
					//add loop structure alignment info here.
					start	= preRightStemPos + 1;
					end		= curLeftStemPos - 1;

					double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, NULL, NULL);
					cout<<"Score for Loop ("<<this->pStem[preStemIdx].charid<<"->"<<this->pStem[stemIdx].charid<<"): "<<dLoopScore<<endl;
				} else {
					cout<<"Something wrong in computing the score"<<endl;
				}
				//add the next stem sequence/alignment info
				cout<<"Score for Stem "<<this->pStem[stemIdx].charid<<" : "<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<endl;
			}
		} 
		else 
		{	
			//even - right_arm
			curLeftStemPos	= candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c;
			curRightStemPos = candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d;
			
			//20081021 check the score consistency
//			cout<<"right:  stem["<<stemIdx<<"]'s "<<stemCandIdxArray[m]<<"-th candidate pos info ["<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].c
//				<<", "<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].d
//				<<"]"<<endl;

            loopid	= this->identifyLoopId(pre_gid, cur_gid);
			preStemIdx = (pre_gid - 1)/2;
			if((preRightStemPos+1) <= (curLeftStemPos-1)+1+this->getNumOfNtsOverlapBetweenStem())
			{
				//add loop structure alignment info here.
				start	= preRightStemPos + 1;
				end		= curLeftStemPos - 1;
				double dLoopScore = vtber.VTBSearchMax(&pLoopModel[loopid], cWindowSequence, start, end, NULL, NULL);
				cout<<"Score for Loop ("<<this->pStem[preStemIdx].charid<<"->"<<this->pStem[stemIdx].charid<<"): "<<dLoopScore<<endl;
				//add the next stem sequence/alignment info
//				cout<<"Score for Right Arm of Stem "<<this->pStem[stemIdx].charid<<":"<<candidateStems[stemIdx].s_locinfo[stemCandIdxArray[m]].p_score<<endl;
			} else {
				cout<<"Something wrong in computing the score"<<endl;
			}
		}
		preLeftStemPos	= curLeftStemPos;
		preRightStemPos	= curRightStemPos;
		
		pre_gid	= cur_gid;
	}//for
}

int DynamicP::getJumpStrategy()
{
	return this->iJumpStrategy;
}
void DynamicP::setJumpStrategy(int strategy)
{
	this->iJumpStrategy = strategy;
}

int DynamicP::getStepSize()
{
	return this->iStepSize;
}
void DynamicP::setStepSize(int stepsize)
{
	this->iStepSize = stepsize;
}
double DynamicP::getCpuTimeAll()
{
	return this->cpu_time_hours_all;
}
double DynamicP::getCpuTimePreprocessing()
{
	return this->cpu_time_hours_preprocessing;
}
double DynamicP::getCpuTimeDP()
{
	return this->cpu_time_hours_dp;
}
double DynamicP::getScoreThresInJump()
{
	return this->dScoreThresInJump;
}
void DynamicP::setScoreThresInJump(double dScoreThres)
{
	this->dScoreThresInJump = dScoreThres;
}
//flag of searching reverse genome sequence
int	DynamicP::getFlagSearchReverseGenome()
{
	return this->iFlagSearchReverseGenome;
}

void DynamicP::setFlagSearchReverseGenome(int flag)
{
	this->iFlagSearchReverseGenome = flag;
//	if(flag >= 1)
//		this->iFlagSearchReverseGenome = 1;
//	else
//		this->iFlagSearchReverseGenome = 0;
}

void DynamicP::setSearchParams(int		top_k,
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
							   int		iFlagPrintScoreInfo,
							   int		iFlagPrintDebugInfo,
							   int		iFlagPrintDebugInputInfo,
							   int		iFlagPrintDebugTreeInfo,
							   int		iFlagPrintDebugDPWinInfo,
							   int		iFlagPrintDebugTDInfo,
							   int		iFlagPrintDebugDPSearchInfo)
{
	this->setLeadingNum(top_k);	//[k value]
	this->setThreshold(threshold);	//[threshold]
	this->setNumOfNtsOverlapBetweenStem(num_nts_overlap_between_stem);	//[number of nts in stem when overlap is allowed]
	this->setFlagMergeCandInPreprocess(iFlagMergeCandInPreprocess);	//whether taking the merge-candidate strategy in preprocessing
	this->setFlagCandwithShortestLength(iFlagCandwithShortestLength);
	this->setShiftNumMergeCand(iShiftNumMergeCand);		//[num of shift when merge candidate in preprocessing]
	this->setAllowedNullLoopInsNum(iAllowedNullLoopInsNum);//[number of insertion allowed in null loop]
	this->setPcoeff(pcoeff);		//pcoeff
	this->setJumpStrategy(iJumpStrategy);	//[skip-and-jump strategy]
	this->setStepSize(iStepSize);			//[stepsize in skip-and-jump strategy]
	this->setScoreThresInJump(dScoreThresholdInJump);			//[score_filtering threshold in jump strategy]
	this->setFlagPrintStrAlignInfo(iFlagPrintStrAlignInfo);	//[structure alignment info -s|-n]
	this->setFlagPrintScoreInfo(iFlagPrintScoreInfo);	//
	this->setFlagPrintDebugInfo(iFlagPrintDebugInfo);			//[debug info -n|-d]
	this->setFlagPrintDebugInputInfo(iFlagPrintDebugInputInfo);
	this->setFlagPrintDebugTreeInfo(iFlagPrintDebugTreeInfo);
	this->setFlagPrintDebugDPWinInfo(iFlagPrintDebugDPWinInfo);
	this->setFlagPrintDebugTDInfo(iFlagPrintDebugTDInfo);
	this->setFlagPrintDebugDPSearchInfo(iFlagPrintDebugDPSearchInfo);
}

int	DynamicP::getNumHit()
{
	return this->numhit;
}

int	DynamicP::getPreHitNum()
{
	return this->pre_hit_num;
}

void DynamicP::setPreHitNum(int num)
{
	this->pre_hit_num = num;
}

int	DynamicP::getNumPastaLine()
{
	return this->num_pastaline;
}
void DynamicP::setNumPastaLine(int num)
{
	this->num_pastaline = num;
}
