#pragma once 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#pragma warning (disable:4996)

#define MAP_PAGE_SIZE (1<<5)
#define MMTagger_Model 0
#define MultiTag_Model 1

extern float gdmin_maf;
extern int gnwindow_len;
extern float gpmin_r2[3];
extern int gnmax_len;
extern int gnmax_merge_window_len;
extern int gnmax_covered_times;
extern int gnmem_size;
extern double gdmax_tagSNP_num;
extern int gnmin_bin_size;
extern bool gbno_equv_tag_SNPs;
extern int gnmodel;

extern double gdavg_cand_exts;
extern double gdavg_cand_RHS_exts;
extern int gnmax_cand_exts;
extern int gnmax_cand_RHS_exts;
extern double gdmerge_avg_cand_exts;
extern double gdmerge_avg_cand_RHS_exts;
extern int gnmerge_max_cand_exts;
extern int gnmerge_max_cand_RHS_exts;


void FastTagger(char* szmaf_filename, char* szmatrix_filename, char *szoutput_filename);


struct RHS_INFO_NODE 
{
	int nSNPid;
	float dR2;
	int ncomb_bitmap;
	RHS_INFO_NODE *pnext;
};


struct RHS_NODE  //for rule generation
{
	int nSNPid;
	RHS_NODE *pnext;
};

struct LHS_NODE
{
	int *pLHS_set;
	LHS_NODE *pnext;
};

struct MAP_LHS_NODE_PAGE
{
	LHS_NODE pnode_array[MAP_PAGE_SIZE];
	MAP_LHS_NODE_PAGE *pnext;
};

struct MAP_LHS_NODE_BUF
{
	MAP_LHS_NODE_PAGE *phead;
	MAP_LHS_NODE_PAGE *pcur_page;
	int num_of_pages;
	int ncur_pos;
};

struct LHS_NODE_MAP
{
	LHS_NODE **pLHS_nodes;
	MAP_LHS_NODE_BUF the_node_buf;
	LHS_NODE_MAP *pnext;
};

struct RULE
{
	int nLHS_len;
	int *pLHS_set;
	int nRHS_len;
	union
	{
		int *pRHS_set;
		RULE* pnext;
	};
};

struct RULE_NODE
{
	RULE *prule;
	RULE_NODE *pnext;
};


struct SNP
{
	int norder;
	int npos;
	int nwindow_start_pos;
	int nwindow_end_pos;

	union
	{
		char *pbitmap;
		char *pcovered_flags;
	};
	union
	{
		int n1;
		int nflag; //bit0: 1: tag, 0, non-tag; bit1: 1: covered, 0: not covered
	};
	
	union
	{
		int num_of_merged_ids;
		int nSNPid_merged_with;
	};
	int *pmerged_ids;

	union
	{
		int *p1sid_list;
		RULE_NODE *prule_node_head;
	};

	int nbeing_covered_times;
	int num_of_covered_SNPs;

	//int nfile_no;
	bool bhas_rules;
	unsigned int nrule_start_pos;
	unsigned int nrule_end_pos;
	unsigned int nrule_mem_size_start;
	unsigned int nrule_mem_size_end;
};

struct REP_SNP
{
	int nid;
	int nstart_pos;
	int nend_pos;
	int nleft_boundary;
	RHS_NODE *pRHS_node_head;
	LHS_NODE_MAP *pLHS_node_map;
	int nbeing_covered_times;
	
};

//---------------------------------------------------------------------
extern int gnhapmap_num;
extern int gnmax_n1;

extern int *gporig_window_start_pos;
extern int *gporig_window_end_pos;
extern int gnwindow_size_sum;

extern int gnrep_SNP_num;
extern int gnremoved_SNP_num;

extern int gnnonb_SNP_num;

extern int gnmax_merged_ids;

int LoadSNPs(char* szmaf_filename, double dmin_maf, SNP* &pSNPs, char* szoutput_name);
int GetHapMapNum(char *szhapmap_filename);
void CalcWindowSize(SNP *pSNPs, int num_of_SNPs, int nwindow_len);
int MergeSNPs(char* szmatrix_filename, SNP *pSNPs, int num_of_SNPs, REP_SNP* &pRep_SNPs, char* szmerged_ids_filename);
void LoadBitmaps(FILE* fp, int &nline_no, SNP *pSNPs, REP_SNP *pRep_SNPs, int ncur_pos, int &ncur_load_pos, int &ncur_null_pos);
void EndShift(SNP* pSNPs, REP_SNP *pRep_SNPs, int num_of_rep_SNPs, int nleft_boundary);

inline bool IsIdentical(char* pbitmap1, char* pbitmap2)
{
	int i;

	for(i=0;i<gnhapmap_num;i++)
	{
		if(pbitmap1[i]!=pbitmap2[i])
			return false;
	}
	return true;
}

inline bool IsReverseIdentical(char* pbitmap1, char* pbitmap2)
{
	int i;
	for(i=0;i<gnhapmap_num;i++)
	{
		if(pbitmap1[i]==pbitmap2[i])
			return false;
	}
	return true;
}

//---------------------------------------------------------------------


extern int *gpmerged_bins;

extern int *gpmerged_SNP_ids;
extern int gnmerged_id_buf_pos;

extern int gnmax_RHS_len;

extern int *gpLHS_num;
extern int *gprule_num;
extern int *gpmerged_bin_rule_num;
extern int gnmerged_bin_rule_num;
extern int gnsingleton_rules;
extern int gnmultiSNP_rules;
extern int gnmultiSNP_LHSs;

extern int gntotal_calls;
extern double gdsubset_check_times;
extern double gdsubset_check_depth_sum;

extern FILE *gfp_rule;
extern unsigned int gnrule_file_pos;
extern int gntotal_rule_size;
extern int gnfile_no;

extern RULE gorule;

extern int gnum_of_partitions;
extern int gnum_of_iterations;
extern int gnorig_num_of_tag_SNPs;
extern int gnnot_covered_tag_SNPs;

extern int gnum_of_tag_SNPs;
extern int gncovered_SNPs;

extern int gninvalid_rule_num;

extern double gdverify_time;

void GenRules(char* szmatrix_filename, SNP *pSNPs, REP_SNP* prepresent_SNPs, int num_of_represent_SNPs, char* szoutput_name);

void OutputRule(SNP *pSNPs, int nSNP_id1, int nSNP_id2, float dR2, int ncomb_bitmap);
void OutputRule(SNP *pSNPs, int nLHS_len, int *pLHS_set, int nRHS_len, RHS_INFO_NODE *pRHS_set);
RULE* ReadInOneRule(unsigned int nrule_end_pos);
void GetRuleMemSize(SNP *pSNPs, int num_of_SNPs);

void VerifyRule(REP_SNP* pRep_SNPs, int nindex1, int nindex2);
void VerifyRule(REP_SNP *pRep_SNPs, int *pLHS_set, int nLHS_len, int *pRHS_set, int nRHS_len);

void OutputTagSNPs(SNP* pSNPs, int num_of_SNPs, int *ptagSNPs, int num_of_tagSNPs, char* szoutput_name);
void SelectTagSNPs(SNP *pSNPs, int num_of_SNPs, int nmem_size, char* szoutput_name);
int RemoveDupTagSNPs(SNP* pSNPs, int *ptagSNPs, int num_of_tagSNPs);

void Verify(char *szmaf_filename, char *szmatrix_filename, char *szoutput_name);

inline float CalcR2(int n11, int n1x, int nx1, int nhapmap_num, int &ncomb_bitmap)
{
	double dR2;

	dR2 = n11*nhapmap_num-n1x*nx1;
	if(dR2>0)
		ncomb_bitmap = 1;
	else 
		ncomb_bitmap = 0;

	if(dR2==0)
		return 0;
	else 
	{
		dR2 = dR2*dR2;
		dR2 = dR2/(n1x*(gnhapmap_num-n1x)*nx1*(gnhapmap_num-nx1));

		return (float)dR2;
	}	
}
float CalcR2(char* pbitmap1, int n1x, int *p1sid_list, int nx1, int &ncomb_bitmap);
float CalcMultiMarkerR2(char* pbitmap1, int nSNP_num, int* p1sid_list, int nx1, int *pprefix_counter, int &ncomb_bitmap);


//==========================================================
extern int gnused_mem_size;
extern int gnmax_used_mem_size;

inline void IncMemSize(int nsize)
{
	gnused_mem_size += nsize;
	if(gnmax_used_mem_size<gnused_mem_size)
		gnmax_used_mem_size = gnused_mem_size;
}

inline void DecMemSize(int nsize)
{
	gnused_mem_size -= nsize;
}
//===========================================================================




#define PAGE_SIZE (1<<10)


//=======================================================================
inline char* NewBitmap()
{
	IncMemSize(sizeof(char)*gnhapmap_num);
	return new char[gnhapmap_num];
}
inline void DelBitmap(char* pbitmap)
{
	DecMemSize(sizeof(char)*gnhapmap_num);
	delete []pbitmap;
}
inline int* New1SidList()
{
	IncMemSize(sizeof(int)*gnmax_n1);
	return new int[gnmax_n1];
}
inline void Del1SidList(int *p1sid_list)
{
	DecMemSize(sizeof(int)*gnmax_n1);
	delete []p1sid_list;
}


//=======================================================================

//----------------------------------------------------------------------------
//used during rule generation

struct RHS_NODE_PAGE
{
	int nmax_SNP_id;
	RHS_NODE pnode_array[PAGE_SIZE];
	RHS_NODE_PAGE *pnext;
};

struct RHS_NODE_BUF
{
	RHS_NODE_PAGE *phead;
	RHS_NODE_PAGE *pcur_page;
	int num_of_pages;
	int ncur_pos;
};

extern RHS_NODE_BUF goRHSnode_buf;

inline RHS_NODE* NewRHSNode(int nRHS_SNPid)
{
	RHS_NODE_PAGE *pnew_page;
	RHS_NODE *pRHS_node;

	if(goRHSnode_buf.ncur_pos==PAGE_SIZE)
	{
		if(goRHSnode_buf.pcur_page->pnext!=NULL)
			goRHSnode_buf.pcur_page = goRHSnode_buf.pcur_page->pnext;
		else
		{
			pnew_page = new RHS_NODE_PAGE;
			IncMemSize(sizeof(RHS_NODE_PAGE));
			pnew_page->pnext = NULL;
			goRHSnode_buf.pcur_page->pnext = pnew_page;
			goRHSnode_buf.pcur_page = pnew_page;
			goRHSnode_buf.num_of_pages++;
		}
		goRHSnode_buf.ncur_pos = 0;
		goRHSnode_buf.pcur_page->nmax_SNP_id = nRHS_SNPid;
	}
	else if(goRHSnode_buf.pcur_page->nmax_SNP_id<nRHS_SNPid)
		goRHSnode_buf.pcur_page->nmax_SNP_id = nRHS_SNPid;

	pRHS_node = &(goRHSnode_buf.pcur_page->pnode_array[goRHSnode_buf.ncur_pos]);
	goRHSnode_buf.ncur_pos++;

	pRHS_node->nSNPid = nRHS_SNPid;

	return pRHS_node;
}


inline void InitRHSNodeBuf()
{
	goRHSnode_buf.phead = new RHS_NODE_PAGE;
	IncMemSize(sizeof(RHS_NODE_PAGE));
	goRHSnode_buf.phead->pnext = NULL;
	goRHSnode_buf.pcur_page = goRHSnode_buf.phead;
	goRHSnode_buf.num_of_pages = 1;
	goRHSnode_buf.ncur_pos = 0;
	goRHSnode_buf.pcur_page->nmax_SNP_id = 0;
}

inline void DelRHSNodeBuf()
{
	RHS_NODE_PAGE *pRHS_page;

	while(goRHSnode_buf.phead!=NULL)
	{
		pRHS_page = goRHSnode_buf.phead;
		goRHSnode_buf.phead = goRHSnode_buf.phead->pnext;
		delete pRHS_page;
		DecMemSize(sizeof(RHS_NODE_PAGE));
		goRHSnode_buf.num_of_pages--;
	}
	if(goRHSnode_buf.num_of_pages!=0)
		printf("Error: the number of pages should be 0\n");

}
//---------------------------------------


//=======================================
struct SNPSET_PAGE
{
	int nmax_SNP_id;
	int parray[PAGE_SIZE];
	SNPSET_PAGE *pnext;
};

struct SNPSET_BUF
{
	SNPSET_PAGE *phead;
	SNPSET_PAGE *pcur_page;
	int num_of_pages;
	int ncur_pos;
};

extern SNPSET_BUF goSNPset_buf;

inline int* NewSNPSet(int nlen, int nmax_SNP_id)
{
	SNPSET_PAGE *pnew_page;
	int *pset;

	if(goSNPset_buf.ncur_pos+nlen>PAGE_SIZE)
	{
		if(goSNPset_buf.pcur_page->pnext!=NULL)
			goSNPset_buf.pcur_page = goSNPset_buf.pcur_page->pnext;
		else
		{
			pnew_page = new SNPSET_PAGE;
			IncMemSize(sizeof(SNPSET_PAGE));
			pnew_page->pnext = NULL;
			goSNPset_buf.pcur_page->pnext = pnew_page;
			goSNPset_buf.pcur_page = pnew_page;
			goSNPset_buf.num_of_pages++;
		}
		goSNPset_buf.ncur_pos = 0;
		goSNPset_buf.pcur_page->nmax_SNP_id = nmax_SNP_id;
	}
	else if(goSNPset_buf.pcur_page->nmax_SNP_id<nmax_SNP_id)
		goSNPset_buf.pcur_page->nmax_SNP_id = nmax_SNP_id;

	pset = &goSNPset_buf.pcur_page->parray[goSNPset_buf.ncur_pos];
	goSNPset_buf.ncur_pos += nlen;

	return pset;
}

inline void InitSNPsetBuf()
{
	goSNPset_buf.phead = new SNPSET_PAGE;
	IncMemSize(sizeof(SNPSET_PAGE));
	goSNPset_buf.phead->pnext = NULL;
	goSNPset_buf.pcur_page = goSNPset_buf.phead;
	goSNPset_buf.num_of_pages = 1;
	goSNPset_buf.ncur_pos = 0;
	goSNPset_buf.pcur_page->nmax_SNP_id = 0;
}

inline void DelSNPsetBuf()
{
	SNPSET_PAGE *pSNPset_page;

	while(goSNPset_buf.phead!=NULL)
	{
		pSNPset_page = goSNPset_buf.phead;
		goSNPset_buf.phead = goSNPset_buf.phead->pnext;
		delete pSNPset_page;
		DecMemSize(sizeof(SNPSET_PAGE));
		goSNPset_buf.num_of_pages--;
	}
	if(goSNPset_buf.num_of_pages!=0)
		printf("Error: the number of pages should be 0\n");
}


//=======================================================


//=====================================================
//for subset checking
extern LHS_NODE_MAP *gpfree_LHS_node_maps;

inline void DelLHSNodeMap(LHS_NODE_MAP* pLHS_node_map)
{
	pLHS_node_map->pnext = gpfree_LHS_node_maps;
	gpfree_LHS_node_maps = pLHS_node_map;
}

inline LHS_NODE_MAP* NewLHSNodeMap()
{
	LHS_NODE_MAP *pLHS_node_map;

	if(gpfree_LHS_node_maps!=NULL)
	{
		pLHS_node_map = gpfree_LHS_node_maps;
		gpfree_LHS_node_maps = gpfree_LHS_node_maps->pnext;
		memset(pLHS_node_map->pLHS_nodes, 0, sizeof(LHS_NODE*)*gnmerge_max_cand_RHS_exts);
		pLHS_node_map->the_node_buf.pcur_page = pLHS_node_map->the_node_buf.phead;
		pLHS_node_map->the_node_buf.ncur_pos = 0;
	}
	else
	{
		pLHS_node_map = new LHS_NODE_MAP;
		IncMemSize(sizeof(LHS_NODE_MAP));
		pLHS_node_map->pLHS_nodes = new LHS_NODE*[gnmerge_max_cand_RHS_exts];
		IncMemSize(sizeof(LHS_NODE*)*gnmerge_max_cand_RHS_exts);
		memset(pLHS_node_map->pLHS_nodes, 0, sizeof(LHS_NODE*)*gnmerge_max_cand_RHS_exts);
		pLHS_node_map->the_node_buf.phead = new MAP_LHS_NODE_PAGE;
		IncMemSize(sizeof(MAP_LHS_NODE_PAGE));
		pLHS_node_map->the_node_buf.phead->pnext = NULL;
		pLHS_node_map->the_node_buf.pcur_page = pLHS_node_map->the_node_buf.phead;
		pLHS_node_map->the_node_buf.num_of_pages = 1;
		pLHS_node_map->the_node_buf.ncur_pos = 0;
	}
	pLHS_node_map->pnext = NULL;

	return pLHS_node_map;
}

inline void DelLHSNodeMaps()
{
	LHS_NODE_MAP *pLHS_node_map;
	MAP_LHS_NODE_PAGE *pLHS_node_page;

	while(gpfree_LHS_node_maps!=NULL)
	{
		pLHS_node_map = gpfree_LHS_node_maps;
		gpfree_LHS_node_maps = gpfree_LHS_node_maps->pnext;
		delete []pLHS_node_map->pLHS_nodes;
		DecMemSize(sizeof(RULE_NODE*)*gnmerge_max_cand_RHS_exts);
		
		while(pLHS_node_map->the_node_buf.phead!=NULL)
		{
			pLHS_node_page = pLHS_node_map->the_node_buf.phead;
			pLHS_node_map->the_node_buf.phead = pLHS_node_map->the_node_buf.phead->pnext;
			delete pLHS_node_page;
			DecMemSize(sizeof(MAP_LHS_NODE_PAGE));
			pLHS_node_map->the_node_buf.num_of_pages--;
		}
		if(pLHS_node_map->the_node_buf.num_of_pages!=0)
			printf("Error: the number of pages should be 0\n");

		delete pLHS_node_map;
		DecMemSize(sizeof(LHS_NODE_MAP));
	}
}

inline LHS_NODE* NewLHSNode(MAP_LHS_NODE_BUF *pLHS_node_buf)
{
	MAP_LHS_NODE_PAGE *pnew_page;
	LHS_NODE *pLHS_node;

	if(pLHS_node_buf->ncur_pos==MAP_PAGE_SIZE)
	{
		if(pLHS_node_buf->pcur_page->pnext==NULL)
		{
			pnew_page = new MAP_LHS_NODE_PAGE;
			IncMemSize(sizeof(MAP_LHS_NODE_PAGE));
			pnew_page->pnext = NULL;
			pLHS_node_buf->pcur_page->pnext = pnew_page;
			pLHS_node_buf->pcur_page = pnew_page;
			pLHS_node_buf->num_of_pages++;
		}
		else
			pLHS_node_buf->pcur_page = pLHS_node_buf->pcur_page->pnext;
		pLHS_node_buf->ncur_pos = 0;
	}

	pLHS_node = &(pLHS_node_buf->pcur_page->pnode_array[pLHS_node_buf->ncur_pos]);
	pLHS_node_buf->ncur_pos++;

	return pLHS_node;
}

inline void InsertLHSNode(REP_SNP *pRep_SNPs, int nRHS_index, int nmap_start_pos, int *pLHS_set)
{
	LHS_NODE_MAP *pLHS_node_map;
	LHS_NODE *pLHS_node;
	int nlast_index;

	if(pRep_SNPs[nRHS_index].pLHS_node_map==NULL)
		pRep_SNPs[nRHS_index].pLHS_node_map = NewLHSNodeMap();

	pLHS_node_map = pRep_SNPs[nRHS_index].pLHS_node_map;

	pLHS_node = NewLHSNode(&pLHS_node_map->the_node_buf);
	pLHS_node->pLHS_set = pLHS_set;

	nlast_index = pLHS_set[pLHS_set[0]];

	pLHS_node->pnext = pLHS_node_map->pLHS_nodes[nlast_index-nmap_start_pos];
	pLHS_node_map->pLHS_nodes[nlast_index-nmap_start_pos] = pLHS_node;
}

//=====================================================




//=============================================================
//for tag SNP selection
struct SNP_ID_NODE
{
	int nid;
	int nflag;
	SNP_ID_NODE *pprev;
	SNP_ID_NODE *pnext;
};

extern SNP_ID_NODE *gpSNP_ID_nodes;
extern SNP_ID_NODE **gpSNP_ID_heads;

inline void MoveSNPIDNode(SNP *pSNPs, int nid, int noffset, int ninc, int ncovered_id)
{
	SNP_ID_NODE *pSNP_ID_node;


	pSNP_ID_node = &gpSNP_ID_nodes[nid-noffset];

	if(pSNP_ID_node->nflag==-1)
		return;

	if(pSNP_ID_node->pprev!=NULL)
		pSNP_ID_node->pprev->pnext = pSNP_ID_node->pnext;
	else 
		gpSNP_ID_heads[pSNPs[nid].num_of_covered_SNPs] = pSNP_ID_node->pnext;
	if(pSNP_ID_node->pnext!=NULL)
		pSNP_ID_node->pnext->pprev = pSNP_ID_node->pprev;

	pSNPs[nid].num_of_covered_SNPs += ninc;

	pSNP_ID_node->pprev = NULL;
	pSNP_ID_node->pnext = gpSNP_ID_heads[pSNPs[nid].num_of_covered_SNPs];
	if(gpSNP_ID_heads[pSNPs[nid].num_of_covered_SNPs]!=NULL)
		gpSNP_ID_heads[pSNPs[nid].num_of_covered_SNPs]->pprev = pSNP_ID_node;
	gpSNP_ID_heads[pSNPs[nid].num_of_covered_SNPs] = pSNP_ID_node;

}

inline void RemoveSNPIDNode(SNP* pSNPs, int nid, int noffset)
{
	SNP_ID_NODE *pSNP_ID_node;


	pSNP_ID_node = &gpSNP_ID_nodes[nid-noffset];

	if(pSNP_ID_node->nflag==-1)
		return;

	if(pSNP_ID_node->pprev!=NULL)
		pSNP_ID_node->pprev->pnext = pSNP_ID_node->pnext;
	else 
	{
		if(gpSNP_ID_heads[pSNPs[nid].num_of_covered_SNPs]!=pSNP_ID_node)
			printf("Error: the head shuold point to this SNP_ID_NODE\n");
		gpSNP_ID_heads[pSNPs[nid].num_of_covered_SNPs] = pSNP_ID_node->pnext;
	}
	if(pSNP_ID_node->pnext!=NULL)
		pSNP_ID_node->pnext->pprev = pSNP_ID_node->pprev;

	pSNP_ID_node->pprev = NULL;
	pSNP_ID_node->pnext = NULL;
	pSNP_ID_node->nflag = -1;

}

//=============================================================



//==================================================
struct RULE_PAGE
{
	int nmax_SNP_id;
	RULE prule_array[PAGE_SIZE];
	RULE_PAGE *pnext;
};

struct RULE_BUF
{
	RULE_PAGE *phead;
	RULE_PAGE *pcur_page;
	int num_of_pages;
	int ncur_pos;
};

extern RULE_BUF gorule_buf;
extern RULE *gpfree_rules;

inline RULE* NewRule(int nLHS_len, int *pLHS_set, int nRHS_len, int *pRHS_set) //called during tag SNP selection
{
	RULE_PAGE *pnew_page;
	RULE *prule;

	if(gpfree_rules!=NULL)
	{
		prule = gpfree_rules;
		gpfree_rules = gpfree_rules->pnext;		
	}
	else
	{
		if(gorule_buf.ncur_pos==PAGE_SIZE)
		{
			pnew_page = new RULE_PAGE;
			IncMemSize(sizeof(RULE_PAGE));
			pnew_page->pnext = NULL;
			gorule_buf.pcur_page->pnext = pnew_page;
			gorule_buf.pcur_page = pnew_page;
			gorule_buf.num_of_pages++;
			gorule_buf.ncur_pos = 0;
			gorule_buf.pcur_page->nmax_SNP_id = 0;
		}
		prule = &(gorule_buf.pcur_page->prule_array[gorule_buf.ncur_pos]);
		gorule_buf.ncur_pos++;
	}

	prule->nLHS_len = nLHS_len;
	prule->pLHS_set = NewSNPSet(nLHS_len, 0);
	memcpy(prule->pLHS_set, pLHS_set, sizeof(int)*nLHS_len);

	prule->nRHS_len = nRHS_len;
	prule->pRHS_set = NewSNPSet(nRHS_len, 0);
	memcpy(prule->pRHS_set, pRHS_set, sizeof(int)*nRHS_len);

	return prule;
}

inline void DelRule(RULE *prule)
{
	prule->pnext = gpfree_rules;
	gpfree_rules = prule;
}

inline void InitRuleBuf()
{
	gorule_buf.phead = new RULE_PAGE;
	IncMemSize(sizeof(RULE_PAGE));
	gorule_buf.phead->pnext = NULL;
	gorule_buf.pcur_page = gorule_buf.phead;
	gorule_buf.num_of_pages = 1;
	gorule_buf.ncur_pos = 0;
	gorule_buf.pcur_page->nmax_SNP_id = 0;

	goSNPset_buf.phead = new SNPSET_PAGE;
	IncMemSize(sizeof(SNPSET_PAGE));
	goSNPset_buf.phead->pnext = NULL;
	goSNPset_buf.pcur_page = goSNPset_buf.phead;
	goSNPset_buf.num_of_pages = 1;
	goSNPset_buf.ncur_pos = 0;
	goSNPset_buf.pcur_page->nmax_SNP_id = 0;

	gpfree_rules = NULL;
}

inline void DelRuleBuf()
{
	RULE_PAGE *prule_page;
	SNPSET_PAGE *pSNPset_page;

	while(gorule_buf.phead!=NULL)
	{
		prule_page = gorule_buf.phead;
		gorule_buf.phead = gorule_buf.phead->pnext;
		delete prule_page;
		DecMemSize(sizeof(RULE_PAGE));
		gorule_buf.num_of_pages--;
	}
	if(gorule_buf.num_of_pages!=0)
		printf("Error: the number of pages should be 0\n");

	while(goSNPset_buf.phead!=NULL)
	{
		pSNPset_page = goSNPset_buf.phead;
		goSNPset_buf.phead = goSNPset_buf.phead->pnext;
		delete pSNPset_page;
		DecMemSize(sizeof(SNPSET_PAGE));
		goSNPset_buf.num_of_pages--;
	}
	if(goSNPset_buf.num_of_pages!=0)
		printf("Error: the number of pages should be 0\n");
}
//=========================================================


//------------------------------------------------------
struct RULE_NODE_PAGE
{
	RULE_NODE pnode_array[PAGE_SIZE];
	RULE_NODE_PAGE *pnext;
};

struct RULE_NODE_BUF
{
	RULE_NODE_PAGE *phead;
	RULE_NODE_PAGE *pcur_page;
	int num_of_pages;
	int ncur_pos;
};

extern RULE_NODE_BUF gorule_node_buf;
extern RULE_NODE *gpfree_rule_nodes;

inline RULE_NODE* NewRuleNode()
{
	RULE_NODE_PAGE *pnew_page;
	RULE_NODE *prule_node;

	if(gpfree_rule_nodes!=NULL)
	{
		prule_node = gpfree_rule_nodes;
		gpfree_rule_nodes = gpfree_rule_nodes->pnext;
	}
	else
	{
		if(gorule_node_buf.ncur_pos==PAGE_SIZE)
		{
			pnew_page = new RULE_NODE_PAGE;
			IncMemSize(sizeof(RULE_NODE_PAGE));
			pnew_page->pnext = NULL;
			gorule_node_buf.pcur_page->pnext = pnew_page;
			gorule_node_buf.pcur_page = pnew_page;
			gorule_node_buf.num_of_pages++;
			gorule_node_buf.ncur_pos = 0;
		}

		prule_node = &(gorule_node_buf.pcur_page->pnode_array[gorule_node_buf.ncur_pos]);
		gorule_node_buf.ncur_pos++;
	}
	return prule_node;
}

inline void DelRuleNode(RULE_NODE *prule_node)
{
	prule_node->pnext = gpfree_rule_nodes;
	gpfree_rule_nodes = prule_node;
}


inline void InitRuleNodeBuf()
{
	gorule_node_buf.phead = new RULE_NODE_PAGE;
	IncMemSize(sizeof(RULE_NODE_PAGE));
	gorule_node_buf.phead->pnext = NULL;
	gorule_node_buf.pcur_page = gorule_node_buf.phead;
	gorule_node_buf.num_of_pages = 1;
	gorule_node_buf.ncur_pos = 0;

	gpfree_rule_nodes = NULL;
}

inline void DelRuleNodeBuf()
{
	RULE_NODE_PAGE *prule_node_page;

	while(gorule_node_buf.phead!=NULL)
	{
		prule_node_page = gorule_node_buf.phead;
		gorule_node_buf.phead = gorule_node_buf.phead->pnext;
		delete prule_node_page;
		DecMemSize(sizeof(RULE_NODE_PAGE));
		gorule_node_buf.num_of_pages--;
	}
	if(gorule_node_buf.num_of_pages!=0)
		printf("Error: the number of pages should be 0\n");
}
//------------------------------------------------------


