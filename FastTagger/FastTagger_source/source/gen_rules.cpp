#include "global.h"
#include <time.h>
#include <sys/timeb.h>


int gnprefix_len;
int *gpprefix;
int* gpDFS_cand_RHSs;
RHS_INFO_NODE *gprule_RHS_nodes;

int *gpLHS_indices;
int *gpRHS_indices;

char *gpDFS_bitmap_buf;

int *gpmerged_bins;
int *gpmerged_bin_rule_num;
int gnmerged_bin_rule_num;

int *gpmerged_SNP_ids;
int gnmerged_id_buf_pos;

RHS_NODE_BUF goRHSnode_buf;
SNPSET_BUF goSNPset_buf;
LHS_NODE_MAP *gpfree_LHS_node_maps;

int *gpLHS_num;
int *gprule_num;
int gnmax_RHS_len;

int gntotal_calls;

double gdsubset_check_times;
double gdsubset_check_depth_sum;


float CalcR2(char* pbitmap1, int n1x, int *p1sid_list, int nx1, int &ncomb_bitmap)
{
	int i, n11;
	float dR2;

	n11 = 0;
	for(i=0;i<nx1;i++)
	{
		if(pbitmap1[p1sid_list[i]]==1)
			n11++;
	}
	dR2 = CalcR2(n11, n1x, nx1, gnhapmap_num, ncomb_bitmap);

	return dR2;
}

float CalcMultiMarkerR2(char* pbitmap1, int nSNP_num, int* p1sid_list, int nx1, int *pprefix_counters, int &ncomb_bitmap)
{
	int ncounter_num, i, n11, n1x, nflag, num_of_0s, nmerged_1s, nmerged_0s, nmerged_bin_num;
	int nmaxR2_bin_no, nmaxR2_flag, pcounters[8];
	float dR2, dmax_R2;

	ncounter_num = 1<<nSNP_num;

	memset(pcounters, 0, sizeof(int)*ncounter_num);
	for(i=0;i<nx1;i++)
		pcounters[pbitmap1[p1sid_list[i]]]++;

	if(gnmodel==MMTagger_Model)
	{
		nmerged_1s = 0;
		nmerged_0s = 0;
		nmerged_bin_num = 0;
		n11 = 0;
		n1x = 0;
		ncomb_bitmap = 0;
		for(i=0;i<ncounter_num;i++)
		{
			if(pprefix_counters[i]>=gnmin_bin_size)
			{
				num_of_0s = pprefix_counters[i]-pcounters[i];
				if(pcounters[i]>num_of_0s || pcounters[i]==num_of_0s && nx1<=gnhapmap_num-nx1)
				{
					n11 += pcounters[i];
					n1x += pprefix_counters[i];
					ncomb_bitmap |= (1<<i);
				}
			}
			else
			{
				nmerged_1s += pcounters[i];
				nmerged_0s += pprefix_counters[i]-pcounters[i];
				gpmerged_bins[nmerged_bin_num++] = i;
			}
		}
		if(nmerged_bin_num>0)
		{
			if(nmerged_1s>nmerged_0s || nmerged_1s==nmerged_0s && nx1<=gnhapmap_num-nx1)
			{
				n11 += nmerged_1s;
				n1x += nmerged_1s+nmerged_0s;
				for(i=0;i<nmerged_bin_num;i++)
					ncomb_bitmap |= (1<<gpmerged_bins[i]);
			}
		}
		dR2 = CalcR2(n11, n1x, nx1, gnhapmap_num, nflag);
		if(dR2>0 && nflag==0)
			printf("Error: the correlation should be positive\n");

		if(dR2>=gpmin_r2[nSNP_num-1] && nmerged_bin_num>1)
			gpmerged_bin_rule_num[nSNP_num]++;
	}
	else
	{
		dmax_R2 = 0;
		for(i=0;i<ncounter_num;i++)
		{
			dR2 = CalcR2(pcounters[i], pprefix_counters[i], nx1, gnhapmap_num, nflag);
			if(dR2>dmax_R2)
			{
				dmax_R2 = dR2;
				nmaxR2_bin_no = i;
				nmaxR2_flag = nflag;
			}
		}

		dR2 = dmax_R2;
		ncomb_bitmap = 0;
		if(dmax_R2>=gpmin_r2[nSNP_num-1])
		{
			if(nmaxR2_flag==1) //positive association
				ncomb_bitmap |= (1<<nmaxR2_bin_no);
			else 
			{
				for(i=0;i<ncounter_num;i++)
				{
					if(i!=nmaxR2_bin_no)
						ncomb_bitmap |= (1<<i);
				}
			}
		}
	}

	return dR2;
}

void AddOneRule(REP_SNP *pRep_SNPs, int nindex1, int nindex2)
{
	RHS_NODE *pRHS_node;

	//VerifyRule(pRep_SNPs, nindex1, nindex2);

	pRHS_node = NewRHSNode(nindex1);
	pRHS_node->pnext = pRep_SNPs[nindex2].pRHS_node_head;
	pRep_SNPs[nindex2].pRHS_node_head = pRHS_node;

	pRHS_node = NewRHSNode(nindex2);
	pRHS_node->pnext = pRep_SNPs[nindex1].pRHS_node_head;
	pRep_SNPs[nindex1].pRHS_node_head = pRHS_node;

}

void AddOneRule(REP_SNP *pRep_SNPs, int *pLHS_set, int nLHS_len, int *pRHS_set, int nRHS_len)
{
	int *pmap_LHS_set;
	int i, nRHS_index, nmax_index;

	//VerifyRule(pRep_SNPs, pLHS_set, nLHS_len, pRHS_set, nRHS_len);

	nmax_index = pLHS_set[0];
	if(nmax_index<pRHS_set[nRHS_len-1])
		nmax_index = pRHS_set[nRHS_len-1];

	pmap_LHS_set = NewSNPSet(nLHS_len+1, nmax_index);
	pmap_LHS_set[0] = nLHS_len;
	memcpy(&pmap_LHS_set[1], pLHS_set, sizeof(int)*nLHS_len);
	for(i=0;i<nRHS_len;i++)
	{
		nRHS_index = pRHS_set[i];
		pRep_SNPs[nRHS_index].nbeing_covered_times++;
		if(gnmax_covered_times==0 || pRep_SNPs[nRHS_index].nbeing_covered_times<gnmax_covered_times)
			InsertLHSNode(pRep_SNPs, nRHS_index, pRep_SNPs[nRHS_index].nstart_pos, pmap_LHS_set);
	}
}

//true if set 1 is a subset of set2
bool IsSubset(int *pset1, int nlen1, int* pset2, int nlen2)
{
	int i, j;

	i = 0;
	j = 0;
	while(i<nlen1 && j<nlen2)
	{
		if(pset1[i]==pset2[j])
		{
			i++;
			j++;
		}
		else if(pset1[i] > pset2[j])
			return false;
		else 
			j++;
	}
	if(i<nlen1)
		return false;
	else
		return true;
}

bool IsCoveredBySubset(REP_SNP *pRep_SNPs, int nRHS_index, int *pLHS_set, int nLHS_len)
{
	LHS_NODE_MAP *pLHS_node_map;
	LHS_NODE *pLHS_node;
	int *pmap_LHS_set;
	bool bis_covered;
	int nlast_LHS_index;

	if(gnmax_covered_times==-1)
		return false;
	if(gnmax_covered_times>0 && pRep_SNPs[nRHS_index].nbeing_covered_times>=gnmax_covered_times)
		return true;

	gdsubset_check_times++;

	bis_covered = false;

	pLHS_node_map = pRep_SNPs[nRHS_index].pLHS_node_map;
	nlast_LHS_index = pLHS_set[nLHS_len-1];
	pLHS_node = pLHS_node_map->pLHS_nodes[nlast_LHS_index-pRep_SNPs[nRHS_index].nstart_pos];

	while(pLHS_node!=NULL)
	{
		gdsubset_check_depth_sum++;
		pmap_LHS_set = pLHS_node->pLHS_set;
		if(pmap_LHS_set[pmap_LHS_set[0]]!=nlast_LHS_index)
			printf("Error: the last LHS id should be %d\n", nlast_LHS_index);
		if(pmap_LHS_set[1]<pLHS_set[nLHS_len-2])
			break;
		else if(pmap_LHS_set[0]<nLHS_len)
		{
			if(IsSubset(&pmap_LHS_set[1], pmap_LHS_set[0], pLHS_set, nLHS_len))
			{
				bis_covered = true;
				break;
			}
		}
		pLHS_node = pLHS_node->pnext;
	}

	return bis_covered;
}


void DFSGenRules(char* pdfs_bitmap, int *pcand_RHSs, int num_of_cand_exts, int num_of_cand_RHSs, SNP *pSNPs, REP_SNP* pRep_SNPs, int nstart_pos, int nend_pos, int *pprefix_counters)
{
	int i, j, nnew_end_pos, nnew_start_pos, nRHS_len, ncomb_bitmap, ncand_ext_SNP_id, nid, nlast_RHS_pos;
	int *pnew_cand_RHSs, num_of_remain_exts, num_of_remain_RHSs, ncand_ext_pos, num_of_new_cand_exts, num_of_new_cand_RHSs;
	char *pbitmap;
	RHS_NODE *pRHS_node;
	float dR2;
	bool bis_covered;
	int p1counters[4], pcounters[8], ncounter_num;

	gntotal_calls++;

	pbitmap = &(gpDFS_bitmap_buf[gnprefix_len*gnhapmap_num]);

	nlast_RHS_pos = 0;

	for(i=0;i<num_of_cand_exts;i++)
	{
//if(gnprefix_len>=2 && gpprefix[0]==1041 && gpprefix[1]==1040 && i==3)
//printf("stop\n");

		pnew_cand_RHSs = &(gpDFS_cand_RHSs[gnprefix_len*gnmerge_max_cand_RHS_exts]);

		gpLHS_indices[gnprefix_len] = pcand_RHSs[i];
		ncand_ext_SNP_id = pRep_SNPs[pcand_RHSs[i]].nid;
		if(ncand_ext_SNP_id>gpprefix[gnprefix_len-1])
			printf("Error: the SNPs should be in descending order\n");
		gpprefix[gnprefix_len++] = ncand_ext_SNP_id;

		nnew_end_pos = nend_pos;
		if(nnew_end_pos>pRep_SNPs[pcand_RHSs[i]].nend_pos)
			nnew_end_pos = pRep_SNPs[pcand_RHSs[i]].nend_pos;
		nnew_start_pos = nstart_pos;
		if(nnew_start_pos<pRep_SNPs[pcand_RHSs[i]].nstart_pos)
			nnew_start_pos = pRep_SNPs[pcand_RHSs[i]].nstart_pos;

		num_of_remain_RHSs = num_of_cand_RHSs-1;
		ncand_ext_pos = -1;

		//j = num_of_cand_RHSs-1;
		//while(j>=0 && pcand_RHSs[j]>=nnew_end_pos)
		//	j--;

		while(nlast_RHS_pos<num_of_cand_RHSs && pcand_RHSs[nlast_RHS_pos]<nnew_end_pos)
			nlast_RHS_pos++;
		nlast_RHS_pos--;
		while(nlast_RHS_pos>=0 && pcand_RHSs[nlast_RHS_pos]>=nnew_end_pos)
			nlast_RHS_pos--;
		j = nlast_RHS_pos;
		if(nlast_RHS_pos<0)
			nlast_RHS_pos = 0;
		pRHS_node = pRep_SNPs[pcand_RHSs[i]].pRHS_node_head;
		while(pRHS_node!=NULL && pRHS_node->nSNPid>=nnew_start_pos && j>=0)
		{
			if(pRHS_node->nSNPid==pcand_RHSs[j])
			{
				pRHS_node = pRHS_node->pnext;
				j--;
			}
			else if(j==i)
			{
				ncand_ext_pos = num_of_remain_RHSs;
				j--;
			}
			else if(pRHS_node->nSNPid < pcand_RHSs[j])
			{
				pnew_cand_RHSs[num_of_remain_RHSs--] = pcand_RHSs[j];
				j--;
			}
			else 
				pRHS_node = pRHS_node->pnext;
		}
		while(j>=0 && pcand_RHSs[j]>=nnew_start_pos)
		{
			if(j==i)
				ncand_ext_pos = num_of_remain_RHSs;
			else
				pnew_cand_RHSs[num_of_remain_RHSs--] = pcand_RHSs[j];
			j--;
		}
		if(ncand_ext_pos==-1)
			printf("Error: the value of ncand_ext_pos should not be -1\n");
		num_of_remain_exts = ncand_ext_pos-num_of_remain_RHSs;
		num_of_remain_RHSs++;
		pnew_cand_RHSs = &pnew_cand_RHSs[num_of_remain_RHSs];
		num_of_remain_RHSs = num_of_cand_RHSs-num_of_remain_RHSs;

		if(gnprefix_len>=3 && num_of_remain_RHSs>0)
		{
			num_of_new_cand_RHSs = 0;
			num_of_new_cand_exts = 0;
			for(j=0;j<num_of_remain_RHSs;j++)
			{
				if(pRep_SNPs[pnew_cand_RHSs[j]].pLHS_node_map!=NULL)
					bis_covered = IsCoveredBySubset(pRep_SNPs, pnew_cand_RHSs[j], gpLHS_indices, gnprefix_len);
				else
					bis_covered = false;
				if(!bis_covered)
					pnew_cand_RHSs[num_of_new_cand_RHSs++] = pnew_cand_RHSs[j];
				if(j==num_of_remain_exts-1)
					num_of_new_cand_exts = num_of_new_cand_RHSs;
			}
			num_of_remain_RHSs = num_of_new_cand_RHSs;
			num_of_remain_exts = num_of_new_cand_exts;
		}

		if(num_of_remain_RHSs>0)
		{
			ncounter_num = 1<<(gnprefix_len-1);
			memset(p1counters, 0, sizeof(int)*ncounter_num);
			for(j=0;j<pSNPs[ncand_ext_SNP_id].n1;j++)
				p1counters[pdfs_bitmap[pSNPs[ncand_ext_SNP_id].p1sid_list[j]]]++;
			memset(pcounters, 0, sizeof(int)*ncounter_num*2);
			for(j=0;j<ncounter_num;j++)
			{
				pcounters[j<<1] = pprefix_counters[j]-p1counters[j];
				pcounters[(j<<1)+1] = p1counters[j];
			}

			for(j=0;j<gnhapmap_num;j++)
				pbitmap[j] = (pdfs_bitmap[j]<<1) | pSNPs[ncand_ext_SNP_id].pbitmap[j];

			nRHS_len = 0;
			num_of_new_cand_exts = 0;
			for(j=0;j<num_of_remain_exts;j++)
			{
				nid = pRep_SNPs[pnew_cand_RHSs[j]].nid;
				dR2 = CalcMultiMarkerR2(pbitmap, gnprefix_len, pSNPs[nid].p1sid_list, pSNPs[nid].n1, pcounters, ncomb_bitmap);
				if(dR2>=gpmin_r2[gnprefix_len-1])
				{
					gpRHS_indices[nRHS_len] = pnew_cand_RHSs[j];
					gprule_RHS_nodes[nRHS_len].nSNPid = nid;
					gprule_RHS_nodes[nRHS_len].dR2 = dR2;
					gprule_RHS_nodes[nRHS_len].ncomb_bitmap = ncomb_bitmap;
					nRHS_len++;
				}
				else
					pnew_cand_RHSs[num_of_new_cand_exts++] = pnew_cand_RHSs[j];
			}
			num_of_new_cand_RHSs = num_of_new_cand_exts;
			for(j=num_of_remain_exts;j<num_of_remain_RHSs;j++)
			{
				nid = pRep_SNPs[pnew_cand_RHSs[j]].nid;
				dR2 = CalcMultiMarkerR2(pbitmap, gnprefix_len, pSNPs[nid].p1sid_list, pSNPs[nid].n1, pcounters, ncomb_bitmap);
				if(dR2>=gpmin_r2[gnprefix_len-1])
				{
					gpRHS_indices[nRHS_len] = pnew_cand_RHSs[j];
					gprule_RHS_nodes[nRHS_len].nSNPid = nid;
					gprule_RHS_nodes[nRHS_len].dR2 = dR2;
					gprule_RHS_nodes[nRHS_len].ncomb_bitmap = ncomb_bitmap;
					nRHS_len++;
				}
				else
					pnew_cand_RHSs[num_of_new_cand_RHSs++] = pnew_cand_RHSs[j];
			}

			if(nRHS_len>0)
			{
				OutputRule(pSNPs, gnprefix_len, gpprefix, nRHS_len, gprule_RHS_nodes);
				if(gnprefix_len<gnmax_len && gnmax_covered_times>=0)
					AddOneRule(pRep_SNPs, gpLHS_indices, gnprefix_len, gpRHS_indices, nRHS_len);
			}

			if(gnprefix_len<gnmax_len && num_of_new_cand_exts>0 && num_of_new_cand_RHSs>1)
				DFSGenRules(pbitmap, pnew_cand_RHSs, num_of_new_cand_exts, num_of_new_cand_RHSs, pSNPs, pRep_SNPs, nnew_start_pos, nnew_end_pos, pcounters);
		}

		gnprefix_len--;
	}
}

void ShiftWindow(SNP *pSNPs, REP_SNP *pRep_SNPs, int ncur_pos, int num_of_SNPs, int &nleft_boundary, int &ncur_null_pos)
{
	RHS_NODE_PAGE *pRHS_node_page;
	SNPSET_PAGE *pSNPset_page;
	int i, nid, nnull_pos_id;

	for(i=nleft_boundary;i<pRep_SNPs[ncur_pos].nleft_boundary;i++)
	{
		nid = pRep_SNPs[i].nid;

		if(ncur_null_pos<num_of_SNPs)
		{
			nnull_pos_id = pRep_SNPs[ncur_null_pos].nid;
			pSNPs[nnull_pos_id].pbitmap = pSNPs[nid].pbitmap;
			pSNPs[nnull_pos_id].p1sid_list = pSNPs[nid].p1sid_list;
			ncur_null_pos++;
		}
		else
		{
			DelBitmap(pSNPs[nid].pbitmap);
			Del1SidList(pSNPs[nid].p1sid_list);
		}

		if(pRep_SNPs[i].pLHS_node_map!=NULL)
		{
			DelLHSNodeMap(pRep_SNPs[i].pLHS_node_map);
			pRep_SNPs[i].pLHS_node_map = NULL;
		}
		if(pRep_SNPs[i].pRHS_node_head!=NULL)
			pRep_SNPs[i].pRHS_node_head = NULL;

		while(goRHSnode_buf.phead->nmax_SNP_id<=i && goRHSnode_buf.phead!=goRHSnode_buf.pcur_page)
		{
			pRHS_node_page = goRHSnode_buf.phead;
			goRHSnode_buf.phead = goRHSnode_buf.phead->pnext;
			pRHS_node_page->pnext = goRHSnode_buf.pcur_page->pnext;
			goRHSnode_buf.pcur_page->pnext = pRHS_node_page;
		}
		while(goSNPset_buf.phead->nmax_SNP_id<=i && goSNPset_buf.phead!=goSNPset_buf.pcur_page)
		{
			pSNPset_page = goSNPset_buf.phead;
			goSNPset_buf.phead = goSNPset_buf.phead->pnext;
			pSNPset_page->pnext = goSNPset_buf.pcur_page->pnext;
			goSNPset_buf.pcur_page->pnext = pSNPset_page;
		}
	}
	nleft_boundary = pRep_SNPs[ncur_pos].nleft_boundary;
}

void EndShift(SNP* pSNPs, REP_SNP *pRep_SNPs, int num_of_rep_SNPs, int nleft_boundary)
{
	int i, nid; 

	for(i=nleft_boundary;i<num_of_rep_SNPs;i++)
	{
		nid = pRep_SNPs[i].nid;

		DelBitmap(pSNPs[nid].pbitmap);
		Del1SidList(pSNPs[nid].p1sid_list);

		if(pRep_SNPs[i].pLHS_node_map!=NULL)
		{
			DelLHSNodeMap(pRep_SNPs[i].pLHS_node_map);
			pRep_SNPs[i].pLHS_node_map = NULL;
		}
		if(pRep_SNPs[i].pRHS_node_head!=NULL)
			pRep_SNPs[i].pRHS_node_head = NULL;
	}
	DelRHSNodeBuf();
	DelLHSNodeMaps();
	DelSNPsetBuf();
}


void GenRules(char* szmatrix_filename, SNP *pSNPs, REP_SNP* pRep_SNPs, int num_of_rep_SNPs, char* szoutput_name)
{
	FILE *fpmatrix;
	int i, j, nid1, nid2, nleft_boundary;
	int ncounter_num, *pcand_RHSs, num_of_cand_exts, num_of_cand_RHSs, ncomb_bitmap, pcounters[2];
	char szrule_filename[200];
	RHS_NODE *pRHS_node;
	float dR2;
	int ncur_load_pos, nline_no, ncur_null_pos;

	fpmatrix = fopen(szmatrix_filename, "rt");
	if(fpmatrix==NULL)
	{
		printf("Error: cannot open file %s for read\n", szmatrix_filename);
		return;
	}	

	sprintf(szrule_filename, "%s.rule.bin", szoutput_name);
	gfp_rule = fopen(szrule_filename, "wb");
	if(gfp_rule==NULL)
	{
		printf("Error: cannot open file %s for write\n", szrule_filename);
		return;
	}
	gnrule_file_pos = 0;
	gnmax_RHS_len = 1;

	InitRHSNodeBuf();
	InitSNPsetBuf();

	gpDFS_cand_RHSs = new int[gnmerge_max_cand_RHS_exts*gnmax_len];
	IncMemSize(sizeof(int)*gnmerge_max_cand_RHS_exts*gnmax_len);
	gprule_RHS_nodes = new RHS_INFO_NODE[gnmerge_max_cand_RHS_exts];
	IncMemSize(sizeof(RHS_INFO_NODE)*gnmerge_max_cand_RHS_exts);

	gpDFS_bitmap_buf = new char[gnhapmap_num*gnmax_len];
	IncMemSize(sizeof(char)*gnhapmap_num*gnmax_len);
	ncounter_num = 1<<(gnmax_len+1);
	gpmerged_bins = new int[ncounter_num];
	IncMemSize(sizeof(int)*ncounter_num);

	gpprefix = new int[gnmax_len];
	IncMemSize(sizeof(int)*gnmax_len);
	gnprefix_len = 0;
	gpLHS_indices = new int[gnmax_len];
	IncMemSize(sizeof(int)*gnmax_len);
	gpRHS_indices = new int[gnmerge_max_cand_RHS_exts];
	IncMemSize(sizeof(int)*gnmerge_max_cand_RHS_exts);


	gpfree_LHS_node_maps = NULL;

	gntotal_calls = 0;
	gdsubset_check_times = 0;
	gdsubset_check_depth_sum = 0;

	nleft_boundary = 0;
	ncur_load_pos = 0;
	nline_no = 0;
	ncur_null_pos = 0;

	for(i=0;i<num_of_rep_SNPs;i++)
	{
		nid1 = pRep_SNPs[i].nid;
		ShiftWindow(pSNPs, pRep_SNPs, i, num_of_rep_SNPs, nleft_boundary, ncur_null_pos);
		LoadBitmaps(fpmatrix, nline_no, pSNPs, pRep_SNPs, i, ncur_load_pos, ncur_null_pos);

		if(gnmax_len==1)
		{
			for(j=i+1;j<pRep_SNPs[i].nend_pos;j++)
			{
				if(i>=pRep_SNPs[j].nstart_pos)
				{
					nid2 = pRep_SNPs[j].nid;
					dR2 = CalcR2(pSNPs[nid1].pbitmap, pSNPs[nid1].n1, pSNPs[nid2].p1sid_list, pSNPs[nid2].n1, ncomb_bitmap);
					if(dR2>=gpmin_r2[0])
						OutputRule(pSNPs, nid1, nid2, dR2, ncomb_bitmap);
				}
			}
		}
		else
		{
			pcand_RHSs = gpDFS_cand_RHSs;
			num_of_cand_exts = i-pRep_SNPs[i].nstart_pos-1; //the first position
			
			j = i-1;
			pRHS_node = pRep_SNPs[i].pRHS_node_head;
			while(pRHS_node!=NULL && j>=pRep_SNPs[i].nstart_pos)
			{
				if(pRHS_node->nSNPid==j)
				{
					pRHS_node = pRHS_node->pnext;
					j--;
				}
				else if(pRHS_node->nSNPid<j)
				{
					if(i<pRep_SNPs[j].nend_pos)
						pcand_RHSs[num_of_cand_exts--] = j;
					j--;
				}
				else 
					printf("Error: the RHS should not be larger than %d\n", j);
			}
			while(j>=pRep_SNPs[i].nstart_pos)
			{
				if(i<pRep_SNPs[j].nend_pos)
					pcand_RHSs[num_of_cand_exts--] = j;
				j--;
			}
			num_of_cand_exts++;
			pcand_RHSs = &(pcand_RHSs[num_of_cand_exts]);
			num_of_cand_exts = i-pRep_SNPs[i].nstart_pos-num_of_cand_exts;

			num_of_cand_RHSs = num_of_cand_exts;
			for(j=i+1;j<pRep_SNPs[i].nend_pos;j++)
			{
				if(i>=pRep_SNPs[j].nstart_pos)
				{
					nid2 = pRep_SNPs[j].nid;
					dR2 = CalcR2(pSNPs[nid1].pbitmap, pSNPs[nid1].n1, pSNPs[nid2].p1sid_list, pSNPs[nid2].n1, ncomb_bitmap);
					if(dR2>=gpmin_r2[0])
					{
						OutputRule(pSNPs, nid1, nid2, dR2, ncomb_bitmap);
						AddOneRule(pRep_SNPs, i, j);
					}
					else 
						pcand_RHSs[num_of_cand_RHSs++] = j;
				}
			}				

			if(num_of_cand_exts>0 && num_of_cand_RHSs>1)
			{
				gpLHS_indices[gnprefix_len] = i;
				gpprefix[gnprefix_len++] = nid1;
				
				pcounters[0] = gnhapmap_num-pSNPs[nid1].n1;
				pcounters[1] = pSNPs[nid1].n1;
				DFSGenRules(pSNPs[nid1].pbitmap, pcand_RHSs, num_of_cand_exts, num_of_cand_RHSs, pSNPs, pRep_SNPs, pRep_SNPs[i].nstart_pos, pRep_SNPs[i].nend_pos, pcounters);

				gnprefix_len--;
			}
		}
	}
	EndShift(pSNPs, pRep_SNPs, num_of_rep_SNPs, nleft_boundary);

	fclose(fpmatrix);
	fclose(gfp_rule);
	gntotal_rule_size = gnrule_file_pos;

	delete []gpLHS_indices;
	DecMemSize(sizeof(int)*gnmax_len);
	delete []gpRHS_indices;
	DecMemSize(sizeof(int)*gnmerge_max_cand_RHS_exts);
	delete []gpprefix;
	DecMemSize(sizeof(int)*gnmax_len);
	delete []gprule_RHS_nodes;
	DecMemSize(sizeof(RHS_INFO_NODE)*gnmerge_max_cand_RHS_exts);

	delete []gpDFS_bitmap_buf;
	DecMemSize(sizeof(char)*gnhapmap_num*gnmax_len);
	delete []gpmerged_bins;
	DecMemSize(sizeof(int)*ncounter_num);
	delete []gpDFS_cand_RHSs;
	DecMemSize(sizeof(int)*gnmerge_max_cand_RHS_exts*gnmax_len);
}



//============================================================================================

void VerifyRule(REP_SNP* pRep_SNPs, int nindex1, int nindex2)
{
	if(nindex1==nindex2)
		printf("Error: the two SNPs should not be the same\n");

	if(pRep_SNPs[nindex1].pRHS_node_head!=NULL)
	{
		if(nindex2<=pRep_SNPs[nindex1].pRHS_node_head->nSNPid)
			printf("Error: the RHS should be in desending order\n");
	}
	if(pRep_SNPs[nindex2].pRHS_node_head!=NULL)
	{
		if(nindex1<=pRep_SNPs[nindex2].pRHS_node_head->nSNPid)
			printf("Error: the RHS should be in desending order\n");
	}
}



void VerifyRule(REP_SNP *pRep_SNPs, int *pLHS_set, int nLHS_len, int *pRHS_set, int nRHS_len)
{
	int i;

	for(i=0;i<nLHS_len-1;i++)
	{
		if(pLHS_set[i]<=pLHS_set[i+1])
			printf("Error: the LHS SNPs should be in descending order\n");
	}
	for(i=0;i<nRHS_len-1;i++)
	{
		if(pRHS_set[i]>=pRHS_set[i+1])
			printf("Error: the RHS SNPs should be in ascending order\n");
	}

}

