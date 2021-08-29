#include "global.h"


SNP_ID_NODE *gpSNP_ID_nodes;
SNP_ID_NODE **gpSNP_ID_heads;

int gncur_highest_SCS;

char* gpcovered_flag_buf;
FILE *gfp_cvg;

int gnum_of_partitions;
int gnum_of_iterations;
int gnorig_num_of_tag_SNPs;
int gnnot_covered_tag_SNPs;

int gnum_of_tag_SNPs;
int gncovered_SNPs;

int gninvalid_rule_num;
int *gpunfolded_LHS_set;

int *gptemp_merged_ids;

RULE_BUF gorule_buf;
RULE_NODE_BUF gorule_node_buf;
RULE_NODE *gpfree_rule_nodes;
RULE *gpfree_rules;

void SetRHSCoverFlags(SNP* pSNPs, int *pLHS_set, int nLHS_len, int *pRHS_set, int nRHS_len, int *punfolded_LHS_set, int ncur_pos, int nleft_boundary, int nright_boundary, int &nvalid_rule_num)
{
	int i, j, k, ncur_sid, nmerged_id, nRHS_id, nRHS_merged_id;

	ncur_sid = pLHS_set[ncur_pos];
	if(pSNPs[ncur_sid].pmerged_ids==NULL)
		printf("Error: the pmerged_ids field should not be NULL\n");

	for(i=0;i<pSNPs[ncur_sid].num_of_merged_ids;i++)
	{
		nmerged_id = pSNPs[ncur_sid].pmerged_ids[i];
		if(nmerged_id>=nleft_boundary && nmerged_id<nright_boundary)
		{
			punfolded_LHS_set[ncur_pos] = nmerged_id;
			if(nleft_boundary<gporig_window_start_pos[nmerged_id])
				nleft_boundary = gporig_window_start_pos[nmerged_id];
			if(nright_boundary>gporig_window_end_pos[nmerged_id])
				nright_boundary = gporig_window_end_pos[nmerged_id];
			if(ncur_pos+1==nLHS_len)
			{
				for(j=0;j<nRHS_len;j++)
				{
					nRHS_id = pRHS_set[j];
					if(pSNPs[nRHS_id].pmerged_ids==NULL)
						printf("Error: the pmerged_ids field should not be NULL\n");
					for(k=0;k<pSNPs[nRHS_id].num_of_merged_ids;k++)
					{
						nRHS_merged_id = pSNPs[nRHS_id].pmerged_ids[k];
						if(nRHS_merged_id>=nleft_boundary && nRHS_merged_id<nright_boundary)
						{
							pSNPs[nRHS_merged_id].nbeing_covered_times++;
							nvalid_rule_num++;
						}
						else if(nRHS_merged_id>=nright_boundary)
							break;
					}
				}
			}
			else
				SetRHSCoverFlags(pSNPs, pLHS_set, nLHS_len, pRHS_set, nRHS_len, punfolded_LHS_set, ncur_pos+1, nleft_boundary, nright_boundary, nvalid_rule_num);
		}
		else if(nmerged_id>=nright_boundary)
			break;
	}
}

void SetCoverFlags(SNP *pSNPs, int nstart_pos, int nend_pos, int nrule_end_pos)
{
	int i, j, nSNP_id_merged_with, nSNP_id1, nSNP_id2, nLHS_id, nRHS_id, num_of_rules, nvalid_rule_num;
	RULE *prule, *pnew_rule;
	RULE_NODE *prule_node;

	for(i=nstart_pos;i<nend_pos;i++)
	{
		pSNPs[i].nflag = 0;
		pSNPs[i].num_of_covered_SNPs = 1;
		pSNPs[i].nbeing_covered_times = 0;
		pSNPs[i].prule_node_head = NULL;
	}

	for(i=nstart_pos;i<nend_pos;i++)
	{
		if(pSNPs[i].nwindow_start_pos>=0)
			nSNP_id_merged_with = i;
		else
			nSNP_id_merged_with = pSNPs[i].nSNPid_merged_with;
		for(j=0;j<pSNPs[nSNP_id_merged_with].num_of_merged_ids;j++)
		{
			nSNP_id1 = pSNPs[nSNP_id_merged_with].pmerged_ids[j];
			if(nSNP_id1!=i && nSNP_id1>=nstart_pos && nSNP_id1<nend_pos && 
				nSNP_id1>=gporig_window_start_pos[i] && nSNP_id1<gporig_window_end_pos[i])
			{
				pSNPs[i].pcovered_flags[nSNP_id1-gporig_window_start_pos[i]] = 3;
				pSNPs[i].num_of_covered_SNPs++;
				//pSNPs[i].nbeing_covered_times++;
			}
		}
	}

	num_of_rules = 0;

	prule = ReadInOneRule(nrule_end_pos);
	while(prule!=NULL)
	{
		if(prule->nLHS_len>gnmax_len)
			printf("Error: the length of the rule cannot be larger than %d\n", gnmax_len);
		else if(prule->nLHS_len==1)
		{
			nvalid_rule_num = 0;
			nSNP_id1 = prule->pLHS_set[0];
			nSNP_id2 = prule->pRHS_set[0];
			i = 0;
			for(i=0;i<pSNPs[nSNP_id1].num_of_merged_ids;i++)
			{
				nLHS_id = pSNPs[nSNP_id1].pmerged_ids[i];
				if(nLHS_id>=nstart_pos && nLHS_id<nend_pos)
				{
					for(j=0;j<pSNPs[nSNP_id2].num_of_merged_ids;j++)
					{
						nRHS_id = pSNPs[nSNP_id2].pmerged_ids[j];
						if(nRHS_id>=nstart_pos && nRHS_id<nend_pos && 
							nRHS_id>=gporig_window_start_pos[nLHS_id] && nRHS_id<gporig_window_end_pos[nLHS_id])
						{
							nvalid_rule_num++;
							if(nLHS_id<gporig_window_start_pos[nRHS_id] || nRHS_id>=gporig_window_end_pos[nRHS_id])
								printf("Error: nLHS_id should be in the window of nRHS_id\n");

							pSNPs[nLHS_id].pcovered_flags[nRHS_id-gporig_window_start_pos[nLHS_id]] = 3;
							pSNPs[nLHS_id].num_of_covered_SNPs++;
							pSNPs[nLHS_id].nbeing_covered_times++;

							pSNPs[nRHS_id].pcovered_flags[nLHS_id-gporig_window_start_pos[nRHS_id]] = 3;
							pSNPs[nRHS_id].num_of_covered_SNPs++;
							pSNPs[nRHS_id].nbeing_covered_times++;
						}
						else if(nRHS_id>=nend_pos)
							break;
					}
				}
				else if(nLHS_id>=nend_pos)
					break;
			}
			if(nvalid_rule_num==0)
				gninvalid_rule_num++;
		}
		else
		{
/*
			nvalid_rule_num = 0;
			SetRHSCoverFlags(pSNPs, prule->pLHS_set, prule->nLHS_len, prule->pRHS_set, prule->nRHS_len, gpunfolded_LHS_set, 0, nstart_pos, nend_pos, nvalid_rule_num);
			if(nvalid_rule_num==0)
				gninvalid_rule_num++;
			else
			{
				pnew_rule = NewRule(prule->nLHS_len, prule->pLHS_set, prule->nRHS_len, prule->pRHS_set);
				for(i=0;i<prule->nLHS_len;i++)
				{
					nLHS_id = prule->pLHS_set[i];
					prule_node = NewRuleNode();
					prule_node->prule = pnew_rule;
					prule_node->pnext = pSNPs[nLHS_id].prule_node_head;
					pSNPs[nLHS_id].prule_node_head = prule_node;
				}
			}
*/

			for(i=0;i<prule->nLHS_len;i++)
			{
				if(prule->pLHS_set[i]>=nend_pos)
					break;
			}
			if(i>=prule->nLHS_len)
			{
				pnew_rule = NewRule(prule->nLHS_len, prule->pLHS_set, prule->nRHS_len, prule->pRHS_set);
				for(i=0;i<prule->nLHS_len;i++)
				{
					nLHS_id = prule->pLHS_set[i];
					prule_node = NewRuleNode();
					prule_node->prule = pnew_rule;
					prule_node->pnext = pSNPs[nLHS_id].prule_node_head;
					pSNPs[nLHS_id].prule_node_head = prule_node;
				}
				for(i=0;i<prule->nRHS_len;i++)
				{
					nRHS_id = prule->pRHS_set[i];
					if(nRHS_id<nend_pos)
					{
						for(j=0;j<pSNPs[nRHS_id].num_of_merged_ids;j++)
						{
							if(pSNPs[nRHS_id].pmerged_ids[j]<nend_pos)
								pSNPs[pSNPs[nRHS_id].pmerged_ids[j]].nbeing_covered_times++;
						}
					}
				}
			}
		}
		num_of_rules++;
		prule = ReadInOneRule(nrule_end_pos);
	}
	//printf("#rules: %d, #invalid rules: %d\n", num_of_rules, gninvalid_rule_num);

}

void InitilizeSNPIDHead(SNP *pSNPs, int nstart_pos, int nend_pos)
{
	int i, nSCS, num_of_SNPs;
	SNP_ID_NODE *pSNP_ID_node;

	num_of_SNPs = nend_pos-nstart_pos;

	gpSNP_ID_nodes = new SNP_ID_NODE[num_of_SNPs];
	IncMemSize(sizeof(SNP_ID_NODE)*num_of_SNPs);
	memset(gpSNP_ID_heads, 0, sizeof(SNP_ID_NODE*)*(gnmax_cand_RHS_exts+1));

	for(i=nend_pos-1;i>=nstart_pos;i--)
	{
		nSCS = pSNPs[i].num_of_covered_SNPs;
		pSNP_ID_node = &gpSNP_ID_nodes[i-nstart_pos];
		pSNP_ID_node->nid = i;
		pSNP_ID_node->nflag = 0;
		pSNP_ID_node->pprev = NULL;
		pSNP_ID_node->pnext = gpSNP_ID_heads[nSCS];
		if(gpSNP_ID_heads[nSCS]!=NULL)
			gpSNP_ID_heads[nSCS]->pprev = pSNP_ID_node;
		gpSNP_ID_heads[nSCS] = pSNP_ID_node;
	}
	gncur_highest_SCS = gnmax_cand_RHS_exts;
}

RULE_NODE* SplitRules(SNP *pSNPs, int nTag_SNP_id, int nSNPid_merged_with)
{
	int i, *pmerged_ids, num_of_merged_ids, nindex;
	RULE_NODE *prule_node, *pnew_rule_node, *pnew_head;
	RULE *prule, *pnew_rule;
	bool bsplit;

	num_of_merged_ids = pSNPs[nSNPid_merged_with].num_of_merged_ids;
	if(num_of_merged_ids==1)
	{
		if(nTag_SNP_id!=nSNPid_merged_with)
			printf("Error: the two ids should be the same\n");
		return NULL;
	}

	pnew_head = NULL;

	pmerged_ids = pSNPs[nSNPid_merged_with].pmerged_ids;
	if(gbno_equv_tag_SNPs && pmerged_ids[0]>=gporig_window_start_pos[nTag_SNP_id] && pmerged_ids[num_of_merged_ids-1]<gporig_window_end_pos[nTag_SNP_id])
		bsplit = false;
	else 
		bsplit = true;

	if(bsplit)
	{
		prule_node = pSNPs[nSNPid_merged_with].prule_node_head;
		pSNPs[nTag_SNP_id].prule_node_head = NULL;

		pnew_head = prule_node;

		while(prule_node!=NULL)
		{
			prule = prule_node->prule;
			if(prule->nLHS_len>0)
			{
				nindex = -1;
				for(i=0;i<prule->nLHS_len;i++)
				{
					if(prule->pLHS_set[i]==nSNPid_merged_with)
					{
						nindex = i;
						break;
					}
				}
				if(nindex==-1)
					printf("Error: cannot find %d in LHS\n", nSNPid_merged_with);
				else if(prule->nLHS_len>0)
				{
					pnew_rule = NewRule(prule->nLHS_len, prule->pLHS_set, prule->nRHS_len, prule->pRHS_set);
					pnew_rule->pLHS_set[nindex] = nTag_SNP_id;

					for(i=0;i<pnew_rule->nLHS_len;i++)
					{
						pnew_rule_node = NewRuleNode();
						pnew_rule_node->prule = pnew_rule;
						pnew_rule_node->pnext = pSNPs[pnew_rule->pLHS_set[i]].prule_node_head;
						pSNPs[pnew_rule->pLHS_set[i]].prule_node_head = pnew_rule_node;				
					}
				}
			}

			prule_node = prule_node->pnext;
		}
	}
	else if(nTag_SNP_id!=nSNPid_merged_with)
	{
		pSNPs[nTag_SNP_id].prule_node_head = pSNPs[nSNPid_merged_with].prule_node_head;
	}

	return pnew_head;
}


void UpdateCoverage(SNP *pSNPs, int nstart_pos, int nend_pos, int nTagSNP_id, int &nuncovered_SNPs)
{
	int num_of_covered, i, j, k;
	int num_of_untaged_SNPs, nindex, nuntaged_SNP_id, nLHS_id, nRHS_id;
	RULE_NODE *prule_node, *psplited_rule_node_head;
	RULE *prule;
	int *pmerged_ids, num_of_merged_ids, nSNPid_merged_with, nRHS_merged_id;
	int ntag_window_start_pos, ni_window_start_pos, nleft_boundary, nright_boundary;

	if((pSNPs[nTagSNP_id].nflag & 1) || (pSNPs[nTagSNP_id].nflag & 4))
		printf("Error: the SNP or its equivalent SNP has been selected as tag before\n");

	ntag_window_start_pos = gporig_window_start_pos[nTagSNP_id];
	if(ntag_window_start_pos<nstart_pos)
		ntag_window_start_pos = nstart_pos;

	if(pSNPs[nTagSNP_id].nflag & 2) // already been coverd
		num_of_covered = 0;
	else
	{
		nuncovered_SNPs--;
		num_of_covered = 1;

		//decrease SCS of those SNPs that covering this Tag SNP
		for(i=ntag_window_start_pos;i<gporig_window_end_pos[nTagSNP_id] && i<nend_pos;i++)
		{
			if(!(pSNPs[i].nflag & 1) && !(pSNPs[i].nflag & 4) && i!=nTagSNP_id) 
			{
				if((pSNPs[nTagSNP_id].pcovered_flags[i-gporig_window_start_pos[nTagSNP_id]] & 2)) // nTagSNP_id is covered by i
				{
					if(pSNPs[i].pcovered_flags[nTagSNP_id-gporig_window_start_pos[i]] & 1) //i covers nTagSNP_id
					{
						MoveSNPIDNode(pSNPs, i, nstart_pos, -1, i);
						pSNPs[i].pcovered_flags[nTagSNP_id-gporig_window_start_pos[i]] -= 1;
					}
					else 
						printf("Error: the flag should be 1\n");
				}
			}
		}
	}

	pSNPs[nTagSNP_id].nflag = 3;
	
	//decrease SCS of those SNPs that cover a SNP that is covered by this tag SNP
	for(i=ntag_window_start_pos;i<gporig_window_end_pos[nTagSNP_id] && i<nend_pos;i++)
	{
		if(!(pSNPs[i].nflag & 1) && !(pSNPs[i].nflag & 4)) 
		{
			if(pSNPs[nTagSNP_id].pcovered_flags[i-gporig_window_start_pos[nTagSNP_id]] & 1) //i is covered by nTagSNP_id
			{
				if(pSNPs[i].pcovered_flags[nTagSNP_id-gporig_window_start_pos[i]] & 2) //nTagSNP_id covers i
					pSNPs[i].pcovered_flags[nTagSNP_id-gporig_window_start_pos[i]] -= 2;
				else
					printf("Error: the flag should be 2\n");

				if(pSNPs[i].nflag & 2)
					printf("Error: the SNP should not be covered before\n");
				else 
					pSNPs[i].nflag |= 2;

				num_of_covered++;
				nuncovered_SNPs--;

				MoveSNPIDNode(pSNPs, i, nstart_pos, -1, i); //it covers itself, so minus 1
				
				ni_window_start_pos = gporig_window_start_pos[i];
				if(ni_window_start_pos<nstart_pos)
					ni_window_start_pos = nstart_pos;
				for(j=ni_window_start_pos;j<gporig_window_end_pos[i] && j<nend_pos;j++)
				{
					if(!(pSNPs[j].nflag & 1) && !(pSNPs[j].nflag & 4) )
					{
						if(pSNPs[i].pcovered_flags[j-gporig_window_start_pos[i]] & 2) // SNP j covers SNP i
						{
							if(pSNPs[j].pcovered_flags[i-gporig_window_start_pos[j]] & 1) //i is covered by j
							{
								MoveSNPIDNode(pSNPs, j, nstart_pos, -1, i);
								pSNPs[j].pcovered_flags[i-gporig_window_start_pos[j]] -= 1;
							}
							else 
								printf("Error: the flag should be 1\n");
						}
					}
				}
			}
		}
	}
	if(num_of_covered!=pSNPs[nTagSNP_id].num_of_covered_SNPs)
		printf("Erorr: inconsistent number of covered SNPs\n");

	if(pSNPs[nTagSNP_id].pmerged_ids!=NULL)
		nSNPid_merged_with = nTagSNP_id;
	else
		nSNPid_merged_with = pSNPs[nTagSNP_id].nSNPid_merged_with;

	//exclude those SNPs that are identical to the tag SNP selected and are within maximum window distance as candidate tag SNPs
	if(gbno_equv_tag_SNPs)
	{
		pmerged_ids = pSNPs[nSNPid_merged_with].pmerged_ids;
		num_of_merged_ids = pSNPs[nSNPid_merged_with].num_of_merged_ids;
		for(i=0;i<num_of_merged_ids;i++)
		{
			if(pmerged_ids[i]!=nTagSNP_id && pmerged_ids[i]<nend_pos && 
				pmerged_ids[i]>=ntag_window_start_pos && pmerged_ids[i]<gporig_window_end_pos[nTagSNP_id])
			{
				if(!(pSNPs[pmerged_ids[i]].nflag & 2))
					printf("Error: this SNP should be covered already\n");
				RemoveSNPIDNode(pSNPs, pmerged_ids[i], nstart_pos);
				pSNPs[pmerged_ids[i]].nflag |= 4;
			}
			else if(pmerged_ids[i]>=nend_pos)
				break;
		}
	}

	//split rules if necessary
	psplited_rule_node_head = SplitRules(pSNPs, nTagSNP_id, nSNPid_merged_with);

	nuntaged_SNP_id = -1;
	num_of_untaged_SNPs = 0;
	//update multi SNPs 
	prule_node = pSNPs[nTagSNP_id].prule_node_head;
	while(prule_node!=NULL)
	{
		pSNPs[nTagSNP_id].prule_node_head = prule_node->pnext;
		prule = prule_node->prule;

		if(prule->nLHS_len>0)
		{
			nindex = -1;
			num_of_untaged_SNPs = 0;
			nuntaged_SNP_id = -1;
			nleft_boundary = nstart_pos;
			nright_boundary = nend_pos;
			for(i=0;i<prule->nLHS_len;i++)
			{
				if(prule->pLHS_set[i]==nSNPid_merged_with || prule->pLHS_set[i]==nTagSNP_id)
				{
					nindex = i;
					if(nleft_boundary<gporig_window_start_pos[nTagSNP_id])
						nleft_boundary = gporig_window_start_pos[nTagSNP_id];
					if(nright_boundary>gporig_window_end_pos[nTagSNP_id])
						nright_boundary = gporig_window_end_pos[nTagSNP_id];
				}
				else if(pSNPs[prule->pLHS_set[i]].nflag & 1)
				{
					if(nleft_boundary<gporig_window_start_pos[prule->pLHS_set[i]])
						nleft_boundary = gporig_window_start_pos[prule->pLHS_set[i]];
					if(nright_boundary>gporig_window_end_pos[prule->pLHS_set[i]])
						nright_boundary = gporig_window_end_pos[prule->pLHS_set[i]];
				}
				else
				{
					num_of_untaged_SNPs++;
					nuntaged_SNP_id = prule->pLHS_set[i];
					if(pSNPs[prule->pLHS_set[i]].nwindow_start_pos==-1)
						printf("Error: window start position should not be -1\n");
					if(nleft_boundary<pSNPs[prule->pLHS_set[i]].nwindow_start_pos)
						nleft_boundary = pSNPs[prule->pLHS_set[i]].nwindow_start_pos;
					if(nright_boundary>pSNPs[prule->pLHS_set[i]].nwindow_end_pos)
						nright_boundary = pSNPs[prule->pLHS_set[i]].nwindow_end_pos;
				}
			}
			if(nleft_boundary>=nright_boundary || nTagSNP_id<nleft_boundary || nTagSNP_id>=nright_boundary)
				prule->nLHS_len = 0;
			if(nindex==-1)
			{
				printf("Error: the tag SNP is not found in LHS\n");
				prule->nLHS_len = 0;
			}
		}

		if(prule->nLHS_len>0 && num_of_untaged_SNPs>=1)
		{
			if(nTagSNP_id<gporig_window_start_pos[prule->pLHS_set[nindex]] || 
				nTagSNP_id>=gporig_window_end_pos[prule->pLHS_set[nindex]])
				printf("Error: tag SNP selected should be within the window\n");
			
			prule->pLHS_set[nindex] = nTagSNP_id;

			if(num_of_untaged_SNPs==1 && nuntaged_SNP_id>=nstart_pos && nuntaged_SNP_id<nright_boundary)
			{
				for(k=0;k<pSNPs[nuntaged_SNP_id].num_of_merged_ids;k++)
				{
					nLHS_id = pSNPs[nuntaged_SNP_id].pmerged_ids[k];
					if(!(pSNPs[nLHS_id].nflag & 1) && !(pSNPs[nLHS_id].nflag & 4) && 
						nLHS_id>=nleft_boundary && nLHS_id<nright_boundary)
					{
						for(i=0;i<prule->nRHS_len;i++)
						{
							nRHS_merged_id = prule->pRHS_set[i];
							for(j=0;j<pSNPs[nRHS_merged_id].num_of_merged_ids;j++)
							{
								nRHS_id = pSNPs[nRHS_merged_id].pmerged_ids[j]; 
								if(!(pSNPs[nRHS_id].nflag & 2) && 
									nRHS_id>=gporig_window_start_pos[nLHS_id] && nRHS_id<gporig_window_end_pos[nLHS_id] &&
									nRHS_id>=nleft_boundary && nRHS_id<nright_boundary)
								{
									if(!(pSNPs[nLHS_id].pcovered_flags[nRHS_id-gporig_window_start_pos[nLHS_id]] & 1))
									{
										if(pSNPs[nRHS_id].pcovered_flags[nLHS_id-gporig_window_start_pos[nRHS_id]] & 2)
											printf("Error: RHS_id should not be covered by LHS_id before\n");
										pSNPs[nLHS_id].pcovered_flags[nRHS_id-gporig_window_start_pos[nLHS_id]] |= 1;
										pSNPs[nRHS_id].pcovered_flags[nLHS_id-gporig_window_start_pos[nRHS_id]] |= 2;

										MoveSNPIDNode(pSNPs, nLHS_id, nstart_pos, 1, nRHS_id);
										if(gncur_highest_SCS < pSNPs[nLHS_id].num_of_covered_SNPs)
											gncur_highest_SCS = pSNPs[nLHS_id].num_of_covered_SNPs;
									}
								}
							}
						}
					}
				}
			}
		}
		else if(prule->nLHS_len>0 && num_of_untaged_SNPs==0)
			DelRule(prule);

		DelRuleNode(prule_node);
		prule_node = pSNPs[nTagSNP_id].prule_node_head;
	}

	if(psplited_rule_node_head!=NULL)
		pSNPs[nSNPid_merged_with].prule_node_head = psplited_rule_node_head;

}


int AddNoneCoveredSNPs(SNP *pSNPs, int nstart_pos, int nend_pos, int &nuncovered_SNPs, int* ptagSNPs)
{
	int i, num_of_tag_SNPs;

	if(gdmax_tagSNP_num!=0 && gdmax_tagSNP_num!=1)
		return 0;

	num_of_tag_SNPs = 0;
	for(i=nstart_pos;i<nend_pos;i++)
	{
		if(pSNPs[i].nbeing_covered_times==0 && pSNPs[i].nflag==0)
		{
			RemoveSNPIDNode(pSNPs, i, nstart_pos);
			ptagSNPs[num_of_tag_SNPs++] = i;

			UpdateCoverage(pSNPs, nstart_pos, nend_pos, i, nuncovered_SNPs);
		}
	}
	//printf("#Tag SNPs not covered by any other SNPs: %d\n", num_of_tag_SNPs);
	return num_of_tag_SNPs;
}

void AddExistTagSNPs(SNP *pSNPs, int nstart_pos, int nend_pos, int* ptagSNPs, int num_of_tagSNPs, int &nuncovered_SNPs)
{
	int i;

	for(i=0;i<num_of_tagSNPs;i++)
	{
		RemoveSNPIDNode(pSNPs, ptagSNPs[i], nstart_pos);
		UpdateCoverage(pSNPs, nstart_pos, nend_pos, ptagSNPs[i], nuncovered_SNPs);
	}
}


int SelectTagSNPs(SNP *pSNPs, int nstart_pos, int nend_pos, int nrule_end_pos, int *pexist_tagSNPs, int nexist_tagSNP_num, int *ptagSNPs)
{
	int nuncovered_SNPs, nTagSNP_id, i, nflag_buf_pos;
	double dcoverage;
	SNP_ID_NODE *pSNP_ID_node;
	int num_of_tag_SNPs, nSNP_num;

	InitRuleBuf();
	InitRuleNodeBuf();

	nSNP_num = nend_pos-nstart_pos;

	gpcovered_flag_buf = new char[gnwindow_size_sum];
	IncMemSize(sizeof(char)*gnwindow_size_sum);
	memset(gpcovered_flag_buf, 0, sizeof(char)*gnwindow_size_sum);
	nflag_buf_pos = 0;
	for(i=nstart_pos;i<nend_pos;i++)
	{
		pSNPs[i].pcovered_flags = &(gpcovered_flag_buf[nflag_buf_pos]);
		nflag_buf_pos += gporig_window_end_pos[i]-gporig_window_start_pos[i];
	}
	if(nflag_buf_pos!=gnwindow_size_sum)
		printf("Error: inconsistent flag buffer size\n");


	SetCoverFlags(pSNPs, nstart_pos, nend_pos, nrule_end_pos);
	InitilizeSNPIDHead(pSNPs, nstart_pos, nend_pos);

	dcoverage = 0.05;

	nuncovered_SNPs = nend_pos-nstart_pos;

	AddExistTagSNPs(pSNPs, nstart_pos, nend_pos, pexist_tagSNPs, nexist_tagSNP_num, nuncovered_SNPs);

	//add those SNPs that are not covered by any other SNPs as Tag SNP
	num_of_tag_SNPs = AddNoneCoveredSNPs(pSNPs, nstart_pos, nend_pos, nuncovered_SNPs, ptagSNPs);
	gnnot_covered_tag_SNPs += num_of_tag_SNPs;

	fprintf(gfp_cvg, "%d  %d  %d  %.2f  1\n", num_of_tag_SNPs, nSNP_num-nuncovered_SNPs, nuncovered_SNPs, (double)(nSNP_num-nuncovered_SNPs)/nSNP_num);

	while(nuncovered_SNPs>0 && (gdmax_tagSNP_num==0 || gdmax_tagSNP_num==1 || 
		gdmax_tagSNP_num>1 && num_of_tag_SNPs<gdmax_tagSNP_num || 
		gdmax_tagSNP_num<1 && (double)(nSNP_num-nuncovered_SNPs)/nSNP_num<gdmax_tagSNP_num))
	{
		while(gncur_highest_SCS>0 && gpSNP_ID_heads[gncur_highest_SCS]==NULL)
			gncur_highest_SCS--;

		if(gncur_highest_SCS==0)
			printf("Error: all the SNPs should be covered already\n");


		pSNP_ID_node = gpSNP_ID_heads[gncur_highest_SCS];
		nTagSNP_id = pSNP_ID_node->nid;

		if(gncur_highest_SCS!=pSNPs[nTagSNP_id].num_of_covered_SNPs)
			printf("Error: inconsistent number of covered SNPs\n");
		if(pSNPs[nTagSNP_id].nflag & 1)
			printf("Error: the SNP has been selected as Tag SNP\n");

		gpSNP_ID_heads[gncur_highest_SCS] = gpSNP_ID_heads[gncur_highest_SCS]->pnext;
		if(gpSNP_ID_heads[gncur_highest_SCS]!=NULL)
			gpSNP_ID_heads[gncur_highest_SCS]->pprev = NULL;
		pSNP_ID_node->pnext = NULL;
		pSNP_ID_node->pprev = NULL;

		ptagSNPs[num_of_tag_SNPs++] = nTagSNP_id;

		UpdateCoverage(pSNPs, nstart_pos, nend_pos, nTagSNP_id, nuncovered_SNPs);

		if(num_of_tag_SNPs>0 && num_of_tag_SNPs%10000==0)
			fprintf(gfp_cvg, "%d  %d  %d  %.2f  %d\n", num_of_tag_SNPs, nSNP_num-nuncovered_SNPs, nuncovered_SNPs, (double)(nSNP_num-nuncovered_SNPs)/nSNP_num, gncur_highest_SCS);
		if(nSNP_num-nuncovered_SNPs>=(int)(nSNP_num*dcoverage))
		{
			fprintf(gfp_cvg, "%d  %d  %d  %.2f  %d\n", num_of_tag_SNPs, nSNP_num-nuncovered_SNPs, nuncovered_SNPs, (double)(nSNP_num-nuncovered_SNPs)/nSNP_num, gncur_highest_SCS);
			dcoverage += 0.05;
		}
	}
	fprintf(gfp_cvg, "%d  %d  %d  %.2f  %d\n", num_of_tag_SNPs, nSNP_num-nuncovered_SNPs, nuncovered_SNPs, (double)(nSNP_num-nuncovered_SNPs)/nSNP_num, gncur_highest_SCS);
	fprintf(gfp_cvg, "\n");

	for(i=nstart_pos;i<nend_pos;i++)
		pSNPs[i].prule_node_head = NULL;

	if(gdmax_tagSNP_num==0)
	{
		for(i=nstart_pos;i<nend_pos;i++)
		{
			if(!(pSNPs[i].nflag & 2))
				printf("Error: SNP %d is not covered\n", i);
		}
		for(i=gnmax_cand_RHS_exts;i>=1;i--)
		{
			if(gpSNP_ID_heads[i]!=NULL)
				printf("Error: the %d-th SNP ID node head is not empty\n", i);
		}
	}
	
	delete []gpSNP_ID_nodes;
	DecMemSize(sizeof(SNP_ID_NODE)*nSNP_num);
	delete []gpcovered_flag_buf;
	DecMemSize(sizeof(char)*gnwindow_size_sum);

	DelRuleBuf();
	DelRuleNodeBuf();

	//printf("#Tag SNPs selected: %d\n", num_of_tag_SNPs);

	return num_of_tag_SNPs;
}

void SelectTagSNPs(SNP *pSNPs, int num_of_SNPs, int nmem_size, char* szoutput_name)
{
	char szrule_filename[200], szcvg_filename[200];
	int nstart_pos, nend_pos, nprev_start_pos, i;
	unsigned int  ncur_file_pos, nfile_start_pos, nfile_end_pos, nrule_mem_size_start;
	int *ptagSNPs, *pcur_tagSNPs, num_of_new_tagSNPs;
	int *pexist_tagSNPs, nexist_tagSNP_num;

	gbno_equv_tag_SNPs = 1;

	sprintf(szcvg_filename, "%s.cvg.txt", szoutput_name);
	gfp_cvg = fopen(szcvg_filename, "wt");
	if(gfp_cvg==NULL)
	{
		printf("Error: cannot open file %s for write\n", szcvg_filename);
		return;
	}
	sprintf(szrule_filename, "%s.rule.bin", szoutput_name);
	gfp_rule = fopen(szrule_filename, "rb");
	if(gfp_rule==NULL)
	{
		printf("Error: cannot open file %s for read\n", szrule_filename);
		return;
	}
	gnrule_file_pos = 0;

	if(gnmax_RHS_len==0)
		gnmax_RHS_len = 100;

	gorule.pLHS_set = new int[gnmax_len];
	IncMemSize(sizeof(int)*gnmax_len);
	gorule.pRHS_set = new int[gnmax_RHS_len];
	IncMemSize(sizeof(int)*gnmax_RHS_len);
	gpunfolded_LHS_set = new int[gnmax_len];
	IncMemSize(sizeof(int)*gnmax_len);

	gptemp_merged_ids = new int[gnmax_merged_ids];
	IncMemSize(sizeof(int)*gnmax_merged_ids);

	gpSNP_ID_heads = new SNP_ID_NODE*[gnmax_cand_RHS_exts+1];
	IncMemSize(sizeof(SNP_ID_NODE*)*(gnmax_cand_RHS_exts+1));

	ptagSNPs = new int[num_of_SNPs];
	IncMemSize(sizeof(int)*num_of_SNPs);

	gnum_of_tag_SNPs = 0;
	gnnot_covered_tag_SNPs = 0;
	gninvalid_rule_num = 0;

	if(gdmax_tagSNP_num>0)
	{
		gnwindow_size_sum = 0;
		for(i=0;i<num_of_SNPs;i++)
			gnwindow_size_sum += gporig_window_end_pos[i]-gporig_window_start_pos[i];
		gnum_of_tag_SNPs = SelectTagSNPs(pSNPs, 0, num_of_SNPs, gntotal_rule_size, NULL, 0, ptagSNPs);
		gnum_of_partitions = 1;
	}
	else 
	{
		pexist_tagSNPs = new int[gnmax_cand_exts];
		IncMemSize(sizeof(int)*gnmax_cand_exts);
		nexist_tagSNP_num = 0;

		GetRuleMemSize(pSNPs, num_of_SNPs);
		nstart_pos = 0;
		nprev_start_pos = -1;
		nend_pos = 0;
		gnum_of_tag_SNPs = 0;
		gnum_of_partitions = 0;
		nfile_start_pos = 0;
		nfile_end_pos = 0;
		while(nend_pos<num_of_SNPs)
		{
			gnwindow_size_sum = 0;
			while(nend_pos<num_of_SNPs && !pSNPs[nend_pos].bhas_rules)
			{
				gnwindow_size_sum += (gporig_window_end_pos[nend_pos]-gporig_window_start_pos[nend_pos]);
				nend_pos++;
			}
			if(nend_pos<num_of_SNPs)
			{
				nfile_start_pos = pSNPs[nend_pos].nrule_start_pos;
				nrule_mem_size_start = pSNPs[nend_pos].nrule_mem_size_start;
			}
			nfile_end_pos = nfile_start_pos;
			
			while(nend_pos<num_of_SNPs && (nend_pos-nstart_pos<1000 || pSNPs[nend_pos].npos-pSNPs[nend_pos-1].npos<=gnwindow_len) && 
				(gnmem_size==0 || nend_pos-nstart_pos<1000 || !pSNPs[nend_pos].bhas_rules || 
				gnused_mem_size+
				(gnwindow_size_sum+gporig_window_end_pos[nend_pos]-gporig_window_start_pos[nend_pos])*sizeof(char)+
				(nend_pos-nstart_pos)*sizeof(SNP_ID_NODE)+
				pSNPs[nend_pos].nrule_mem_size_end-nrule_mem_size_start<(unsigned int)gnmem_size*(1<<20)))
			{
				gnwindow_size_sum += (gporig_window_end_pos[nend_pos]-gporig_window_start_pos[nend_pos]);
				if(pSNPs[nend_pos].bhas_rules && nfile_end_pos<pSNPs[nend_pos].nrule_end_pos)
					nfile_end_pos = pSNPs[nend_pos].nrule_end_pos;
				nend_pos++;
			}
			ncur_file_pos = ftell(gfp_rule);
			if(ncur_file_pos!=nfile_start_pos)
				fseek(gfp_rule, nfile_start_pos-ncur_file_pos, SEEK_CUR);
			gnrule_file_pos = nfile_start_pos;

			pcur_tagSNPs = &ptagSNPs[gnum_of_tag_SNPs]; 
			num_of_new_tagSNPs = SelectTagSNPs(pSNPs, nstart_pos, nend_pos, nfile_end_pos, pexist_tagSNPs, nexist_tagSNP_num, pcur_tagSNPs);
			gnum_of_tag_SNPs += num_of_new_tagSNPs;

			gnum_of_partitions++;

			if(nend_pos>=num_of_SNPs)
				break;

			nprev_start_pos = nstart_pos;
			nstart_pos = gporig_window_start_pos[nend_pos];
			if(nstart_pos==nprev_start_pos)
				nstart_pos = nend_pos;
			nexist_tagSNP_num = 0;
			for(i=nstart_pos;i<nend_pos;i++)
			{
				if(pSNPs[i].nflag & 1)
					pexist_tagSNPs[nexist_tagSNP_num++] = i;
			}
			nend_pos = nstart_pos;
		}
		delete []pexist_tagSNPs;
		DecMemSize(sizeof(int)*gnmax_cand_exts);

		if(gnum_of_partitions>1)
			gnum_of_tag_SNPs = RemoveDupTagSNPs(pSNPs, ptagSNPs, gnum_of_tag_SNPs);
	}
	fclose(gfp_cvg);
	fclose(gfp_rule);

	delete []gpSNP_ID_heads;
	DecMemSize(sizeof(SNP_ID_NODE*)*(gnmax_cand_RHS_exts+1));

	delete []gorule.pLHS_set;
	DecMemSize(sizeof(int)*gnmax_len);
	delete []gorule.pRHS_set;
	DecMemSize(sizeof(int)*gnmax_RHS_len);
	delete []gpunfolded_LHS_set;
	DecMemSize(sizeof(int)*gnmax_len);

	delete []gptemp_merged_ids;
	DecMemSize(sizeof(int)*gnmax_merged_ids);

	gncovered_SNPs = 0;
	for(i=0;i<num_of_SNPs;i++)
	{
		if(pSNPs[i].nflag!=0)
			gncovered_SNPs++;
	}
	printf("#Tag SNPs selected: %d\n", gnum_of_tag_SNPs);
	printf("#covered SNPs: %d, covered percentage: %.3f\n", gncovered_SNPs, (double)gncovered_SNPs/num_of_SNPs);

	OutputTagSNPs(pSNPs, num_of_SNPs, ptagSNPs, gnum_of_tag_SNPs, szoutput_name);

	delete []ptagSNPs;
	DecMemSize(sizeof(int)*num_of_SNPs);
}


int comp_int(const void *e1, const void *e2)
{
	int n1, n2;

	n1 = *(int*)e1;
	n2 = *(int*)e2;

	if(n1<n2)
		return -1;
	else if(n1>n2)
		return 1;
	else
		return 0;
}

int RemoveDupTagSNPs(SNP* pSNPs, int *ptagSNPs, int num_of_tagSNPs)
{
	int i, num_of_remain_tagSNPs;

	qsort(ptagSNPs, num_of_tagSNPs, sizeof(int), comp_int);

	num_of_remain_tagSNPs = 0;
	for(i=0;i<num_of_tagSNPs;i++)
	{
		if(i==0 || ptagSNPs[i]!=ptagSNPs[i-1])
			ptagSNPs[num_of_remain_tagSNPs++] = ptagSNPs[i];
	}

	//printf("#Removed duplicate tag SNPs: %d, #Remaining tag SNPs: %d\n", num_of_tagSNPs-num_of_remain_tagSNPs, num_of_remain_tagSNPs);

	return num_of_remain_tagSNPs;
}


