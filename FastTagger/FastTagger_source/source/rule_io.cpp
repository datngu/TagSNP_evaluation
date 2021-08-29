#include "global.h"


FILE *gfp_rule;
unsigned int gnrule_file_pos;
int gntotal_rule_size;
int gnfile_no;

RULE gorule;


void OutputRule(SNP *pSNPs, int nSNP_id1, int nSNP_id2, float dR2, int ncomb_bitmap)
{
	int nlen;

	gprule_num[1]++;

	nlen = 1;
	fwrite(&nlen, sizeof(int), 1, gfp_rule);
	fwrite(&nSNP_id1, sizeof(int), 1, gfp_rule);
	fwrite(&nSNP_id2, sizeof(int), 1, gfp_rule);
	fwrite(&dR2, sizeof(float), 1, gfp_rule);
	fwrite(&ncomb_bitmap, sizeof(int), 1, gfp_rule);
	gnrule_file_pos += sizeof(int)*4+sizeof(float);


}

void OutputRule(SNP *pSNPs, int nLHS_len, int *pLHS_set, int nRHS_len, RHS_INFO_NODE *pRHS_set)
{
	int i;

	gpLHS_num[nLHS_len]++;
	gprule_num[nLHS_len] += nRHS_len;

	if(gnmax_RHS_len<nRHS_len)
		gnmax_RHS_len = nRHS_len;

	fwrite(&nLHS_len, sizeof(int), 1, gfp_rule);
	fwrite(pLHS_set, sizeof(int), nLHS_len, gfp_rule);
	gnrule_file_pos += sizeof(int)+sizeof(int)*nLHS_len;

	fwrite(&nRHS_len, sizeof(int), 1, gfp_rule);
	gnrule_file_pos += sizeof(int);
	for(i=0;i<nRHS_len;i++)
	{
		fwrite(&pRHS_set[i].nSNPid, sizeof(int), 1, gfp_rule);
		fwrite(&pRHS_set[i].dR2, sizeof(float), 1, gfp_rule);
		fwrite(&pRHS_set[i].ncomb_bitmap, sizeof(int), 1, gfp_rule);
		gnrule_file_pos += sizeof(int)+sizeof(int)+sizeof(float);
	}
}

RULE* ReadInOneRule(unsigned int nrule_end_pos)
{
	float dR2;
	int i, ncomb_bitmap;

	gorule.nLHS_len = 0;

	if(nrule_end_pos>0 && gnrule_file_pos>=nrule_end_pos)
		return NULL;

	fread(&gorule.nLHS_len, sizeof(int), 1, gfp_rule);
	if(feof(gfp_rule))
		return NULL;

	if(gorule.nLHS_len==1)
	{
		fread(&gorule.pLHS_set[0], sizeof(int), 1, gfp_rule);
		gorule.nRHS_len = 1;
		fread(&gorule.pRHS_set[0], sizeof(int), 1, gfp_rule);
		fread(&dR2, sizeof(float), 1, gfp_rule);
		fread(&ncomb_bitmap, sizeof(int), 1, gfp_rule);
		gnrule_file_pos += sizeof(int)*4+sizeof(float);
	}
	else
	{
		fread(gorule.pLHS_set, sizeof(int), gorule.nLHS_len, gfp_rule);
		gnrule_file_pos += sizeof(int)+sizeof(int)*gorule.nLHS_len;

		fread(&gorule.nRHS_len, sizeof(int), 1, gfp_rule);
		gnrule_file_pos += sizeof(int);
		for(i=0;i<gorule.nRHS_len;i++)
		{
			fread(&gorule.pRHS_set[i], sizeof(int), 1, gfp_rule);
			fread(&dR2, sizeof(float), 1, gfp_rule);
			fread(&ncomb_bitmap, sizeof(int), 1, gfp_rule);
			gnrule_file_pos += sizeof(int)+sizeof(int)+sizeof(float);
		}
	}

	if(gorule.nLHS_len==0)
		return NULL;

	return &gorule;
}

void ReadInOneRule(FILE *fp, int &nLHS_len, int *pLHS_set, int &nRHS_len, RHS_INFO_NODE *pRHS_set)
{
	int i;

	if(feof(fp))
		return;

	fread(&nLHS_len, sizeof(int), 1, fp);
	if(feof(fp))
	{
		nLHS_len = 0;
		return;
	}

	if(nLHS_len>gnmax_len)
		printf("Error: the length of LHS cannot be larger than %d\n", gnmax_len);
	if(nLHS_len==1)
	{
		fread(&pLHS_set[0], sizeof(int), 1, fp);
		nRHS_len = 1;
		fread(&pRHS_set[0].nSNPid, sizeof(int), 1, fp);
		fread(&pRHS_set[0].dR2, sizeof(float), 1, fp);
		fread(&pRHS_set[0].ncomb_bitmap, sizeof(int), 1, fp);
	}
	else
	{
		fread(pLHS_set, sizeof(int), nLHS_len, fp);

		fread(&nRHS_len, sizeof(int), 1, fp);
		if(nRHS_len>gnmax_RHS_len)
			printf("Error: the length of LHS cannot be larger than %d\n", gnmax_RHS_len);
		for(i=0;i<nRHS_len;i++)
		{
			fread(&pRHS_set[i].nSNPid, sizeof(int), 1, fp);
			fread(&pRHS_set[i].dR2, sizeof(float), 1, fp);
			fread(&pRHS_set[i].ncomb_bitmap, sizeof(int), 1, fp);
		}
	}
}

void GetRuleMemSize(SNP *pSNPs, int num_of_SNPs)
{
	RULE *prule;
	int i, nrule_mem_size;
	unsigned int nb4read_pos;

	rewind(gfp_rule);

	for(i=0;i<num_of_SNPs;i++)
	{
		pSNPs[i].bhas_rules = false;
		pSNPs[i].nrule_start_pos = -1;
		pSNPs[i].nrule_end_pos = -1;
		pSNPs[i].nrule_mem_size_start = -1;
		pSNPs[i].nrule_mem_size_end = -1;
	}

	gnrule_file_pos = 0;
	nb4read_pos = gnrule_file_pos;
	nrule_mem_size = 0;
	prule = ReadInOneRule(0);
	while(prule!=NULL)
	{
		for(i=0;i<prule->nLHS_len;i++)
		{
			if(!pSNPs[prule->pLHS_set[i]].bhas_rules)
			{
				pSNPs[prule->pLHS_set[i]].bhas_rules = true;
				pSNPs[prule->pLHS_set[i]].nrule_start_pos = nb4read_pos;
				pSNPs[prule->pLHS_set[i]].nrule_mem_size_start = nrule_mem_size;
			}
		}
		for(i=0;i<prule->nRHS_len;i++)
		{
			if(!pSNPs[prule->pRHS_set[i]].bhas_rules)
			{
				pSNPs[prule->pRHS_set[i]].bhas_rules = true;
				pSNPs[prule->pRHS_set[i]].nrule_start_pos = nb4read_pos;
				pSNPs[prule->pRHS_set[i]].nrule_mem_size_start = nrule_mem_size;
			}
		}

		nrule_mem_size += sizeof(RULE)+sizeof(int)*(prule->nLHS_len+prule->nRHS_len);
		nrule_mem_size += sizeof(RULE_NODE)*prule->nLHS_len;

		for(i=0;i<prule->nLHS_len;i++)
		{
			pSNPs[prule->pLHS_set[i]].nrule_end_pos = gnrule_file_pos;
			pSNPs[prule->pLHS_set[i]].nrule_mem_size_end = nrule_mem_size;
		}
		for(i=0;i<prule->nRHS_len;i++)
		{
			pSNPs[prule->pRHS_set[i]].nrule_end_pos = gnrule_file_pos;
			pSNPs[prule->pRHS_set[i]].nrule_mem_size_end = nrule_mem_size;
		}

		nb4read_pos = gnrule_file_pos;
		prule = ReadInOneRule(0);
	}
	rewind(gfp_rule);

	for(i=0;i<num_of_SNPs;i++)
	{
		if(pSNPs[i].nwindow_start_pos==-1 && pSNPs[pSNPs[i].nSNPid_merged_with].bhas_rules)
		{
			pSNPs[i].bhas_rules = true;
			pSNPs[i].nrule_start_pos = pSNPs[pSNPs[i].nSNPid_merged_with].nrule_start_pos;
			pSNPs[i].nrule_end_pos = pSNPs[pSNPs[i].nSNPid_merged_with].nrule_end_pos;
			pSNPs[i].nrule_mem_size_start = pSNPs[pSNPs[i].nSNPid_merged_with].nrule_mem_size_start;
			pSNPs[i].nrule_mem_size_end = pSNPs[pSNPs[i].nSNPid_merged_with].nrule_mem_size_end;
		}
		if(i>0 && pSNPs[i].bhas_rules && pSNPs[i-1].bhas_rules && pSNPs[i-1].nrule_end_pos>pSNPs[i].nrule_end_pos)
		{
			pSNPs[i].nrule_end_pos = pSNPs[i-1].nrule_end_pos;
			pSNPs[i].nrule_mem_size_end = pSNPs[i-1].nrule_mem_size_end;
		}
	}
	for(i=num_of_SNPs-1;i>0;i--)
	{
		if(pSNPs[i].bhas_rules && pSNPs[i-1].bhas_rules && pSNPs[i-1].nrule_start_pos>pSNPs[i].nrule_start_pos)
		{
			pSNPs[i-1].nrule_start_pos = pSNPs[i].nrule_start_pos;
			pSNPs[i-1].nrule_mem_size_start = pSNPs[i].nrule_mem_size_start;
		}
	}

}


void OutputOneTagRule(SNP* pSNPs, int *pLHS_set, int nLHS_len, RHS_INFO_NODE *pRHS_set, int nRHS_len, int *punfolded_LHS_set, int ncur_pos, int nleft_boundary, int nright_boundary, FILE *fpout, int &nrule_num)
{
	int i, j, ncur_sid, nmerged_id;

	ncur_sid = pLHS_set[ncur_pos];
	if(pSNPs[ncur_sid].pmerged_ids==NULL)
		printf("Error: the pmerged_ids field should not be NULL\n");

	for(i=0;i<pSNPs[ncur_sid].num_of_merged_ids;i++)
	{
		nmerged_id = pSNPs[ncur_sid].pmerged_ids[i];
		if(pSNPs[nmerged_id].nflag==1)
		{
			if(nmerged_id>=nleft_boundary && nmerged_id<nright_boundary)
			{
				punfolded_LHS_set[ncur_pos] = nmerged_id;
				if(ncur_pos+1==nLHS_len)
				{
					fprintf(fpout, "%d ", nLHS_len);
					for(j=0;j<nLHS_len;j++)
						fprintf(fpout, "%d ", punfolded_LHS_set[j]);
					fprintf(fpout, "\n");
					fprintf(fpout, "%d ", nRHS_len);
					for(j=0;j<nRHS_len;j++)
						fprintf(fpout, "%d %.3f %d\t", pRHS_set[j].nSNPid, pRHS_set[j].dR2, pRHS_set[j].ncomb_bitmap);
					fprintf(fpout, "\n");
					nrule_num++;
				}
				else
				{
					if(nleft_boundary<gporig_window_start_pos[nmerged_id])
						nleft_boundary = gporig_window_start_pos[nmerged_id];
					if(nright_boundary>gporig_window_end_pos[nmerged_id])
						nright_boundary = gporig_window_end_pos[nmerged_id];
					OutputOneTagRule(pSNPs, pLHS_set, nLHS_len, pRHS_set, nRHS_len, punfolded_LHS_set, ncur_pos+1, nleft_boundary, nright_boundary, fpout, nrule_num);
				}
			}
		}
	}
}


void OutputTagSNPs(SNP* pSNPs, int num_of_SNPs, int *ptagSNPs, int num_of_tagSNPs, char* szoutput_name)
{
	FILE *fp, *fpout;
	char sztagSNP_filename[200], szrulebin_filename[200], sztag_rule_filename[200];
	int i, num_of_rules, nLHS_len, *pLHS_set, nRHS_len, *punfolded_LHS_set;
	RHS_INFO_NODE *pRHS_set;
	int nid, nmerged_id, num_of_tag_rules, nrule_num;

	//-------------------------------------------------------
	// output tag SNP id
	sprintf(sztagSNP_filename, "%s.tagSNP.txt", szoutput_name);
	fpout = fopen(sztagSNP_filename, "wt");
	if(fpout==NULL)
	{
		printf("Error: cannot open file %s for write\n", sztagSNP_filename);
		return;		
	}
	for(i=0;i<num_of_tagSNPs;i++)
		fprintf(fpout, "%d\n", ptagSNPs[i]);
	fclose(fpout);
	//----------------------------------------------------


	//==================================================
	//output rules
	sprintf(sztag_rule_filename, "%s.tagrule.txt", szoutput_name);
	fpout = fopen(sztag_rule_filename, "wt");
	if(fpout==NULL)
	{
		printf("Error: cannot open file %s for write\n", sztag_rule_filename);
		return;		
	}
	sprintf(szrulebin_filename, "%s.rule.bin", szoutput_name);
	fp = fopen(szrulebin_filename, "rb");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szrulebin_filename);
		return;
	}
	gnrule_file_pos = 0;

	pLHS_set = new int[gnmax_len];
	IncMemSize(sizeof(int)*gnmax_len);
	pRHS_set = new RHS_INFO_NODE[gnmax_RHS_len];
	IncMemSize(sizeof(RHS_INFO_NODE)*gnmax_RHS_len);
	punfolded_LHS_set = new int[gnmax_len];
	IncMemSize(sizeof(int)*gnmax_len);

	for(i=0;i<num_of_SNPs;i++)
		pSNPs[i].nflag = 0;
	for(i=0;i<num_of_tagSNPs;i++)
		pSNPs[ptagSNPs[i]].nflag = 1;

	num_of_tag_rules = 0;

	num_of_rules = 0;
	ReadInOneRule(fp, nLHS_len, pLHS_set, nRHS_len, pRHS_set);
	while(nLHS_len>0)
	{
		if(nLHS_len==1)
		{
			if(nRHS_len!=1)
				printf("Error: RHS_len should be 1 too\n");
			nid = pLHS_set[0];
			for(i=0;i<pSNPs[nid].num_of_merged_ids;i++)
			{
				nmerged_id = pSNPs[nid].pmerged_ids[i];
				if(pSNPs[nmerged_id].nflag==1)
				{
					fprintf(fpout, "1 %d\n", nmerged_id);
					fprintf(fpout, "1 ");
					fprintf(fpout, "%d %.3f %d\t", pRHS_set[0].nSNPid, pRHS_set[0].dR2, pRHS_set[0].ncomb_bitmap);
					fprintf(fpout, "\n");
					num_of_tag_rules++;
				}
			}

			nid = pRHS_set[0].nSNPid;
			for(i=0;i<pSNPs[nid].num_of_merged_ids;i++)
			{
				nmerged_id = pSNPs[nid].pmerged_ids[i];
				if(pSNPs[nmerged_id].nflag==1)
				{
					fprintf(fpout, "1 %d\n", nmerged_id);
					fprintf(fpout, "1 ");
					fprintf(fpout, "%d %.3f %d\t", pLHS_set[0], pRHS_set[0].dR2, pRHS_set[0].ncomb_bitmap);
					fprintf(fpout, "\n");
					num_of_tag_rules++;
				}
			}
		}
		else 
		{
			nrule_num = 0;
			OutputOneTagRule(pSNPs, pLHS_set, nLHS_len, pRHS_set, nRHS_len, punfolded_LHS_set, 0, 0, num_of_SNPs, fpout, nrule_num);
			num_of_tag_rules += nrule_num;
		}
		num_of_rules++;
		ReadInOneRule(fp, nLHS_len, pLHS_set, nRHS_len, pRHS_set);
	}
	fclose(fpout);
	fclose(fp);

	delete []pLHS_set;
	DecMemSize(sizeof(int)*gnmax_len);
	delete []pRHS_set;
	DecMemSize(sizeof(RHS_INFO_NODE)*gnmax_RHS_len);
	delete []punfolded_LHS_set;
	DecMemSize(sizeof(int)*gnmax_len);

	printf("#tagging rules: %d\n\n", num_of_tag_rules);

}




