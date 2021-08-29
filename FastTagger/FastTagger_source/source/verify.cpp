#include <time.h>
#include <sys/timeb.h>

#include "global.h"

double gdverify_time;


float CalcR2(char* pbitmap1, int n1x, char* pbitmap2, int nx1, int &ncomb_bitmap)
{
	int i, n11;
	float dR2;

	n11 = 0;
	for(i=0;i<gnhapmap_num;i++)
	{
		if(pbitmap1[i]==1 && pbitmap2[i]==1)
			n11++;
	}
	dR2 = CalcR2(n11, n1x, nx1, gnhapmap_num, ncomb_bitmap);

	return dR2;
}

float CalcMultiMarkerR2(char* pbitmap1, int nSNP_num, char* pbitmap2, int nx1, int &ncomb_bitmap)
{
	int ncounter_num, i, n11, n1x, num_of_1s, num_of_0s, nflag, nmerged_1s, nmerged_0s, nmerged_bin_num;
	int nmaxR2_bin_no, nmaxR2_flag, pcounters[16];
	float dR2, dmax_R2;

	ncounter_num = 1<<nSNP_num;
	memset(pcounters, 0, sizeof(int)*(ncounter_num<<1));

	for(i=0;i<gnhapmap_num;i++)
		pcounters[(pbitmap1[i]<<1)|pbitmap2[i]]++;

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
			num_of_1s = pcounters[(i<<1)|1];
			num_of_0s = pcounters[i<<1];
			if(num_of_1s+num_of_0s>=gnmin_bin_size)
			{
				if(num_of_1s>num_of_0s || num_of_1s==num_of_0s && nx1<=gnhapmap_num-nx1)
				{
					n11 += num_of_1s;
					n1x += num_of_1s+num_of_0s;
					ncomb_bitmap |= (1<<i);
				}
			}
			else
			{
				nmerged_1s += num_of_1s;
				nmerged_0s += num_of_0s;
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
		n11 = 0;
		n1x = 0;
		dmax_R2 = 0;
		for(i=0;i<ncounter_num;i++)
		{
			n11 = pcounters[(i<<1)|1];
			n1x = n11+pcounters[i<<1];
			dR2 = CalcR2(n11, n1x, nx1, gnhapmap_num, nflag);
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


int GetRowNum(char* szfilename)
{
	FILE *fp;
	char ch;
	int num_of_rows;

	num_of_rows = 0;

	fp = fopen(szfilename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szfilename);
		return 0;
	}

	ch = fgetc(fp);
	while(!feof(fp))
	{
		while(!feof(fp) && ch!='\n')
			ch = fgetc(fp);
		num_of_rows++;

		ch = fgetc(fp);
	}
	fclose(fp);

	return num_of_rows;
}

void LoadSNPs(char* szmaf_filename, char* szSNPlineno_filename, double dmin_maf, SNP* pSNPs, int num_of_SNPs)
{
	FILE *fpmaf, *fplineno;
	char ch, allele1, allele2, szSNPid[100], szmaf[50];
	int nlen, nline_no, nSNP_pos, nSNP_no, npicked_lineno;
	double dmaf;

	fpmaf = fopen(szmaf_filename, "rt");
	if(fpmaf==NULL)
	{
		printf("Error: cannot open file %s for read\n", szmaf_filename);
		return;
	}
	fplineno = fopen(szSNPlineno_filename, "rt");
	if(fplineno==NULL)
	{
		printf("Error: cannot open file %s for read\n", szSNPlineno_filename);
		return;
	}

	nline_no = 0;
	nSNP_no = 0;

	fscanf(fplineno, "%d", &npicked_lineno);

	ch = fgetc(fpmaf);
	while(!feof(fpmaf))
	{
		nlen = 0;
		while(!feof(fpmaf) && ch!='\t' && ch!='\n')
		{
			szSNPid[nlen++]= ch;
			ch = fgetc(fpmaf);
		}
		szSNPid[nlen] = 0;

		while(!feof(fpmaf) && (ch=='\t' || ch==' '))
			ch = fgetc(fpmaf);

		if(ch<'0' || ch>'9')
			printf("Error: position expected\n");
		 
		nSNP_pos = 0;
		while(!feof(fpmaf) && ch>='0' && ch<='9')
		{
			nSNP_pos = nSNP_pos*10+ch-'0';
			ch = fgetc(fpmaf);
		}
		
		while(!feof(fpmaf) && (ch=='\t' || ch==' '))
			ch = fgetc(fpmaf);

		allele1 = ch;
		if(allele1!='A' && allele1!='C' && allele1!='G' && allele1!='T')
			printf("Error: 'A', 'C', 'G' or 'T' is expected here\n");

		ch = fgetc(fpmaf);
		while(!feof(fpmaf) && (ch=='\t' || ch==' '))
			ch = fgetc(fpmaf);
		allele2 = ch;
		if(allele2!='A' && allele2!='C' && allele2!='G' && allele2!='T')
			printf("Error: 'A', 'C', 'G' or 'T' is expected here\n");

		ch = fgetc(fpmaf);
		while(!feof(fpmaf) && (ch=='\t' || ch==' '))
			ch = fgetc(fpmaf);

		nlen = 0;
		while(!feof(fpmaf) && ch!='\n')
		{
			szmaf[nlen++] = ch;
			ch = fgetc(fpmaf);
		}
		szmaf[nlen] = 0;

		dmaf = atof(szmaf);
		if(nline_no==npicked_lineno)
		{
			if(dmaf<dmin_maf)
				printf("Error: the MAF of the SNP should be no less than %.3f\n", dmin_maf);
			pSNPs[nSNP_no].norder = nline_no;
			pSNPs[nSNP_no].npos = nSNP_pos;
			nSNP_no++;
			fscanf(fplineno, "%d", &npicked_lineno);
		}
		nline_no++;

		if(ch=='\n')
			ch = fgetc(fpmaf);
	}
	fclose(fpmaf);
	fclose(fplineno);

	if(nSNP_no!=num_of_SNPs)
		printf("Error: inconsistent number of SNPs\n");
}

void LoadBitmaps(char* szmatrix_filename, SNP* pSNPs, int num_of_SNPs)
{
	FILE *fp;
	int i, nline_no, nhapmap_no;
	char ch;

	fp = fopen(szmatrix_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szmatrix_filename);
		return;
	}

	nline_no = 0;
	for(i=0;i<num_of_SNPs;i++)
	{
		while(!feof(fp) && nline_no<pSNPs[i].norder)
		{
			ch = fgetc(fp);
			while(!feof(fp) && ch!='\n')
				ch = fgetc(fp);
			nline_no++;
		}
		if(feof(fp))
			break;

		pSNPs[i].pbitmap = new char[gnhapmap_num];
		IncMemSize(sizeof(char)*gnhapmap_num);

		pSNPs[i].n1 = 0;
		nhapmap_no = 0;
		ch = fgetc(fp);
		while(!feof(fp) && ch!='\n')
		{
			while(!feof(fp) && (ch<'0' || ch>'9') && ch!='\n')
				ch = fgetc(fp);
			if(ch>='0' && ch<='9')
			{
				pSNPs[i].pbitmap[nhapmap_no++] = ch-'0';
				if(ch=='1')
					pSNPs[i].n1++;
				else if(ch!='0')
					printf("Error: wrong allele\n");
				ch = fgetc(fp);
			}
		}
		if(nhapmap_no!=gnhapmap_num)
			printf("Error: inconsistent number of hapmaps at %d-th hapmap\n", i);
		nline_no++;
	}
	fclose(fp);

	if(i<num_of_SNPs)
		printf("Error: inconsistent number of SNPs\n");
}



void LoadMergedIds(char* szmergedSNPs_filename, SNP* pSNPs, int num_of_SNPs, int* pmerged_SNPid_buf)
{
	FILE *fp;
	int i, num_of_merged_ids, nid1, nid2, nbuf_pos;

	for(i=0;i<num_of_SNPs;i++)
	{
		pSNPs[i].num_of_merged_ids = -1;
		pSNPs[i].pmerged_ids = NULL;
	}

	fp = fopen(szmergedSNPs_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szmergedSNPs_filename);
		return;
	}
	nbuf_pos = 0;

	fscanf(fp, "%d", &num_of_merged_ids);
	while(!feof(fp))
	{
		fscanf(fp, "%d", &nid1);
		pSNPs[nid1].num_of_merged_ids = num_of_merged_ids;
		pSNPs[nid1].pmerged_ids = &pmerged_SNPid_buf[nbuf_pos];
		nbuf_pos += num_of_merged_ids;

		pSNPs[nid1].pmerged_ids[0] = nid1;
		for(i=1;i<num_of_merged_ids;i++)
		{
			fscanf(fp, "%d", &nid2);
			pSNPs[nid1].pmerged_ids[i] = nid2;
			if(nid2<0)
				pSNPs[-nid2].nSNPid_merged_with = nid1;
			else
				pSNPs[nid2].nSNPid_merged_with = nid1;
		}
		fscanf(fp, "%d", &num_of_merged_ids);
	}
	fclose(fp);

	if(nbuf_pos!=num_of_SNPs)
		printf("Error: inconsistant number of SNPs\n");

	for(i=0;i<num_of_SNPs;i++)
	{
		if(pSNPs[i].num_of_merged_ids==-1)
			printf("Error: SNP %d should be either merged with others or be itself\n", i);
	}
}


void VerifyMergedSNPs(SNP* pSNPs, int num_of_SNPs)
{
	int i, j, nmerged_id;

	for(i=0;i<num_of_SNPs;i++)
	{
		if(pSNPs[i].pmerged_ids!=NULL)
		{
			for(j=1;j<pSNPs[i].num_of_merged_ids;j++)
			{
				nmerged_id = pSNPs[i].pmerged_ids[j];
				if(nmerged_id>=0)
				{
					if(!IsIdentical(pSNPs[i].pbitmap, pSNPs[nmerged_id].pbitmap))
						printf("Error: the two SNPs should be identical\n");
				}
				else
				{
					nmerged_id = -nmerged_id;
					if(!IsReverseIdentical(pSNPs[i].pbitmap, pSNPs[nmerged_id].pbitmap))
						printf("Error: the two SNPs should be reverse identical\n");
				}
				delete []pSNPs[nmerged_id].pbitmap;
				pSNPs[nmerged_id].pbitmap = NULL;
			}
		}
		else if(pSNPs[i].pbitmap!=NULL)
			printf("Error: the bitmap should have been released\n");
	}
}

void VerifyRules(char* szrule_filename, SNP* pSNPs)
{
	FILE *fp;
	int nLHS_len, *pLHS_set, nRHS_len, nid1, nid2, ncomb_bitmap, nrule_comb_bitmap, i, j, ncounter_num;
	float dR2, drule_r2;
	char *ptemp_bitmap;

	fp = fopen(szrule_filename, "rb");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szrule_filename);
		return;
	}

	ncounter_num = 1<<(gnmax_len+1);
	gpmerged_bins = new int[ncounter_num];
	gpmerged_bin_rule_num = new int[gnmax_len+1];
	memset(gpmerged_bin_rule_num, 0, sizeof(int)*(gnmax_len+1));

	pLHS_set = new int[gnmax_len];
	ptemp_bitmap = new char[gnhapmap_num];

	fread(&nLHS_len, sizeof(int), 1, fp);
	while(!feof(fp))
	{
		if(nLHS_len==1)
		{
			fread(&nid1, sizeof(int), 1, fp);
			fread(&nid2, sizeof(int), 1, fp);
			fread(&dR2, sizeof(float), 1, fp);
			fread(&ncomb_bitmap, sizeof(int), 1, fp);

			drule_r2 = CalcR2(pSNPs[nid1].pbitmap, pSNPs[nid1].n1, pSNPs[nid2].pbitmap, pSNPs[nid2].n1, nrule_comb_bitmap);

			if(drule_r2!=dR2)
				printf("Error: inconsistent R2 value\n");
			if(nrule_comb_bitmap!=ncomb_bitmap)
				printf("Error: inconsistent flag bitmap\n");

		}
		else if(nLHS_len<=gnmax_len)
		{
			fread(pLHS_set, sizeof(int), nLHS_len, fp);
			memcpy(ptemp_bitmap, pSNPs[pLHS_set[0]].pbitmap, sizeof(char)*gnhapmap_num);
			for(i=1;i<nLHS_len;i++)
			{
				for(j=0;j<gnhapmap_num;j++)
					ptemp_bitmap[j] = (ptemp_bitmap[j]<<1)+pSNPs[pLHS_set[i]].pbitmap[j];
			}

			fread(&nRHS_len, sizeof(int), 1, fp);
			for(i=0;i<nRHS_len;i++)
			{
				fread(&nid2, sizeof(int), 1, fp);
				fread(&dR2, sizeof(float), 1, fp);
				fread(&ncomb_bitmap, sizeof(int), 1, fp);

				drule_r2 = CalcMultiMarkerR2(ptemp_bitmap, nLHS_len, pSNPs[nid2].pbitmap, pSNPs[nid2].n1, nrule_comb_bitmap);
				if(drule_r2!=dR2)
					printf("Error: inconsistent R2 value\n");
				if(nrule_comb_bitmap!=ncomb_bitmap)
					printf("Error: inconsistent flag bitmap\n");
			}
		}
		else 
			printf("Error: the length of left hand side cannot be larger than %d\n", gnmax_len);

		fread(&nLHS_len, sizeof(int), 1, fp);
	}
	fclose(fp);

	delete []pLHS_set;
	delete []ptemp_bitmap;
	delete []gpmerged_bins;
	delete []gpmerged_bin_rule_num;

}

int GetTagSNPs(char* sztagSNP_filename, char* pcovered)
{
	FILE *fp;
	int nid, num_of_tag_SNPs;

	fp = fopen(sztagSNP_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", sztagSNP_filename);
		return 0;
	}

	num_of_tag_SNPs = 0;

	fscanf(fp, "%d", &nid);
	while(!feof(fp))
	{
		pcovered[nid] = 3;
		num_of_tag_SNPs++;

		fscanf(fp, "%d", &nid);
	}
	fclose(fp);

	return num_of_tag_SNPs;
}

void GetTagSNPCoverage(char* sztagrule_filename, SNP *pSNPs, char *pcovered)
{
	FILE *fp;
	int nLHS_len, *pLHS_set, nRHS_len, nid1, nid2, i, j, k, nmerged_id, ncomb_bitmap, nrule_comb_bitmap;
	float dR2, drule_r2;
	char *ptemp_bitmap;
	int ncounter_num, nrule_no, nleft_boundary, nright_boundary;
	
	fp = fopen(sztagrule_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", sztagrule_filename);
		return;
	}

	ncounter_num = 1<<(gnmax_len+1);
	gpmerged_bins = new int[ncounter_num];
	gpmerged_bin_rule_num = new int[gnmax_len+1];
	memset(gpmerged_bin_rule_num, 0, sizeof(int)*(gnmax_len+1));

	pLHS_set = new int[gnmax_len];
	ptemp_bitmap = new char[gnhapmap_num];

	nrule_no = 0;

	fscanf(fp, "%d", &nLHS_len);
	while(!feof(fp))
	{
		for(i=0;i<nLHS_len;i++)
		{
			fscanf(fp, "%d", &pLHS_set[i]);
			if(pcovered[pLHS_set[i]]!=3)
				printf("Error: this SNP should be a tag SNP\n");
		}

		if(nLHS_len==1)
		{
			nid1 = pLHS_set[0];
			nleft_boundary = pSNPs[nid1].nwindow_start_pos;
			nright_boundary = pSNPs[nid1].nwindow_end_pos;

			if(pSNPs[nid1].pmerged_ids==NULL)
				nid1 = pSNPs[nid1].nSNPid_merged_with;
		}
		else
		{
			nid1 = pLHS_set[0];
			nleft_boundary = pSNPs[nid1].nwindow_start_pos;
			nright_boundary = pSNPs[nid1].nwindow_end_pos;
			if(pSNPs[nid1].pmerged_ids==NULL)
				nid1 = pSNPs[nid1].nSNPid_merged_with;
			memcpy(ptemp_bitmap, pSNPs[nid1].pbitmap, sizeof(char)*gnhapmap_num);
			for(i=1;i<nLHS_len;i++)
			{
				nid1 = pLHS_set[i];
				if(nleft_boundary<pSNPs[nid1].nwindow_start_pos)
					nleft_boundary = pSNPs[nid1].nwindow_start_pos;
				if(nright_boundary>pSNPs[nid1].nwindow_end_pos)
					nright_boundary = pSNPs[nid1].nwindow_end_pos;

				if(pSNPs[nid1].pmerged_ids==NULL)
					nid1 = pSNPs[nid1].nSNPid_merged_with;
				for(j=0;j<gnhapmap_num;j++)
					ptemp_bitmap[j] = (ptemp_bitmap[j]<<1)+pSNPs[nid1].pbitmap[j];
			}
		}


		fscanf(fp, "%d", &nRHS_len);
		for(i=0;i<nRHS_len;i++)
		{
			fscanf(fp, "%d", &nid2);
			fscanf(fp, "%f", &dR2);
			fscanf(fp, "%d", &ncomb_bitmap);

			if(pSNPs[nid2].pmerged_ids==NULL)
				printf("Error: the pmerged_ids field should not be NULL\n");

			if(nLHS_len==1)
				drule_r2 = CalcR2(pSNPs[nid1].pbitmap, pSNPs[nid1].n1, pSNPs[nid2].pbitmap, pSNPs[nid2].n1, nrule_comb_bitmap);
			else
				drule_r2 = CalcMultiMarkerR2(ptemp_bitmap, nLHS_len, pSNPs[nid2].pbitmap, pSNPs[nid2].n1, nrule_comb_bitmap);
			if(drule_r2-dR2>0.001)
				printf("Error: inconsistent R2 value\n");
			if(nrule_comb_bitmap!=ncomb_bitmap)
				printf("Error: inconsistent flag bitmap\n");

			for(k=0;k<pSNPs[nid2].num_of_merged_ids;k++)
			{
				nmerged_id = pSNPs[nid2].pmerged_ids[k];
				if(nmerged_id<0)
					nmerged_id = -nmerged_id;
				if(nmerged_id>=nleft_boundary && nmerged_id<nright_boundary && pcovered[nmerged_id]==0)
					pcovered[nmerged_id] = 2;
			}
		}
		nrule_no++;

		fscanf(fp, "%d", &nLHS_len);
	}
	fclose(fp);

	delete []pLHS_set;
	delete []ptemp_bitmap;
	delete []gpmerged_bins;
	delete []gpmerged_bin_rule_num;
}

void Verify(char *szmaf_filename, char *szmatrix_filename, char *szoutput_name)
{
	char szSNPlineno_filename[200], szmergedSNPs_filename[200], szrule_filename[200];
	char sztagSNP_filename[200], sztagrule_filename[200];
	int num_of_SNPs, *pmerged_SNPid_buf, i, j, num_of_tag_SNPs, nuncovered_SNP_num, nid, nmerged_id, ndistance;
	int nstart_pos, nend_pos;
	SNP *pSNPs;
	char *pcovered;
	struct timeb start, end;

	ftime(&start);

	printf("Verifying results...\n");

	sprintf(szSNPlineno_filename, "%s.SNPlineno.txt", szoutput_name);
	num_of_SNPs = GetRowNum(szSNPlineno_filename);
	pSNPs = new SNP[num_of_SNPs];

	LoadSNPs(szmaf_filename, szSNPlineno_filename, gdmin_maf, pSNPs, num_of_SNPs);

	nstart_pos = 0;
	nend_pos = 0;
	for(i=0;i<num_of_SNPs;i++)
	{
		while(pSNPs[i].npos-pSNPs[nstart_pos].npos>gnwindow_len)
			nstart_pos++;
		while(nend_pos<num_of_SNPs && pSNPs[nend_pos].npos-pSNPs[i].npos<=gnwindow_len)
			nend_pos++;
		pSNPs[i].nwindow_start_pos = nstart_pos;
		pSNPs[i].nwindow_end_pos = nend_pos;
	}

	sprintf(szmergedSNPs_filename, "%s.merged-ids.txt", szoutput_name);
	pmerged_SNPid_buf = new int[num_of_SNPs];
	LoadMergedIds(szmergedSNPs_filename, pSNPs, num_of_SNPs, pmerged_SNPid_buf);

	gnhapmap_num = GetHapMapNum(szmatrix_filename);
	LoadBitmaps(szmatrix_filename, pSNPs, num_of_SNPs);

	VerifyMergedSNPs(pSNPs, num_of_SNPs);

	sprintf(szrule_filename, "%s.rule.bin", szoutput_name);
	VerifyRules(szrule_filename, pSNPs);

	pcovered = new char[num_of_SNPs];
	memset(pcovered, 0, sizeof(char)*num_of_SNPs);
	sprintf(sztagSNP_filename, "%s.tagSNP.txt", szoutput_name);
	num_of_tag_SNPs = GetTagSNPs(sztagSNP_filename, pcovered);

	for(i=0;i<num_of_SNPs;i++)
	{
		if(pcovered[i]==3)
		{
			nid = i;
			if(pSNPs[nid].pmerged_ids==NULL)
				nid = pSNPs[nid].nSNPid_merged_with;
			for(j=0;j<pSNPs[nid].num_of_merged_ids;j++)
			{
				nmerged_id = pSNPs[nid].pmerged_ids[j];
				if(nmerged_id<0)
					nmerged_id = -nmerged_id;
				if(nmerged_id!=i)
				{
					ndistance = pSNPs[nmerged_id].npos-pSNPs[i].npos;
					if(ndistance<0)
						ndistance = -ndistance;
					if(ndistance<=gnwindow_len && pcovered[nmerged_id]==0)
						pcovered[nmerged_id] = 2;
				}
			}
		}
	}

	sprintf(sztagrule_filename, "%s.tagrule.txt", szoutput_name);
	GetTagSNPCoverage(sztagrule_filename, pSNPs, pcovered);

	nuncovered_SNP_num = 0;
	for(i=0;i<num_of_SNPs;i++)
	{
		if(pcovered[i]==0)
		{
			nuncovered_SNP_num++;
			//printf("Error: SNP %d is not covered\n", i);
		}
	}
	printf("#tag SNPs: %d, #uncovered SNPs: %d\n", num_of_tag_SNPs, nuncovered_SNP_num);
	printf("\n");

	delete []pcovered;

	for(i=0;i<num_of_SNPs;i++)
	{
		if(pSNPs[i].pbitmap!=NULL)
			delete []pSNPs[i].pbitmap;
	}

	delete []pSNPs;
	delete []pmerged_SNPid_buf;

	ftime(&end);
	gdverify_time = end.time-start.time+(double)(end.millitm-start.millitm)/1000;

}

