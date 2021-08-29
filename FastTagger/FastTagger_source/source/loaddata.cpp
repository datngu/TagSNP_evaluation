#include <string>
#include <vector>
#include <iostream>
using namespace std;

#include "global.h"

int gnhapmap_num;
int gnmax_n1;

int *gporig_window_start_pos;
int *gporig_window_end_pos;
int gnwindow_size_sum;

int gnremoved_SNP_num;
int gnrep_SNP_num;

int gnnonb_SNP_num;

int gnmax_merged_ids;

void fix_string(char *s)
{
	for (int i = 0; ; ++i)
	{
		if ((int)s[i] == 13) 
		{
			s[i] = '\0';
		}
		if ((int)s[i] == 0) break;
	}
}

int LoadSNPs(char* szmaf_filename, double dmin_maf, SNP* &pSNPs, char* szoutput_name)
{
	FILE *fp, *fpout, *fpid;
	//string in_file = *szmaf_filename;
	char szSNPlineno_filename[200], szSNPid_filename[200], ch, szSNPid[100], szmaf[50];
	int nlen, nline_no, num_of_SNPs, i, nprev_pos, nmax_dist, nmin_dist, nprev_distance, ndistance;
	double dmaf, davg_distance;
	vector<SNP> vecSNPs;
	SNP oneSNP;
	char allele1, allele2;
	//fix_string(szmaf_filename);
	fp = fopen( szmaf_filename, "rt");
	if(fp==NULL)
	{
		cout << "Error: cannot open file: |" << szmaf_filename << "|" << endl;
		//printf("Error: cannot open file %s for read\n", szmaf_filename);
		return 0;
	}

	sprintf(szSNPlineno_filename, "%s.SNPlineno.txt", szoutput_name);
	fpout = fopen(szSNPlineno_filename, "wt");
	if(fpout==NULL)
	{
		printf("Error: cannot open file %s for write\n", szSNPlineno_filename);
		return 0;
	}
	sprintf(szSNPid_filename, "%s.SNPid.txt", szoutput_name);
	fpid = fopen(szSNPid_filename, "wt");
	if(fpid==NULL)
	{
		printf("Error: cannot open file %s for write\n", szSNPid_filename);
		return 0;
	}

	oneSNP.pbitmap = NULL;
	oneSNP.n1 = 0;
	oneSNP.p1sid_list = NULL;
	oneSNP.nbeing_covered_times = 0;
	oneSNP.num_of_covered_SNPs = 0;
	oneSNP.num_of_merged_ids = 0;
	oneSNP.pmerged_ids = NULL;
	oneSNP.nrule_start_pos = -1;
	oneSNP.nrule_mem_size_start = -1;
	oneSNP.nrule_mem_size_end = -1;

	nline_no = 0;

	davg_distance = 0;
	ndistance = 0;
	nprev_distance = gnwindow_len+1;
	nmax_dist = 0 ;
	nmin_dist = -1;
	nprev_pos = -1;
	gnnonb_SNP_num = 0;

	ch = fgetc(fp);
	while(!feof(fp))
	{
		nlen = 0;
		while(!feof(fp) && ch!='\t' && ch!='\n')
		{
			szSNPid[nlen++]= ch;
			ch = fgetc(fp);
		}
		szSNPid[nlen] = 0;

		while(!feof(fp) && (ch=='\t' || ch==' '))
			ch = fgetc(fp);

		if(ch<'0' || ch>'9')
			printf("Error: position expected\n");
		 
		oneSNP.npos = 0;
		while(!feof(fp) && ch>='0' && ch<='9')
		{
			oneSNP.npos = oneSNP.npos*10+ch-'0';
			ch = fgetc(fp);
		}
		oneSNP.norder = nline_no;
		
		while(!feof(fp) && (ch=='\t' || ch==' '))
			ch = fgetc(fp);
		allele1 = ch;
		if(allele1!='A' && allele1!='C' && allele1!='G' && allele1!='T')
			printf("Error: 'A', 'C', 'G' or 'T' is expected here\n");

		ch = fgetc(fp);
		while(!feof(fp) && (ch=='\t' || ch==' '))
			ch = fgetc(fp);
		allele2 = ch;
		if(allele2!='A' && allele2!='C' && allele2!='G' && allele2!='T')
			printf("Error: 'A', 'C', 'G' or 'T' is expected here\n");

		ch = fgetc(fp);
		while(!feof(fp) && (ch=='\t' || ch==' '))
			ch = fgetc(fp);
		nlen = 0;
		while(!feof(fp) && ch!='\n')
		{
			szmaf[nlen++] = ch;
			ch = fgetc(fp);
		}
		szmaf[nlen] = 0;

		dmaf = atof(szmaf);
		if(dmaf>=dmin_maf)
		{
			vecSNPs.push_back(oneSNP);
			fprintf(fpid, "%s\t%d\t%c\t%c\n", szSNPid, oneSNP.npos, allele1, allele2);
			fprintf(fpout, "%d\n", nline_no);

			if(nprev_pos>=0)
			{
				ndistance = oneSNP.npos-nprev_pos;
				davg_distance += ndistance;
				if(nmin_dist==-1 || ndistance<nmin_dist)
					nmin_dist = ndistance;
				if(ndistance>nmax_dist)
					nmax_dist = ndistance;
				if(nprev_distance>gnwindow_len && ndistance>gnwindow_len)
					gnnonb_SNP_num++;
				nprev_distance = ndistance;
			}
			nprev_pos = oneSNP.npos;
		}

		nline_no++;

		if(ch=='\n')
			ch = fgetc(fp);
	}
	fclose(fp);
	fclose(fpout);
	fclose(fpid);
	if(ndistance>gnwindow_len)
		gnnonb_SNP_num++;

	num_of_SNPs = (int)vecSNPs.size();
	if(num_of_SNPs>0)
	{
		pSNPs = new SNP[num_of_SNPs];
		IncMemSize(num_of_SNPs*sizeof(SNP));
		for(i=0;i<num_of_SNPs;i++)
			pSNPs[i] = vecSNPs[i];
	}

	davg_distance /= (nline_no-1);

	printf("Distance between SNPs: avg=%.2f, min=%d, max=%d\n", davg_distance, nmin_dist, nmax_dist);
	printf("#selected SNPs: %d\n", num_of_SNPs);

	return num_of_SNPs;
}


void CalcWindowSize(SNP *pSNPs, int num_of_SNPs, int nwindow_len)
{
	int i, nstart_pos, nend_pos;

	gnmax_cand_exts = 0;
	gnmax_cand_RHS_exts = 0;
	gdavg_cand_exts = 0;
	gdavg_cand_RHS_exts = 0;
	gnwindow_size_sum = 0;

	nstart_pos = 0;
	nend_pos = 0;
	for(i=0;i<num_of_SNPs;i++)
	{
		while(pSNPs[i].npos-pSNPs[nstart_pos].npos>nwindow_len)
			nstart_pos++;
		while(nend_pos<num_of_SNPs && pSNPs[nend_pos].npos-pSNPs[i].npos<=nwindow_len)
			nend_pos++;
		pSNPs[i].nwindow_start_pos = nstart_pos;
		pSNPs[i].nwindow_end_pos = nend_pos;

		gporig_window_start_pos[i] = nstart_pos;
		gporig_window_end_pos[i] = nend_pos;

		gdavg_cand_exts += i-nstart_pos;
		gdavg_cand_RHS_exts += nend_pos-nstart_pos-1;

		if(gnmax_cand_exts < i-nstart_pos)
			gnmax_cand_exts = i-nstart_pos;
		if(gnmax_cand_RHS_exts < nend_pos-nstart_pos-1)
			gnmax_cand_RHS_exts = nend_pos-nstart_pos-1;
		gnwindow_size_sum += nend_pos-nstart_pos;
	}

	gdavg_cand_exts /= num_of_SNPs;
	gdavg_cand_RHS_exts /= num_of_SNPs;

	//printf("#candidate extensions: avg=%.2f, max=%d\n", gdavg_cand_exts, gnmax_cand_exts);
	//printf("#candidate RHS extensions: avg=%.2f, max=%d\n", gdavg_cand_RHS_exts, gnmax_cand_RHS_exts);
	printf("\n");

}

int GetHapMapNum(char *szhapmap_filename)
{
	FILE *fp;
	char ch;
	int num_of_hapmaps;

	fp = fopen(szhapmap_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szhapmap_filename);
		return 0;
	}

	num_of_hapmaps = 0;

	ch = fgetc(fp);
	while(!feof(fp) && ch!='\n')
	{
		while(!feof(fp) && (ch<'0' || ch>'9') && ch!='\n')
			ch = fgetc(fp);
		if(ch>='0' && ch<='9')
		{
			num_of_hapmaps++;
			ch = fgetc(fp);
		}
	}
	fclose(fp);

	return num_of_hapmaps;
}

void LoadOneBitmap(FILE *fp, int &nline_no, SNP* pSNPs, int ncur_SNP_no)
{
	int nhapmap_no;
	char ch;

	while(!feof(fp) && nline_no<pSNPs[ncur_SNP_no].norder)
	{
		ch = fgetc(fp);
		while(!feof(fp) && ch!='\n')
			ch = fgetc(fp);
		nline_no++;
	}
	if(feof(fp))
	{
		printf("Error: fail to read bitmap for SNP %d\n", ncur_SNP_no);
		return;
	}

	nhapmap_no = 0;
	ch = fgetc(fp);
	while(!feof(fp) && ch!='\n')
	{
		while(!feof(fp) && (ch<'0' || ch>'9') && ch!='\n')
			ch = fgetc(fp);
		if(ch>='0' && ch<='9')
		{
			pSNPs[ncur_SNP_no].pbitmap[nhapmap_no++] = ch-'0';
			if(ch=='1')
				pSNPs[ncur_SNP_no].n1++;
			else if(ch!='0')
				printf("Error: wrong allele\n");
			ch = fgetc(fp);
		}
	}
	if(gnmax_n1<pSNPs[ncur_SNP_no].n1)
		gnmax_n1 = pSNPs[ncur_SNP_no].n1;

	if(pSNPs[ncur_SNP_no].n1==0 || pSNPs[ncur_SNP_no].n1==gnhapmap_num)
		printf("Warning: the number of 1 or 0 of SNP %d is 0\n", pSNPs[ncur_SNP_no].norder);
	if(nhapmap_no!=gnhapmap_num)
		printf("Error: inconsistent number of hapmaps at %d-th hapmap\n", ncur_SNP_no);
	nline_no++;

}

int MergeSNPs(char* szmatrix_filename, SNP *pSNPs, int num_of_SNPs, REP_SNP* &pRep_SNPs, char* szmerged_ids_filename)
{
	FILE *fpout, *fpmatrix;
	int i, j, num_of_merged, nadjust_pos, num_of_cand_exts, num_of_cand_RHS_exts, nleft_boundary;
	int nline_no, ncur_load_pos, ncur_null_pos;

	gnmerge_max_cand_exts = 0;
	gnmerge_max_cand_RHS_exts = 0;
	gdmerge_avg_cand_exts = 0;
	gdmerge_avg_cand_RHS_exts = 0;

	gnremoved_SNP_num = 0;
	gnrep_SNP_num = num_of_SNPs;
	gnmerged_id_buf_pos = 0;

	gnmax_merged_ids = 0;
	gnmax_n1 = 0;

	//fix_string(szmatrix_filename);

	fpmatrix = fopen(szmatrix_filename, "rt");
	if(fpmatrix==NULL)
	{
		printf("Error: cannot open file %s for read\n", szmatrix_filename);
		return 0;
	}	
	
	fpout = fopen(szmerged_ids_filename, "wt");
	if(fpout==NULL)
	{
		printf("Error: cannot open file %s for write\n", szmerged_ids_filename);
		return 0;
	}

	ncur_load_pos = 0;
	nline_no = 0;
	ncur_null_pos = 0;

	for(i=0;i<num_of_SNPs;i++)
	{
		if(i>ncur_load_pos)
			printf("Error: i cannot be larger than ncur_load_pos\n");
		if(i==ncur_load_pos)
		{
			if(pSNPs[i].pbitmap==NULL)
			{
				if(ncur_null_pos!=i)
					printf("Error: ncur_null_pos should be equal to i\n");
				pSNPs[i].pbitmap = NewBitmap();
				ncur_null_pos++;
			}
			LoadOneBitmap(fpmatrix, nline_no, pSNPs, i);
			ncur_load_pos++;
		}
		if(pSNPs[i].pbitmap!=NULL)
		{
			pSNPs[i].pmerged_ids = &gpmerged_SNP_ids[gnmerged_id_buf_pos];
			pSNPs[i].pmerged_ids[0] = i;
			num_of_merged = 1;

			nadjust_pos = pSNPs[i].nwindow_end_pos;

			for(j=i+1;j<num_of_SNPs && pSNPs[j].npos-pSNPs[i].npos<=gnmax_merge_window_len;j++)
			{
				if(j>ncur_load_pos)
					printf("Error: j cannot be larger than ncur_load_pos\n");
				if(j==ncur_load_pos)
				{
					if(pSNPs[j].pbitmap==NULL)
					{
						if(j!=ncur_null_pos)
							printf("Error: ncur_null_pos should be equal to j\n");
						pSNPs[j].pbitmap = NewBitmap();
						ncur_null_pos++;
					}
					LoadOneBitmap(fpmatrix, nline_no, pSNPs, j);
					ncur_load_pos++;
				}
				if(pSNPs[j].pbitmap!=NULL)
				{
					if(IsIdentical(pSNPs[i].pbitmap, pSNPs[j].pbitmap))
					{
						if(nadjust_pos<pSNPs[j].nwindow_start_pos)
							nadjust_pos = pSNPs[j].nwindow_start_pos;
						while(nadjust_pos<pSNPs[j].nwindow_end_pos)
						{
							if(pSNPs[nadjust_pos].nwindow_start_pos>i)
								pSNPs[nadjust_pos].nwindow_start_pos = i;
							nadjust_pos++;
						}
						pSNPs[i].nwindow_end_pos = pSNPs[j].nwindow_end_pos;
						if(pSNPs[i].nwindow_start_pos>pSNPs[j].nwindow_start_pos)
							pSNPs[i].nwindow_start_pos = pSNPs[j].nwindow_start_pos;

						pSNPs[j].nwindow_start_pos = -1;
						
						if(ncur_null_pos<num_of_SNPs)
							pSNPs[ncur_null_pos++].pbitmap = pSNPs[j].pbitmap;
						else
							DelBitmap(pSNPs[j].pbitmap);
						pSNPs[j].pbitmap = NULL;

						pSNPs[j].nSNPid_merged_with = i;
						pSNPs[i].pmerged_ids[num_of_merged++] = j;

						gnremoved_SNP_num++;
					}
				}
			}

			pSNPs[i].num_of_merged_ids = num_of_merged;
			gnmerged_id_buf_pos += num_of_merged;

			fprintf(fpout, "%d ", num_of_merged);
			for(j=0;j<num_of_merged;j++)
				fprintf(fpout, "%d ", pSNPs[i].pmerged_ids[j]);
			fprintf(fpout, "\n");

			if(gnmax_merged_ids<num_of_merged)
				gnmax_merged_ids = num_of_merged;

			if(ncur_null_pos<num_of_SNPs)
				pSNPs[ncur_null_pos++].pbitmap = pSNPs[i].pbitmap;
			else
				DelBitmap(pSNPs[i].pbitmap);
			pSNPs[i].pbitmap = NULL;
		}
		else
			gnrep_SNP_num--;
	}
	fclose(fpout);
	fclose(fpmatrix);

	if(gnremoved_SNP_num+gnrep_SNP_num!=num_of_SNPs)
		printf("Error: inconsistent number of SNPs\n");

	pRep_SNPs = new REP_SNP[gnrep_SNP_num];
	IncMemSize(sizeof(REP_SNP)*gnrep_SNP_num);

	gnrep_SNP_num = 0;
	for(i=0;i<num_of_SNPs;i++)
	{
		if(pSNPs[i].pmerged_ids!=NULL)
			pRep_SNPs[gnrep_SNP_num++].nid = i;
	}

	for(i=0;i<gnrep_SNP_num;i++)
	{
		j = i-1;
		while(j>=0 && pRep_SNPs[j].nid>=pSNPs[pRep_SNPs[i].nid].nwindow_start_pos)
			j--;
		pRep_SNPs[i].nstart_pos = j+1;
		num_of_cand_exts = i-1-j;
		num_of_cand_RHS_exts = num_of_cand_exts;

		j = i+1;
		while(j<gnrep_SNP_num && pRep_SNPs[j].nid<pSNPs[pRep_SNPs[i].nid].nwindow_end_pos)
			j++;
		pRep_SNPs[i].nend_pos = j;
		num_of_cand_RHS_exts += j-1-i;

		pRep_SNPs[i].pLHS_node_map = NULL;
		pRep_SNPs[i].pRHS_node_head = NULL;
		pRep_SNPs[i].nbeing_covered_times = 0;

		if(gnmerge_max_cand_exts<num_of_cand_exts)
			gnmerge_max_cand_exts = num_of_cand_exts;
		if(gnmerge_max_cand_RHS_exts<num_of_cand_RHS_exts)
			gnmerge_max_cand_RHS_exts = num_of_cand_RHS_exts;

		gdmerge_avg_cand_exts += num_of_cand_exts;
		gdmerge_avg_cand_RHS_exts += num_of_cand_RHS_exts;
	}

	gdmerge_avg_cand_exts /= gnrep_SNP_num;
	gdmerge_avg_cand_RHS_exts /= gnrep_SNP_num;

	nleft_boundary = pRep_SNPs[gnrep_SNP_num-1].nstart_pos;
	for(i=gnrep_SNP_num-1;i>=0;i--)
	{
		pRep_SNPs[i].nleft_boundary = pRep_SNPs[i].nstart_pos;
		if(pRep_SNPs[i].nleft_boundary>nleft_boundary)
			pRep_SNPs[i].nleft_boundary = nleft_boundary;
		else
			nleft_boundary = pRep_SNPs[i].nleft_boundary;
	}

	printf("#Representative SNPs: %d, #Removed SNPs: %d\n", gnrep_SNP_num, gnremoved_SNP_num);
	//printf("#candidate extensions after merging: avg=%.2f, max=%d\n", gdmerge_avg_cand_exts, gnmerge_max_cand_exts);
	//printf("#candidate RHS extensions after merging: avg=%.2f, max=%d\n", gdmerge_avg_cand_RHS_exts, gnmerge_max_cand_RHS_exts);
	printf("\n");

	return gnrep_SNP_num;
}

void LoadBitmaps(FILE* fp, int &nline_no, SNP *pSNPs, REP_SNP *pRep_SNPs, int ncur_pos, int &ncur_load_pos, int &ncur_null_pos)
{
	char ch;
	int nid, nhapmap_no, n1;

	while(ncur_load_pos<pRep_SNPs[ncur_pos].nend_pos)
	{
		nid = pRep_SNPs[ncur_load_pos].nid;
		if(pSNPs[nid].pbitmap==NULL)
		{
			pSNPs[nid].pbitmap = NewBitmap();
			pSNPs[nid].p1sid_list = New1SidList();
			if(ncur_null_pos!=ncur_load_pos)
				printf("Error: the NULL position should be the same as the loading position\n");
			ncur_null_pos = ncur_load_pos+1;
		}

		while(!feof(fp) && nline_no<pSNPs[nid].norder)
		{
			ch = fgetc(fp);
			while(!feof(fp) && ch!='\n')
				ch = fgetc(fp);
			nline_no++;
		}

		n1 = 0;
		nhapmap_no = 0;
		ch = fgetc(fp);
		while(!feof(fp) && ch!='\n')
		{
			while(!feof(fp) && (ch<'0' || ch>'9') && ch!='\n')
				ch = fgetc(fp);
			if(ch>='0' && ch<='9')
			{
				pSNPs[nid].pbitmap[nhapmap_no] = ch-'0';
				if(ch=='1')
					pSNPs[nid].p1sid_list[n1++] = nhapmap_no;
				else if(ch!='0')
					printf("Error: wrong allele\n");
				nhapmap_no++;
				ch = fgetc(fp);
			}
		}
		nline_no++;
		ncur_load_pos++;

		if(nhapmap_no!=gnhapmap_num)
			printf("Error: inconsistent number of hapmaps at %d-th hapmap\n", nid);
		if(pSNPs[nid].n1!=n1)
			printf("Error: inconsistent number of 1 for SNP %d\n", nid);
	}
}


