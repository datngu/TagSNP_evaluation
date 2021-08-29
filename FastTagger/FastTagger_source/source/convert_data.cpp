#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "global.h"


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

void LoadMatrix(char* szphased_hapmap_filename, int num_of_SNPs, int num_of_hapmaps, bool *pmatrix)
{
	FILE *fp;
	char ch;
	int nrow_no, nclmn_no;

	fp = fopen(szphased_hapmap_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szphased_hapmap_filename);
		return;
	}

	nrow_no = 0;

	ch = fgetc(fp);
	while(!feof(fp))
	{
		nclmn_no = 0;
		while(!feof(fp) && ch!='\n' && ch!='\r')
		{
			if(ch=='0')
				pmatrix[nclmn_no*num_of_hapmaps+nrow_no] = 0;
			else if(ch=='1')
				pmatrix[nclmn_no*num_of_hapmaps+nrow_no] = 1;
			else
				printf("Error: '0' or '1' expected\n");

			nclmn_no++;

			ch = fgetc(fp);
			while(!feof(fp) && (ch==' ' || ch=='\t'))
				ch = fgetc(fp);
		}
		if(ch=='\r')
			ch = fgetc(fp);
		if(ch=='\n')
			ch = fgetc(fp);
		if(nclmn_no>0)
		{
			if(nclmn_no!=num_of_SNPs)
				printf("Error: inconsistent number of SNPs\n");
			nclmn_no = 0;
			nrow_no++;
		}
	}
	fclose(fp);

	if(nrow_no!=num_of_hapmaps)
		printf("Error: inconsistent number of samples\n");

}

void OutputMatrix(bool *pmatrix, int num_of_SNPs, int num_of_hapmaps, char* szoutput_filename)
{
	FILE *fp;
	int i, j;

	fp = fopen(szoutput_filename, "wt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for write\n", szoutput_filename);
		return;
	}

	for(i=0;i<num_of_SNPs;i++)
	{
		for(j=0;j<num_of_hapmaps-1;j++)
			fprintf(fp, "%d ", pmatrix[i*num_of_hapmaps+j]);
		fprintf(fp, "%d\n", pmatrix[i*num_of_hapmaps+num_of_hapmaps-1]);
	}

	fclose(fp);
}

void OutputMatrixData(char* szlegend_filename, bool *pmatrix, int num_of_SNPs, int num_of_hapmaps, double dmin_maf, int nvalid_SNP_num, char* szoutput_name)
{
	FILE *fp, *fpmaf, *fpmatrix, *fpout;
	char ch, allele1, allele2, szSNPid[100], szposition[100], szmaf_filename[200], szmatrix_filename[200];
	int nlen, nSNP_no, pfreqs[2], i, num_of_clmns, nselected_SNPs;
	double dmin_freq, davg_maf;


	fp = fopen(szlegend_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szlegend_filename);
		return;
	}

	if(nvalid_SNP_num==-1)
	{
		sprintf(szmaf_filename, "%s.maf.txt", szoutput_name);
		fpmaf = fopen(szmaf_filename, "wt");
		if(fpmaf==NULL)
		{
			printf("Error: cannot open file %s for write\n", szmaf_filename);
			return;
		}
		sprintf(szmatrix_filename, "%s.haps", szoutput_name);
		fpmatrix = fopen(szmatrix_filename, "wt");
		if(fpmatrix==NULL)
		{
			printf("Error: cannot open file %s for write\n", szmatrix_filename);
			return;
		}
		fpout = NULL;
		for(i=0;i<num_of_hapmaps/2-1;i++)
			fprintf(fpmatrix, "Sample_%d\t", i+1);
		fprintf(fpmatrix, "Sample_%d\n", num_of_hapmaps/2);
		for(i=0;i<num_of_hapmaps/2-1;i++)
			fprintf(fpmatrix, "T\tU\t");
		fprintf(fpmatrix, "T\tU\n");
	}
	else if(nvalid_SNP_num==0)
	{
		sprintf(szmaf_filename, "%s.maf", szoutput_name);
		fpmaf = fopen(szmaf_filename, "wt");
		if(fpmaf==NULL)
		{
			printf("Error: cannot open file %s for write\n", szmaf_filename);
			return;
		}
		sprintf(szmatrix_filename, "%s.matrix", szoutput_name);
		fpmatrix = fopen(szmatrix_filename, "wt");
		if(fpmatrix==NULL)
		{
			printf("Error: cannot open file %s for write\n", szmatrix_filename);
			return;
		}
		fpout = NULL;
	}
	else
		printf("Error: wrong valid SNP num\n");

	davg_maf = 0;
	nselected_SNPs = 0;

	nSNP_no = 0;
	num_of_clmns = 0;
	ch = fgetc(fp);
	while(!feof(fp) && ch!='\n')
	{
		nlen = 0;
		while(!feof(fp) && ch!='\t' && ch!='\n')
		{
			ch = fgetc(fp);
			nlen++;
		}
		if(nlen>0)
			num_of_clmns++;
		if(ch=='\t')
			ch = fgetc(fp);
	}
	ch = fgetc(fp);

	while(!feof(fp))
	{
		if(num_of_clmns==4)
		{
			nlen = 0;
			while(!feof(fp) && ch!='\t' && ch!='\n')
			{
				szSNPid[nlen++] = ch;
				ch = fgetc(fp);
			}
			szSNPid[nlen] = 0;

			if(ch!='\t')
				printf("Error: tab expected\n");
			else
				ch = fgetc(fp);
		}
		else 
			sprintf(szSNPid, "SNP%d", nSNP_no);

		nlen = 0;
		while(!feof(fp) && ch!='\t' && ch!='\n')
		{
			szposition[nlen++] = ch;
			ch = fgetc(fp);
		}
		szposition[nlen] = 0;

		if(ch!='\t')
			printf("Error: tab expected\n");
		else
			ch = fgetc(fp);

		if(ch!='A' && ch!='C' && ch!='G' && ch!='T')
			printf("Error: 'A' or 'C' or 'G' or 'T' expected\n");
		allele1 = ch;

		if(ch!='\n')
			ch = fgetc(fp);

		if(ch!='\t')
			printf("Error: tab expected\n");
		else
			ch = fgetc(fp);

		if(ch!='A' && ch!='C' && ch!='G' && ch!='T')
			printf("Error: 'A' or 'C' or 'G' or 'T' expected\n");
		allele2 = ch;

		ch = fgetc(fp);
		if(ch=='\r')
			ch = fgetc(fp);
		if(ch!='\n')
			printf("Error: new line expected\n");
		else 
			ch = fgetc(fp);

		pfreqs[0] = 0;
		pfreqs[1] = 0;
		for(i=0;i<num_of_hapmaps;i++)
			pfreqs[pmatrix[nSNP_no*num_of_hapmaps+i]]++;

		if(pfreqs[0]<pfreqs[1])
			dmin_freq = (double)pfreqs[0]/num_of_hapmaps;
		else
			dmin_freq = (double)pfreqs[1]/num_of_hapmaps;

		if(dmin_freq>=dmin_maf)
		{
			if(nvalid_SNP_num==-1)
			{
				fprintf(fpmaf, "%s\t%s\t%.6f\n", szSNPid, szposition, dmin_freq);
				for(i=0;i<num_of_hapmaps-1;i++)
					fprintf(fpmatrix, "%d\t", pmatrix[nSNP_no*num_of_hapmaps+i]+1);
				fprintf(fpmatrix, "%d\n", pmatrix[nSNP_no*num_of_hapmaps+num_of_hapmaps-1]+1);
			}
			else if(nvalid_SNP_num==0) //0: major allele, 1: minor allele
			{
				if(pfreqs[0]>pfreqs[1] || pfreqs[0]==pfreqs[1] && pmatrix[nSNP_no*num_of_hapmaps]==0)
				{
					fprintf(fpmaf, "%s\t%s\t%c\t%c\t%.6f\n", szSNPid, szposition, allele1, allele2, dmin_freq);
					for(i=0;i<num_of_hapmaps-1;i++)
						fprintf(fpmatrix, "%d ", pmatrix[nSNP_no*num_of_hapmaps+i]);
					fprintf(fpmatrix, "%d\n", pmatrix[nSNP_no*num_of_hapmaps+num_of_hapmaps-1]);
				}
				else
				{
					fprintf(fpmaf, "%s\t%s\t%c\t%c\t%.6f\n", szSNPid, szposition, allele2, allele1, dmin_freq);
					for(i=0;i<num_of_hapmaps-1;i++)
						fprintf(fpmatrix, "%d ", 1-pmatrix[nSNP_no*num_of_hapmaps+i]);
					fprintf(fpmatrix, "%d\n", 1-pmatrix[nSNP_no*num_of_hapmaps+num_of_hapmaps-1]);
				}
			}
			davg_maf += dmin_freq;
			nselected_SNPs++;
		}

		nSNP_no++;
	}
	fclose(fp);
	fclose(fpmaf);
	fclose(fpmatrix);

	if(nSNP_no!=num_of_SNPs)
		printf("Error: inconsistent number of SNPs\n");
	davg_maf /= nselected_SNPs;
	printf("Average MAF: %.3f\n", davg_maf);

}

void OutputMMTaggerMatrixData(char* szlegend_filename, bool *pmatrix, int num_of_SNPs, int num_of_hapmaps, double dmin_maf, int nvalid_SNP_num, char* szoutput_name)
{
	FILE *fp, *fpout;
	char ch, allele1, allele2, szSNPid[100], szposition[100];
	int nlen, nSNP_no, pfreqs[2], i, num_of_clmns;
	double dmin_freq;


	fp = fopen(szlegend_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szlegend_filename);
		return;
	}

	fpout = fopen(szoutput_name, "wt");
	if(fpout==NULL)
	{
		printf("Error: cannot open file %s for write\n", szoutput_name);
		return;
	}
	fprintf(fpout, "%d %d\n", num_of_hapmaps, nvalid_SNP_num);


	nSNP_no = 0;
	num_of_clmns = 0;
	ch = fgetc(fp);
	while(!feof(fp) && ch!='\n')
	{
		nlen = 0;
		while(!feof(fp) && ch!='\t' && ch!='\n')
		{
			ch = fgetc(fp);
			nlen++;
		}
		if(nlen>0)
			num_of_clmns++;
		if(ch=='\t')
			ch = fgetc(fp);
	}
	ch = fgetc(fp);

	while(!feof(fp))
	{
		if(num_of_clmns==4)
		{
			nlen = 0;
			while(!feof(fp) && ch!='\t' && ch!='\n')
			{
				szSNPid[nlen++] = ch;
				ch = fgetc(fp);
			}
			szSNPid[nlen] = 0;

			if(ch!='\t')
				printf("Error: tab expected\n");
			else
				ch = fgetc(fp);
		}
		else 
			sprintf(szSNPid, "SNP%d", nSNP_no);

		nlen = 0;
		while(!feof(fp) && ch!='\t' && ch!='\n')
		{
			szposition[nlen++] = ch;
			ch = fgetc(fp);
		}
		szposition[nlen] = 0;

		if(ch!='\t')
			printf("Error: tab expected\n");
		else
			ch = fgetc(fp);

		if(ch!='A' && ch!='C' && ch!='G' && ch!='T')
			printf("Error: 'A' or 'C' or 'G' or 'T' expected\n");
		allele1 = ch;

		if(ch!='\n')
			ch = fgetc(fp);

		if(ch!='\t')
			printf("Error: tab expected\n");
		else
			ch = fgetc(fp);

		if(ch!='A' && ch!='C' && ch!='G' && ch!='T')
			printf("Error: 'A' or 'C' or 'G' or 'T' expected\n");
		allele2 = ch;

		ch = fgetc(fp);
		if(ch=='\r')
			ch = fgetc(fp);
		if(ch!='\n')
			printf("Error: new line expected\n");
		else 
			ch = fgetc(fp);

		pfreqs[0] = 0;
		pfreqs[1] = 0;
		for(i=0;i<num_of_hapmaps;i++)
			pfreqs[pmatrix[nSNP_no*num_of_hapmaps+i]]++;

		if(pfreqs[0]<pfreqs[1])
			dmin_freq = (double)pfreqs[0]/num_of_hapmaps;
		else
			dmin_freq = (double)pfreqs[1]/num_of_hapmaps;

		if(dmin_freq>=dmin_maf)
		{
			if(pfreqs[0]>=pfreqs[1])
			{
				fprintf(fpout, "%c %c %s %.3f ", allele1, allele2, szposition, dmin_freq);
				for(i=0;i<num_of_hapmaps-1;i++)
					fprintf(fpout, "%d ", pmatrix[nSNP_no*num_of_hapmaps+i]);
				fprintf(fpout, "%d\n", pmatrix[nSNP_no*num_of_hapmaps+num_of_hapmaps-1]);
			}
			else
			{
				fprintf(fpout, "%c %c %s %.3f ", allele2, allele1, szposition, dmin_freq);
				for(i=0;i<num_of_hapmaps-1;i++)
					fprintf(fpout, "%d ", 1-pmatrix[nSNP_no*num_of_hapmaps+i]);
				fprintf(fpout, "%d\n", 1-pmatrix[nSNP_no*num_of_hapmaps+num_of_hapmaps-1]);
			}
		}

		nSNP_no++;
	}
	fclose(fp);
	fclose(fpout);

	if(nSNP_no!=num_of_SNPs)
		printf("Error: inconsistent number of SNPs\n");
}

void Phased4ToMatrix(char* szphased_hapmap_filename, char *szdata_prefix, double dmin_maf, char *szoutput_name)
{
	char szphased_legend_filename[200], szphased_samples_filename[200];
	int num_of_SNPs, num_of_samples, num_of_hapmaps;
	bool *pmatrix;

	sprintf(szphased_legend_filename, "%s_legend.txt", szdata_prefix);
	sprintf(szphased_samples_filename, "%s_sample.txt", szdata_prefix);

	num_of_SNPs = GetRowNum(szphased_legend_filename)-1;
	num_of_samples = GetRowNum(szphased_samples_filename);
	num_of_hapmaps = num_of_samples/3*4;

	pmatrix = new bool[num_of_SNPs*num_of_hapmaps];
	LoadMatrix(szphased_hapmap_filename, num_of_SNPs, num_of_hapmaps, pmatrix);

	OutputMatrixData(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf, 0, szoutput_name);

	delete []pmatrix;
}



void Phased2ToMatrix(char* szphased_hapmap_filename, char *szdata_prefix, double dmin_maf, char *szoutput_name)
{
	char szphased_legend_filename[200], szphased_samples_filename[200];

	int num_of_SNPs, num_of_samples, num_of_hapmaps;
	bool *pmatrix;

	sprintf(szphased_legend_filename, "%s_legend.txt", szdata_prefix);
	sprintf(szphased_samples_filename, "%s_sample.txt", szdata_prefix);

	num_of_SNPs = GetRowNum(szphased_legend_filename)-1;
	num_of_samples = GetRowNum(szphased_samples_filename);
	num_of_hapmaps = num_of_samples*2;

	pmatrix = new bool[num_of_SNPs*num_of_hapmaps];
	LoadMatrix(szphased_hapmap_filename, num_of_SNPs, num_of_hapmaps, pmatrix);

	OutputMatrixData(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf, 0, szoutput_name);

	delete []pmatrix;
}


//----------------------------------------------------------------------------------------

int GetValidSNPNum(char* szlegend_filename, bool *pmatrix, int num_of_SNPs, int num_of_hapmaps, double dmin_maf)
{
	FILE *fp;
	char ch, allele1, allele2, szSNPid[100], szposition[100];
	int nlen, nSNP_no, pfreqs[2], i, num_of_clmns, num_of_valid_SNPs;
	double dmin_freq;


	fp = fopen(szlegend_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szlegend_filename);
		return 0;
	}

	num_of_valid_SNPs = 0;

	nSNP_no = 0;
	num_of_clmns = 0;
	ch = fgetc(fp);
	while(!feof(fp) && ch!='\n')
	{
		nlen = 0;
		while(!feof(fp) && ch!='\t' && ch!='\n')
		{
			ch = fgetc(fp);
			nlen++;
		}
		if(nlen>0)
			num_of_clmns++;
		if(ch=='\t')
			ch = fgetc(fp);
	}
	ch = fgetc(fp);

	while(!feof(fp))
	{
		if(num_of_clmns==4)
		{
			nlen = 0;
			while(!feof(fp) && ch!='\t' && ch!='\n')
			{
				szSNPid[nlen++] = ch;
				ch = fgetc(fp);
			}
			szSNPid[nlen] = 0;

			if(ch!='\t')
				printf("Error: tab expected\n");
			else
				ch = fgetc(fp);
		}
		else 
			sprintf(szSNPid, "SNP%d", nSNP_no);

		nlen = 0;
		while(!feof(fp) && ch!='\t' && ch!='\n')
		{
			szposition[nlen++] = ch;
			ch = fgetc(fp);
		}
		szposition[nlen] = 0;

		if(ch!='\t')
			printf("Error: tab expected\n");
		else
			ch = fgetc(fp);

		if(ch!='A' && ch!='C' && ch!='G' && ch!='T')
			printf("Error: 'A' or 'C' or 'G' or 'T' expected\n");
		allele1 = ch;

		if(ch!='\n')
			ch = fgetc(fp);

		if(ch!='\t')
			printf("Error: tab expected\n");
		else
			ch = fgetc(fp);

		if(ch!='A' && ch!='C' && ch!='G' && ch!='T')
			printf("Error: 'A' or 'C' or 'G' or 'T' expected\n");
		allele2 = ch;

		ch = fgetc(fp);
		if(ch=='\r')
			ch = fgetc(fp);
		if(ch!='\n')
			printf("Error: new line expected\n");
		else 
			ch = fgetc(fp);

		pfreqs[0] = 0;
		pfreqs[1] = 0;
		for(i=0;i<num_of_hapmaps;i++)
			pfreqs[pmatrix[nSNP_no*num_of_hapmaps+i]]++;

		if(pfreqs[0]<pfreqs[1])
			dmin_freq = (double)pfreqs[0]/num_of_hapmaps;
		else
			dmin_freq = (double)pfreqs[1]/num_of_hapmaps;

		if(dmin_freq>=dmin_maf)
			num_of_valid_SNPs++;

		nSNP_no++;
	}
	fclose(fp);

	if(nSNP_no!=num_of_SNPs)
		printf("Error: inconsistent number of SNPs\n");

	return num_of_valid_SNPs;

}

void Phased4ToMMTagger(char* szphased_hapmap_filename, char *szdata_prefix, double dmin_maf, char *szoutput_name)
{
	char szphased_legend_filename[200], szphased_samples_filename[200];
	int num_of_SNPs, num_of_samples, num_of_hapmaps, num_of_valid_SNPs;
	bool *pmatrix;

	sprintf(szphased_legend_filename, "%s_legend.txt", szdata_prefix);
	sprintf(szphased_samples_filename, "%s_sample.txt", szdata_prefix);

	num_of_SNPs = GetRowNum(szphased_legend_filename)-1;
	num_of_samples = GetRowNum(szphased_samples_filename);
	num_of_hapmaps = num_of_samples/3*4;

	pmatrix = new bool[num_of_SNPs*num_of_hapmaps];
	LoadMatrix(szphased_hapmap_filename, num_of_SNPs, num_of_hapmaps, pmatrix);

	num_of_valid_SNPs = GetValidSNPNum(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf);

	OutputMMTaggerMatrixData(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf, num_of_valid_SNPs, szoutput_name);

	delete []pmatrix;
}

void Phased2ToMMTagger(char* szphased_hapmap_filename, char *szdata_prefix, double dmin_maf, char *szoutput_name)
{
	char szphased_legend_filename[200], szphased_samples_filename[200];
	int num_of_SNPs, num_of_samples, num_of_hapmaps, num_of_valid_SNPs;
	bool *pmatrix;

	sprintf(szphased_legend_filename, "%s_legend.txt", szdata_prefix);
	sprintf(szphased_samples_filename, "%s_sample.txt", szdata_prefix);

	num_of_SNPs = GetRowNum(szphased_legend_filename)-1;
	num_of_samples = GetRowNum(szphased_samples_filename);
	num_of_hapmaps = num_of_samples*2;

	pmatrix = new bool[num_of_SNPs*num_of_hapmaps];
	LoadMatrix(szphased_hapmap_filename, num_of_SNPs, num_of_hapmaps, pmatrix);

	num_of_valid_SNPs = GetValidSNPNum(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf);

	OutputMMTaggerMatrixData(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf, num_of_valid_SNPs, szoutput_name);

	delete []pmatrix;
}

//-------------------------------------------------------------------------------------------

void Phased4ToMultiTag(char* szphased_hapmap_filename, char *szdata_prefix, double dmin_maf, char *szoutput_name)
{
	char szphased_legend_filename[200], szphased_samples_filename[200];
	int num_of_SNPs, num_of_samples, num_of_hapmaps;
	bool *pmatrix;

	sprintf(szphased_legend_filename, "%s_legend.txt", szdata_prefix);
	sprintf(szphased_samples_filename, "%s_sample.txt", szdata_prefix);

	num_of_SNPs = GetRowNum(szphased_legend_filename)-1;
	num_of_samples = GetRowNum(szphased_samples_filename);
	num_of_hapmaps = num_of_samples/3*4;

	pmatrix = new bool[num_of_SNPs*num_of_hapmaps];
	LoadMatrix(szphased_hapmap_filename, num_of_SNPs, num_of_hapmaps, pmatrix);

	OutputMatrixData(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf, -1, szoutput_name);

	delete []pmatrix;
}

void Phased2ToMultiTag(char* szphased_hapmap_filename, char *szdata_prefix, double dmin_maf, char *szoutput_name)
{
	char szphased_legend_filename[200], szphased_samples_filename[200];
	int num_of_SNPs, num_of_samples, num_of_hapmaps;
	bool *pmatrix;

	sprintf(szphased_legend_filename, "%s_legend.txt", szdata_prefix);
	sprintf(szphased_samples_filename, "%s_sample.txt", szdata_prefix);

	num_of_SNPs = GetRowNum(szphased_legend_filename)-1;
	num_of_samples = GetRowNum(szphased_samples_filename);
	num_of_hapmaps = num_of_samples*2;

	pmatrix = new bool[num_of_SNPs*num_of_hapmaps];
	LoadMatrix(szphased_hapmap_filename, num_of_SNPs, num_of_hapmaps, pmatrix);

	OutputMatrixData(szphased_legend_filename, pmatrix, num_of_SNPs, num_of_hapmaps, dmin_maf, -1, szoutput_name);

	delete []pmatrix;
}


//============================================================================================

void VerifyMAF(char* szhaptype_filename, int noffset, char* szmaf_filename)
{
	FILE *fp, *fpmaf;
	char szSNPid[100], szposition[100];
	double dmaf, dmin_freq;
	int pfreqs[2], num_of_hapmaps, nhapmap_no, nSNP_no, nallele_no;
	char ch;

	fp = fopen(szhaptype_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szhaptype_filename);
		return;
	}
	fpmaf = fopen(szmaf_filename, "rt");
	if(fpmaf==NULL)
	{
		printf("Error: cannot open file %s for read\n", szhaptype_filename);
		return;
	}

	nSNP_no = 0;
	num_of_hapmaps = 0;

	ch = fgetc(fp);
	while(!feof(fp))
	{
		memset(pfreqs, 0, sizeof(int)*2);
		nhapmap_no = 0;

		while(!feof(fp) && ch!='\n')
		{
			if(ch<'0' || ch>'2')
			{
				printf("Error: number 0, 1, 2 expected\n");
				nallele_no = -1;
			}
			else
			{
				nallele_no = ch - '0'-noffset;
				if(nallele_no<0 || nallele_no>1)
				{
					printf("Error: allele value should be either 0 or 1\n");
					nallele_no = -1;
				}
			}
			nhapmap_no++;

			if(nallele_no>=0)
				pfreqs[nallele_no]++;

			ch = fgetc(fp);
			while(!feof(fp) && (ch==' ' || ch=='\t'))
				ch = fgetc(fp);
		}
		if(nSNP_no==0)
			num_of_hapmaps = nhapmap_no;
		else if(nhapmap_no!=num_of_hapmaps)
			printf("Error: inconsistent number of hapmaps\n");

		if(pfreqs[0]<pfreqs[1])
			dmin_freq = (double)pfreqs[0]/num_of_hapmaps;
		else 
			dmin_freq = (double)pfreqs[1]/num_of_hapmaps;

		fscanf(fpmaf, "%s", szSNPid);
		fscanf(fpmaf, "%s", szposition);
		fscanf(fpmaf, "%lf", &dmaf);

		if(dmin_freq-dmaf>0.00000001 || dmaf-dmin_freq>0.00000001)
			printf("Error: inconsistent maf\n");

		nSNP_no++;
		if(ch=='\n')
			ch = fgetc(fp);
	}

	fclose(fpmaf);
	fclose(fp);
}

