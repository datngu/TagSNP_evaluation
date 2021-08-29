//#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>
#include <iostream>

#include "global.h"

using namespace std;

char gszmaf_filename[200];
char gszmatrix_filename[200];
float gdmin_maf;
int gnwindow_len;
float gpmin_r2[3];
int gnmax_len;
int gnmax_merge_window_len;
int gnmax_covered_times;
int gnmem_size;
double gdmax_tagSNP_num;
int gnmin_bin_size;
char gszoutput_name[200];
int gnmodel;

bool gbno_equv_tag_SNPs;

int gnum_of_SNPs;

double gdavg_cand_exts;
double gdavg_cand_RHS_exts;
int gnmax_cand_exts;
int gnmax_cand_RHS_exts;
double gdmerge_avg_cand_exts;
double gdmerge_avg_cand_RHS_exts;
int gnmerge_max_cand_exts;
int gnmerge_max_cand_RHS_exts;

int gnsingleton_rules;
int gnmultiSNP_rules;
int gnmultiSNP_LHSs;

int gnused_mem_size;
int gnmax_used_mem_size;
int gnmax_genrule_used_mem_size;

double gdtotal_run_time;
double gdgen_rule_time;
double gdpre_process_time;

char* maf_path;
char* matrix_path;
char* out_tag;

void GetParameters(char* szpara_filename);
void PrintSum(char* szdataset_name);

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		printf("FastTagger in_maf in_matrix out_tag para_file\n");
		return 0;
	}

	maf_path = argv[1];
	matrix_path = argv[2];
	out_tag = argv[3];
	GetParameters(argv[4]);
	// these arguments had been moved out to command line!
	/*
	if(gszmaf_filename[0]==0)
	{
		printf("Please specify the .maf file \n");
		return 0;
	}
	if(gszmatrix_filename[0]==0)
	{
		printf("Please specify the .matrix file \n");
		return 0;
	}
	if(gszoutput_name[0]==0)
	{
		printf("Please specify the output name\n");
		return 0;
	}
	*/

	FastTagger(maf_path, matrix_path, out_tag);

	PrintSum(gszmatrix_filename);

	//Verify(gszmaf_filename, gszmatrix_filename, gszoutput_name);

//	_CrtDumpMemoryLeaks();

	return 0;
}


void GetParameters(char* szpara_filename)
{
	FILE *fp;
	char ch, szpara[100], szpara_value[200];
	int nname_len, nvalue_len, i;
	float dmin_r2;

	gszmaf_filename[0] = 0;
	gszmatrix_filename[0] = 0;
	gdmin_maf = (float)0.05;
	gnwindow_len = 100000;
	gpmin_r2[0] = (float)0.8;
	gpmin_r2[1] = (float)0.9;
	gpmin_r2[2] = (float)0.95;
	gnmax_len = 3;
	gnmax_merge_window_len = gnwindow_len;
	gnmin_bin_size = 0;
	gnmax_covered_times = 0;
	gnmem_size = 0;
	gdmax_tagSNP_num = 0;
	gnmodel = 0;

	fp = fopen(szpara_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szpara_filename);
		return;
	}

	ch = fgetc(fp);
	while(!feof(fp))
	{
		nname_len = 0;
		while(!feof(fp) && ch!='=' && ch!='\n')
		{
			szpara[nname_len++] = ch;
			ch = fgetc(fp);
		}
		szpara[nname_len] = 0;
		if(ch=='=')
			ch = fgetc(fp);
		nvalue_len = 0;
		while(!feof(fp) && ch!='\n')
		{
			szpara_value[nvalue_len++] = ch;
			ch = fgetc(fp);
		}
		szpara_value[nvalue_len] = 0;

		if(nname_len>0 && nvalue_len>0)
		{
			if(strcmp(szpara, "data_maf")==0)
				strcpy(gszmaf_filename, szpara_value);
			else if(strcmp(szpara, "data_matrix")==0)
				strcpy(gszmatrix_filename, szpara_value);
			else if(strcmp(szpara, "output")==0)
				strcpy(gszoutput_name, szpara_value);
			else if(strcmp(szpara, "min_maf")==0)
			{
				gdmin_maf = (float)atof(szpara_value);
				if(gdmin_maf<0 || gdmin_maf>1)
				{
					printf("The value of parameter min_maf should be between 0 and 1. Use default value 0.05\n");
					gdmin_maf = (float)0.05;
				}
			}
			else if(strcmp(szpara, "window_len")==0)
			{
				gnwindow_len = atoi(szpara_value);
				if(gnwindow_len<=0)
				{
					printf("The value of parameter window_len should be larger than 0. Use default value 100k\n");
					gnwindow_len = 100000;
				}
			}
			else if(strstr(szpara, "min_r2"))
			{
				dmin_r2 = (float)atof(szpara_value);
				if(dmin_r2<=0 && dmin_r2>1)
					printf("The value of min_r2 parameter should be between 0 and 1.\n");
				else if(strcmp(szpara, "min_r2_1")==0)
					gpmin_r2[0] = dmin_r2;
				else if(strcmp(szpara, "min_r2_2")==0)
					gpmin_r2[1] = dmin_r2;
				else if(strcmp(szpara, "min_r2_3")==0)
					gpmin_r2[2] = dmin_r2;
				else 
					printf("Invalid parameter: %s\n", szpara);
			}
			else if(strcmp(szpara, "max_len")==0)
			{
				gnmax_len = atoi(szpara_value);
				if(gnmax_len<=0 && gnmax_len>3)
				{
					printf("The value of max_len parameter should be between 1 and 3. Reset to default value 3.\n");
					gnmax_len = 3;
				}				
			}
			else if(strcmp(szpara, "max_merge_window_len")==0)
			{
				gnmax_merge_window_len = atoi(szpara_value);
				if(gnmax_merge_window_len<0)
					gnmax_merge_window_len = gnwindow_len;
			}
			else if(strcmp(szpara, "min_bin_size")==0)
			{
				gnmin_bin_size = atoi(szpara_value);
				if(gnmin_bin_size<0)
					gnmin_bin_size = 5;
			}
			else if(strcmp(szpara, "max_covered_times")==0)
			{
				gnmax_covered_times = atoi(szpara_value);
				if(gnmax_covered_times<-1)
					gnmax_covered_times = 0;
			}
			else if(strcmp(szpara, "mem_size")==0)
			{
				gnmem_size = atoi(szpara_value);
				if(gnmem_size<0)
					gnmem_size = 0;
			}
			else if(strcmp(szpara, "max_tagSNP_num")==0)
			{
				gdmax_tagSNP_num = atof(szpara_value);
				if(gdmax_tagSNP_num<0)
					gdmax_tagSNP_num = 0;
			}
			else if(strcmp(szpara, "model")==0)
			{
				for(i=0;i<nvalue_len;i++)
				{
					if(szpara_value[i]>='a' && szpara_value[i]<='z')
						szpara_value[i] += 'A'-'a';
				}
				if(strcmp(szpara_value, "MULTITAG")==0)
					gnmodel = 1;
			}
			else
				printf("Invalid parameter: %s\n", szpara);
		}

		if(ch=='\n')
			ch = fgetc(fp);
	}
	fclose(fp);

	printf("Parameter settings:\n");
	printf("-------------------\n");
	printf("data_maf=%s\n", gszmaf_filename);
	printf("data_matrix=%s\n", gszmatrix_filename);
	printf("min_maf=%.2f\n", gdmin_maf);
	printf("window_len=%d\n", gnwindow_len);
	printf("max_len=%d\n", gnmax_len);
	for(int i=0;i<gnmax_len;i++)
		printf("  len=%d, min_r2=%.2f\n", i+1, gpmin_r2[i]);
	printf("nmax_merge_window_len=%d\n", gnmax_merge_window_len);
	printf("min_bin_size=%d\n", gnmin_bin_size);
	printf("nmax_covered_times=%d\n", gnmax_covered_times);
	printf("mem_size=%d MB\n", gnmem_size);
	printf("max_tagSNP_num=%.2f\n", gdmax_tagSNP_num);
	if(gnmodel==MMTagger_Model)
		printf("model=MMTagger\n");
	else if(gnmodel==MultiTag_Model)
		printf("model=MultiTag\n");
	else
		printf("model=unknow model\n");
	printf("\n");
}

void FastTagger(char* szmaf_filename, char* szmatrix_filename, char *szoutput_name)
{
	//cout << "loading SNP" << endl;
	SNP *pSNPs;
	char szmerged_ids_filename[200];
	int num_of_SNPs, i, num_of_Rep_SNPs;
	REP_SNP *pRep_SNPs;
	struct timeb start, end;

	ftime(&start);

	gnused_mem_size = 0;
	gnmax_used_mem_size = 0;

	gprule_num = new int[gnmax_len+1];
	IncMemSize(sizeof(int)*(gnmax_len+1));
	memset(gprule_num, 0, sizeof(int)*(gnmax_len+1));
	gpLHS_num = new int[gnmax_len+1];
	IncMemSize(sizeof(int)*(gnmax_len+1));
	memset(gpLHS_num, 0, sizeof(int)*(gnmax_len+1));
	gpmerged_bin_rule_num = new int[gnmax_len+1];
	IncMemSize(sizeof(int)*(gnmax_len+1));
	memset(gpmerged_bin_rule_num, 0, sizeof(int)*(gnmax_len+1));

	gnsingleton_rules = 0;
	gnmultiSNP_rules = 0;
	gnmultiSNP_LHSs = 0;

	//cout << "loading SNP ... check point 2" << endl;
	num_of_SNPs = LoadSNPs(szmaf_filename, gdmin_maf, pSNPs, szoutput_name);
	//cout << szmaf_filename << endl;
	if(num_of_SNPs==0)
		return;
	gnum_of_SNPs = num_of_SNPs;

	gporig_window_start_pos = new int[num_of_SNPs];
	IncMemSize(sizeof(int)*num_of_SNPs);
	gporig_window_end_pos = new int[num_of_SNPs];
	IncMemSize(sizeof(int)*num_of_SNPs);
	CalcWindowSize(pSNPs, num_of_SNPs, gnwindow_len);

	gnhapmap_num = GetHapMapNum(szmatrix_filename);

	gpmerged_SNP_ids = new int[num_of_SNPs];
	IncMemSize(sizeof(int)*num_of_SNPs);
	sprintf(szmerged_ids_filename, "%s.merged-ids.txt", szoutput_name);
	pRep_SNPs = NULL;
	num_of_Rep_SNPs = MergeSNPs(szmatrix_filename, pSNPs, num_of_SNPs, pRep_SNPs, szmerged_ids_filename);

	ftime(&end);
	gdpre_process_time = end.time-start.time+(double)(end.millitm-start.millitm)/1000;
	printf("Merging SNP time: %.3f\n", gdpre_process_time);

	gntotal_rule_size = 0;

	if(num_of_Rep_SNPs>0)
	{
		GenRules(szmatrix_filename, pSNPs, pRep_SNPs, num_of_Rep_SNPs, szoutput_name);
		//EndShift(pSNPs, pRep_SNPs, num_of_Rep_SNPs, 0);

		delete []pRep_SNPs;
		DecMemSize(sizeof(REP_SNP)*num_of_Rep_SNPs);

		ftime(&end);
		gdgen_rule_time = end.time-start.time+(double)(end.millitm-start.millitm)/1000;
		gnmax_genrule_used_mem_size = gnmax_used_mem_size;

		gnsingleton_rules = gprule_num[1];
		gnmultiSNP_rules = 0;
		gnmultiSNP_LHSs = 0;
		gnmerged_bin_rule_num = 0;
		for(i=2;i<=gnmax_len;i++)
		{
			printf("len=%d\t%d\t%d\t%d\n", i, gpLHS_num[i], gprule_num[i], gpmerged_bin_rule_num[i]);
			gnmultiSNP_LHSs += gpLHS_num[i];
			gnmultiSNP_rules += gprule_num[i];
			gnmerged_bin_rule_num += gpmerged_bin_rule_num[i];
		}

		printf("#Singleton rules: %d\n", gnsingleton_rules);
		printf("#multi-SNP LHSs: %d\n", gnmultiSNP_LHSs);
		printf("#multi-SNP rules: %d\n", gnmultiSNP_rules);
		printf("#multi-SNP rules with merged bins: %d\n", gnmerged_bin_rule_num);
		printf("\n");

		if(gdsubset_check_times>0)
			gdsubset_check_depth_sum /= gdsubset_check_times;
		else 
			gdsubset_check_depth_sum = 0;
		printf("subset check times: %.f, average check depth: %.2f\n", gdsubset_check_times, gdsubset_check_depth_sum);
		printf("rule generation time: %.2f\n", gdgen_rule_time);
		printf("maximal memory consumption: %.2f MB\n\n", (double)gnmax_genrule_used_mem_size/(1<<20));

		//--------------------------------------------------------------------------

		gnmax_used_mem_size = gnused_mem_size;
		SelectTagSNPs(pSNPs, num_of_SNPs, gnmem_size, szoutput_name);

	}

	delete []gporig_window_start_pos;
	DecMemSize(sizeof(int)*num_of_SNPs);
	delete []gporig_window_end_pos;
	DecMemSize(sizeof(int)*num_of_SNPs);
	delete []gpmerged_SNP_ids;
	DecMemSize(sizeof(int)*num_of_SNPs);

	delete []pSNPs;
	DecMemSize(num_of_SNPs*sizeof(SNP));

	delete []gpmerged_bin_rule_num;
	DecMemSize(sizeof(int)*(gnmax_len+1));
	delete []gprule_num;
	DecMemSize(sizeof(int)*(gnmax_len+1));
	delete []gpLHS_num;
	DecMemSize(sizeof(int)*(gnmax_len+1));

	ftime(&end);
	gdtotal_run_time = end.time-start.time+(double)(end.millitm-start.millitm)/1000;

	printf("#recursive calls: %d\n", gntotal_calls);
	printf("total running time: %.2f\n", gdtotal_run_time);
	printf("maximal memory consumption: %.2f MB\n\n", (double)gnmax_used_mem_size/(1<<20));

	if(gnused_mem_size!=0)
		printf("Error: %d bytes are not released\n", gnused_mem_size);

}

void PrintSum(char* szdataset_name)
{
	FILE *fpsum;

	fpsum = fopen("FastTagger.sum", "a+");
	if(fpsum==NULL)
	{
		printf("Error: cannot open file FastTagger.sum for appending\n");
		return;
	}

	fprintf(fpsum, "%s ", szdataset_name);
	if(gnmodel==0)
		fprintf(fpsum, "MMTagger ");
	else
		fprintf(fpsum, "MultiTag ");
	fprintf(fpsum, "%.2f %d %.2f %.2f %.2f %d\t", gdmin_maf, gnwindow_len, gpmin_r2[0], gpmin_r2[1], gpmin_r2[2], gnmax_len);
	fprintf(fpsum, "%d %d %d %dMB %d %.2f\t", gnmax_merge_window_len, gnmin_bin_size, gnmax_covered_times, gnmem_size, gbno_equv_tag_SNPs, gdmax_tagSNP_num);
	fprintf(fpsum, "%d %d %d %d %d %d\t", gnum_of_SNPs, gnrep_SNP_num, gnum_of_tag_SNPs, gncovered_SNPs, gnnot_covered_tag_SNPs, gnnonb_SNP_num);
	fprintf(fpsum, "%.3f %.3f %d\t", gdgen_rule_time, gdtotal_run_time, gntotal_calls);
	fprintf(fpsum, "%.2fMB %.2fMB %d\t", (double)gnmax_genrule_used_mem_size/(1<<20), (double)gnmax_used_mem_size/(1<<20), gnum_of_partitions);
	fprintf(fpsum, "%d %d %d %d\t", gnsingleton_rules, gnmultiSNP_LHSs, gnmultiSNP_rules, gnmerged_bin_rule_num);
	fprintf(fpsum, "%.f %.2f\t", gdsubset_check_times, gdsubset_check_depth_sum);
	fprintf(fpsum, "%.2f %d %.2f %d\t", gdavg_cand_exts, gnmax_cand_exts, gdavg_cand_RHS_exts, gnmax_cand_RHS_exts);
	fprintf(fpsum, "%.2f %d %.2f %d\t", gdmerge_avg_cand_exts, gnmerge_max_cand_exts, gdmerge_avg_cand_RHS_exts, gnmerge_max_cand_RHS_exts);
	fprintf(fpsum, "\n");

	fclose(fpsum);
}


