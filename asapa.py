#!/usr/bin/python
# -*- coding: UTF-8 -*-
# This pipeline was created by JZX in September, 2022.
# To establish the connection among AS (Alternative splicing), ATI (alternative transcriptional initiation) and APA (Alternative polyadenylation), 
# polyA-removed ccs reads were mapped to reference genome thus obtained the transcription start site (TSS) and polyA site (PAS), and were
# also collapsed into transcripts and then divided by suppa2 according AS event.
# Version = 1.0.0

import subprocess
import csv
import os
import os.path
import rpy2.robjects as robjects
import timeit 
import time
import sys
from scipy.stats import ks_2samp       
from scipy.stats import chi2_contingency
from scipy.stats import spearmanr
from scipy.stats import fisher_exact
import numpy as np

help_txt="""
python asapa.py-----help:
Version = 1.0.0
\t
Preparation: the subreads files were processed using pbCCS, Lima, and SUPPA2
\tUsage: python asapa.py build ref.fa  subreads.bam/qry_dir(contain *.bam)/all(current dir)
\tOptional parameters:
\t\t-n                  \tdefault=15, CPU thread number
\t\t-log                \tdefault=no, yes will write terminal print in output.log
\t\t-max_fuzzy_TSS      \tdefault=5,  Max fuzzy TSS dist
\t\t-max_fuzzy_PAS      \tdefault=5,  Max fuzzy PAS dist
\t\t-max_fuzzy_junction \tdefault=5,  Max fuzzy junction dist(from cDNA_cupcake)

Function1: AS vs AS
\tUsage: python asapa.py AS_AS 
\tOptional parameters: 
\t\t-min_ccsnum         \tdefault=10, min ccs number in AS1(2)form1(2)
\t\t-min_dSegmentlen    \tdefault=10, min AS differential segment length
\t\t-min_ccs_usage      \tdefault=0, min ccs usage(used/geneccs)

Function2: AS vs ATI 
\tUsage: python asapa.py AS_ATI
\t\t-min_ccsnum         \tdefault=10, min ccs number in ASform1(2)
\t\t-min_dSegmentlen    \tdefault=10, Min AS differential segment length
\t\t-min_ccs_usage      \tdefault=0, min ccs usage(used/geneccs)
\t\t-min_KS_statistic   \tdefault=0.2, Min KS_statistic

Function3: AS vs APA 
\tUsage: python asapa.py AS_APA
\tOptional parameters:
\t\t-min_ccsnum         \tdefault=10, min ccs number in ASform1(2)
\t\t-min_dSegmentlen    \tdefault=10, Min AS differential segment length
\t\t-min_ccs_usage      \tdefault=0, min ccs usage(used/geneccs)
\t\t-min_KS_statistic   \tdefault=0.2, Min KS_statistic

Function4: ATI vs APA
\tUsage: python asapa.py ATI_APA
\tOptional parameters:
\t\t-min_ccsnum         \tdefault=10, min ccs number of gene/transcript
\t\t-min_geneccs_usage  \tdefault=0, min ccs usage(used/geneccs)
\t\t-min_TSSPASccs_usage\tdefault=0.5, min ccs usage(used/TSSPASccs)
\t\t-min_correlation    \tdefault=0.5, min spearman correlation
\t\t-max_bin_extent     \tdefault=1000,  Max distance of binning extension

The output folder will be created in the current path:
\toutput0_preparation (preparation: subreads to ccs, lima, minimap2, cDNA_cupcake and SUPPA2)
\toutput1_ASAS\t(function1: coupling bewteen AS and AS)
\toutput2_ASATI\t(function2: coupling bewteen AS and ATI)
\toutput3_ASAPA\t(function3: coupling bewteen AS and APA)
\toutput4_ATIAPA\t(function4: coupling bewteen ATI and APA)

Dependency:
Conda is recommended
\t conda install pbbam
\t conda install pbccs==6.4
\t conda install lima
\t conda install minimap2
\t conda install samtools
\t conda install bedtools
\t conda install perl
\t conda install R
\t conda install rpy2
\t conda install r-dplyr
\t conda install -c bioconda trim_isoseq_polya
\t conda install -c bioconda suppa

Independent installation required
\t cDNA_Cupcake

"""
#Sequencing primers of IsoSeq should be specified correctly. Two common primers are provided.
isoseq_primer1=">primer_5p\nAAGCAGTGGTATCAACGCAGAGTACATGGGG\n>primer_3p\nAAGCAGTGGTATCAACGCAGAGTAC\n"
isoseq_primer2=">primer_5p\nATGTAATACGACTCACTATAGGGC\n>primer_3p\nAAAAAAAAAACGCCTGAGA\n"
#help
if (len(sys.argv)==1) or sys.argv[1] in ("h","-h","help","-help"):print(help_txt);sys.exit()
if sys.argv[1] not in ("build","AS_AS","AS_ATI","AS_APA","ATI_APA"):print(help_txt);sys.exit()
if   sys.argv[1] =="build":	outputfile="output0_preparation"
elif sys.argv[1] =="AS_AS":	outputfile="toutput1_ASAS"
elif sys.argv[1] =="AS_ATI":	outputfile="toutput2_ASATI"
elif sys.argv[1] =="AS_APA":	outputfile="toutput3_ASAPA"
elif sys.argv[1] =="ATI_APA":	outputfile="toutput4_ATIAPA"
argument_name_list=[]
argument_index_list=[]  
if sys.argv[1]=="build":
    if len(sys.argv)==4:
        thread="15"
        log="no"
        max_fuzzy_TSS="5"
        max_fuzzy_PAS="5"
        max_fuzzy_junction="5"
        ref=sys.argv[2]
        qry=sys.argv[3]
    elif len(sys.argv)>4:
        if (len(sys.argv)-3)%2!=0:	print("ERROR, the number of parameters is incorrect.");exit()
        if sys.argv[-1][0]=="-":	print("ERROR, the value of "+sys.argv[-1]+" was not sepecified.");exit()
        thread_index=""
        max_fuzzy_TSS_index=""
        max_fuzzy_PAS_index=""
        max_fuzzy_junction_index=""
        for x in sys.argv:
            if "-"==x[0]:
                if x[1:] not in ["n","log","max_fuzzy_TSS","max_fuzzy_PAS","max_fuzzy_junction"]:print("Error, unrecognized parameter: "+x);exit()  
                elif x[1:] in argument_name_list:print("ERROR: duplicated parameter: "+x);exit()
                argument_name_list.append(x[1:])
                argument_index_list.append(i);argument_index_list.append(i+1)
                if	x[1:]=="n":			thread_index=i
                if	x[1:]=="log":			log_index=i
                if	x[1:]=="max_fuzzy_TSS":		max_fuzzy_TSS_index=i
                if	x[1:]=="max_fuzzy_PAS":		max_fuzzy_PAS_index=i
                if	x[1:]=="max_fuzzy_junction":	max_fuzzy_junction_index=i
            i+=1    
        if thread_index=="": 			thread	="15"
        else: 					thread	=sys.argv[thread_index+1]
        if log_index=="":	        	log	="no"
        else: 					log	=sys.argv[log_index+1]
        if max_fuzzy_TSS_index=="":		max_fuzzy_TSS="5"
        else: 					max_fuzzy_TSS=sys.argv[max_fuzzy_TSS_index+1]
        if max_fuzzy_PAS_index=="":		max_fuzzy_PAS="5"
        else: 					max_fuzzy_PAS=sys.argv[max_fuzzy_PAS_index+1]
        if max_fuzzy_junction_index=="":	max_fuzzy_junction="5"
        else: 					max_fuzzy_junction=sys.argv[max_fuzzy_junction_index+1]     
        if int(thread)<=0: print("ERROR, thread must more than 0.");exit()
        if log not in ("yes","no"): print("ERROR, -log should be yes or no.");exit()
        if int(max_fuzzy_TSS)<0: print("ERROR, max_fuzzy_TSS has to be at least 0.");exit()
        if int(max_fuzzy_PAS)<0: print("ERROR, max_fuzzy_PAS has to be at least 0.");exit()
        if int(max_fuzzy_junction)<0: print("ERROR, max_fuzzy_junction has to be at least 0.");exit()
        inputoutput_index_list=[]       
        i=0
        while i<len(sys.argv):
            if i>0 and i not in argument_index_list:
                inputoutput_index_list.append(i)
            i+=1    
        if len(inputoutput_index_list)!=2:
            print("ERROR, Ref.fa and subreads.bam must be specified.");exit()
        else:
            ref=sys.argv[inputoutput_index_list[0]]
            qry=sys.argv[inputoutput_index_list[1]]
    if thread=="15":
        print("\t-n                 \tCPU thread num:         \t15 (default)")
    else:
        print("\t-n                 \tCPU thread num:         \t"+thread+" (default = 15)")
    if max_fuzzy_TSS=="5":
        print("\t-max_fuzzy_TSS     \tMax fuzzy TSS dist:     \t5 (default)")
    else:
        print("\t-max_fuzzy_TSS     \tMax fuzzy TSS dist:     \t"+max_fuzzy_TSS+" (default = 5)")
    if max_fuzzy_PAS=="5":
        print("\t-max_fuzzy_PAS     \tMax fuzzy PAS dist:     \t5 (default)")
    else:
        print("\t-max_fuzzy_PAS     \tMax fuzzy PAS dist:     \t"+max_fuzzy_TSS+" (default = 5)")
    if max_fuzzy_junction=="5":
        print("\t-max_fuzzy_junction\tmax_fuzzy junction dist:\t5 (default)")
    else:
        print("\t-max_fuzzy_junction\tmax_fuzzy_junction dist:\t"+max_fuzzy_junction+" (default = 5)")

if sys.argv[1]=="AS_AS":
    if len(sys.argv)==2:
        thread="15"
        log="no"
        min_ccsnum="10"
        min_dSegmentlen="10"
        min_ccs_usage="0"
    elif len(sys.argv)>2:
        if (len(sys.argv)-2)%2!=0:	print("ERROR, the number of parameters is incorrect.");exit()
        if sys.argv[-1][0]=="-":	print("ERROR, the value of "+sys.argv[-1]+" was not sepecified.");exit()
        i=0
        thread_index=""
        log_index=""
        min_ccsnum_index=""
        min_dSegmentlen_index=""
        min_ccs_usage_index=""
        for x in sys.argv:
            if "-"==x[0]:
                if x[1:] not in ["n","log","min_ccsnum","min_dSegmentlen","min_ccs_usage"]:print("Error, unrecognized parameter: "+x);exit()  
                elif x[1:] in argument_name_list:print("ERROR: duplicated parameter: "+x);exit()
                argument_name_list.append(x[1:])
                argument_index_list.append(i);argument_index_list.append(i+1)
                if	x[1:]=="n":			thread_index=i
                if	x[1:]=="log":			log_index=i
                if	x[1:]=="min_ccsnum":		min_ccsnum_index=i
                if	x[1:]=="min_dSegmentlen":	min_dSegmentlen_index=i
                if	x[1:]=="min_ccs_usage":		min_ccs_usage_index=i
            i+=1    
        if thread_index=="": 			thread		="15"
        else: 					thread		=sys.argv[thread_index+1]
        if log_index=="":	        	log		="no"
        else: 					log		=sys.argv[log_index+1]
        if min_ccsnum_index=="":		min_ccsnum	="10"
        else: 					min_ccsnum	=sys.argv[min_ccsnum_index+1]
        if min_dSegmentlen_index=="":		min_dSegmentlen	="10"
        else: 					min_dSegmentlen	=sys.argv[min_dSegmentlen_index+1]  
        if min_ccs_usage_index=="":		min_ccs_usage	="0"
        else: 					min_ccs_usage	=sys.argv[min_ccs_usage_index+1]  
        if min_ccs_usage_index=="":		min_ccs_usage	="0"
        else: 					min_ccs_usage	=sys.argv[min_ccs_usage_index+1]  
        if int(thread)<=0: print("ERROR, thread must more than 0.");exit()
        if log not in ("yes","no"): print("ERROR, -log should be yes or no.");exit()
        if int(min_ccsnum)<1: print("ERROR, min_ccsnum has to be at least 1.");exit()
        if int(min_dSegmentlen)<0: print("ERROR, min_dSegmentlen has to be at least 0.");exit()
        if float(min_ccs_usage)<0 or float(min_ccs_usage)>1 : print("ERROR, 0<=min_ccs_usage<=1.");exit()
    if min_ccsnum=="10":
        print("\t-min_ccsnum        \t\t10 (default)")
    else:
        print("\t-min_ccsnum        \t\t"+min_ccsnum+" (default = 10)")
    if min_ccs_usage=="0":
        print("\t-min_ccs_usage     \t\t0 (default)")
    else:
        print("\t-min_ccs_usage     \t\t"+min_ccs_usage+" (default = 0)")     
    if min_dSegmentlen=="15":
        print("\t-min_dSegmentlen   \t\t10 (default)")
    else:
        print("\t-min_dSegmentlen   \t\t"+min_dSegmentlen+" (default = 10)")   
  
if sys.argv[1]=="AS_APA":
    if len(sys.argv)==2:
        thread="15"
        log="no"
        min_ccsnum="10"
        min_dSegmentlen="10"
        min_KS_statistic="0.2"
        min_ccs_usage="0"
    elif len(sys.argv)>2:
        if (len(sys.argv)-2)%2!=0:	print("ERROR, the number of parameters is incorrect.");exit()
        if sys.argv[-1][0]=="-":	print("ERROR, the value of "+sys.argv[-1]+" was not sepecified.");exit()
        i=0
        thread_index=""
        log_index=""
        min_ccsnum_index=""
        min_dSegmentlen_index=""
        min_KS_statistic_index=""
        min_ccs_usage_index=""
        for x in sys.argv:
            if "-"==x[0]:
                if x[1:] not in ["n","log","min_ccsnum","min_dSegmentlen","min_KS_statistic","min_ccs_usage"]:print("Error, unrecognized parameter: "+x);exit()  
                elif x[1:] in argument_name_list:print("ERROR: duplicated parameter: "+x);exit()
                argument_name_list.append(x[1:])
                argument_index_list.append(i);argument_index_list.append(i+1)
                if	x[1:]=="n":			thread_index=i
                if	x[1:]=="log":			log_index=i
                if	x[1:]=="min_ccsnum":		min_ccsnum_index=i
                if	x[1:]=="min_dSegmentlen":	min_dSegmentlen_index=i
                if	x[1:]=="min_KS_statistic":	min_KS_statistic_index=i
                if	x[1:]=="min_ccs_usage":		min_ccs_usage_index=i
            i+=1    
        if thread_index=="": 			thread		="15"
        else: 					thread		=sys.argv[thread_index+1]
        if log_index=="":	        	log		="no"
        else: 					log		=sys.argv[log_index+1]
        if min_ccsnum_index=="":		min_ccsnum	="10"
        else: 					min_ccsnum	=sys.argv[min_ccsnum_index+1]
        if min_dSegmentlen_index=="":		min_dSegmentlen	="10"
        else: 					min_dSegmentlen	=sys.argv[min_dSegmentlen_index+1]
        if min_KS_statistic_index=="":		min_KS_statistic="0.2"
        else: 					min_KS_statistic=sys.argv[min_KS_statistic_index+1] 
        if min_ccs_usage_index=="":		min_ccs_usage	="0"
        else: 					min_ccs_usage	=sys.argv[min_ccs_usage_index+1]    
        if int(thread)<=0: print("ERROR, thread must more than 0.");exit()
        if log not in ("yes","no"): print("ERROR, -log should be yes or no.");exit()
        if int(min_ccsnum)<1: print("ERROR, min_ccsnum has to be at least 1.");exit()
        if int(min_dSegmentlen)<0: print("ERROR, min_dSegmentlen has to be at least 0.");exit()
        if float(min_KS_statistic)<0 or float(min_KS_statistic)>1:print("ERROR, 0<=min_KS_statistic<=1.");exit()
        if float(min_ccs_usage)<0 or float(min_ccs_usage)>1 : print("ERROR, 0<=min_ccs_usage<=1.");exit()
    if min_ccsnum=="10":
        print("\t-min_ccsnum        \t\t10 (default)")
    else:
        print("\t-min_ccsnum        \t\t"+min_ccsnum+" (default = 10)")
    if min_ccs_usage=="0":
        print("\t-min_ccs_usage     \t\t0 (default)")
    else:
        print("\t-min_ccs_usage     \t\t"+min_ccs_usage+" (default = 0)")  
    if min_dSegmentlen=="10":
        print("\t-min_dSegmentlen   \t\t10 (default)")
    else:
        print("\t-min_dSegmentlen   \t\t"+min_dSegmentlen+" (default = 10)")       
    if min_KS_statistic=="0.2":
        print("\t-min_KS_statistic  \t\t0.2 (default)")
    else:
        print("\t-min_KS_statistic  \t\t"+min_KS_statistic+" (default = 0.2)")   
  
if sys.argv[1]=="AS_ATI":
    if len(sys.argv)==2:
        thread="15"
        log="no"
        min_ccsnum="10"
        min_dSegmentlen="10"
        min_KS_statistic="0.2"
        min_ccs_usage="0"
    elif len(sys.argv)>2:
        if (len(sys.argv)-2)%2!=0:	print("ERROR, the number of parameters is incorrect.");exit()
        if sys.argv[-1][0]=="-":	print("ERROR, the value of "+sys.argv[-1]+" was not sepecified.");exit()
        i=0
        thread_index=""
        log_index=""
        min_ccsnum_index=""
        min_dSegmentlen_index=""
        min_KS_statistic_index=""
        min_ccs_usage_index=""
        for x in sys.argv:
            if "-"==x[0]:
                if x[1:] not in ["n","log","min_ccsnum_index","min_dSegmentlen","min_KS_statistic","min_ccs_usage"]:print("Error, unrecognized parameter: "+x);exit()  
                elif x[1:] in argument_name_list:print("ERROR: duplicated parameter: "+x);exit()
                argument_name_list.append(x[1:])
                argument_index_list.append(i);argument_index_list.append(i+1)
                if	x[1:]=="n":			thread_index=i
                if	x[1:]=="log":			log_index=i
                if	x[1:]=="min_ccsnum":		min_ccsnum_index=i
                if	x[1:]=="min_dSegmentlen":	min_dSegmentlen_index=i
                if	x[1:]=="min_KS_statistic":	min_KS_statistic_index=i
                if	x[1:]=="min_ccs_usage":		min_ccs_usage_index=i
            i+=1    
        if thread_index=="": 			thread		="15"
        else: 					thread		=sys.argv[thread_index+1]
        if log_index=="":	        	log		="no"
        else: 					log		=sys.argv[log_index+1]
        if min_ccsnum_index=="":		min_ccsnum	="10"
        else: 					min_ccsnum	=sys.argv[min_ccsnum_index+1]
        if min_dSegmentlen_index=="":		min_dSegmentlen	="10"
        else: 					min_dSegmentlen	=sys.argv[min_dSegmentlen_index+1]
        if min_KS_statistic_index=="":		min_KS_statistic="0.2"
        else: 					min_KS_statistic=sys.argv[min_KS_statistic_index+1]   
        if min_ccs_usage_index=="":		min_ccs_usage	="0"
        else: 					min_ccs_usage	=sys.argv[min_ccs_usage_index+1]   
        if int(thread)<=0: print("ERROR, thread must more than 0.");exit()
        if log not in ("yes","no"): print("ERROR, -log should be yes or no.");exit()
        if int(min_ccsnum)<1: print("ERROR, min_ccsnum has to be at least 1.");exit()
        if int(min_dSegmentlen)<0: print("ERROR, min_dSegmentlen has to be at least 0.");exit()
        if float(min_KS_statistic)<0 or float(min_KS_statistic)>1:print("ERROR, 0<=min_KS_statistic<=1.");exit()
        if float(min_ccs_usage)<0 or float(min_ccs_usage)>1 : print("ERROR, 0<=min_ccs_usage<=1.");exit()
    if min_ccsnum=="15":
        print("\t-min_ccsnum        \t\t10 (default)")
    else:
        print("\t-min_ccsnum        \t\t"+min_ccsnum+" (default = 10)")
    if min_ccs_usage=="0":
        print("\t-min_ccs_usage     \t\t0 (default)")
    else:
        print("\t-min_ccs_usage     \t\t"+min_ccs_usage+" (default = 0)")    
    if min_dSegmentlen=="10":
        print("\t-min_dSegmentlen   \t\t10 (default)")
    else:
        print("\t-min_dSegmentlen   \t\t"+min_dSegmentlen+" (default = 10)")       
    if min_KS_statistic=="0.2":
        print("\t-min_KS_statistic  \t\t0.2 (default)")
    else:
        print("\t-min_KS_statistic  \t\t"+min_KS_statistic+" (default = 0.2)")   

if sys.argv[1]=="ATI_APA":
    if len(sys.argv)==2:
        thread="15"
        log="no"
        min_ccsnum="10"
        min_correlation="0.5"
        max_bin_extent="1000"
        min_geneccs_usage="0"
        min_TSSPASccs_usage="0.5"
    elif len(sys.argv)>2:
        if (len(sys.argv)-2)%2!=0:	print("ERROR, the number of parameters is incorrect.");exit()
        if sys.argv[-1][0]=="-":	print("ERROR, the value of "+sys.argv[-1]+" was not sepecified.");exit()
        i=0
        thread_index=""
        log_index=""
        min_ccsnum_index=""
        min_correlation_index=""
        min_geneccs_usage_index=""
        min_TSSPASccs_usage_index=""
        for x in sys.argv:
            if "-"==x[0]:
                if x[1:] not in ["n","log","min_ccsnum","min_correlation","max_bin_extent","min_geneccs_usage","min_TSSPASccs_usage"]:print("Error, unrecognized parameter: "+x);exit()  
                elif x[1:] in argument_name_list:print("ERROR: duplicated parameter: "+x);exit()
                argument_name_list.append(x[1:])
                argument_index_list.append(i);argument_index_list.append(i+1)
                if	x[1:]=="n":			thread_index=i
                if	x[1:]=="log":			log_index=i
                if	x[1:]=="min_ccsnum":		min_ccsnum_index=i
                if	x[1:]=="min_correlation":	min_correlation_index=i
                if	x[1:]=="max_bin_extent":	max_bin_extent_index=i
                if	x[1:]=="min_geneccs_usage":		min_geneccs_usage_index=i
                if	x[1:]=="min_TSSPASccs_usage":		min_TSSPASccs_usage_index=i
            i+=1    
        if thread_index=="": 			thread		="15"
        else: 					thread		=sys.argv[thread_index+1]
        if log_index=="":	        	log		="no"
        else: 					log		=sys.argv[log_index+1]
        if min_ccsnum_index=="":		min_ccsnum	="10"
        else: 					min_ccsnum	=sys.argv[min_ccsnum_index+1]
        if min_correlation_index=="":		min_correlation	="0.5"
        else: 					min_correlation	=sys.argv[min_correlation_index+1]   
        if max_bin_extent_index=="":		max_bin_extent	="1000"
        else: 					max_bin_extent	=sys.argv[min_correlation_index+1]  
        if min_geneccs_usage_index=="":		min_geneccs_usage	="0"
        else: 					min_geneccs_usage	=sys.argv[min_geneccs_usage_index+1] 
        if min_TSSPASccs_usage_index=="":	min_TSSPASccs_usage	="0.5"
        else: 					min_TSSPASccs_usage	=sys.argv[min_TSSPASccs_usage_index+1] 
        if int(thread)<=0: print("ERROR, thread must more than 0.");exit()
        if log not in ("yes","no"): print("ERROR, -log should be yes or no.");exit()
        if int(min_ccsnum)<1: print("ERROR, min_ccsnum has to be at least 1.");exit()
        if float(min_correlation)<0 or float(min_correlation)>1:print("ERROR, 0<=min_KS_statistic<=1.");exit()
        if int(max_bin_extent)<=0: print("ERROR, max_bin_extent must more than 0.");exit()
        if float(min_geneccs_usage_index)<0 or float(min_geneccs_usage_index)>1 : print("ERROR, 0<=min_ccs_usage<=1.");exit()
        if float(min_TSSPASccs_usage_index)<0 or float(min_TSSPASccs_usage_index)>1 : print("ERROR, 0<=min_ccs_usage<=1.");exit()
    if min_ccsnum=="10":
        print("\t-min_ccsnum         \t\t10 (default)")
    else:
        print("\t-min_ccsnum         \t\t"+min_ccsnum+" (default = 10)")
    if min_geneccs_usage=="":
        print("\t-min_geneccs_usage  \t\t0 (default)")
    else:
        print("\t-min_geneccs_usage  \t\t"+min_geneccs_usage+" (default = 0)")   
    if min_geneccs_usage=="":
        print("\t-min_TSSPASccs_usage\t\t0.5 (default)")
    else:
        print("\t-min_TSSPASccs_usage\t\t"+min_TSSPASccs_usage+" (default = 0.5)")  
    if min_correlation=="0.5":
        print("\t-max_bin_extent     \t\t1000 (default)")
    else:
        print("\t-max_bin_extent     \t\t"+max_bin_extent+" (default = 1000)")  
    if min_correlation=="0.5":
        print("\t-min_correlation    \t\t0.5 (default)")
    else:
        print("\t-min_correlation    \t\t"+min_correlation+" (default = 0.5)")   



#pivot.py
pivot_script="""
#!/usr/bin/python
# -*- coding: UTF-8 -*-
import subprocess
import csv
import sys
import timeit
import os
if (len(sys.argv)<2) or sys.argv[1] in ("h","-h","help","-help"):
    print ("pivot.py-----help:")
    print ("pivot.py pivot/unpivot comma/Semicolon head/nohead replace/noreplace input output")     
    print (" ")
    sys.exit()
elif len(sys.argv)!=7:
    print ("pivot.py-----ERROR, argument number must be 6:")
    print ("pivot.py pivot/unpivot comma/Semicolon head/nohead replace/noreplace input output")   
    print (" ")
    sys.exit()
runtype		=sys.argv[1]
delimiter	=sys.argv[2]
inputhead	=sys.argv[3]
inputreplace	=sys.argv[4]
inputfile	=sys.argv[5]
outputfile	=sys.argv[6]
tmpfile0="pivot_tmp0_"+str(timeit.default_timer()).replace(",","")
tmpfile1="pivot_tmp1_"+str(timeit.default_timer()).replace(",","")
tmpfile2="pivot_tmp2_"+str(timeit.default_timer()).replace(",","")
tmpfile3="pivot_tmp3_"+str(timeit.default_timer()).replace(",","")
if os.path.exists(outputfile):
    subprocess.run(["rm "+outputfile],shell=True)  
if runtype not in ["pivot","unpivot"]:
    print("error, argument1 must be pivot/unpivot")
    exit()
if delimiter not in ["comma","semicolon"]:
    print("error, argument2 must be comma/semicolon")
    exit()
if inputhead not in ["head","nohead"]:
    print("error, argument3 must be head/nohead")
    exit()
if inputreplace not in ["replace","noreplace"]:
    print("error, argument4 must be replace/noreplace")
    exit()
if inputhead=="head":
    subprocess.run(["head -n 1 "+inputfile+"  > "+tmpfile0],shell=True) 
    subprocess.run(["sed '1d' "+inputfile+" |cut -f 1,2 > "+tmpfile1],shell=True) 
else:
    subprocess.run(["cut -f 1,2  "+inputfile+" > "+tmpfile1],shell=True) 
if inputreplace=="replace":
    subprocess.run(["sort -n  "+tmpfile1+"  > "+tmpfile2],shell=True) 

else:
    subprocess.run(["sort -n  "+tmpfile1+" |uniq > "+tmpfile2],shell=True) 
if runtype=="unpivot":  
    with open (tmpfile2,"r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            eachline=line.strip()
            eachline_arr=eachline.replace("\\t","_||_").split("_||_")
            col1=eachline_arr[0]
            col2=eachline_arr[1]
            if delimiter=="comma":
                col2_arr=col2.split(",")
            elif delimiter=="semicolon":
                col2_arr=col2.split(";")
            for x in col2_arr:
                with open (tmpfile3,"a",encoding="utf-8") as f2:
                    f2.write(col1+"\\t"+x+"\\n")
                    f2.close()
elif runtype=="pivot":     
    with open (tmpfile2,"r",encoding="ISO-8859-1") as f:
        line_num=len(f.readlines())
    i=0
    col2_arr=[]
    with open (tmpfile2,"r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            eachline=line.strip()
            eachline_arr=eachline.replace("\\t","_||_").split("_||_")
            col1=eachline_arr[0]
            col2=eachline_arr[1];
            if i==1:
                col1_last=col1
                col2_arr.append(col2)
            if i>1 and col1!=col1_last:
                if delimiter=="comma":		
                    col2_str=",".join(col2_arr)
                    if inputreplace=="noreplace":
                        col2_arr=list(set(col2_str.split(",")))
                        col2_str=",".join(col2_arr)
                elif delimiter=="semicolon":	
                    col2_str=";".join(col2_arr)
                    if inputreplace=="noreplace":
                        col2_arr=list(set(col2_str.split(";")))
                        col2_str=";".join(col2_arr)
                with open (tmpfile3,"a",encoding="utf-8") as f2:
                    f2.write(col1_last+"\\t"+col2_str+"\\n")
                    f2.close()
                col1_last=col1
                col2_arr=[]
                col2_arr.append(col2)
            if inputreplace=="replace":
                if col1==col1_last :	
                    col2_arr.append(col2)
            else:
                if col1==col1_last and col2 not in col2_arr:
                    col2_arr.append(col2)
            if i==line_num:
                if delimiter=="comma":		col2_str=",".join(col2_arr)
                elif delimiter=="semicolon":	col2_str=";".join(col2_arr)
                with open (tmpfile3,"a",encoding="utf-8") as f2:
                    f2.write(col1_last+"\\t"+col2_str+"\\n")
                    f2.close()  
if inputhead=="head": 
    subprocess.run(["cat  "+tmpfile0+" "+tmpfile3+" > "+outputfile],shell=True) 
    subprocess.run(["rm  "+tmpfile0],shell=True)
else:
    subprocess.run(["cat  "+tmpfile3+" > "+outputfile],shell=True) 
subprocess.run(["rm  "+tmpfile1+" "+tmpfile2+" "+tmpfile3],shell=True) 
"""
pivot_script=pivot_script.replace('"','\"')
#getfastabylist.pl
getfastabylist_script = """
#! /usr/bin/perl -w
use strict;
die "perl $0 <lst><fa>\\n" unless  @ARGV==2;
my ($lst,$fa)=@ARGV;
open IN,$lst||die;
my %ha;
map{chomp;$ha{(split)[0]}=1}<IN>;
close IN;

$fa=~/gz$/?(open IN,"gzip -cd $fa|"||die):(open IN,$fa||die);
$/=">";<IN>;$/="\\n";
my %out;
while(<IN>){
	my $info=$1 if(/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	$/="\\n";
	$seq=~s/>|\\r|\*//g;
	print ">$info\\n$seq" if(exists $ha{$info} && ! exists $out{$info});
	$out{$info}=1;
}
close IN;
"""
getfastabylist_script=getfastabylist_script.replace('"','\"') 
#Calculate the median
def calcMedian(data):
    if len(data)%2==0:
        median=(sorted(data)[int(len(data)/2)]+sorted(data)[int(len(data)/2)+1])/2
    else:
        median=sorted(data)[int(len(data)/2)]
    return median

print ()
##Write a log.
base_path=os.path.abspath("./")  
if log=="yes":
    class Logger(object):
        logfile =""
        def __init__(self, filename=""):
            self.logfile = filename
            self.terminal = sys.stdout
            return
        def write(self, message):
            self.terminal.write(message)
            if self.logfile != "":
                try:
                    self.log = open(self.logfile, "a")
                    self.log.write(message)
                    self.log.close()
                except:
                    pass
        def flush(self):
            pass
    sys.stdout = Logger(outputfile+"/output.log")
    sys.stderr = Logger(outputfile+"/output.log") 

print()
time_start=timeit.default_timer()
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))  

#################################################################################################################################################################
#################################################################################################################################################################
#Start preparation step.
if sys.argv[1] =="build":
    
    print("Start preparation step.")
    print()
    if "output0_preparation" in os.listdir("./"):
        print("Note: ./output0_preparation exist")
    else: 
        subprocess.run(["mkdir ./output0_preparation"],shell=True)

    ##	
    if "qry_complete" in os.listdir("./output0_preparation/"):
        print("./output0_preparation/qry_complete exist, the completed QRY file will be moved to this folder")
    else: 
        print("Create ./output0_preparation/qry_complete, the completed QRY.bam file will be moved to this folder")
        subprocess.run(["mkdir ./output0_preparation/qry_complete"],shell=True)
    
    ##Specify the input files.
    if qry in ("all","ALL"):
        qryfile_dir = "./"
        qryfolderfile_list=os.listdir(qryfile_dir)
        qry_list=[]  
        for onefile in qryfolderfile_list:
            if onefile[-4:]=='.bam':
                qry_list.append(onefile)
        qry_len=len(qry_list)         
    elif qry not in ("all","ALL") and os.path.isdir(qry)==True:
        if qry[-1]=="/":
            qryfile_dir = qry
        else:
            qryfile_dir = qry+"/"
        qryfolderfile_list=os.listdir(qryfile_dir)
        qry_list=[]  
        for onefile in qryfolderfile_list:
            if onefile[-4:]=='.bam':
                qry_list.append(onefile)
        qry_len=len(qry_list)  
    else:
        if "/" in qry:
            qryfile_strarr = qry.split('/')
            qryfile_dir_arr=qryfile_strarr[:-1]
            qryfile_dir="/".join(qryfile_dir_arr)+"/"
            qryfile = qryfile_strarr[-1]
        else:
            qryfile_dir = "./"
            qryfile=qry   
        qry_list=[qryfile]  
        qry_len=len(qry_list) 
    print ("qry_path:\t"+str(qryfile_dir))   
    print ("qry_num:\t"+str(qry_len))      
    if qry_len==1:
        print ("qry_name:\t"+str(qry_list[0]))  
    if qry_len>1:
        qry_list_str="\n".join(qry_list)
        with open (base_path+"/output0_preparation/qry_list.log","w",encoding='utf-8') as f:
            f.write(qry_list_str+"\n")
            f.close() 
    ##Specify the ref file.
    if os.path.isfile(ref)!=True:
        print ("Error, the ref.fa is not found.")
        exit()
    if "/" in ref:    
        reffile_strarr = ref.split('/')
        reffile_dir_arr=reffile_strarr[:-1]
        reffile_dir="/".join(reffile_dir_arr)+"/"
        reffile = reffile_strarr[-1]
    else:
        reffile_dir = "./"
        reffile=ref  
    print ("ref_path:\t"+str(reffile_dir))   
    if ".fasta" == reffile[-6:]:
        ref_name=reffile[:-6]
    elif ".fas" == reffile[-4:]:
        ref_name=reffile[:-4]
    elif ".fa" == reffile[-3:]:
        ref_name=reffile[:-3]    
    else:
        print ("ERROR: the name of ref file must end with .fasta/.fas/.fa")
        exit()
    print ("ref_name:\t"+str(ref_name))   
    ##step1step2-Start the main program.
    
    print("Steps in /output0_preparation/")
    os.chdir("./output0_preparation")  
    os.system("pwd")
    subprocess.run(["mkdir ./script"],shell=True)
    with open ("./script/pivot.py","w",encoding="utf-8") as f:
        f.write(pivot_script)
        f.close()
    with open ("./script/getfastabylist.pl","w",encoding="utf-8") as f:
        f.write(getfastabylist_script)
        f.close()
    
    if "./1ccs_2lima" in os.listdir("./"):
        print("Note: ./1ccs_2lima exitst")
    else: 
        subprocess.run(["mkdir ./1ccs_2lima"],shell=True)
    
    os.chdir("./1ccs_2lima") 
    
    print()
    i=0
    while i<qry_len:
        qry_file=qry_list[i]
        qry_name=qry_list[i][:-4]
        middle_time1=timeit.default_timer()
        print("     Subreads file progress "+str(i+1)+"/"+str(qry_len)+": "+qry_name)
        print("     "+time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))) 
        subprocess.run(["mkdir out_"+qry_name],shell=True)
        os.chdir("./out_"+qry_name)      
        print("     1-ccs")
        subprocess.run(["mkdir 1-ccs"],shell=True)
        #For RSII data, conda install pbccs=3.4
        #cmd="ccs ../../"+qryfile_dir+qry_file+" --minPasses 1 ./1-ccs/ROI.bam";subprocess.run([cmd],shell=True)      #pbccs=3.4
        #cmd="ccs ../../../"+qryfile_dir+qry_file+" --minPasses 1 --min-rq 0.9 ./1-ccs/ROI.bam";subprocess.run([cmd],shell=True) 	#pbccs=6.4	
        print("     2-lima")
        subprocess.run(["mkdir 2-lima"],shell=True)
        with open ("./2-lima/primers.fasta","w",encoding="ISO-8859-1") as f:
            f.write(isoseq_primer1)
            f.close()    
        cmd="lima  ./1-ccs/ROI.bam ./2-lima/primers.fasta ./2-lima/FLNC.bam --isoseq ";subprocess.run([cmd],shell=True) 
        if "FLNC.primer_5p--primer_3p.bam" not in os.listdir("./2-lima/"):
            with open ("./2-lima/primers.fasta","w",encoding="ISO-8859-1") as f:
                f.write(isoseq_primer2)
                f.close()    
            cmd="lima ./1-ccs/ROI.bam ./2-lima/primers.fasta ./2-lima/FLNC.bam --isoseq";subprocess.run([cmd],shell=True) 
        if "FLNC.primer_5p--primer_3p.bam" not in os.listdir("./2-lima/"):
            print("ERROR, premade IsoSeq primer might not suitable, please manually change it to the correct primer sequence.")  
        if "FLNC.lima.summary" in os.listdir("./2-lima/"):
            subprocess.run(["sed -n '2p' ./2-lima/FLNC.lima.summary > ./2-lima/summary_2p"],shell=True)
            with open("./2-lima/summary_2p","r",encoding="ISO-8859-1") as f:
                summary_2p=f.read()
            subprocess.run(["rm ./2-lima/summary_2p"],shell=True)
            ZMWs_above_thresholds=summary_2p.split("(")[-1].split(")")[0].strip()[:-1]
            if float(ZMWs_above_thresholds)<50:
                print("ZMWs_above_thresholds:"+ZMWs_above_thresholds)
                print("ERROR, premade IsoSeq primer might not suitable, please manually change it to the correct primer sequence.")  
        subprocess.run(["mv ../../../"+qryfile_dir+qry_file+" "+base_path+"/output0_preparation/qry_complete/"],shell=True)
        middle_time2=timeit.default_timer()
        print('     The current subreads file is processed: %.0f Seconds'%(middle_time2-middle_time1))
        
        os.chdir("../") 
        print()
        i+=1   
    os.chdir("../") 
    print ("     Step1 and step2 in all samples were completed.")
    print()
    ##step3-Collect files from step1step2.
    print("     3-all_FLNC")
    subprocess.run(["mkdir ./3-all_FLNC"],shell=True)
    for root, dirs, files in os.walk("./", topdown=False):
        current_folder_arr=dirs
    if "1ccs_2lima" not in current_folder_arr:
        print ("ERROR, folder of ./1ccs_2lima missing.")
        exit()
    for root2, dirs2, files2 in os.walk("./1ccs_2lima/", topdown=False):
        sample_num=len(dirs2)
        outdir_arr=dirs2
    sample_arr=[]
    i=0
    for onesample_outdir in outdir_arr:
        i+=1
        onesample=onesample_outdir[4:]
        sample_arr.append(onesample)
        FLNC_dir="./1ccs_2lima/"+onesample_outdir+"/2-lima/"
        if "FLNC.primer_5p--primer_3p.bam" not in os.listdir(FLNC_dir): 
            print("ERROR, "+FLNC_dir+" doesn't have FLNC.primer_5p--primer_3p.bam");exit()
        qry_size = os.path.getsize(FLNC_dir+"/FLNC.primer_5p--primer_3p.bam")
        print("     File collection:"+str(i)+"/"+str(sample_num)+": "+onesample+'_FLNC.bam\t%.3f' % (qry_size / 1024 / 1024)+' Mbytes')
        subprocess.run(["cp "+FLNC_dir+"/FLNC.primer_5p--primer_3p.bam ./"+onesample+"_FLNC.bam"],shell=True)
        with open (FLNC_dir+"/FLNC.lima.report","r",encoding="ISO-8859-1") as f:
            for line in f.readlines():
                eachline=line.strip()
                eachline_arr=eachline.split("\t")
                newline=eachline_arr[0]+"/ccs\t"+onesample
                with open ("./3-all_FLNC/ccs_sample","a",encoding="utf-8") as f2:
                    f2.write(newline+"\n")
                    f2.close()
    if sample_num==1:
        subprocess.run(["mv *.bam ./3-all_FLNC/all_FLNC.bam"],shell=True)
    else:
        print("     Merge 5p-3p bam files")
        cmd="samtools merge ./3-all_FLNC/all_FLNC.bam *.bam";subprocess.run([cmd],shell=True) 
        subprocess.run(["rm  *.bam"],shell=True)
    cmd="samtools sort -@  "+thread+" -o ./3-all_FLNC/all_FLNC_sort.bam ./3-all_FLNC/all_FLNC.bam 1>/dev/null 2>&1";subprocess.run([cmd],shell=True) 
    subprocess.run(["rm ./3-all_FLNC/all_FLNC.bam"],shell=True)
    cmd="pbindex ./3-all_FLNC/all_FLNC_sort.bam ";subprocess.run([cmd],shell=True) 
    cmd="bedtools bamtofastq -i ./3-all_FLNC/all_FLNC_sort.bam -fq ./3-all_FLNC/all_FLNC.fq";subprocess.run([cmd],shell=True) 
    cmd="awk '{if(NR%4 == 1){print \">\" substr($0, 2)}}{if(NR%4 == 2){print}}' ./3-all_FLNC/all_FLNC.fq > ./3-all_FLNC/all_FLNC.fa";subprocess.run([cmd],shell=True) 
    subprocess.run(["rm ./3-all_FLNC/all_FLNC.fq"],shell=True)
    with open("./3-all_FLNC/all_FLNC.fa",'r') as r:
        with open("./3-all_FLNC/all_FLNC.fa_polyA10","w") as w:
            for eachline in r:
                if eachline.startswith(">"):oneline=eachline.strip()+'\n'
                else:oneline=eachline.strip()+"AAAAAAAAAAAAAAA"+'\n'
                w.write(oneline)
    subprocess.run(["rm ./3-all_FLNC/all_FLNC.fa"],shell=True)
    cmd="trim_isoseq_polyA -i ./3-all_FLNC/all_FLNC.fa_polyA10  -t "+thread+" -G > ./3-all_FLNC/all_FLNC_nopolyA.fa 2> ./3-all_FLNC/all_FLNC_nopolyA.log";subprocess.run([cmd],shell=True) 
    subprocess.run(["rm ./3-all_FLNC/all_FLNC.fa_polyA10"],shell=True)
    subprocess.run(["seqkit rmdup -n ./3-all_FLNC/all_FLNC_nopolyA.fa > ./3-all_FLNC/all_FLNC_nopolyA.fa.uniq 2>/dev/null"],shell=True)
    subprocess.run(["rm ./3-all_FLNC/all_FLNC_nopolyA.fa"],shell=True)
    subprocess.run(["mv ./3-all_FLNC/all_FLNC_nopolyA.fa.uniq ./3-all_FLNC/all_FLNC_nopolyA.fa"],shell=True)

    print ("     4-all_FLNC_minimap2ref")
    subprocess.run(["mkdir ./4-all_FLNC_minimap2ref"],shell=True)
    print ("          Map and sort")
    subprocess.run(["cp "+ref+" ./4-all_FLNC_minimap2ref/ref.fa"],shell=True)
    cmd="minimap2 -ax splice -uf -k 14 -t "+thread+" --secondary=no ./4-all_FLNC_minimap2ref/ref.fa ./3-all_FLNC/all_FLNC_nopolyA.fa > ./4-all_FLNC_minimap2ref/minimap.sam 2>/dev/null";subprocess.run([cmd],shell=True) 
    cmd="samtools view -bS ./4-all_FLNC_minimap2ref/minimap.sam > ./4-all_FLNC_minimap2ref/minimap.bam -@ "+thread;subprocess.run([cmd],shell=True)
    cmd="samtools sort ./4-all_FLNC_minimap2ref/minimap.bam -@ "+thread+" -o ./4-all_FLNC_minimap2ref/minimap.sort.bam 1>/dev/null 2>&1";subprocess.run([cmd],shell=True) 
    cmd="samtools view -h ./4-all_FLNC_minimap2ref/minimap.sort.bam > ./4-all_FLNC_minimap2ref/minimap.sort.sam -@ "+thread;subprocess.run([cmd],shell=True) 
    cmd="samtools index  ./4-all_FLNC_minimap2ref/minimap.sort.bam ";subprocess.run([cmd],shell=True) 
    subprocess.run(["rm ./4-all_FLNC_minimap2ref/minimap.sam ./4-all_FLNC_minimap2ref/minimap.bam"],shell=True) 
    cmd="paftools.js sam2paf ./4-all_FLNC_minimap2ref/minimap.sort.sam > ./4-all_FLNC_minimap2ref/minimap.sort.paf";subprocess.run([cmd],shell=True) 
    
    print ("          Get information of ccs alignment(align_start,align_end,TSS,PAS)")
    with open ("./4-all_FLNC_minimap2ref/FLNC_inform","w",encoding="utf-8") as f2:
        f2.write("ccs_name"+"\t"+"align_chr"+"\t"+"strand"+"\t"+"align_start"+"\t"+"align_end"+"\t"+"TSS"+"\t"+"PAS"+"\t"+"TSS_PAS_mark"+"\t"+"intron_mark"+"\n")
        f2.close() 
    with open  ("./4-all_FLNC_minimap2ref/minimap.sort.paf","r",encoding="ISO-8859-1") as f:
        paf_line_num=len(f.readlines())
    i=0
    with open  ("./4-all_FLNC_minimap2ref/minimap.sort.paf","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            print("          Get ccs information: "+str(i-1)+"/"+str(paf_line_num-1), end="\r")
            eachline= line.strip()
            eachline_arr=eachline.split("\t")
            ccs_name	=eachline_arr[0]
            ccs_len	        =eachline_arr[1]
            ccs_align_start	=eachline_arr[2]    
            ccs_align_end	=eachline_arr[3]
            strand		=eachline_arr[4]
            align_chr	=eachline_arr[5]
            align_start	=str(int(eachline_arr[7])+1)   
            align_end	=eachline_arr[8]
            if strand=="+":	PAS=align_end;		TSS=align_start
            else: 		PAS=align_start;	TSS=align_end
            if    int(ccs_align_start)<=int(max_fuzzy_TSS) and int(ccs_len)-int(ccs_align_end)<=int(max_fuzzy_PAS): TSS_PAS_mark="TSS_PAS"
            elif  int(ccs_align_start)>int(max_fuzzy_TSS)  and int(ccs_len)-int(ccs_align_end)<=int(max_fuzzy_PAS): TSS_PAS_mark="noTSS_PAS"
            elif  int(ccs_align_start)<=int(max_fuzzy_TSS) and int(ccs_len)-int(ccs_align_end)>int(max_fuzzy_PAS):  TSS_PAS_mark="TSS_noPAS"
            else: TSS_PAS_mark="noTSS_noPAS"
            if abs((int(align_end)-int(align_start))-(int(ccs_align_end)-int(ccs_align_start)))<40: intron_mark="nointron"
            else:intron_mark="normal"
            new_line=ccs_name+"\t"+align_chr+"\t"+strand+"\t"+align_start+"\t"+align_end+"\t"+TSS+"\t"+PAS+"\t"+TSS_PAS_mark+"\t"+intron_mark
            with open ("./4-all_FLNC_minimap2ref/FLNC_inform","a",encoding="utf-8") as f2:
                f2.write(new_line+"\n") 
                f2.close() 
    print()  
    #Delete multiple alignment
    subprocess.run(["sort -n ./4-all_FLNC_minimap2ref/FLNC_inform | uniq > ./4-all_FLNC_minimap2ref/FLNC_inform.uniq"],shell=True)
    subprocess.run(["mv ./4-all_FLNC_minimap2ref/FLNC_inform.uniq ./4-all_FLNC_minimap2ref/FLNC_inform"],shell=True)
    subprocess.run(["cut -f 1 ./4-all_FLNC_minimap2ref/FLNC_inform > ./4-all_FLNC_minimap2ref/FLNC_inform.f1"],shell=True)
    subprocess.run(["sort  ./4-all_FLNC_minimap2ref/FLNC_inform.f1|uniq -d > ./4-all_FLNC_minimap2ref/FLNC_inform.f1.non-uniq"],shell=True) 
    
    with open ("./4-all_FLNC_minimap2ref/FLNC_inform.f1.non-uniq","r",encoding="ISO-8859-1") as f:
        delete_list=[col[0] for col in csv.reader(f,delimiter='\t')]
    with open ("./4-all_FLNC_minimap2ref/FLNC_inform","r",encoding="ISO-8859-1") as f:
        FLNA_PAS_len=len(f.readlines())
    i=0;last_ccs="";last_deleted_ccs=""
    with open ("./4-all_FLNC_minimap2ref/FLNC_inform","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            print("          Delete Non-specific alignment: "+str(i-1)+"/"+str(FLNA_PAS_len-1),end="\r")
            eachline=line.strip()     
            eachline_arr=eachline.split("\t")
            one_ccs=eachline_arr[0]
            if last_deleted_ccs!="" and one_ccs!=last_ccs:
                delete_list.remove(last_deleted_ccs)
                last_deleted_ccs=""
            if one_ccs  not in delete_list:
                last_deleted_ccs==""
                with open ("./4-all_FLNC_minimap2ref/FLNC_inform.uniq","a",encoding="utf-8") as f2:
                    f2.write(eachline+"\n")
                    f2.close()
            else: last_deleted_ccs= one_ccs
            last_ccs=one_ccs
    print()
    subprocess.run(["rm ./4-all_FLNC_minimap2ref/FLNC_inform ./4-all_FLNC_minimap2ref/FLNC_inform.f1 ./4-all_FLNC_minimap2ref/FLNC_inform.f1.non-uniq"],shell=True)
    exit()
    print ("     5-cDNA_cupcake")
    subprocess.run(["mkdir ./5-cDNA_cupcake"],shell=True)
    cmd="collapse_isoforms_by_sam.py -c 0.95 -i 0.85 --max_5_diff  10000 --max_3_diff  10000 --max_fuzzy_junction "+str(max_fuzzy_junction)+" --input ./3-all_FLNC/all_FLNC_nopolyA.fa -s ./4-all_FLNC_minimap2ref/minimap.sort.sam -o ./5-cDNA_cupcake/cDNA_cupcake --cpus "+thread+" 1> ./5-cDNA_cupcake/cDNA_cupcake.log1 2> ./5-cDNA_cupcake/cDNA_cupcake.log2"
    subprocess.run([cmd],shell=True) 
    dict_gene_info={};gene_list=[]
    print("          gff to geneid——transcript——transcript.num")
    with open ("./5-cDNA_cupcake/cDNA_cupcake.collapsed.gff","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            if eachline_arr[2]=="transcript":
                chromosome	=eachline_arr[0]
                strand		=eachline_arr[6]
                start		=eachline_arr[3]
                end		=eachline_arr[4]
                gffcol9_arr=eachline_arr[8].split(";")                                          
                transcript_id	=gffcol9_arr[0].strip()[15:-1]
                gene_id		=gffcol9_arr[1].strip()[9:-1]
                if gene_id not in gene_list:	gene_list.append(gene_id)
                with open("./5-cDNA_cupcake/gene_transcript","a",encoding="utf-8") as f2:
                    f2.write(gene_id+"\t"+transcript_id+"\n")
                    f2.close()
                with open("./5-cDNA_cupcake/transcript_info","a",encoding="utf-8") as f3:
                    f3.write(transcript_id+"\t"+chromosome+"\t"+strand+"\t"+start+"\t"+end+"\n")
                    f3.close()                
                if gene_id not in dict_gene_info:		dict_gene_info[gene_id]=[chromosome,strand,start,end]
                else:
                    if int(start)<int(dict_gene_info[gene_id][2]):	dict_gene_info[gene_id][2]=start
                    if int(end)>int(dict_gene_info[gene_id][3]):	dict_gene_info[gene_id][3]=end
    for onegene in gene_list:
        with open ("./5-cDNA_cupcake/gene_info","a",encoding="utf-8") as f4:
            onegene_info	=dict_gene_info[onegene]
            newline		=onegene+"\t"+onegene_info[0]+"\t"+onegene_info[1]+"\t"+onegene_info[2]+"\t"+onegene_info[3]
            f4.write(newline+"\n")
            f4.close()   
    subprocess.run(["python ./script/pivot.py  pivot  comma  nohead  noreplace ./5-cDNA_cupcake/gene_transcript ./5-cDNA_cupcake/gene_transcript.pivot"],shell=True)
    with open ("./5-cDNA_cupcake/gene_transcript.pivot","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            transcript_arr=eachline_arr[1].split(",")
            transcript_num=len(transcript_arr)
            with open("./5-cDNA_cupcake/gene_transcript_transcript.num","a",encoding="utf-8") as f2:
                f2.write(eachline+"\t"+str(transcript_num)+"\n")
                f2.close()        
    print("          cDNA_cupcake.collapsed.group.txt to geneid——ccs——ccs.num")
    r_script = '''
                library(dplyr,warn.conflicts = F)
                file<- read.table("./5-cDNA_cupcake/gene_transcript",header=F,stringsAsFactors = F,check.names = F, sep="")
                names(file) <- c("gene_name","transcript_name")
                index<- read.table("./5-cDNA_cupcake/cDNA_cupcake.collapsed.group.txt",header=F,stringsAsFactors = F,check.names = F, sep="")
                names(index) <- c("transcript_name","ccs_name")
                result <- merge(file,index,by="transcript_name",all.x=FALSE,all.y=FALSE);
                result2 <- result[,c("gene_name","ccs_name")]
                write.table (result2,file ="./5-cDNA_cupcake/gene_ccs0", row.names = FALSE, col.names =TRUE, quote = FALSE,sep="\t");  
    '''
    robjects.r(r_script)
    subprocess.run(["python ./script/pivot.py  pivot   comma    nohead  noreplace    ./5-cDNA_cupcake/gene_ccs0  ./5-cDNA_cupcake/gene_ccs.pivot"],shell=True)
    subprocess.run(["rm   ./5-cDNA_cupcake/gene_ccs0 "],shell=True)
    with open ("./5-cDNA_cupcake/gene_ccs.pivot","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            ccs_arr=eachline_arr[1].split(",")
            ccs_num=len(ccs_arr)
            with open("./5-cDNA_cupcake/gene_ccs_ccs.num","a",encoding="utf-8") as f2:
                f2.write(eachline+"\t"+str(ccs_num)+"\n")
                f2.close()  
    subprocess.run(["python ./script/pivot.py  unpivot   comma    nohead  noreplace    ./5-cDNA_cupcake/gene_ccs.pivot ./5-cDNA_cupcake/gene_ccs  "],shell=True)
    subprocess.run(["rm ./5-cDNA_cupcake/gene_ccs.pivot"],shell=True)
    r_script = '''
                library(dplyr,warn.conflicts = F)
                file<- read.table("./5-cDNA_cupcake/gene_transcript_transcript.num",header=F,stringsAsFactors = F,check.names = F, sep="")
                names(file) <- c("gene_name","transcript_arr","transcript_num")
                index<- read.table("./5-cDNA_cupcake/gene_ccs_ccs.num",header=F,stringsAsFactors = F,check.names = F, sep="")
                names(index) <- c("gene_name","ccs_arr","ccs_num")
                result <- right_join(index,file,by="gene_name");
                result2 <- result[,c("gene_name","transcript_arr","transcript_num","ccs_arr","ccs_num")]
                write.table (result2,file ="./5-cDNA_cupcake/gene_transcript_num_ccs_num", row.names = FALSE, col.names =TRUE, quote = FALSE,sep="\t");  
    '''
    robjects.r(r_script)
    subprocess.run(["rm ./5-cDNA_cupcake/gene_ccs_ccs.num"],shell=True)
    
    print ("     6-suppa") 
    subprocess.run(["mkdir ./6-suppa"],shell=True)
    cmd="suppa.py generateEvents -i ./5-cDNA_cupcake/cDNA_cupcake.collapsed.gff -o ./6-suppa/suppa -e RI SS SE MX -f ioe  2>/dev/null";subprocess.run([cmd],shell=True) 
    
    cmd="cat ./6-suppa/*.ioe > ./6-suppa/AS_all.ioe";subprocess.run([cmd],shell=True) 
    with open ("./6-suppa/AS_all.ioe","r",encoding="ISO-8859-1") as f:
        ioe_line_num=len(f.readlines())
    i=0
    with open ("./6-suppa/AS_All.ioe.simple","w",encoding="utf-8") as f2:
        f2.write("event_id"+"\t"+"transcript1"+"\t"+"transcript2"+"\n")
        f2.close()
    with open ("./6-suppa/AS_all.ioe","r",encoding="ISO-8859-1") as f:
       for line in f.readlines():
            i+=1;transcript2_arr=[]
            print("          Process ioe file: "+str(i-1)+"/"+str(ioe_line_num-1),end="\r")
            eachline=line.strip()
            if "alternative_transcripts" not in eachline:   
                eachline_arr=eachline.split("\t")
                event_id			=eachline_arr[2]
                transcript1_str	=eachline_arr[3]
                transcript1_arr=transcript1_str.split(",")
                total_transcript_str	=eachline_arr[4]
                total_transcript_arr=total_transcript_str.split(",")
                for k in total_transcript_arr:
                    if k not in transcript1_arr:
                        transcript2_arr.append(k)
                transcript2_str=",".join(transcript2_arr)
                newline=event_id+"\t"+transcript1_str+"\t"+transcript2_str
                with open ("./6-suppa/AS_All.ioe.simple","a",encoding="utf-8") as f2:
                    f2.write(newline+"\n")
                    f2.close()
    print() 
    subprocess.run(["rm ./6-suppa/AS_all.ioe"],shell=True) 
    
    print ("     7-ccs_inform")  
    subprocess.run(["mkdir ./7-ccs_inform"],shell=True)
    print("          AS to transcript")
    subprocess.run(["cut -f 1,2 ./6-suppa/AS_All.ioe.simple > ./7-ccs_inform/1-ioe_simple.transcript1"],shell=True)
    subprocess.run(["cut -f 1,3 ./6-suppa/AS_All.ioe.simple > ./7-ccs_inform/1-ioe_simple.transcript2"],shell=True)
    subprocess.run(["sed '1d' ./6-suppa/AS_All.ioe.simple| cut  -f 1 | cut -d ';' -f 1  | sort -n | uniq > ./7-ccs_inform/ASgene.list"],shell=True)
    r_script = '''
                library(dplyr,warn.conflicts = F)
                file<- read.table("./7-ccs_inform/ASgene.list",header=F,stringsAsFactors = F,check.names = F, sep="\t")
                names(file) <- c("gene_name")
                index<- read.table("./5-cDNA_cupcake/gene_transcript_num_ccs_num",header=T,stringsAsFactors = F,check.names = F, sep="\t")
                names(index) <- c("gene_name","transcript_arr","transcript_num","ccs_arr","ccs_num")
                index2 <- index[,c("gene_name","ccs_num")]
                result <- merge(file,index2,by="gene_name",all.x=FALSE,all.y=FALSE);
                result2 <- result[,c("gene_name","ccs_num")]
                write.table (result2,file ="./7-ccs_inform/ASgene_ccsnum", row.names = FALSE, col.names =F, quote = FALSE,sep="\t");  
    '''
    robjects.r(r_script)
    subprocess.run(["python ./script/pivot.py  unpivot   comma    head  noreplace    ./7-ccs_inform/1-ioe_simple.transcript1   ./7-ccs_inform/1-ioe_simple.transcript1.unpivot"],shell=True)
    subprocess.run(["python ./script/pivot.py  unpivot   comma    head  noreplace    ./7-ccs_inform/1-ioe_simple.transcript2   ./7-ccs_inform/1-ioe_simple.transcript2.unpivot"],shell=True)
    print("          AS to transcript to ccs 1/2")
    r_script = '''
                library(dplyr,warn.conflicts = F)
                file<- read.table("./7-ccs_inform/1-ioe_simple.transcript1.unpivot",header=T,stringsAsFactors = F,check.names = F, sep="")
                names(file) <- c("AS_name","PB_name")
                index<- read.table("./5-cDNA_cupcake/cDNA_cupcake.collapsed.group.txt",header=F,stringsAsFactors = F,check.names = F, sep="")
                names(index) <- c("PB_name","ccs_arr")
                result <- merge(file,index,by="PB_name",all.x=FALSE,all.y=FALSE);
                result2 <- result[,c("AS_name","ccs_arr")]
                write.table (result2,file ="./7-ccs_inform/2-AS_ccs.transcript1", row.names = FALSE, col.names =TRUE, quote = FALSE,sep="\t");  
    '''
    robjects.r(r_script)
    print("          AS to transcript to ccs 2/2")
    r_script = '''
                library(dplyr,warn.conflicts = F)
                file<- read.table("./7-ccs_inform/1-ioe_simple.transcript2.unpivot",header=T,stringsAsFactors = F,check.names = F, sep="")
                names(file) <- c("AS_name","PB_name")
                index<- read.table("./5-cDNA_cupcake/cDNA_cupcake.collapsed.group.txt",header=F,stringsAsFactors = F,check.names = F, sep="")
                names(index) <- c("PB_name","ccs_arr")
                result <- merge(file,index,by="PB_name",all.x=FALSE,all.y=FALSE);
                result2 <- result[,c("AS_name","ccs_arr")]
                write.table (result2,file ="./7-ccs_inform/2-AS_ccs.transcript2", row.names = FALSE, col.names =TRUE, quote = FALSE,sep="\t");  
    '''
    robjects.r(r_script)
    subprocess.run(["rm ./7-ccs_inform/1-ioe_simple.transcript1.unpivot ./7-ccs_inform/1-ioe_simple.transcript2.unpivot"],shell=True)
    subprocess.run(["python ./script/pivot.py  unpivot   comma    head  noreplace    ./7-ccs_inform/2-AS_ccs.transcript1   ./7-ccs_inform/2-AS_ccs.transcript1.unpivot"],shell=True)
    subprocess.run(["python ./script/pivot.py  unpivot   comma    head  noreplace    ./7-ccs_inform/2-AS_ccs.transcript2   ./7-ccs_inform/2-AS_ccs.transcript2.unpivot"],shell=True)
    print("          AS to ccs to ccs_inform 1/2")
    r_script = '''
                library(dplyr,warn.conflicts = F)
                file<- read.table("./7-ccs_inform/2-AS_ccs.transcript1.unpivot",header=T,stringsAsFactors = F,check.names = F, sep="")
                names(file) <- c("AS_name","ccs_name")
                index<- read.table("./4-all_FLNC_minimap2ref/FLNC_inform.uniq",header=T,stringsAsFactors = F,check.names = F, sep="")
                names(index) <- c("ccs_name","align_chr","strand","align_start","align_end","TSS","PAS","TSS_PAS_mark")
                result <- merge(file,index,by="ccs_name",all.x=FALSE,all.y=FALSE);
                result2 <- result[,c("AS_name","ccs_name","align_chr","strand","align_start","align_end","TSS","PAS","TSS_PAS_mark")]
                result3 <- arrange(result2,AS_name)
                write.table (result3,file ="./7-ccs_inform/3-AS_ccs_info.transcript1", row.names = FALSE, col.names =TRUE, quote = FALSE,sep="\t");
    '''
    robjects.r(r_script)
    print("          AS to ccs to ccs_inform 2/2")
    r_script = '''
                library(dplyr,warn.conflicts = F)
                file<- read.table("./7-ccs_inform/2-AS_ccs.transcript2.unpivot",header=T,stringsAsFactors = F,check.names = F, sep="")
                names(file) <- c("AS_name","ccs_name")
                index<- read.table("./4-all_FLNC_minimap2ref/FLNC_inform.uniq",header=T,stringsAsFactors = F,check.names = F, sep="")
                names(index) <- c("ccs_name","align_chr","strand","align_start","align_end","TSS","PAS","TSS_PAS_mark")
                result <- merge(file,index,by="ccs_name",all.x=FALSE,all.y=FALSE);
                result2 <- result[,c("AS_name","ccs_name","align_chr","strand","align_start","align_end","TSS","PAS","TSS_PAS_mark")]
                result3 <- arrange(result2,AS_name)
                write.table (result3,file ="./7-ccs_inform/3-AS_ccs_info.transcript2", row.names = FALSE, col.names =TRUE, quote = FALSE,sep="\t");  
    '''
    robjects.r(r_script)
    subprocess.run(["rm ./7-ccs_inform/2-AS_ccs.transcript1.unpivot ./7-ccs_inform/2-AS_ccs.transcript2.unpivot"],shell=True)
    print("          split AS_ccs_info by gene")
    subprocess.run(["mkdir ./7-ccs_inform/ASccs_split"],shell=True)
    with open ("./7-ccs_inform/3-AS_ccs_info.transcript1","r",encoding="ISO-8859-1") as f:
        transcript1_file_num=len(f.readlines())
    i=0
    with open ("./7-ccs_inform/3-AS_ccs_info.transcript1","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i>1:
                print("          Processing transcript1: "+str(i)+"/"+str(transcript1_file_num),end="\r")
                eachline	=line.strip() 
                geneid		=eachline.split(";")[0]
                with open ("./7-ccs_inform/ASccs_split/"+geneid+"_1","a",encoding="utf-8") as f2:
                    f2.write(eachline+"\n")
                    f2.close()
    print()   
    with open ("./7-ccs_inform/3-AS_ccs_info.transcript2","r",encoding="ISO-8859-1") as f:
        transcript2_file_num=len(f.readlines())
    i=0
    with open ("./7-ccs_inform/3-AS_ccs_info.transcript2","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i>1:
                print("          Processing transcript2: "+str(i)+"/"+str(transcript2_file_num),end="\r")
                eachline	=line.strip() 
                geneid		=eachline.split(";")[0]
                with open ("./7-ccs_inform/ASccs_split/"+geneid+"_2","a",encoding="utf-8") as f2:
                    f2.write(eachline+"\n")
                    f2.close()
    print()
#################################################################################################################################################################
#################################################################################################################################################################
#Start AS-AS analysis.
if sys.argv[1] =="AS_AS":
    print("Start AS-AS analysis")
    print("Analysis in /output1_ASAS/")
    if "output1_ASAS" in os.listdir("./"):
        print("Note: ./output1_ASAS/ exist")
    else: 
        subprocess.run(["mkdir ./output1_ASAS"],shell=True)
    os.chdir("./output1_ASAS")
    os.system("pwd")
    print()
    with open ("../output0_preparation/7-ccs_inform/ASgene_ccsnum","r",encoding="ISO-8859-1") as f:
        gene_list=[col[0] for col in csv.reader(f,delimiter='\t')] 
    with open ("../output0_preparation/7-ccs_inform/ASgene_ccsnum","r",encoding="ISO-8859-1") as f:
        ccsnum_list=[col[1] for col in csv.reader(f,delimiter='\t')] 
    gene_num=len(gene_list)
    head_str="eventid\tgene\tgene_reads.num\tAS1\tAS2\tdSegment1\tdSegment2\tdSegmentlen1\tdSegmentlen2\t"+	\
         "AS1form1_AS2form1.reads\tAS1form1_AS2form2.reads\tAS1form2_AS2form1.reads\tAS1form2_AS2form2.reads\t"+	\
         "AS1form1_AS2form1.num\tAS1form1_AS2form2.num\tAS1form2_AS2form1.num\tAS1form2_AS2form2.num\t"+	\
         "read_usage\tfisher_oddsratio\tfisher_pvalue\tchisquare_value\tchisquare_pvalue"
    with open ("AS2AS_fisherchi2","w",encoding="utf-8") as f:
        f.write(head_str+"\n")
        f.close()
    i=0
    while i< gene_num:#< gene_num
        one_gene	=gene_list[i]
        gene_ccs_num	=ccsnum_list[i]
        i+=1
        print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene,end="\r")
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            transcript1_ASlist=[col[0] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            transcript1_ccslist=[col[1] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            transcript2_ASlist=[col[0] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            transcript2_ccslist=[col[1] for col in csv.reader(f,delimiter='\t')]
        transcript1_ASlist_len=len(transcript1_ASlist)
        transcript2_ASlist_len=len(transcript2_ASlist)
        transcript_ASlist_len=transcript1_ASlist_len+transcript2_ASlist_len
        gene_AS_list	=list(set(transcript1_ASlist))
        gene_AS_num	=len(gene_AS_list)
        ccs_list	=list(set(transcript1_ccslist+transcript2_ccslist))
        ccs_num		=len(ccs_list)
        #Get dict of AS infomation: dict[AS]=[oneAS_type/chr/strand/min/max/start/end/dSegment_pos/dSegment_length]
        dict_ASinfo={}
        for oneAS in gene_AS_list:
            oneAS_arr=oneAS.replace(":-",";minus").replace("-",";").replace(":",";").split(";")
            oneAS_type				=oneAS_arr[1]
            oneAS_chr				=oneAS_arr[2]
            oneAS_strand			=oneAS_arr[-1]
            oneAS_pos_arr			=list(map(int,oneAS_arr[3:-1]))   
            oneAS_min				=min(oneAS_pos_arr)
            oneAS_max				=max(oneAS_pos_arr) 
            if oneAS_type=="RI":		
                oneAS_start			=oneAS_pos_arr[1]
                oneAS_end			=oneAS_pos_arr[2]
            else:		
                oneAS_start			=oneAS_min
                oneAS_end			=oneAS_max
            if oneAS_type=="RI":		
                oneAS_dSegment_start		=oneAS_pos_arr[1]+1
                oneAS_dSegment_end		=oneAS_pos_arr[2]-1 
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)       
            if oneAS_type=="A3":
                if oneAS_strand=="+":	
                    oneAS_dSegment_start	=oneAS_pos_arr[1]
                    oneAS_dSegment_end		=oneAS_pos_arr[3]-1
                else:		
                    oneAS_dSegment_start	=oneAS_pos_arr[2]+1
                    oneAS_dSegment_end		=oneAS_pos_arr[0]	
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   
            if oneAS_type=="A5":
                if oneAS_strand=="+":
                    oneAS_dSegment_start	=oneAS_pos_arr[2]+1
                    oneAS_dSegment_end		=oneAS_pos_arr[0]
                else:		
                    oneAS_dSegment_start	=oneAS_pos_arr[1]
                    oneAS_dSegment_end		=oneAS_pos_arr[3]-1
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   	
            if oneAS_type=="SE":
                oneAS_dSegment_start		=oneAS_pos_arr[1]
                oneAS_dSegment_end		=oneAS_pos_arr[2]
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   
            if oneAS_type=="MX":
                oneAS_dSegment_start1		=oneAS_pos_arr[1]
                oneAS_dSegment_end1		=oneAS_pos_arr[2]
                oneAS_dSegment_start2		=oneAS_pos_arr[5]
                oneAS_dSegment_end2		=oneAS_pos_arr[6]
                oneAS_dSegment_start		=oneAS_dSegment_start1
                oneAS_dSegment_end		=oneAS_dSegment_end2
                oneAS_dSegment_length		=oneAS_dSegment_end1-oneAS_dSegment_start1+1+oneAS_dSegment_end2-oneAS_dSegment_start2+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start1)+"-"+str(oneAS_dSegment_end1)+";"+str(oneAS_dSegment_start2)+"-"+str(oneAS_dSegment_end2)
            dict_ASinfo[oneAS]=[oneAS_type,oneAS_chr,oneAS_strand,oneAS_min,oneAS_max,oneAS_start,oneAS_end,oneAS_dSegment_pos,oneAS_dSegment_length]
        ##Get of dict_ASccs[AS,ccs]=["0/1/2","map_start","map_end"]
        dict_ASccs={}
        for oneccs in ccs_list:
            for oneAS in gene_AS_list:
                dict_ASccs[oneAS,oneccs]=["0","",""]
        k=0
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            for line in f.readlines():
                k+=1
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len),end="\r")
                eachline_arr=line.strip().split("\t")
                dict_ASccs[eachline_arr[0],eachline_arr[1]]=["1",eachline_arr[4],eachline_arr[5]]
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            for line in f.readlines():
                k+=1
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len),end="\r")
                eachline_arr=line.strip().split("\t")
                dict_ASccs[eachline_arr[0],eachline_arr[1]]=["2",eachline_arr[4],eachline_arr[5]]        
        #Get paired_AS
        if gene_AS_num>1 : #and gene_AS_num<300
            AS_pairs=[]
            RI_pairs=[]
            for paired_AS1 in gene_AS_list:
                for paired_AS2 in gene_AS_list:
                    paired_AS1_info=dict_ASinfo[paired_AS1]
                    paired_AS2_info=dict_ASinfo[paired_AS2]
                    AS1_type			=paired_AS1_info[0]
                    AS2_type			=paired_AS2_info[0]
                    strand			=paired_AS1_info[2]
                    paired_AS1_min = paired_AS1_info[3];	        paired_AS1_max = paired_AS1_info[4] 
                    paired_AS2_min = paired_AS2_info[3];		paired_AS2_max = paired_AS2_info[4] 
                    mark="no"
                    if paired_AS1_max<=paired_AS2_min:	mark="yes"
                    if paired_AS1_max==paired_AS2_min:
                        if AS1_type=="RI":
                            if AS2_type=="A5" and strand=="+":	mark="no"
                            if AS2_type=="A3" and strand=="-":	mark="no"
                        if AS2_type=="RI":
                            if AS1_type=="A3" and strand=="+":	mark="no"
                            if AS1_type=="A5" and strand=="-":	mark="no"
                    if mark=="yes":
                        AS_pairs.append(paired_AS1+"_||_"+paired_AS2) 
                    if AS1_type=="RI" and AS2_type=="RI" and paired_AS1_min<paired_AS2_max and paired_AS2_min<paired_AS1_max:
                        if paired_AS1_info[5]==paired_AS2_info[5] and paired_AS1_info[6]==paired_AS2_info[6]: 
                            if paired_AS1_min==paired_AS2_min and paired_AS1_max<paired_AS2_max:RI_pairs.append(paired_AS1+"_||_"+paired_AS2)
                            if paired_AS1_min>paired_AS2_min and paired_AS1_max==paired_AS2_max:RI_pairs.append(paired_AS1+"_||_"+paired_AS2)
            #Get_adjacent RI-AS
            RI1234_arr=[]
            RI123_arr=[]
            for oneRIpair1 in RI_pairs:
                for oneRIpair2 in RI_pairs:
                    RI1=oneRIpair1.split("_||_")[0]
                    RI2=oneRIpair1.split("_||_")[1]
                    RI3=oneRIpair2.split("_||_")[0]
                    RI4=oneRIpair2.split("_||_")[1]
                    RI1_info=dict_ASinfo[RI1]
                    RI2_info=dict_ASinfo[RI2]
                    RI3_info=dict_ASinfo[RI3]
                    RI4_info=dict_ASinfo[RI4] 
                    strand			=RI1_info[2]
                    RI1_min = RI1_info[3];	        	RI1_max = RI1_info[4] 
                    RI2_min = RI2_info[3];	        	RI2_max = RI2_info[4] 
                    RI3_min = RI3_info[3];	       	 	RI3_max = RI3_info[4] 
                    RI4_min = RI4_info[3];	       		RI4_max = RI4_info[4]
                    RI1_intron_start=RI1_info[5];		RI1_intron_end=RI1_info[6]
                    RI2_intron_start=RI2_info[5];		RI2_intron_end=RI2_info[6]
                    RI3_intron_start=RI3_info[5];		RI3_intron_end=RI3_info[6]
                    RI4_intron_start=RI4_info[5];		RI4_intron_end=RI4_info[6] 
                    if RI1_min==RI2_min and RI3_max==RI4_max and RI1_intron_end==RI3_min and RI3_intron_start==RI1_max:
                        RI1234_arr.append(RI1+"_||_"+RI2+"_||_"+RI3+"_||_"+RI4)
                for oneAS in gene_AS_list:
                    oneAS_info=dict_ASinfo[oneAS]
                    oneAS_type	=oneAS_info[0]
                    oneAS_start	=oneAS_info[5];	        
                    oneAS_end	=oneAS_info[6];	 
                    if strand=="+" and oneAS_type=="A3" and oneAS_end==RI1_min:	RI123_arr.append(RI1+"_||_"+RI2+"_||_"+oneAS)
                    if strand=="+" and oneAS_type=="A5" and oneAS_start==RI1_max:	RI123_arr.append(RI1+"_||_"+RI2+"_||_"+oneAS)
                    if strand=="-" and oneAS_type=="A3" and oneAS_start==RI1_max:	RI123_arr.append(RI1+"_||_"+RI2+"_||_"+oneAS)
                    if strand=="-" and oneAS_type=="A5" and oneAS_end==RI1_min:	RI123_arr.append(RI1+"_||_"+RI2+"_||_"+oneAS)
            #Get Fisherchi2
            j=0;j_all=len(AS_pairs)
            for onepair in AS_pairs:
                j+=1
                AS1=onepair.split("_||_")[0]
                AS2=onepair.split("_||_")[1]
                AS1_info	=dict_ASinfo[AS1]
                AS2_info	=dict_ASinfo[AS2]
                AS1_start=AS1_info[5];	AS2_start=AS2_info[5];	all_AS_start=min(AS1_start,AS2_start)
                AS1_end=AS1_info[6];	AS2_end=AS2_info[6];	all_AS_end=max(AS1_end,AS2_end)
                AS1form1_ccsnum=AS1form2_ccsnum=AS2form1_ccsnum=AS2form2_ccsnum=0
                AS1form1_AS2form1_ccsnum=AS1form1_AS2form2_ccsnum=AS1form2_AS2form1_ccsnum=AS1form2_AS2form2_ccsnum=0
                AS1form1_AS2form1_ccsarr=[]
                AS1form1_AS2form2_ccsarr=[]
                AS1form2_AS2form1_ccsarr=[]
                AS1form2_AS2form2_ccsarr=[]            
                for oneccs in ccs_list:
                    oneccs_AS1_info			=dict_ASccs[AS1,oneccs]
                    oneccs_AS2_info			=dict_ASccs[AS2,oneccs]
                    if oneccs_AS2_info!=["0","",""]:
                        oneccs_mapstart		=int(oneccs_AS2_info[1])
                        oneccs_mapend		=int(oneccs_AS2_info[2])
                        if all_AS_start>=oneccs_mapstart and all_AS_end<=oneccs_mapend:
                            if     oneccs_AS1_info[0]=="1" and oneccs_AS2_info[0]=="1":	AS1form1_AS2form1_ccsnum+=1;AS1form1_AS2form1_ccsarr.append(oneccs)
                            elif   oneccs_AS1_info[0]=="1" and oneccs_AS2_info[0]=="2":	AS1form1_AS2form2_ccsnum+=1;AS1form1_AS2form2_ccsarr.append(oneccs)
                            elif   oneccs_AS1_info[0]=="2" and oneccs_AS2_info[0]=="1":	AS1form2_AS2form1_ccsnum+=1;AS1form2_AS2form1_ccsarr.append(oneccs)
                            elif   oneccs_AS1_info[0]=="2" and oneccs_AS2_info[0]=="2":	AS1form2_AS2form2_ccsnum+=1;AS1form2_AS2form2_ccsarr.append(oneccs)
                AS1form1_ccsnum	=AS1form1_AS2form1_ccsnum + AS1form1_AS2form2_ccsnum
                AS1form2_ccsnum	=AS1form2_AS2form1_ccsnum + AS1form2_AS2form2_ccsnum
                AS2form1_ccsnum	=AS1form1_AS2form1_ccsnum + AS1form2_AS2form1_ccsnum 
                AS2form2_ccsnum	=AS1form1_AS2form2_ccsnum + AS1form2_AS2form2_ccsnum 
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len)+"\tAS_pair:"+str(j)+"/"+str(j_all)+"\t",end="\r")
                if AS1form1_ccsnum>int(min_ccsnum) and AS1form2_ccsnum>int(min_ccsnum) and AS2form1_ccsnum>int(min_ccsnum) and AS2form2_ccsnum>int(min_ccsnum):
                    obs = np.array([[AS1form1_AS2form1_ccsnum, AS1form1_AS2form2_ccsnum], [AS1form2_AS2form1_ccsnum, AS1form2_AS2form2_ccsnum]])
                    chi2_raw=str(chi2_contingency(obs))
                    chi2_arr=chi2_raw.split(",")
                    chi_square_value	        =chi2_arr[0][1:]
                    chi_square_pvalue		=chi2_arr[1]
                    F_value			=chi2_arr[2]
                    oddsr, fisher_pvalue = fisher_exact(obs, alternative='two-sided')
                    if float(chi_square_pvalue)<1:
                        dSegment1_pos	=AS1_info[7];		dSegment1_length	=AS1_info[8]
                        dSegment2_pos	=AS2_info[7];		dSegment2_length	=AS2_info[8]
                        AS1form1_AS2form1_ccsstr=",".join(AS1form1_AS2form1_ccsarr)
                        AS1form1_AS2form2_ccsstr=",".join(AS1form1_AS2form2_ccsarr)
                        AS1form2_AS2form1_ccsstr=",".join(AS1form2_AS2form1_ccsarr)
                        AS1form2_AS2form2_ccsstr=",".join(AS1form2_AS2form2_ccsarr)
                        eventid=one_gene+"_"+AS1_info[0]+AS2_info[0]+"_"+dSegment1_pos+"_"+dSegment2_pos
                        ccs_usage=(AS1form1_AS2form1_ccsnum+AS1form1_AS2form2_ccsnum+AS1form2_AS2form1_ccsnum+AS1form2_AS2form2_ccsnum)/int(gene_ccs_num)
                        newline=eventid+"\t"+one_gene+"\t"+gene_ccs_num+"\t"+AS1+"\t"+AS2+"\t"+										\
                            dSegment1_pos+"\t"+dSegment2_pos+"\t"+str(dSegment1_length)+"\t"+str(dSegment2_length)+"\t"+							\
                            str(AS1form1_AS2form1_ccsstr)+"\t"+str(AS1form1_AS2form2_ccsstr)+"\t"+str(AS1form2_AS2form1_ccsstr)+"\t"+str(AS1form2_AS2form2_ccsstr)+"\t"+	\
                            str(AS1form1_AS2form1_ccsnum)+"\t"+str(AS1form1_AS2form2_ccsnum)+"\t"+str(AS1form2_AS2form1_ccsnum)+"\t"+str(AS1form2_AS2form2_ccsnum)+"\t"+	\
                            str(ccs_usage)+"\t"+str(oddsr)+"\t"+str(fisher_pvalue)+"\t"+chi_square_value+"\t"+chi_square_pvalue
                        with open ("AS2AS_fisherchi2","a",encoding="utf-8") as f:
                            f.write(newline+"\n")
                            f.close()
            #Get Fisherchi2.RI1234
            j=0;j_all=len(RI1234_arr)
            if RI1234_arr!=[]:print()
            for one4AS in RI1234_arr:
                j+=1
                AS1=one4AS.split("_||_")[0]
                AS2=one4AS.split("_||_")[1]
                AS3=one4AS.split("_||_")[2]
                AS4=one4AS.split("_||_")[3]
                AS1_info=dict_ASinfo[AS1];	AS2_info=dict_ASinfo[AS2];	AS3_info=dict_ASinfo[AS3];	AS4_info=dict_ASinfo[AS4]	
                AS1_start=AS1_info[5];	AS2_start=AS2_info[5];		AS3_start=AS3_info[5];		AS4_start=AS4_info[5];		all_AS_start=min(AS1_start,AS2_start,AS3_start,AS4_start)	
                AS1_end=AS1_info[6];	AS2_end=AS2_info[6];		AS3_end=AS3_info[6];		AS4_end=AS4_info[6];		all_AS_end=max(AS1_end,AS2_end,AS3_end,AS4_end)
                RIform1_ccsnum=RIform2_ccsnum=ASform1_ccsnum=ASform2_ccsnum=0
                RIform1_ASform1_ccsnum=RIform1_ASform2_ccsnum=RIform2_ASform1_ccsnum=RIform2_ASform2_ccsnum=0
                RIform1_ASform1_ccsarr=[]
                RIform1_ASform2_ccsarr=[]
                RIform2_ASform1_ccsarr=[]
                RIform2_ASform2_ccsarr=[] 
                for oneccs in ccs_list:
                    oneccs_AS1_info			=dict_ASccs[AS1,oneccs]
                    oneccs_AS2_info			=dict_ASccs[AS2,oneccs]
                    oneccs_AS3_info			=dict_ASccs[AS3,oneccs]
                    oneccs_AS4_info			=dict_ASccs[AS4,oneccs]
                    if oneccs_AS3_info!=["0","",""] or oneccs_AS4_info!=["0","",""]:
                        if oneccs_AS3_info!=["0","",""]: 	oneccs_mapstart=int(oneccs_AS3_info[1]);	oneccs_mapend=int(oneccs_AS3_info[2])
                        else:				oneccs_mapstart=int(oneccs_AS4_info[1]);	oneccs_mapend=int(oneccs_AS4_info[2])
                        if all_AS_start>=oneccs_mapstart and all_AS_end<=oneccs_mapend:    
                            if   	 oneccs_AS2_info[0]	=="1" and oneccs_AS4_info[0]	=="1":	RIform1_ASform1_ccsnum+=1;RIform1_ASform1_ccsarr.append(oneccs)
                            elif     oneccs_AS1_info[0]	=="1" and oneccs_AS4_info[0]	=="2":	RIform1_ASform2_ccsnum+=1;RIform1_ASform2_ccsarr.append(oneccs)
                            elif     oneccs_AS2_info[0]	=="2" and oneccs_AS3_info[0]	=="1":	RIform2_ASform1_ccsnum+=1;RIform2_ASform1_ccsarr.append(oneccs)
                            elif     oneccs_AS1_info[0]	=="2" and oneccs_AS3_info[0]	=="2":	RIform2_ASform2_ccsnum+=1;RIform2_ASform2_ccsarr.append(oneccs)
                RIform1_ccsnum	=RIform1_ASform1_ccsnum + RIform1_ASform2_ccsnum
                RIform2_ccsnum	=RIform2_ASform1_ccsnum + RIform2_ASform2_ccsnum
                ASform1_ccsnum	=RIform1_ASform1_ccsnum + RIform2_ASform1_ccsnum
                ASform2_ccsnum	=RIform1_ASform2_ccsnum + RIform2_ASform2_ccsnum
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len)+"\tAS_pair:"+str(j)+"/"+str(j_all)+"\tRI12-RI34*",end="\r")
                if RIform1_ccsnum>0 and RIform2_ccsnum>0 and ASform1_ccsnum>0 and ASform2_ccsnum>0:
                    obs = np.array([[RIform1_ASform1_ccsnum, RIform1_ASform2_ccsnum], [RIform2_ASform1_ccsnum, RIform2_ASform2_ccsnum]])
                    chi2_raw=str(chi2_contingency(obs))
                    chi2_arr=chi2_raw.split(",")
                    chi_square_value	        =chi2_arr[0][1:]
                    chi_square_pvalue				=chi2_arr[1]
                    F_value				=chi2_arr[2]
                    oddsr, fisher_pvalue = fisher_exact(obs, alternative='two-sided')
                    if float(chi_square_pvalue)<1:
                        dSegment1_pos	=AS1_info[7];		dSegment1_length	=AS1_info[8]
                        dSegment2_pos	=AS3_info[7];		dSegment2_length	=AS3_info[8]
                        RIform1_ASform1_ccsstr=",".join(RIform1_ASform1_ccsarr)
                        RIform1_ASform2_ccsstr=",".join(RIform1_ASform2_ccsarr)
                        RIform2_ASform1_ccsstr=",".join(RIform2_ASform1_ccsarr)
                        RIform2_ASform2_ccsstr=",".join(RIform2_ASform2_ccsarr)
                        eventid=one_gene+"_"+AS1_info[0]+AS2_info[0]+"_"+dSegment1_pos+"_"+dSegment2_pos
                        ccs_usage=(RIform1_ASform1_ccsnum+RIform1_ASform2_ccsnum+RIform2_ASform1_ccsnum+RIform2_ASform2_ccsnum)/int(gene_ccs_num)
                        newline=eventid+"\t"+one_gene+"\t"+gene_ccs_num+"\t"+AS1+","+AS2+"\t"+AS3+","+AS4+"\t"+									\
                            dSegment1_pos+"\t"+dSegment2_pos+"\t"+str(dSegment1_length)+"\t"+str(dSegment2_length)+"\t"+						\
                            str(RIform1_ASform1_ccsstr)+"\t"+str(RIform1_ASform2_ccsstr)+"\t"+str(RIform2_ASform1_ccsstr)+"\t"+str(RIform2_ASform2_ccsstr)+"\t"+	\
                            str(RIform1_ASform1_ccsnum)+"\t"+str(RIform1_ASform2_ccsnum)+"\t"+str(RIform2_ASform1_ccsnum)+"\t"+str(RIform2_ASform2_ccsnum)+"\t"+	\
                            str(ccs_usage)+"\t"+str(oddsr)+"\t"+str(fisher_pvalue)+"\t"+chi_square_value+"\t"+chi_square_pvalue
                        with open ("AS2AS_fisherchi2","a",encoding="utf-8") as f:
                            f.write(newline+"\n")
                            f.close()
            #Get Fisherchi2.RI123
            j=0;j_all=len(RI123_arr)
            if RI123_arr!=[]:print()
            for one3AS in RI123_arr:
                j+=1
                AS1=one3AS.split("_||_")[0]
                AS2=one3AS.split("_||_")[1]
                AS3=one3AS.split("_||_")[2]
                AS1_info=dict_ASinfo[AS1];	AS2_info=dict_ASinfo[AS2];	AS3_info=dict_ASinfo[AS3]		
                AS1_start=AS1_info[5];	AS2_start=AS2_info[5];		AS3_start=AS3_info[5];		all_AS_start=min(AS1_start,AS2_start,AS3_start)	
                AS1_end=AS1_info[6];	AS2_end=AS2_info[6];		AS3_end	=AS3_info[6];		all_AS_end=max(AS1_end,AS2_end,AS3_end)
                dSegment1_pos	=AS1_info[7];		dSegment1_length	=AS1_info[8]
                dSegment2_pos	=AS3_info[7];		dSegment2_length	=AS3_info[8]
                RIform1_ccsnum=RIform2_ccsnum=ASform1_ccsnum=ASform2_ccsnum=0
                RIform1_ASform1_ccsnum=RIform1_ASform2_ccsnum=RIform2_ASform1_ccsnum=RIform2_ASform2_ccsnum=0
                RIform1_ASform1_ccsarr=[]
                RIform1_ASform2_ccsarr=[]
                RIform2_ASform1_ccsarr=[]
                RIform2_ASform2_ccsarr=[] 
                for oneccs in ccs_list:
                    oneccs_AS1_info			=dict_ASccs[AS1,oneccs]
                    oneccs_AS2_info			=dict_ASccs[AS2,oneccs]
                    oneccs_AS3_info			=dict_ASccs[AS3,oneccs]
                    if oneccs_AS3_info!=["0","",""]:
                        oneccs_mapstart		=int(oneccs_AS3_info[1])
                        oneccs_mapend		=int(oneccs_AS3_info[2])
                        if all_AS_start>=oneccs_mapstart and all_AS_end<=oneccs_mapend: 
                            if   	 oneccs_AS2_info[0]	=="1" and oneccs_AS3_info[0]	=="1":	RIform1_ASform1_ccsnum+=1;RIform1_ASform1_ccsarr.append(oneccs)
                            elif     oneccs_AS1_info[0]	=="1" and oneccs_AS3_info[0]	=="2":	RIform1_ASform2_ccsnum+=1;RIform1_ASform2_ccsarr.append(oneccs)
                            elif     oneccs_AS2_info[0]	=="2" and oneccs_AS3_info[0]	=="1":	RIform2_ASform1_ccsnum+=1;RIform2_ASform1_ccsarr.append(oneccs)
                            elif     oneccs_AS1_info[0]	=="2" and oneccs_AS3_info[0]	=="2":	RIform2_ASform2_ccsnum+=1;RIform2_ASform2_ccsarr.append(oneccs)
                RIform1_ccsnum	=RIform1_ASform1_ccsnum + RIform1_ASform2_ccsnum
                RIform2_ccsnum	=RIform2_ASform1_ccsnum + RIform2_ASform2_ccsnum
                ASform1_ccsnum	=RIform1_ASform1_ccsnum + RIform2_ASform1_ccsnum
                ASform2_ccsnum	=RIform1_ASform2_ccsnum + RIform2_ASform2_ccsnum
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len)+"\tAS_pair:"+str(j)+"/"+str(j_all)+"\tRI12-SS3*",end="\r")
                if RIform1_ccsnum>0 and RIform2_ccsnum>0 and ASform1_ccsnum>0 and ASform2_ccsnum>0:
                    obs = np.array([[RIform1_ASform1_ccsnum, RIform1_ASform2_ccsnum], [RIform2_ASform1_ccsnum, RIform2_ASform2_ccsnum]])
                    chi2_raw=str(chi2_contingency(obs))
                    chi2_arr=chi2_raw.split(",")
                    chi_square_value	        =chi2_arr[0][1:]
                    chi_square_pvalue		=chi2_arr[1]
                    F_value			=chi2_arr[2]
                    oddsr, fisher_pvalue = fisher_exact(obs, alternative='two-sided')
                    if float(chi_square_pvalue):
                        RIform1_ASform1_ccsstr=",".join(RIform1_ASform1_ccsarr)
                        RIform1_ASform2_ccsstr=",".join(RIform1_ASform2_ccsarr)
                        RIform2_ASform1_ccsstr=",".join(RIform2_ASform1_ccsarr)
                        RIform2_ASform2_ccsstr=",".join(RIform2_ASform2_ccsarr)
                        ccs_usage=(RIform1_ASform1_ccsnum+RIform1_ASform2_ccsnum+RIform2_ASform1_ccsnum+RIform2_ASform2_ccsnum)/int(gene_ccs_num)
                        eventid=one_gene+"_"+AS1_info[0]+AS3_info[0]+"_"+dSegment1_pos+"_"+dSegment2_pos
                        newline=eventid+"\t"+one_gene+"\t"+gene_ccs_num+"\t"+AS1+","+AS2+"\t"+AS3+"\t"+										\
                            dSegment1_pos+"\t"+dSegment2_pos+"\t"+str(dSegment1_length)+"\t"+str(dSegment2_length)+"\t"+						\
                            str(RIform1_ASform1_ccsstr)+"\t"+str(RIform1_ASform2_ccsstr)+"\t"+str(RIform2_ASform1_ccsstr)+"\t"+str(RIform2_ASform2_ccsstr)+"\t"+	\
                            str(RIform1_ASform1_ccsnum)+"\t"+str(RIform1_ASform2_ccsnum)+"\t"+str(RIform2_ASform1_ccsnum)+"\t"+str(RIform2_ASform2_ccsnum)+"\t"+	\
                            str(ccs_usage)+"\t"+str(oddsr)+"\t"+str(fisher_pvalue)+"\t"+chi_square_value+"\t"+chi_square_pvalue
                        with open ("AS2AS_fisherchi2","a",encoding="utf-8") as f:
                            f.write(newline+"\n")
                            f.close()
        print()      
    #AS1form1_AS2form1.num\tAS1form1_AS2form2.num\tAS1form2_AS2form1.num\tAS1form2_AS2form2.num
    print("          Filter by pvalue<0.05,min_ccsnum and min_dSegmentlen")
    with open ("AS2AS_fisherchi2.pvalue0.05","w",encoding="utf-8") as f2:
        f2.write(head_str+"\n")
        f2.close()
    i=0
    with open ("AS2AS_fisherchi2","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1:  continue
            eachline		=line.strip()
            eachline_arr		=eachline.split("\t")
            dSegmentlen1		=int(eachline_arr[7])
            dSegmentlen2		=int(eachline_arr[8])
            AS1form1_AS2form1_ccsnum	=int(eachline_arr[13])
            AS1form1_AS2form2_ccsnum	=int(eachline_arr[14])
            AS1form2_AS2form1_ccsnum	=int(eachline_arr[15])
            AS1form2_AS2form2_ccsnum	=int(eachline_arr[16])
            ccs_usage			=float(eachline_arr[17])
            fisher_pvalue		=float(eachline_arr[19])
            chi_square_pvalue		=float(eachline_arr[21])
            AS1form1_ccsnum	=AS1form1_AS2form1_ccsnum + AS1form1_AS2form2_ccsnum
            AS1form2_ccsnum	=AS1form2_AS2form1_ccsnum + AS1form2_AS2form2_ccsnum
            AS2form1_ccsnum	=AS1form1_AS2form1_ccsnum + AS1form2_AS2form1_ccsnum
            AS2form2_ccsnum	=AS1form1_AS2form2_ccsnum + AS1form2_AS2form2_ccsnum
            if AS1form1_ccsnum>=int(min_ccsnum) and AS1form2_ccsnum>=int(min_ccsnum) and AS2form1_ccsnum>=int(min_ccsnum) and AS2form2_ccsnum>=int(min_ccsnum):
                if dSegmentlen1>=int(min_dSegmentlen) and dSegmentlen2>=int(min_dSegmentlen)  and ccs_usage>=float(min_ccs_usage):
                    if fisher_pvalue<0.05 or chi_square_pvalue<0.05:
                        with open ("AS2AS_fisherchi2.pvalue0.05","a",encoding="utf-8") as f2:
                            f2.write(eachline+"\n")
                            f2.close()                      
    subprocess.run(["sort -g -k 22  AS2AS_fisherchi2.pvalue0.05|cut -f 1-9,14-22   > AS2AS_fisherchi2.pvalue0.05.simple"],shell=True)
    with open("AS2AS_fisherchi2.pvalue0.05.simple","r",encoding="ISO-8859-1") as f:
        result_num=len(f.readlines())-1
    print()
    print("     Complete!")
    print("     Possible results number is "+str(result_num))
    print()
    print("     part_ccs2ref. This section gives the BAM files which can be used as reference for credibility.")
    print("          All ccs were split into two parts: AS1form1 + AS1form2")
    print("          Get ccs(pvalue0.05) list from AS2AS_fisherchi2.pvalue0.05(simple)")
    print("          Get fasta from transcript1 and 2 in each AS and then minimap2ref")
    if "part_ccs2ref" in os.listdir("./"):
        print("          Note: ./part_ccs2ref/ exist and will be deleted.") 
        subprocess.run(["rm -r ./part_ccs2ref"],shell=True) 
    subprocess.run(["mkdir ./part_ccs2ref"],shell=True)
    i=0
    with open("AS2AS_fisherchi2.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i>1:
                eachline=line.strip()
                eachline_arr=eachline.split("\t")
                ccs1_str=eachline_arr[9]
                ccs2_str=eachline_arr[10]
                ccs3_str=eachline_arr[11]
                ccs4_str=eachline_arr[12]
                ccs_list=ccs1_str.split(",")+ccs2_str.split(",")+ccs3_str.split(",")+ccs4_str.split(",")
                for x in ccs_list:
                    with open("./part_ccs2ref/ccs.list0","a",encoding="utf-8") as f2:
                        f2.write(x+"\n")   
                        f2.close()
    subprocess.run(["sort -n ./part_ccs2ref/ccs.list0 | uniq > ./part_ccs2ref/1-ccs.pvalue0.05.list"],shell=True)
    subprocess.run(["rm ./part_ccs2ref/ccs.list0"],shell=True)
    subprocess.run(["pwd"],shell=True)
    print("          Extract fasta sequence")
    subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/1-ccs.pvalue0.05.list ../output0_preparation/3-all_FLNC/all_FLNC_nopolyA.fa > ./part_ccs2ref/2-ccs.pvalue0.05.fa"],shell=True)
    i=0

    if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"):
        subprocess.run(["rm -r ./part_ccs2ref/ccs_pvalue0.05_eachAS/"],shell=True)  
    subprocess.run(["mkdir ./part_ccs2ref/ccs_pvalue0.05_eachAS"],shell=True)
    i=0
    with open("AS2AS_fisherchi2.pvalue0.05","r",encoding="ISO-8859-1") as f:
        file_line_num=len(f.readlines())
    with open("AS2AS_fisherchi2.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            time1=timeit.default_timer()
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            eventid	=eachline_arr[0]
            ccs1_str	=eachline_arr[9]+","+eachline_arr[10]
            ccs2_str	=eachline_arr[11]+","+eachline_arr[12]
            ccs1_list=ccs1_str.split(",")
            ccs2_list=ccs2_str.split(",")
            if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+".list1"):subprocess.run(["rm  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+".list1"],shell=True)  
            if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+".list2"):subprocess.run(["rm  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+".list2"],shell=True)  
            for x in ccs1_list:
                with open("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+".list1","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()            
            for x in ccs2_list:
                with open("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+".list2","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()  
            part_arr=["1","2"]
            for part in part_arr:
                subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+".list"+part+" ./part_ccs2ref/2-ccs.pvalue0.05.fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".fa"],shell=True)
                subprocess.run(["minimap2 -ax splice -uf -k 14 -t "+thread+" --secondary=no ../output0_preparation/4-all_FLNC_minimap2ref/ref.fa ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".sam 2>/dev/null"],shell=True) 
                subprocess.run(["samtools view -bS ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".sam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".bam -@ "+thread],shell=True) 
                subprocess.run(["samtools sort ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".bam -@ "+thread+" -o ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"sort_"+part+".bam 1>/dev/null 2>&1"],shell=True) 
                subprocess.run(["samtools view -h ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"sort_"+part+".bam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"sort_"+part+".sam -@ "+thread],shell=True) 
                subprocess.run(["samtools index  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"sort_"+part+".bam"],shell=True) 
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".sam ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_"+part+".bam "],shell=True) 
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"sort_"+part+".sam "],shell=True)
            subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_1.fa  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid+"_2.fa"],shell=True)
            time2=timeit.default_timer()
            print("          -----Processing "+str(i-1)+"/"+str(file_line_num-1)+': %.0f Seconds'%(time2-time1),end="\r")
    print ()

#################################################################################################################################################################
#################################################################################################################################################################
#Start AS-APA analysis.
if sys.argv[1] =="AS_APA":
    print("Start AS-APA analysis")
    print("Analysis in /output3_ASAPA/")
    if "output3_ASAPA" in os.listdir("./"):
        print("Note: ./output3_ASAPA/ exist")
    else: 
        subprocess.run(["mkdir ./output3_ASAPA"],shell=True)
    os.chdir("./output3_ASAPA")
    os.system("pwd")
    
    print()
    with open ("../output0_preparation/7-ccs_inform/ASgene_ccsnum","r",encoding="ISO-8859-1") as f:
        gene_list=[col[0] for col in csv.reader(f,delimiter='\t')] 
    with open ("../output0_preparation/7-ccs_inform/ASgene_ccsnum","r",encoding="ISO-8859-1") as f:
        ccsnum_list=[col[1] for col in csv.reader(f,delimiter='\t')] 
    gene_num=len(gene_list)

    head_str="gene\tgene_reads.num\tAS\tAS_pos\tdSegment\tdSegmentlen\tPAS_zero\t"+	\
         "PAS_all_median\tAS_PAS_distance\tPAS1_median\tPAS2_median\t"+		\
         "ASform1.reads\tASform2.reads\tASform1.PAS\tASform2.PAS\t"+		\
         "ASform1.dPAS_raw\tASform2.dPAS_raw\tASform1.dPAS\tASform2.dPAS\t"+	\
         "ASform1.reads.num\tASform2.reads.num\t"+					\
         "read_usage\tKS_statistic\tp_value"
    with open ("AS2APA_KS","w",encoding="utf-8") as f:
        f.write(head_str+"\n")
        f.close()
    i=0
    while i< gene_num:#< gene_num
        one_gene	=gene_list[i]
        gene_ccs_num	=ccsnum_list[i]
        i+=1
        print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene,end="\r")
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            transcript1_ASlist=[col[0] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            transcript1_ccslist=[col[1] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            transcript2_ASlist=[col[0] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            transcript2_ccslist=[col[1] for col in csv.reader(f,delimiter='\t')]
        transcript1_ASlist_len=len(transcript1_ASlist)
        transcript2_ASlist_len=len(transcript2_ASlist)
        transcript_ASlist_len=transcript1_ASlist_len+transcript2_ASlist_len
        gene_AS_list	=list(set(transcript1_ASlist))
        gene_AS_num	=len(gene_AS_list)
        ccs_list	=list(set(transcript1_ccslist+transcript2_ccslist))
        ccs_num		=len(ccs_list)
        #Get dict of AS infomation: dict[AS]=[oneAS_type/chr/strand/min/max/start/end/dSegment_pos/dSegment_length]
        dict_ASinfo={}
        for oneAS in gene_AS_list:
            oneAS_arr=oneAS.replace(":-",";minus").replace("-",";").replace(":",";").split(";")
            oneAS_type				=oneAS_arr[1]
            oneAS_chr				=oneAS_arr[2]
            oneAS_strand			=oneAS_arr[-1]
            oneAS_pos_arr			=list(map(int,oneAS_arr[3:-1]))   
            oneAS_min				=min(oneAS_pos_arr)
            oneAS_max				=max(oneAS_pos_arr) 
            if oneAS_type=="RI":		
                oneAS_start			=oneAS_pos_arr[1]
                oneAS_end			=oneAS_pos_arr[2]
            else:		
                oneAS_start			=oneAS_min
                oneAS_end			=oneAS_max
            if oneAS_type=="RI":		
                oneAS_dSegment_start		=oneAS_pos_arr[1]+1
                oneAS_dSegment_end		=oneAS_pos_arr[2]-1 
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)       
            if oneAS_type=="A3":
                if oneAS_strand=="+":	
                    oneAS_dSegment_start	=oneAS_pos_arr[1]
                    oneAS_dSegment_end		=oneAS_pos_arr[3]-1
                else:		
                    oneAS_dSegment_start	=oneAS_pos_arr[2]+1
                    oneAS_dSegment_end		=oneAS_pos_arr[0]	
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   
            if oneAS_type=="A5":
                if oneAS_strand=="+":
                    oneAS_dSegment_start	=oneAS_pos_arr[2]+1
                    oneAS_dSegment_end		=oneAS_pos_arr[0]
                else:		
                    oneAS_dSegment_start	=oneAS_pos_arr[1]
                    oneAS_dSegment_end		=oneAS_pos_arr[3]-1
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   	
            if oneAS_type=="SE":
                oneAS_dSegment_start		=oneAS_pos_arr[1]
                oneAS_dSegment_end		=oneAS_pos_arr[2]
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   
            if oneAS_type=="MX":
                oneAS_dSegment_start1		=oneAS_pos_arr[1]
                oneAS_dSegment_end1		=oneAS_pos_arr[2]
                oneAS_dSegment_start2		=oneAS_pos_arr[5]
                oneAS_dSegment_end2		=oneAS_pos_arr[6]
                oneAS_dSegment_start		=oneAS_dSegment_start1
                oneAS_dSegment_end		=oneAS_dSegment_end2
                oneAS_dSegment_length		=oneAS_dSegment_end1-oneAS_dSegment_start1+1+oneAS_dSegment_end2-oneAS_dSegment_start2+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start1)+"-"+str(oneAS_dSegment_end1)+";"+str(oneAS_dSegment_start2)+"-"+str(oneAS_dSegment_end2)
            dict_ASinfo[oneAS]=[oneAS_type,oneAS_chr,oneAS_strand,oneAS_min,oneAS_max,oneAS_start,oneAS_end,oneAS_dSegment_pos,oneAS_dSegment_length]
        ##Get of dict_ASccs[AS,ccs]=["0/1/2","map_start","map_end"]
        dict_ASccs={}
        for oneccs in ccs_list:
            for oneAS in gene_AS_list:
                dict_ASccs[oneAS,oneccs]=["0","","",""]
        k=0
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            for line in f.readlines():
                k+=1
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len),end="\r")
                eachline_arr=line.strip().split("\t")
                if eachline_arr[8] in ("noTSS_PAS","TSS_PAS"):
                    dict_ASccs[eachline_arr[0],eachline_arr[1]]=["1",eachline_arr[4],eachline_arr[5],eachline_arr[7]]
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            for line in f.readlines():
                k+=1
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len),end="\r")
                eachline_arr=line.strip().split("\t")
                if eachline_arr[8] in ("noTSS_PAS","TSS_PAS"):
                    dict_ASccs[eachline_arr[0],eachline_arr[1]]=["2",eachline_arr[4],eachline_arr[5],eachline_arr[7]]  
        print()
        for oneAS in gene_AS_list:
            oneAS_info		=dict_ASinfo[oneAS];
            oneAS_type		=oneAS_info[0]	
            oneAS_chr		=oneAS_info[1]	
            oneAS_strand	=oneAS_info[2]
            oneAS_start		=oneAS_info[5]
            oneAS_end		=oneAS_info[6]
            oneAS_pos		=oneAS_chr+":"+str(oneAS_start)+"-"+str(oneAS_end)
            oneAS_dSegment	=oneAS_info[7] 
            oneAS_dSegment_len	=oneAS_info[8] 
            ccs1_arr=[];ccs2_arr=[]
            PAS1_arr=[];PAS2_arr=[]
            #Process each ccs
            for oneccs in ccs_list:
                oneccs_AS_info		=dict_ASccs[oneAS,oneccs]
                if oneccs_AS_info!=["0","","",""]:
                    oneccs_mapstart	=oneccs_AS_info[1]
                    oneccs_mapend	=oneccs_AS_info[2]
                    oneccs_PAS		=oneccs_AS_info[3]
                    if int(oneccs_mapstart)<=int(oneAS_start) and int(oneAS_end)<=int(oneccs_mapend):
                        if   oneccs_AS_info[0]=="1": ccs1_arr.append(oneccs);	PAS1_arr.append(oneccs_PAS)
                        elif oneccs_AS_info[0]=="2": ccs2_arr.append(oneccs);	PAS2_arr.append(oneccs_PAS)
            ccs1_num		=len(ccs1_arr)
            ccs2_num		=len(ccs2_arr)
            if ccs1_num>2 and ccs2_num>2:
                ccs1_arr_str	=",".join(ccs1_arr)
                ccs2_arr_str	=",".join(ccs2_arr)
                PAS1_arr_str	=",".join(PAS1_arr)
                PAS2_arr_str	=",".join(PAS2_arr)
                PAS_arr		=list(set(PAS1_arr+PAS2_arr));		
                PAS_int_arr	=list(map(int,PAS_arr));	PAS_int_arr.sort()
                PAS1_int_arr	=list(map(int,PAS1_arr));	PAS1_int_arr.sort();	
                PAS2_int_arr	=list(map(int,PAS2_arr));       PAS2_int_arr.sort();
                PAS1_median	=calcMedian(PAS1_int_arr)
                PAS2_median	=calcMedian(PAS2_int_arr)
                PAS_all_median	=calcMedian(PAS1_int_arr+PAS2_int_arr)
                AS_PAS_distance	=min(abs(PAS_all_median-oneAS_end),abs(PAS_all_median-oneAS_start))
                #Get PAS_zero
                if oneAS_strand=="+": 	PAS_zero=max(PAS_int_arr)
                else:			PAS_zero=min(PAS_int_arr)
                #Standardization by PAS_zero
                dPAS1_int_arr=[] 
                dPAS2_int_arr=[]
                dPAS1_str_arr=[]
                dPAS2_str_arr=[]
                for x in PAS1_int_arr:	
                     dPAS1_int_arr.append(abs(x-PAS_zero))
                     dPAS1_str_arr.append(str(abs(x-PAS_zero)))
                for y in PAS2_int_arr:	
                     dPAS2_int_arr.append(abs(y-PAS_zero))
                     dPAS2_str_arr.append(str(abs(y-PAS_zero)))                 
                standard_PAS1_str	=",".join(dPAS1_str_arr)
                standard_PAS2_str	=",".join(dPAS2_str_arr)
                dPAS1_int_arr_noreplace=list(set(dPAS1_int_arr));	dPAS1_int_arr_noreplace.sort()
                dPAS2_int_arr_noreplace=list(set(dPAS2_int_arr));	dPAS2_int_arr_noreplace.sort()
                standard_PAS1_arr_good=[]
                standard_PAS2_arr_good=[]                    
                name_num_dict = {}
                for key in dPAS1_int_arr:
                    name_num_dict[key] 		= name_num_dict.get(key, 0) + 1
                for x in dPAS1_int_arr_noreplace:
                    each_dPAS_num			=name_num_dict[x]
                    each_dPAS			=str(x)+"("+str(each_dPAS_num)+")"
                    standard_PAS1_arr_good.append(each_dPAS)
                name_num_dict = {}
                for key in dPAS2_int_arr:
                    name_num_dict[key] = name_num_dict.get(key, 0) + 1
                for x in dPAS2_int_arr_noreplace:
                    each_dPAS_num			=name_num_dict[x]
                    each_dPAS			=str(x)+"("+str(each_dPAS_num)+")"
                    standard_PAS2_arr_good.append(each_dPAS)
                standard_PAS1_str_good		=",".join(standard_PAS1_arr_good)
                standard_PAS2_str_good		=",".join(standard_PAS2_arr_good)
                ccs_usage				=int(ccs1_num+ccs2_num)/int(gene_ccs_num)
                KS_str=str(ks_2samp(dPAS1_int_arr,dPAS2_int_arr))[13:][:-1]
                KS_str_arr=KS_str.split(", ")
                KS_statistic	=KS_str_arr[0][10:]
                pvalue		=KS_str_arr[1][7:]
                newline=one_gene+"\t"+gene_ccs_num+"\t"+oneAS+"\t"+oneAS_pos+"\t"+oneAS_dSegment+"\t"+str(oneAS_dSegment_len)+"\t"+str(PAS_zero)+"\t"+		\
                        str(PAS_all_median)+"\t"+str(AS_PAS_distance)+"\t"+str(PAS1_median)+"\t"+str(PAS2_median)+"\t"+			\
                        ccs1_arr_str+"\t"+ccs2_arr_str+"\t"+PAS1_arr_str+"\t"+PAS2_arr_str+"\t"+					\
                        standard_PAS1_str+"\t"+standard_PAS2_str+"\t"+standard_PAS1_str_good+"\t"+standard_PAS2_str_good+"\t"+		\
                        str(ccs1_num)+"\t"+str(ccs2_num)+"\t"+										\
                        str(ccs_usage)+"\t"+KS_statistic+"\t"+pvalue
                with open ("AS2APA_KS","a",encoding="utf-8") as f2:
                     f2.write(newline+"\n")
                     f2.close()
    print("          Filter by pvalue<0.05, min_ccsnum and min_dSegmentlen")
    with open ("AS2APA_KS.pvalue0.05","w",encoding="utf-8") as f2:
        f2.write(head_str+"\n")
        f2.close()
    i=0
    with open ("AS2APA_KS","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            eachline		=line.strip()
            eachline_arr	=eachline.split("\t")
            dSegment_length	=int(eachline_arr[5])
            ccs1_num		=int(eachline_arr[19])
            ccs2_num		=int(eachline_arr[20])
            ccs_usage		=float(eachline_arr[21])
            KS_statistic	=float(eachline_arr[22])
            pvalue		=float(eachline_arr[23])
            if ccs1_num>int(min_ccsnum) and ccs2_num>int(min_ccsnum):
                if dSegment_length>int(min_dSegmentlen) and KS_statistic>=float(min_KS_statistic) and pvalue<0.05 and ccs_usage>=float(min_ccs_usage):
                    with open ("AS2APA_KS.pvalue0.05","a",encoding="utf-8") as f2:
                        f2.write(eachline+"\n")
                        f2.close()                       
    subprocess.run(["sort -g -k 24 AS2APA_KS.pvalue0.05 | cut -f 1-11,18-24   > AS2APA_KS.pvalue0.05.simple"],shell=True)
    with open("AS2APA_KS.pvalue0.05.simple","r",encoding="ISO-8859-1") as f:
        result_num=len(f.readlines())-1
    print()
    print("     Complete!")
    print("     Possible results number is "+str(result_num))
    print()
    print("     part_ccs2ref. This section gives the BAM files which can be used as reference for credibility.")
    print("          All ccs were split into two parts: ASform1 + ASform2")
    print("          Get ccs(pvalue0.05) list from AS2APA_KS.pvalue0.05(simple)")
    print("          Get fasta from transcript1 and 2 in each AS and then minimap2ref")
    if "part_ccs2ref" in os.listdir("./"):
        print("          Note: ./part_ccs2ref/ exist and will be deleted.") 
        subprocess.run(["rm -r ./part_ccs2ref"],shell=True) 
    subprocess.run(["mkdir ./part_ccs2ref"],shell=True)
    i=0
    with open("AS2APA_KS.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            ccs1_str=eachline_arr[11]
            ccs2_str=eachline_arr[12]
            ccs_list=ccs1_str.split(",")+ccs2_str.split(",")
            for x in ccs_list:
                with open("./part_ccs2ref/ccs.list0","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()
    subprocess.run(["sort -n ./part_ccs2ref/ccs.list0 | uniq > ./part_ccs2ref/1-ccs.pvalue0.05.list"],shell=True)
    subprocess.run(["rm ./part_ccs2ref/ccs.list0"],shell=True)

    print("          Extract fasta sequence")
    subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/1-ccs.pvalue0.05.list ../output0_preparation/3-all_FLNC/all_FLNC_nopolyA.fa > ./part_ccs2ref/2-ccs.pvalue0.05.fa"],shell=True)
    i=0
    if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"):
        subprocess.run(["rm -r ./part_ccs2ref/ccs_pvalue0.05_eachAS/"],shell=True)  
    subprocess.run(["mkdir ./part_ccs2ref/ccs_pvalue0.05_eachAS"],shell=True)
    i=0
    with open("AS2APA_KS.pvalue0.05","r",encoding="ISO-8859-1") as f:
        file_line_num=len(f.readlines())
    with open("AS2APA_KS.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            time1=timeit.default_timer()
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            eventid1	=eachline_arr[2]
            eventid2	=eachline_arr[2].replace(";","\;")
            ccs1_str	=eachline_arr[11]
            ccs2_str	=eachline_arr[12]
            ccs1_list=ccs1_str.split(",")
            ccs2_list=ccs2_str.split(",")
            if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list1"):subprocess.run(["rm  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list1"],shell=True)  
            if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list2"):subprocess.run(["rm  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list2"],shell=True)  
            for x in ccs1_list:
                with open("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid1+".list1","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()            
            for x in ccs2_list:
                with open("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid1+".list2","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()  
            part_arr=["1","2"]
            for part in part_arr:
                subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list"+part+" ./part_ccs2ref/2-ccs.pvalue0.05.fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".fa"],shell=True)
                subprocess.run(["minimap2 -ax splice -uf -k 14 -t "+thread+" --secondary=no ../output0_preparation/4-all_FLNC_minimap2ref/ref.fa ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".sam 2>/dev/null"],shell=True) 
                subprocess.run(["samtools view -bS ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".sam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".bam -@ "+thread],shell=True) 
                subprocess.run(["samtools sort ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".bam -@ "+thread+" -o ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".bam 1>/dev/null 2>&1"],shell=True) 
                subprocess.run(["samtools view -h ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".bam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".sam -@ "+thread],shell=True) 
                subprocess.run(["samtools index  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".bam"],shell=True) 
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".sam ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".bam "],shell=True) 
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".sam "],shell=True)
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".fa"],shell=True)
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list"+part],shell=True)
            time2=timeit.default_timer()
            print("          -----Processing "+str(i-1)+"/"+str(file_line_num-1)+': %.0f Seconds'%(time2-time1),end="\r")
    print ()

#################################################################################################################################################################
#################################################################################################################################################################
#Start AS-ATI analysis.                        
if sys.argv[1] =="AS_ATI":
    print("Start AS-ATI analysis")
    print("Analysis in /output2_ASATI/")
    if "output2_ASATI" in os.listdir("./"):
        print("Note: ./output2_ASATI/ exist")
    else: 
        subprocess.run(["mkdir ./output2_ASATI"],shell=True)
    os.chdir("./output2_ASATI")
    os.system("pwd")
    print()
    with open ("../output0_preparation/7-ccs_inform/ASgene_ccsnum","r",encoding="ISO-8859-1") as f:
        gene_list=[col[0] for col in csv.reader(f,delimiter='\t')] 
    with open ("../output0_preparation/7-ccs_inform/ASgene_ccsnum","r",encoding="ISO-8859-1") as f:
        ccsnum_list=[col[1] for col in csv.reader(f,delimiter='\t')] 
    gene_num=len(gene_list)

    head_str="gene\tgene_reads.num\tAS\tAS_pos\tdSegment\tdSegmentlen\tTSS_zero\t"+	\
         "TSS_all_median\tAS_TSS_distance\tTSS1_median\tTSS2_median\t"+		\
         "ASform1.reads\tASform2.reads\tASform1.TSS\tASform2.TSS\t"+		\
         "ASform1.dTSS_raw\tASform2.dTSS_raw\tASform1.dTSS\tASform2.dTSS\t"+	\
         "ASform1.reads.num\tASform2.reads.num\t"+					\
         "read_usage\tKS_statistic\tp_value"
    with open ("AS2ATI_KS","w",encoding="utf-8") as f:
        f.write(head_str+"\n")
        f.close()
    i=0
    while i< gene_num:#< gene_num
        one_gene	=gene_list[i]
        gene_ccs_num	=ccsnum_list[i]
        i+=1
        print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene,end="\r")
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            transcript1_ASlist=[col[0] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            transcript1_ccslist=[col[1] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            transcript2_ASlist=[col[0] for col in csv.reader(f,delimiter='\t')] 
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            transcript2_ccslist=[col[1] for col in csv.reader(f,delimiter='\t')]
        transcript1_ASlist_len=len(transcript1_ASlist)
        transcript2_ASlist_len=len(transcript2_ASlist)
        transcript_ASlist_len=transcript1_ASlist_len+transcript2_ASlist_len
        gene_AS_list	=list(set(transcript1_ASlist))
        gene_AS_num	=len(gene_AS_list)
        ccs_list	=list(set(transcript1_ccslist+transcript2_ccslist))
        ccs_num		=len(ccs_list)
        #Get dict of AS infomation: dict[AS]=[oneAS_type/chr/strand/min/max/start/end/dSegment_pos/dSegment_length]
        dict_ASinfo={}
        for oneAS in gene_AS_list:
            oneAS_arr=oneAS.replace(":-",";minus").replace("-",";").replace(":",";").split(";")
            oneAS_type				=oneAS_arr[1]
            oneAS_chr				=oneAS_arr[2]
            oneAS_strand			=oneAS_arr[-1]
            oneAS_pos_arr			=list(map(int,oneAS_arr[3:-1]))   
            oneAS_min				=min(oneAS_pos_arr)
            oneAS_max				=max(oneAS_pos_arr) 
            if oneAS_type=="RI":		
                oneAS_start			=oneAS_pos_arr[1]
                oneAS_end			=oneAS_pos_arr[2]
            else:		
                oneAS_start			=oneAS_min
                oneAS_end			=oneAS_max
            if oneAS_type=="RI":		
                oneAS_dSegment_start		=oneAS_pos_arr[1]+1
                oneAS_dSegment_end		=oneAS_pos_arr[2]-1 
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)       
            if oneAS_type=="A3":
                if oneAS_strand=="+":	
                    oneAS_dSegment_start	=oneAS_pos_arr[1]
                    oneAS_dSegment_end		=oneAS_pos_arr[3]-1
                else:		
                    oneAS_dSegment_start	=oneAS_pos_arr[2]+1
                    oneAS_dSegment_end		=oneAS_pos_arr[0]	
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   
            if oneAS_type=="A5":
                if oneAS_strand=="+":
                    oneAS_dSegment_start	=oneAS_pos_arr[2]+1
                    oneAS_dSegment_end		=oneAS_pos_arr[0]
                else:		
                    oneAS_dSegment_start	=oneAS_pos_arr[1]
                    oneAS_dSegment_end		=oneAS_pos_arr[3]-1
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   	
            if oneAS_type=="SE":
                oneAS_dSegment_start		=oneAS_pos_arr[1]
                oneAS_dSegment_end		=oneAS_pos_arr[2]
                oneAS_dSegment_length		=oneAS_dSegment_end-oneAS_dSegment_start+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start)+"-"+str(oneAS_dSegment_end)   
            if oneAS_type=="MX":
                oneAS_dSegment_start1		=oneAS_pos_arr[1]
                oneAS_dSegment_end1		=oneAS_pos_arr[2]
                oneAS_dSegment_start2		=oneAS_pos_arr[5]
                oneAS_dSegment_end2		=oneAS_pos_arr[6]
                oneAS_dSegment_start		=oneAS_dSegment_start1
                oneAS_dSegment_end		=oneAS_dSegment_end2
                oneAS_dSegment_length		=oneAS_dSegment_end1-oneAS_dSegment_start1+1+oneAS_dSegment_end2-oneAS_dSegment_start2+1
                oneAS_dSegment_pos		=oneAS_chr+":"+str(oneAS_dSegment_start1)+"-"+str(oneAS_dSegment_end1)+";"+str(oneAS_dSegment_start2)+"-"+str(oneAS_dSegment_end2)
            dict_ASinfo[oneAS]=[oneAS_type,oneAS_chr,oneAS_strand,oneAS_min,oneAS_max,oneAS_start,oneAS_end,oneAS_dSegment_pos,oneAS_dSegment_length]
        ##Get of dict_ASccs[AS,ccs]=["0/1/2","map_start","map_end"]
        dict_ASccs={}
        for oneccs in ccs_list:
            for oneAS in gene_AS_list:
                dict_ASccs[oneAS,oneccs]=["0","","",""]
        k=0
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_1","r",encoding="ISO-8859-1") as f:
            for line in f.readlines():
                k+=1
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len),end="\r")
                eachline_arr=line.strip().split("\t")
                if eachline_arr[8] in ("TSS_PAS","TSS_noPAS"):
                    dict_ASccs[eachline_arr[0],eachline_arr[1]]=["1",eachline_arr[4],eachline_arr[5],eachline_arr[6]]
        with open ("../output0_preparation/7-ccs_inform/ASccs_split/"+one_gene+"_2","r",encoding="ISO-8859-1") as f:
            for line in f.readlines():
                k+=1
                print("          Processing "+str(i)+"/"+str(gene_num)+":\t"+one_gene+"\tGet_ASccs:"+str(k)+"/"+str(transcript_ASlist_len),end="\r")
                eachline_arr=line.strip().split("\t")
                if eachline_arr[8] in ("TSS_PAS","TSS_noPAS"):
                    dict_ASccs[eachline_arr[0],eachline_arr[1]]=["2",eachline_arr[4],eachline_arr[5],eachline_arr[6]]  
        print()
        for oneAS in gene_AS_list:
            oneAS_info		=dict_ASinfo[oneAS];
            oneAS_type		=oneAS_info[0]	
            oneAS_chr		=oneAS_info[1]	
            oneAS_strand	=oneAS_info[2]
            oneAS_start		=oneAS_info[5]
            oneAS_end		=oneAS_info[6] 
            oneAS_pos		=oneAS_chr+":"+str(oneAS_start)+"-"+str(oneAS_end)
            oneAS_dSegment	=oneAS_info[7] 
            oneAS_dSegment_len	=oneAS_info[8] 
            ccs1_arr=[];ccs2_arr=[]
            TSS1_arr=[];TSS2_arr=[]
            #Process each ccs
            for oneccs in ccs_list:
                oneccs_AS_info		=dict_ASccs[oneAS,oneccs]
                if oneccs_AS_info!=["0","","",""]:
                    oneccs_mapstart	=oneccs_AS_info[1]
                    oneccs_mapend	=oneccs_AS_info[2]
                    oneccs_TSS		=oneccs_AS_info[3]
                    if int(oneccs_mapstart)<=int(oneAS_start) and int(oneAS_end)<=int(oneccs_mapend):
                        if   oneccs_AS_info[0]=="1": ccs1_arr.append(oneccs);	TSS1_arr.append(oneccs_TSS)
                        elif oneccs_AS_info[0]=="2": ccs2_arr.append(oneccs);	TSS2_arr.append(oneccs_TSS)
            ccs1_num		=len(ccs1_arr)
            ccs2_num		=len(ccs2_arr)
            if ccs1_num>2 and ccs2_num>2:
                ccs1_arr_str	=",".join(ccs1_arr)
                ccs2_arr_str	=",".join(ccs2_arr)
                TSS1_arr_str	=",".join(TSS1_arr)
                TSS2_arr_str	=",".join(TSS2_arr)
                TSS_arr		=list(set(TSS1_arr+TSS2_arr));		
                TSS_int_arr	=list(map(int,TSS_arr));	TSS_int_arr.sort()
                TSS1_int_arr	=list(map(int,TSS1_arr));	TSS1_int_arr.sort();	
                TSS2_int_arr	=list(map(int,TSS2_arr));       TSS2_int_arr.sort();
                TSS1_median	=calcMedian(TSS1_int_arr)
                TSS2_median	=calcMedian(TSS2_int_arr)
                TSS_all_median	=calcMedian(TSS1_int_arr+TSS2_int_arr)
                AS_TSS_distance	=min(abs(TSS_all_median-oneAS_end),abs(TSS_all_median-oneAS_start))	
                #Get TSS_zero
                if oneAS_strand=="+": 	TSS_zero=min(TSS_int_arr)
                else:			TSS_zero=max(TSS_int_arr)
                #Standardization by TSS_zero
                dTSS1_int_arr=[] 
                dTSS2_int_arr=[]
                dTSS1_str_arr=[]
                dTSS2_str_arr=[]
                for x in TSS1_int_arr:	
                     dTSS1_int_arr.append(abs(x-TSS_zero))
                     dTSS1_str_arr.append(str(abs(x-TSS_zero)))
                for y in TSS2_int_arr:	
                     dTSS2_int_arr.append(abs(y-TSS_zero))
                     dTSS2_str_arr.append(str(abs(y-TSS_zero)))                 
                standard_TSS1_str	=",".join(dTSS1_str_arr)
                standard_TSS2_str	=",".join(dTSS2_str_arr)
                dTSS1_int_arr_noreplace=list(set(dTSS1_int_arr));	dTSS1_int_arr_noreplace.sort()
                dTSS2_int_arr_noreplace=list(set(dTSS2_int_arr));	dTSS2_int_arr_noreplace.sort()
                standard_TSS1_arr_good=[]
                standard_TSS2_arr_good=[]                    
                name_num_dict = {}
                for key in dTSS1_int_arr:
                    name_num_dict[key] 		= name_num_dict.get(key, 0) + 1
                for x in dTSS1_int_arr_noreplace:
                    each_dTSS_num			=name_num_dict[x]
                    each_dTSS			=str(x)+"("+str(each_dTSS_num)+")"
                    standard_TSS1_arr_good.append(each_dTSS)
                name_num_dict = {}
                for key in dTSS2_int_arr:
                    name_num_dict[key] = name_num_dict.get(key, 0) + 1
                for x in dTSS2_int_arr_noreplace:
                    each_dTSS_num			=name_num_dict[x]
                    each_dTSS			=str(x)+"("+str(each_dTSS_num)+")"
                    standard_TSS2_arr_good.append(each_dTSS)
                standard_TSS1_str_good		=",".join(standard_TSS1_arr_good)
                standard_TSS2_str_good		=",".join(standard_TSS2_arr_good)
                ccs_usage				=int(ccs1_num+ccs2_num)/int(gene_ccs_num)
                KS_str=str(ks_2samp(dTSS1_int_arr,dTSS2_int_arr))[13:][:-1]
                KS_str_arr=KS_str.split(", ")
                KS_statistic	=KS_str_arr[0][10:]
                pvalue		=KS_str_arr[1][7:]
                newline=one_gene+"\t"+gene_ccs_num+"\t"+oneAS+"\t"+oneAS_pos+"\t"+oneAS_dSegment+"\t"+str(oneAS_dSegment_len)+"\t"+str(TSS_zero)+"\t"+		\
                        str(TSS_all_median)+"\t"+str(AS_TSS_distance)+"\t"+str(TSS1_median)+"\t"+str(TSS2_median)+"\t"+			\
                        ccs1_arr_str+"\t"+ccs2_arr_str+"\t"+TSS1_arr_str+"\t"+TSS2_arr_str+"\t"+					\
                        standard_TSS1_str+"\t"+standard_TSS2_str+"\t"+standard_TSS1_str_good+"\t"+standard_TSS2_str_good+"\t"+		\
                        str(ccs1_num)+"\t"+str(ccs2_num)+"\t"+										\
                        str(ccs_usage)+"\t"+KS_statistic+"\t"+pvalue
                with open ("AS2ATI_KS","a",encoding="utf-8") as f2:
                     f2.write(newline+"\n")
                     f2.close()
    print("          Filter by pvalue<0.05, min_ccsnum and min_dSegmentlen")
    with open ("AS2ATI_KS.pvalue0.05","w",encoding="utf-8") as f2:
        f2.write(head_str+"\n")
        f2.close()
    i=0
    with open ("AS2ATI_KS","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            eachline		=line.strip()
            eachline_arr	=eachline.split("\t")
            dSegment_length	=int(eachline_arr[5])
            ccs1_num		=int(eachline_arr[19])
            ccs2_num		=int(eachline_arr[20])
            ccs_usage		=float(eachline_arr[21])
            KS_statistic	=float(eachline_arr[22])
            pvalue		=float(eachline_arr[23])
            if ccs1_num>int(min_ccsnum) and ccs2_num>int(min_ccsnum):
                if dSegment_length>int(min_dSegmentlen) and KS_statistic>=float(min_KS_statistic) and pvalue<0.05 and ccs_usage>=float(min_ccs_usage):
                    with open ("AS2ATI_KS.pvalue0.05","a",encoding="utf-8") as f2:
                        f2.write(eachline+"\n")
                        f2.close()                       
    subprocess.run(["sort -g -k 24 AS2ATI_KS.pvalue0.05 | cut -f 1-11,18-24   > AS2ATI_KS.pvalue0.05.simple"],shell=True)
    with open("AS2ATI_KS.pvalue0.05.simple","r",encoding="ISO-8859-1") as f:
        result_num=len(f.readlines())-1
    print()
    print("     Complete!")
    print("     Possible results number is "+str(result_num))
    print()
    print("     part_ccs2ref. This section gives the BAM files which can be used as reference for credibility.")
    print("          All ccs were split into two parts: ASform1 + ASform2")
    print("          Get ccs(pvalue0.05) list from AS2ATI_KS.pvalue0.05(simple)")
    print("          Get fasta from transcript1 and 2 in each AS and then minimap2ref")   
    if "part_ccs2ref" in os.listdir("./"):
        print("          Note: ./part_ccs2ref/ exist and will be deleted.") 
        subprocess.run(["rm -r ./part_ccs2ref"],shell=True) 
    subprocess.run(["mkdir ./part_ccs2ref"],shell=True)
    i=0
    with open("AS2ATI_KS.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            ccs1_str=eachline_arr[11]
            ccs2_str=eachline_arr[12]
            ccs_list=ccs1_str.split(",")+ccs2_str.split(",")
            for x in ccs_list:
                with open("./part_ccs2ref/ccs.list0","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()
    subprocess.run(["sort -n ./part_ccs2ref/ccs.list0 | uniq > ./part_ccs2ref/1-ccs.pvalue0.05.list"],shell=True)
    subprocess.run(["rm ./part_ccs2ref/ccs.list0"],shell=True)
    print("          Extract fasta sequence")
    subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/1-ccs.pvalue0.05.list ../output0_preparation/3-all_FLNC/all_FLNC_nopolyA.fa > ./part_ccs2ref/2-ccs.pvalue0.05.fa"],shell=True)   
    i=0
    if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"):
        subprocess.run(["rm -r ./part_ccs2ref/ccs_pvalue0.05_eachAS/"],shell=True)  
    subprocess.run(["mkdir ./part_ccs2ref/ccs_pvalue0.05_eachAS"],shell=True)
    i=0
    with open("AS2ATI_KS.pvalue0.05","r",encoding="ISO-8859-1") as f:
        file_line_num=len(f.readlines())
    with open("AS2ATI_KS.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            time1=timeit.default_timer()
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            eventid1	=eachline_arr[2]
            eventid2	=eachline_arr[2].replace(";","\;")
            ccs1_str	=eachline_arr[11]
            ccs2_str	=eachline_arr[12]
            ccs1_list=ccs1_str.split(",")
            ccs2_list=ccs2_str.split(",")
            if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list1"):subprocess.run(["rm  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list1"],shell=True)  
            if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list2"):subprocess.run(["rm  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list2"],shell=True)  
            for x in ccs1_list:
                with open("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid1+".list1","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()            
            for x in ccs2_list:
                with open("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid1+".list2","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()  
            part_arr=["1","2"]
            for part in part_arr:
                subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list"+part+" ./part_ccs2ref/2-ccs.pvalue0.05.fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".fa"],shell=True)
                subprocess.run(["minimap2 -ax splice -uf -k 14 -t "+thread+" --secondary=no ../output0_preparation/4-all_FLNC_minimap2ref/ref.fa ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".sam 2>/dev/null"],shell=True) 
                subprocess.run(["samtools view -bS ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".sam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".bam -@ "+thread],shell=True) 
                subprocess.run(["samtools sort ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".bam -@ "+thread+" -o ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".bam 1>/dev/null 2>&1"],shell=True) 
                subprocess.run(["samtools view -h ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".bam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".sam -@ "+thread],shell=True) 
                subprocess.run(["samtools index  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".bam"],shell=True) 
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".sam ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".bam "],shell=True) 
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"sort_"+part+".sam "],shell=True)
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+"_"+part+".fa"],shell=True)
                subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+eventid2+".list"+part],shell=True)
            time2=timeit.default_timer()
            print("          -----Processing "+str(i-1)+"/"+str(file_line_num-1)+': %.0f Seconds'%(time2-time1),end="\r")
    print ()
#################################################################################################################################################################
#################################################################################################################################################################
if sys.argv[1] =="ATI_APA":
    print("Start ATI_APA analysis")
    print("Analysis in /output4_ATIAPA/")
    if "output4_ATIAPA" in os.listdir("./"):
        print("Note: ./output4_ATIAPA/ exist")
    else: 
        subprocess.run(["mkdir ./output4_ATIAPA"],shell=True)
    os.chdir("./output4_ATIAPA")
    os.system("pwd")
    #Get dict_gene_info
    print("     Get dict_gene_info")
    dict_gene_info={}
    with open ("../output0_preparation/5-cDNA_cupcake/gene_info","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            gene		=eachline_arr[0]
            chromosome		=eachline_arr[1]
            strand		=eachline_arr[2]
            start		=eachline_arr[3]
            end			=eachline_arr[4]
            dict_gene_info[gene]=[chromosome,strand,start,end]
    #Get dict_ccs_info
    print("     Get dict_ccs_info")
    dict_ccs_info={}
    i=0
    with open ("../output0_preparation/4-all_FLNC_minimap2ref/FLNC_inform.uniq","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1:continue
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            ccs_name		=eachline_arr[0]     
            align_chr       	=eachline_arr[1] 
            strand  		=eachline_arr[2] 
            align_start     	=eachline_arr[3] 
            align_end       	=eachline_arr[4] 
            TSS     		=eachline_arr[5] 
            PAS     		=eachline_arr[6] 
            length		=abs(int(PAS)-int(TSS))
            TSS_PAS_mark	=eachline_arr[7] 
            intron_mark 	=eachline_arr[8]
            if TSS_PAS_mark=="TSS_PAS" and intron_mark=="normal":
                dict_ccs_info[ccs_name]=[align_chr,strand,align_start,align_end,TSS,PAS,length]
    #Start analysis
    print("     Start spearman analysis")
    head_str="gene\tgene_pos\tgene_reads.num\tTSSPAS_reads.num\t"+			\
         "subclass\tsubclass_pos\tsubclass_reads\tsubclass_reads.num\t"+		\
         "TSS.raw\tPAS.raw\t"+							\
         "TSS_zero\tPAS_zero\tdTSS\tdPAS\t"+					\
         "generead_usage\tTSSPASread_usage\tspearman_correlation\tp_value"
    with open ("ATI2APA_spearman","w",encoding="utf-8") as f:
        f.write(head_str+"\n")
        f.close()
    with open ("../output0_preparation/5-cDNA_cupcake/gene_transcript_num_ccs_num","r",encoding="ISO-8859-1") as f:
        gene_num=len(f.readlines())-1
    dict_gene_ccs={}
    i=0
    with open("../output0_preparation/5-cDNA_cupcake/gene_transcript_num_ccs_num","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1: continue     
            print("          Processing:\t"+str(i-1)+"/"+str(gene_num),end="\r") 
            eachline=line.strip()
            eachline_arr		=eachline.split("\t")   
            gene_name			=eachline_arr[0]
            gene_info			=dict_gene_info[gene_name]			
            gene_pos			=gene_info[0]+":"+gene_info[2]+"-"+gene_info[3]+"("+gene_info[1]+")"
            gene_ccs_num		=eachline_arr[4]
            ccs_arr			=eachline_arr[3].split(",")
            list_ccs_ccsinfo=[]
            for oneccs in ccs_arr:
                if oneccs in dict_ccs_info:
                    oneccs_info		=dict_ccs_info[oneccs]
                    oneccs_info.insert(0,oneccs)
                    list_ccs_ccsinfo.append(oneccs_info)
            list_ccs_ccsinfo.sort(key=lambda x: x[7],reverse=True)
            TSSPAS_ccs_num 		=len(list_ccs_ccsinfo)
            k=0;index_arr=[];dict_subclass={};dict_subclass_info={}
            for oneccs_info_arr in list_ccs_ccsinfo:
                oneccs_name		=oneccs_info_arr[0]
                onegene_strand		=oneccs_info_arr[2]
                oneccs_start		=int(oneccs_info_arr[3])
                oneccs_end		=int(oneccs_info_arr[4])
                oneccs_TSS		=oneccs_info_arr[5]
                oneccs_PAS		=oneccs_info_arr[6]
                if dict_subclass=={}:
                    k+=1;index_arr.append(k)
                    dict_subclass[k]		=[oneccs_info_arr]
                    dict_subclass_info[k]		=[int(oneccs_start),int(oneccs_start),int(oneccs_end),int(oneccs_end)]
                else:
                    split_mark="NO"
                    for oneindex in index_arr:     
                        subclass_info		=dict_subclass_info[oneindex]
                        subclass_info_S1		=subclass_info[0]
                        subclass_info_S2		=subclass_info[1]
                        subclass_info_E1		=subclass_info[2]
                        subclass_info_E2		=subclass_info[3]
                        if (abs(oneccs_start-subclass_info_S1)<int(max_bin_extent) or abs(oneccs_start-subclass_info_S2)<int(max_bin_extent)):
                            if (abs(oneccs_end-subclass_info_E1)<int(max_bin_extent) or abs(oneccs_end-subclass_info_E2)<int(max_bin_extent)):
                                if  oneccs_start<subclass_info_E1 and oneccs_end>subclass_info_S2:
                                    if oneccs_start<subclass_info_S1:		dict_subclass_info[oneindex][0]=oneccs_start
                                    if oneccs_start>subclass_info_S2:		dict_subclass_info[oneindex][1]=oneccs_start
                                    if oneccs_end<subclass_info_E1:   		dict_subclass_info[oneindex][2]=oneccs_end
                                    if oneccs_end>subclass_info_E2:		dict_subclass_info[oneindex][3]=oneccs_end
                                    dict_subclass[oneindex].append(oneccs_info_arr)
                                    split_mark="YES";break
                    if split_mark=="NO":
                        k+=1;index_arr.append(k)
                        dict_subclass[k]			=[oneccs_info_arr]
                        dict_subclass_info[k]		=[int(oneccs_start),int(oneccs_start),int(oneccs_end),int(oneccs_end)]
            for oneindex in index_arr:
                subclass_info		=dict_subclass_info[oneindex]
                subclass_info_str	=str(subclass_info[0])+"-"+str(subclass_info[1])+"-"+str(subclass_info[2])+"-"+str(subclass_info[3])
                subclass_pos	=gene_info[0]+":"+str(subclass_info[0])+"-"+str(subclass_info[3])
                subclass_content	=dict_subclass[oneindex]
                ccs_arr=[];TSS_arr=[];PAS_arr=[];TSS_int_arr=[];PAS_int_arr=[]
                for oneccs_arr in subclass_content:
                    #print(subclass_content)
                    ccs_arr.append(oneccs_arr[0])
                    TSS_arr.append(oneccs_arr[5]);	TSS_int_arr.append(int(oneccs_arr[5]))	
                    PAS_arr.append(oneccs_arr[6]);	PAS_int_arr.append(int(oneccs_arr[6]))	
                ccs_arr_str=",".join(ccs_arr)
                TSS_arr_str=",".join(TSS_arr)
                PAS_arr_str=",".join(PAS_arr)
                goodccs_num=len(ccs_arr)
                if goodccs_num > 2:
                    geneccs_usage	=str(goodccs_num/int(gene_ccs_num))
                    TSSPASccs_usage	=str(goodccs_num/int(TSSPAS_ccs_num))
                    if onegene_strand=="+": 	TSS_zero=min(TSS_int_arr);	PAS_zero=max(PAS_int_arr)
                    else:			TSS_zero=max(TSS_int_arr);	PAS_zero=min(PAS_int_arr)
                    dTSS_str_arr=[] 
                    dPAS_str_arr=[]
                    for x in TSS_int_arr:	
                        dTSS_str_arr.append(str(abs(x-TSS_zero)))
                    for x in PAS_int_arr:
                        dPAS_str_arr.append(str(abs(x-PAS_zero)))
                    dTSS_arr_str=",".join(dTSS_str_arr)
                    dPAS_arr_str=",".join(dPAS_str_arr) 
                    if np.var(TSS_int_arr)>0 and np.var(PAS_int_arr)>0 :        
                        spearman_correlation,pvalue=spearmanr(TSS_int_arr,PAS_int_arr)
                        #Start analysis   
                        with open("ATI2APA_spearman","a",encoding="utf-8") as f2:
                            newline=gene_name+"\t"+gene_pos+"\t"+str(gene_ccs_num)+"\t"+str(TSSPAS_ccs_num)+"\t"+			\
                                subclass_info_str+"\t"+subclass_pos+"\t"+ccs_arr_str+"\t"+str(goodccs_num)+"\t"+				\
                                TSS_arr_str+"\t"+PAS_arr_str+"\t"+									\
                                str(TSS_zero)+"\t"+str(PAS_zero)+"\t"+dTSS_arr_str+"\t"+dPAS_arr_str+"\t"+				\
                                geneccs_usage+"\t"+TSSPASccs_usage+"\t"+str(spearman_correlation)+"\t"+str(pvalue) 
                            f2.write(newline+"\n")
                            f2.close()                                           
    print()
    print("          Filter by pvalue<0.05, min_ccsnum and min_correlation")
    with open ("ATI2APA_spearman.pvalue0.05","w",encoding="utf-8") as f2:
        f2.write(head_str+"\n")
        f2.close()
    i=0
    with open ("ATI2APA_spearman","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            eachline		=line.strip()
            eachline_arr	=eachline.split("\t")
            ccs_num		=int(eachline_arr[7])
            geneccs_usage       =float(eachline_arr[14])
            TSSPASccs_usage     =float(eachline_arr[15])    
            spearman_correlation=float(eachline_arr[16])
            pvalue		=float(eachline_arr[17])
            if ccs_num>int(min_ccsnum) and abs(spearman_correlation)>float(min_correlation) and pvalue<0.05 and geneccs_usage>=float(min_geneccs_usage) and TSSPASccs_usage>=float(min_TSSPASccs_usage):
                with open ("ATI2APA_spearman.pvalue0.05","a",encoding="utf-8") as f2:
                    f2.write(eachline+"\n")
                    f2.close()         
                 
    subprocess.run(["sort -g -k 18 ATI2APA_spearman.pvalue0.05 | cut -f 1-6,8,11-12,15-18   > ATI2APA_spearman.pvalue0.05.simple"],shell=True)
    with open("ATI2APA_spearman.pvalue0.05.simple","r",encoding="ISO-8859-1") as f:
        result_num=len(f.readlines())-1
    print()
    print("     Complete!")
    print("     Possible results number is "+str(result_num))
    print()
    print("     part_ccs2ref. This section gives the BAM files which can be used as reference for credibility.")
    print("          Get ccs(pvalue0.05) list from ATI2APA_spearman.pvalue0.05(simple)")
    print("          Get ccs.fasta in each gene and then minimap2ref")  

    if "part_ccs2ref" in os.listdir("./"):
        print("          Note: ./part_ccs2ref/ exist and will be deleted.") 
        subprocess.run(["rm -r ./part_ccs2ref"],shell=True) 
    subprocess.run(["mkdir ./part_ccs2ref"],shell=True)
    i=0
    with open("ATI2APA_spearman.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            ccs_str=eachline_arr[6]
            ccs_list=ccs_str.split(",")
            for x in ccs_list:
                with open("./part_ccs2ref/ccs.list0","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()
    subprocess.run(["sort -n ./part_ccs2ref/ccs.list0 | uniq > ./part_ccs2ref/1-ccs.pvalue0.05.list"],shell=True)
    subprocess.run(["rm ./part_ccs2ref/ccs.list0"],shell=True)
    print("          Extract fasta sequence")
    subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/1-ccs.pvalue0.05.list ../output0_preparation/3-all_FLNC/all_FLNC_nopolyA.fa > ./part_ccs2ref/2-ccs.pvalue0.05.fa"],shell=True)   
    i=0
    if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"):
        subprocess.run(["rm -r ./part_ccs2ref/ccs_pvalue0.05_eachAS/"],shell=True)  
    subprocess.run(["mkdir ./part_ccs2ref/ccs_pvalue0.05_eachAS"],shell=True)
    i=0
    with open("ATI2APA_spearman.pvalue0.05","r",encoding="ISO-8859-1") as f:
        file_line_num=len(f.readlines())
    with open("ATI2APA_spearman.pvalue0.05","r",encoding="ISO-8859-1") as f:
        for line in f.readlines():
            i+=1
            if i==1 : continue
            time1=timeit.default_timer()
            eachline=line.strip()
            eachline_arr=eachline.split("\t")
            geneid	=eachline_arr[0]+"_"+eachline_arr[4]
            ccs_str	=eachline_arr[6]
            ccs_list=ccs_str.split(",")
            if os.path.exists("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".list"):subprocess.run(["rm  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".list"],shell=True)             
            for x in ccs_list:
                with open("./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".list","a",encoding="utf-8") as f2:
                    f2.write(x+"\n")   
                    f2.close()            
            subprocess.run(["perl ../output0_preparation/script/getfastabylist.pl  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".list ./part_ccs2ref/2-ccs.pvalue0.05.fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".fa"],shell=True)
            subprocess.run(["minimap2 -ax splice -uf -k 14 -t "+thread+" --secondary=no ../output0_preparation/4-all_FLNC_minimap2ref/ref.fa ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".fa > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".sam 2>/dev/null"],shell=True) 
            subprocess.run(["samtools view -bS ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".sam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".bam -@ "+thread],shell=True) 
            subprocess.run(["samtools sort ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".bam -@ "+thread+" -o ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+"sort.bam 1>/dev/null 2>&1"],shell=True) 
            subprocess.run(["samtools view -h ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+"sort.bam > ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+"sort.sam -@ "+thread],shell=True) 
            subprocess.run(["samtools index  ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+"sort.bam"],shell=True) 
            subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".sam ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".bam "],shell=True) 
            subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+"sort.sam "],shell=True)
            subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".fa"],shell=True)
            subprocess.run(["rm ./part_ccs2ref/ccs_pvalue0.05_eachAS/"+geneid+".list"],shell=True)
            time2=timeit.default_timer()
            print("          -----Processing "+str(i-1)+"/"+str(file_line_num-1)+': %.0f Seconds'%(time2-time1),end="\r")
    print ()
print()
#################################################################################################################################################################
#################################################################################################################################################################
print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))      
time_end=timeit.default_timer()
print('All the running time: %.0f Seconds'%(time_end-time_start))












