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



#Sequencing primers of IsoSeq should be specified correctly. Two common primers are provided.
primer1_5p：AAGCAGTGGTATCAACGCAGAGTACATGGGG     primer1_3p：AAGCAGTGGTATCAACGCAGAGTAC
primer2_5p：ATGTAATACGACTCACTATAGGGC            primer2_3p：CGCCTGAGA
Edit the script file if the primer needs to be changed.




#Sequencing primers of IsoSeq should be specified correctly. Two common primers are provided.
primer1_5p：AAGCAGTGGTATCAACGCAGAGTACATGGGG     primer1_3p：AAGCAGTGGTATCAACGCAGAGTAC
primer2_5p：ATGTAATACGACTCACTATAGGGC            primer2_3p：CGCCTGAGA
Edit the script file if the primer needs to be changed.
