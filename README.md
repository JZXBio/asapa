ASAPA
====
python asapa.py-----help:<br>
Version = 1.0.0

Preparation: the subreads files were processed using pbCCS, Lima, and SUPPA2
Usage: **python asapa.py build ref.fa  subreads.bam/qry_dir(contain *.bam)/all(current dir)**
```
    Optional parameters:
        -n                      default=15, CPU thread number
        -log                    default=no, yes will write terminal print in output.log
        -max_fuzzy_TSS          default=5,  Max fuzzy TSS dist
        -max_fuzzy_PAS          default=5,  Max fuzzy PAS dist
        -max_fuzzy_junction     default=5,  Max fuzzy junction dist(from cDNA_cupcake)
```
Function1: AS vs AS
    Usage: python asapa.py AS_AS 
    Optional parameters: 
        -min_ccsnum             default=10, min ccs number in AS1(2)form1(2)
        -min_dSegmentlen        default=10, min AS differential segment length
        -min_ccs_usage          default=0, min ccs usage(used/geneccs)

Function2: AS vs ATI 
    Usage: python asapa.py AS_ATI
        -min_ccsnum             default=10, min ccs number in ASform1(2)
        -min_dSegmentlen        default=10, Min AS differential segment length
        -min_ccs_usage          default=0, min ccs usage(used/geneccs)
        -min_KS_statistic       default=0.2, Min KS_statistic

Function3: AS vs APA 
    Usage: python asapa.py AS_APA
    Optional parameters:
        -min_ccsnum             default=10, min ccs number in ASform1(2)
        -min_dSegmentlen        default=10, Min AS differential segment length
        -min_ccs_usage          default=0, min ccs usage(used/geneccs)
        -min_KS_statistic       default=0.2, Min KS_statistic

Function4: ATI vs APA
    Usage: python asapa.py ATI_APA
    Optional parameters:
        -min_ccsnum             default=10, min ccs number of gene/transcript
        -min_geneccs_usage      default=0, min ccs usage(used/geneccs)
        -min_TSSPASccs_usage    default=0.5, min ccs usage(used/TSSPASccs)
        -min_correlation        default=0.5, min spearman correlation
        -max_bin_extent         default=1000,  Max distance of binning extension

The output folder will be created in the current path:
    output0_preparation (preparation: subreads to ccs, lima, minimap2, cDNA_cupcake and SUPPA2)
    output1_ASAS    (function1: coupling bewteen AS and AS)
    output2_ASATI    (function2: coupling bewteen AS and ATI)
    output3_ASAPA    (function3: coupling bewteen AS and APA)
    output4_ATIAPA    (function4: coupling bewteen ATI and APA)

Dependency:
Conda is recommended
     conda install pbbam
     conda install pbccs==6.4
     conda install lima
     conda install minimap2
     conda install samtools
     conda install bedtools
     conda install perl
     conda install R
     conda install rpy2
     conda install r-dplyr
     conda install -c bioconda trim_isoseq_polya
     conda install -c bioconda suppa

Independent installation required
     cDNA_Cupcake

<br>
#Sequencing primers of IsoSeq should be specified correctly. Two common primers are provided.<br>
primer1_5p：AAGCAGTGGTATCAACGCAGAGTACATGGGG     primer1_3p：AAGCAGTGGTATCAACGCAGAGTAC<br>
primer2_5p：ATGTAATACGACTCACTATAGGGC            primer2_3p：CGCCTGAGA<br>
Edit the script file if the primer needs to be changed.<br>
<br>
