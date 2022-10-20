python asapa.py-----help:<br>
Version = 1.0.0<br>
<br>
Preparation: the subreads files were processed using pbCCS, Lima, and SUPPA2<br>
Usage: python asapa.py build ref.fa  subreads.bam/qry_dir(contain *.bam)/all(current dir)<br>
    Optional parameters:<br>
        -n                       default=15, CPU thread number<br>
        -log                     default=no, yes will write terminal print in output.log<br>
        -max_fuzzy_TSS           default=5,  Max fuzzy TSS dist<br>
        -max_fuzzy_PAS           default=5,  Max fuzzy PAS dist<br>
        -max_fuzzy_junction      default=5,  Max fuzzy junction dist(from cDNA_cupcake)<br>
<br>
Function1: AS vs AS<br>
    Usage: python asapa.py AS_AS <br>
    Optional parameters: <br>
        -min_ccsnum              default=10, min ccs number in AS1(2)form1(2)<br>
        -min_dSegmentlen         default=10, min AS differential segment length<br>
        -min_ccs_usage           default=0, min ccs usage(used/geneccs)<br>
<br>
Function2: AS vs ATI <br>
    Usage: python asapa.py AS_ATI<br>
        -min_ccsnum              default=10, min ccs number in ASform1(2)<br>
        -min_dSegmentlen         default=10, Min AS differential segment length<br>
        -min_ccs_usage           default=0, min ccs usage(used/geneccs)<br>
        -min_KS_statistic        default=0.2, Min KS_statistic<br>
<br>
Function3: AS vs APA <br>
    Usage: python asapa.py AS_APA<br>
    Optional parameters:<br>
        -min_ccsnum              default=10, min ccs number in ASform1(2)<br>
        -min_dSegmentlen         default=10, Min AS differential segment length<br>
        -min_ccs_usage           default=0, min ccs usage(used/geneccs)<br>
        -min_KS_statistic        default=0.2, Min KS_statistic<br>
<br>
Function4: ATI vs APA<br>
    Usage: python asapa.py ATI_APA<br>
    Optional parameters:<br>
        -min_ccsnum              default=10, min ccs number of gene/transcript<br>
        -min_geneccs_usage       default=0, min ccs usage(used/geneccs)<br>
        -min_TSSPASccs_usage     default=0.5, min ccs usage(used/TSSPASccs)<br>
        -min_correlation         default=0.5, min spearman correlation<br>
        -max_bin_extent          default=1000,  Max distance of binning extension<br>
<br>
The output folder will be created in the current path:<br>
    output0_preparation (preparation: subreads to ccs, lima, minimap2, cDNA_cupcake and SUPPA2)<br>
    output1_ASAS              (function1: coupling bewteen AS and AS)<br>
    output2_ASATI             (function2: coupling bewteen AS and ATI)<br>
    output3_ASAPA             (function3: coupling bewteen AS and APA)<br>
    output4_ATIAPA            (function4: coupling bewteen ATI and APA)<br>
<br>
Dependency:<br>
Conda is recommended<br>
    conda install pbbam<br>
    conda install pbccs==6.4<br>
    conda install lima<br>
    conda install minimap2<br>
    conda install samtools<br>
    conda install bedtools<br>
    conda install perl<br>
    conda install R<br>
    conda install rpy2<br>
    conda install r-dplyr<br>
    conda install -c bioconda trim_isoseq_polya<br>
    conda install -c bioconda suppa<br>
<br>
Independent installation required<br>
    cDNA_Cupcake<br>
<br>

<br>
#Sequencing primers of IsoSeq should be specified correctly. Two common primers are provided.<br>
primer1_5p：AAGCAGTGGTATCAACGCAGAGTACATGGGG     primer1_3p：AAGCAGTGGTATCAACGCAGAGTAC<br>
primer2_5p：ATGTAATACGACTCACTATAGGGC            primer2_3p：CGCCTGAGA<br>
Edit the script file if the primer needs to be changed.<br>
<br>
