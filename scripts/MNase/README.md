## Content
+ 01_script_pVII_mapping.sh: bash script to align and preprocess MNase-seq reads 
+ 02_script_pVII_positioning: bash script to generate pVII occupancy maps and call pVII positions
+ 03_accessibility_score
  + 01_bamSieve.sh: bash script to filter low MNase reads for fragment length < 140bp
  + 02_runNF_lowMNase.sh: bash script to run nextflow pipeline *nf_low_MNase*, which calculates accessibility scores over time
  + nextflow: contains nextflow pipeline *nf_low_MNase*
+ 04_nucleosome_calling
  + 01_prepare_input.sh: bash script to prepare the input for nucleR analysis
  + 02_nucleosomeCalling.R: R script to call sites of nucleosome assembly using nucleR
  + 03_filter_nucleosomes.R: R script to filter sites of nucleosome assembly by H3.3 abundance at 4 hpi

