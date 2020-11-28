# WholeExomeSequencing
Takes BAM files and returns SNVs and CNVs. Additional scripts allow for analysis and plotting. 

All shell scripts are designed for use on the Stanford Computing Cluster, and therefore have some dependencies that are not universally accessible, however most dependencies listed in the main shell script are publically available. 

The data used in these examples is used by in-progress research, and will be released upon publication. These graphs are just examples of some of the analysis capabilities of my scripts. These examples come from paired tumor/normal tissue biopsies. 

## Example Outputs

### VEP Mutation Impact Graph

![VEP](ExamplePlots/VEP_screenshot.png)

### Copy Number Variation of Cosmic Genes 

![CNV](ExamplePlots/cnv_screenshot.png)

The following are results from a similar pipeline which uses multiple variant callers using proprietary scripts that are not mine. I just wrote the analysis for this.

### Consensus Calling Results from ExomePlottingHuman.R

![Consensus](ExamplePlots/human_screenshot.png)

### Allele Frequency Distribution

![AF_freq](ExamplePlots/af_freq_screenshot.png)

Big thanks to jhchabon for providing critical dependency scripts and the overall framework design of the MousePipelineDeduped.sh Script
