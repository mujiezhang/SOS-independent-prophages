# SOS-independent-prophages
This is the script repository for the following manuscript:

Yali Hao#, Mujie Zhang#, Xinjuan Lei, Chengrui Zhu, Xiang Xiao, Huahua Jian*, SOS-independent prophages prevail in the bacterial genomes. xxx (2025)

# Documents

## For analysis:

### 1. Prophage categorization
   
- We use the pipeline we developed, PSOSP, to classify prophages into SOS-dependent prophages (SdP), SOS-uncertain prophages (SuPs), and SOS-independent prophages (SiPs).
- The PSOSP scripts and usage instructions are available at https://github.com/mujiezhang/PSOSP.
- The PSOSP web server is available at: https://vee-lab.sjtu.edu.cn/PSOSP/index.html.
   
### 2. vOTU clustering based on ANI
   
- We clustered vOTUs using the [CheckV pipeline](https://bitbucket.org/berkeleylab/checkv/src/master/), based all-versus-all BLASTn search and Leiden algorithm，following MIUViG guidelines (95% average nucleotide identity (ANI); 85% aligned fraction (AF)
  - step1: all-vs-all blastn
    ```
    makeblastdb -in all_virus.fna -dbtype nucl -out all_virus
    blastn -query all_virus.fna -db all_virus -outfmt '6 std qlen slen' -max_target_seqs 100000 -out my_blast.tsv -num_threads 64 -task megablast -evalue 1e-5
    ```
  - step2: calculate ANI using script `anicalc.py`
    ```
    python anicalc.py -i my_blast.tsv -o my_ani.tsv
    ```
  - step3: vOTU clustering using script `aniclust.py`
    ```
    python aniclust.py --fna all_virus.fna --ani my_ani.tsv --out my_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0
    ```
     
### 3.  Genus and Family level clustering based on AAI

- We performed genus/family clustering using the [MGV pipeline](https://github.com/snayfach/MGV/tree/master/aai_cluster) based on all-vs-all BLASTp search and MCL
  - step1: all-vs-all blastp
    ```
    prodigal -a all_votu.faa  -i all_otu.fna   -p meta
    diamond makedb --in all_votu.faa --db viral_proteins --threads 10
    diamond blastp --query all_votu.faa --db viral_proteins.dmnd --out blastp.tsv --outfmt 6 --evalue 1e-5 --max-target-seqs 1000000 --query-cover 50 --subject-cover 50
    ```
  - step2: calculate AAI
    ```
    python amino_acid_identity.py --in_faa query all_votu.faa --in_blast blastp.tsv --out_tsv aai.tsv
    ```
    Note: Modified script `amino_acid_identity.py` for Python3 compatibility: line21:`print "parse"`→`print("parse")`; line38:`print "compute"`→`print("compute")`; line52:`print "write"`→`print("write")`
  - step3: Filter edges and prepare MCL input
    ```
    python filter_aai.py --in_aai aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 50 --out_tsv genus_edges.tsv
    python filter_aai.py --in_aai aai.tsv --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv family_edges.tsv
    ```
  - step4: Genus and family level clustering based on MCL
    ```
    mcl genus_edges.tsv -te 8 -I 2.0 --abc -o genus_clusters.txt
    mcl family_edges.tsv -te 8 -I 1.2 --abc -o family_clusters.txt
    ```
    Note: Adjusted genus filtering to `--min_aai 50` following the parameters in their [paper](https://www.nature.com/articles/s41564-021-00928-6)

### 4. Protein sharing network analysis of viral populations was performed by vConTACT2
   ```
   vcontact2_gene2genome -p all_votu_and_ICTV_phages.faa -o gene2genome.csv -s Prodigal-FAA
   vcontact2 -r all_votu_and_ICTV_phages.faa -p gene2genome.csv --db None -o votu_ICTV_vcontact2 -t 64
   ```

### 5. Calculate the weighted Gene Repertoire Relatedness (wGRR)
- Calculated viral wGRR using the formula from [J. A. M. de Sousa et al. (2023)](https://academic.oup.com/nar/article/51/6/2759/7068371?login=true)
   ```
   diamond blastp --threads 50 --db viral_proteins.dmnd --out diamond.tsv --evalue 0.0001 --max-target-seqs 100000 --query  all_votu.faa --id  35 --query-cover 50 --subject-cover 50
   python calculate_wGRR.py diamond.tsv all_votu.faa result_file
   ```
### 6. Enrichment analysis
   
- We use the script `enrichment_analysis.py` to perform environment enrichment analysis. The fold enrichment was calculated by dividing the proportion of SdPs from a specific environment by the proportion of that environment within all environments.

## For Figures
- [Figure1](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Figure%201):
   `Figure1-c-violin-plot.R`绘制Figure 1 c小提琴图，`ecoli-k12.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt`包含Ecoli K12中所有PBSs的HI, `ecoli-k12_point.tsv`包含24个实验验证的Ecoli K12的PBSs.

- [Figure2](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Figure%202):
   `Figure2-b-bac_cut_genus_30.tree`是Figture 2 b的进化树数据，`Figure2-b-annotation_for_tree.xlsx`是进化树的注释数据，进化树的可视化由[TVBOT](https://www.chiplot.online/tvbot.html)实现. `Figure2-e-violin-plot.R`绘制Figure2 e小提琴图，对应数据文件为``. Figure2 c网络数据文件超过github文件大小限制，如有需要，可通过邮件发送。

- [Supplemental Figure S1](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S1):

- [Supplemental Figure S2](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S2)

- [Supplemental Figure S3](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S3)
  
- [Supplemental Figure S4](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S4)
  
- [Supplemental Figure S5](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S5)
  
- [Supplemental Figure S6](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S6)
  
- [Supplemental Figure S7](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S7)
  
- [Supplemental Figure S8](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S8)
  
- [Supplemental Figure S9](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S9)
  
- [Supplemental Figure S10](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S10)
  
- [Supplemental Figure S11](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S11)
  
- [Supplemental Figure S12](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S12)
  
- [Supplemental Figure S13](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S13)
  
- [Supplemental Figure S14](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S14)
  
- [Supplemental Figure S15](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S15)
  

  


