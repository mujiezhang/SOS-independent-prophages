# SOS-independent-prophages
This is the script repository for the following manuscript:

Yali Hao#, Mujie Zhang#, Xinjuan Lei, Chengrui Zhu, Xiang Xiao, Huahua Jian*, SOS-independent prophages prevail in the bacterial genomes. xxx (2025)

# Documents

## For analysis:

### 1. Prophage categorization
   
- We use the pipeline we developed, PSOSP, to classify prophages into SOS-dependent prophages (SdP), SOS-uncertain prophages (SuPs), and SOS-independent prophages (SiPs).
- The PSOSP scripts and usage instructions are available at **https://github.com/mujiezhang/PSOSP**.
- The PSOSP web server is available at: **https://vee-lab.sjtu.edu.cn/PSOSP/index.html**.
   
### 2. vOTU clustering based on ANI
   
- We clustered vOTUs using the [**CheckV pipeline**](https://bitbucket.org/berkeleylab/checkv/src/master/), based all-versus-all BLASTn search and Leiden algorithm，following MIUViG guidelines (95% average nucleotide identity (ANI); 85% aligned fraction (AF)
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

- We performed genus/family clustering using the [**MGV pipeline**](https://github.com/snayfach/MGV/tree/master/aai_cluster) based on all-vs-all BLASTp search and MCL
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
    Note: Adjusted genus filtering to `--min_aai 50` following the parameters in their [**paper**](https://www.nature.com/articles/s41564-021-00928-6)

### 4. Protein sharing network analysis of viral populations was performed by vConTACT2
   ```
   # prepare input
   vcontact2_gene2genome -p all_votu_and_ICTV_phages.faa -o gene2genome.csv -s Prodigal-FAA
   # run vcontact2
   vcontact2 -r all_votu_and_ICTV_phages.faa -p gene2genome.csv --db None -o votu_ICTV_vcontact2 -t 64
   ```

### 5. Calculate the weighted Gene Repertoire Relatedness (wGRR)
- Calculated viral wGRR using the formula from [**J. A. M. de Sousa et al. (2023)**](https://academic.oup.com/nar/article/51/6/2759/7068371?login=true)
   ```
   # make blastdb
   diamond blastp --threads 50 --db viral_proteins.dmnd --out diamond.tsv --evalue 0.0001 --max-target-seqs 100000 --query  all_votu.faa --id  35 --query-cover 50 --subject-cover 50
   # calculate wGRR
   python calculate_wGRR.py diamond.tsv all_votu.faa result_file
   ```
### 6. Enrichment analysis
   
- We use the script `enrichment_analysis.py` to perform environment enrichment analysis. The fold enrichment was calculated by dividing the proportion of SdPs from a specific environment by the proportion of that environment within all environments.

### 7. Assess nucleotide-level divergence between viral and host genomes

- We calculated nucleotide-level divergence between viral and host genomes using [**VirHostMatcher**](https://github.com/jessieren/VirHostMatcher)
  ```
   python vhm.py -v prophage -b host -o output
  ```

### 8. Host growth features calculation

- We computed host growth characteristics using [**genomeSPOT**](https://github.com/cultivarium/GenomeSPOT) and [**Tome**](https://github.com/EngqvistLab/Tome) 
  ```
  # genomeSPOT
  python -m genome_spot.genome_spot --models models --contigs host_fna  --proteins host_faa --output host_spot
  # Tome
  tome predOGT --indir host_faa -o predicted_ogt.tsv > batch_Tome.log
  ```
### 9. Simulating genomes

- Based on the methods of [**CheckM**](https://genome.cshlp.org/content/25/7/1043.full) and [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7), we simulated bacterial genomes and viral genomes with varying completeness and contamination levels. Script for simulating bacterial genomes: `generate_mock_host_genome.py`; Script for simulating viral genomes: `generate_mock_provirus.py`
   ```
   # to simulate a host genome with target genome completeness and contamination
   python generate_mock_host_genome.py bacterial_genome -c target_completeness -cont target_comtamination  -cg dir-for-contaminant-genomes -o output_file
   # to simulate a viral genome with target genome completeness and contamination
   python generate_mock_provirus.py host_genome -p contig:viral_start-virla_end -c target_completeness -cont target_comtamination -o output_file
   ```

## For Figures
- [**Figure1**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Figure%201):
   - Figure 1c violin plot: Generated by `Figure1-c-violin-plot.R`.
      - `ecoli-k12.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt`: HI scores for all PBSs in E. coli K12.
      - `ecoli-k12_point.tsv`: 24 experimentally validated PBSs in E. coli K12
   
- [**Figure2**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Figure%202):
   - Figure 2b phylogenetic tree:
      - Tree data: `Figure2-b-bac_cut_genus_30.tree`
      - Annotation:`Figure2-b-annotation_for_tree.xlsx
      - Visualization tool: [**TVBOT**](https://www.chiplot.online/tvbot.html)
   - Figure 2c network: Data files exceed GitHub size limits (available upon request)
   - Figure 2e violin plot: Generated by `Figure2-e-violin.R`
     - Data files: `(1~5)_*.txt` (prophage genomic features).
   - Figure 2f phage-host dissimilarity: Generated by `Figure2-f-distance.R`
     - Data files:`6_all_d2star_k6.txt` and `7_phage_codon_distance.txt`.

- [**Supplemental Figure S1**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S1):
   - Generated by `Figure S1.R`
     - Data files: `ecoli-k12-HI-predict-all.tsv` and `25point-all.tsv`.

- [**Supplemental Figure S2**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S2)
   - Fig S2b violin plot: Generated by `Figure S2-b-violin-plot.R`
      - Data files:`wp3.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt` and `wp3_point.tsv`。
   - Fig S2f qPCR data: `qPCR_data.tsv`.
     
- [**Supplemental Figure S3**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S3)
   - Genome map: Generated by [**clinker**](https://github.com/gamcil/clinker)
   - Violin plots: Generated by`Figure S3.R`
     - Data files:`*.tsv` and `*.txt`
     
- [**Supplemental Figure S4**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S4)
   - Fig S4c WP2 violin plot: Generated by `Figure S4-b-violin-plot.R`
     - Data files:`wp2.faa_intergenic_region.fna_pattern.txt_whole_genome_HI.txt` and `wp2_point.tsv`
   - Fig S4c [**IGV visualization**](https://igv.org/):
     - Control group:：`C-1_sorted.tdf` and `C-2_sorted.tdf`
     - MIMIC group：`M-1_sorted.tdf`,`M-2_sorted.tdf`,and `M-3_sorted.tdf`
   - Fig S4d VIP-Seq analysis data: `Figure S4-d. bar-plot-data.tsv`
  
- [**Supplemental Figure S5**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S5)
   - Fig S5b/d/e: Generated by`Figure S5-b-d-e.R`
     - Data files：`class-distribution.tsv` (panel b), `1-6-cluster.txt` (panel d), `3-type-hiscore.tsv` (panel e)
  
- [**Supplemental Figure S6**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S6)
   - Fig S6a/b: Generated by `Figure S6-a-b.R`
     - Data files: All `.tsv` files in directory.
   - Fig S6c/d/e heatmaps: Generated by `Figure S6-c-d-e-heatmap.R`
     - Data files: `host-HI-infor-top10genus-with-family.tsv` and `host-HI-infor-top10genus-with-genus.tsv`
  
- [**Supplemental Figure S7**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S7)
   - Fig S7b/c: Generated by `Figure S7-b-c.R`
     - Data files：`multi-phage-for-r.tsv`(panel b) and `regulation-consistence-for-r.tsv`(panel c)
  
- [**Supplemental Figure S8**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S8)
   - Fig S8b: Generated by `Figure S8.R`.
   - Note: Data files exceed display limits (available upon request)
  
- [**Supplemental Figure S9**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S9)
   - Fig S9a/b/c: Generated by `Figure S9.R`
      - Data files: `phage-distri-by-env-percent.tsv`(all environments), `human-distri-by-env-percent.tsv` (human) and `animal-distri-by-env-percent.tsv` (animal)
  
- [**Supplemental Figure S10**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S10)
   - Fig S10a/b/f: Generated by `FigureS10-a-b-f.R`
     - Data files:`sdp-sip.txt`(genome size/completeness) and `key-gene-data.tsv`(functional genes)
   - Fig S10c functional categories:
     - Data:`functional_category_distribution.tsv`
     - Annotation:`functional_category_annotation.tsv`
   - Fig S10d environmental distribution:
     - Data:`env-distribution.tsv`
     - Annotation:`env_annotation.tsv`
   - Fig S10e host distribution:
     - Data: host-distribution-top10-genus.tsv
     - Annotation: host-distribution-annotation.tsv
   - Visualization: [**Chiplot**](https://www.chiplot.online/#)
   
- [**Supplemental Figure S11**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S11)
   - Fig S11a host features: Generated by `FigureS11-a_host_features.R`
     - Data files: `1-5*.txt/tsv`
   - Fig S11b FAST analysis: Generated by `FigureS11-b_fast.R`
     - Data files: `6_growth.tsv`
   - Fig S11c host growth: Generated by `FigureS11-c_host_growth.R`
     - Data files: `7-10*.txt`
  
- [**Supplemental Figure S12**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S12)
   - Generated by `FigureS12-distance.R`
     - Data files: `3-10*.txt`
  
- **Supplemental Figure S13**: The figure is generated by [**clinker**](https://github.com/gamcil/clinker)
  
- [**Supplemental Figure S14**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S14)
  - Fig S14a/b point plots: Generated by `Fig S14-point-plot.R`
     - Data files: `calculate_rate_distribution2.tsv`
  
- [**Supplemental Figure S15**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S15)
   - Fig S15a/b/c/d: Generated by `Figure S15.R`
     - Data files: `prophagedb-qualify_vs-unqualified.tsv`
  
- [**Supplemental Figure S16**](https://github.com/mujiezhang/SOS-independent-prophages/tree/main/scripts%20and%20data%20for%20figure/Supplemental%20Figure%20S16)
   - Fig S16a/b: Generated by `FigureS16.R`
     - Data files: `retained-top-20-for-r.tsv` and `filtered-out-top-20-for-r.tsv`
  


