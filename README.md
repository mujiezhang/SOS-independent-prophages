# SOS-independent-prophages
This is the script repository for the following manuscript:

Yali Hao#, Mujie Zhang#, Xinjuan Lei, Chengrui Zhu, Xiang Xiao, Huahua Jian*, SOS-independent prophages prevail in the bacterial genomes. xxx (2025)

# Document
## For analysis
1. vOTU clustering based on ANI
   
- 我们使用CheckV提供的pipeline（https://bitbucket.org/berkeleylab/checkv/src/master/）进行了本文中vOTU的聚类，该pipeline基于all-versus-all BLASTn search和Leiden algorithm，following MIUViG guidelines (95% average nucleotide identity (ANI); 85% aligned fraction (AF)
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
     
2.  Genus and Family level clustering based on AAI

- 我们使用snayfach等人提供的pipeline （https://github.com/snayfach/MGV/tree/master/aai_cluster）进行了本文中genus和family水平的聚类，该pipeline基于all-versus-all BLASTp search和MCL
  - step1: all-vs-all blastp
    ```
    prodigal -a all_votu.faa  -i all_otu.fna   -p meta
    diamond makedb --in all_votu.faa --db viral_proteins --threads 10
    diamond blastp --query all_votu.faa --db viral_proteins --out blastp.tsv --outfmt 6 --evalue 1e-5 --max-target-seqs 1000000 --query-cover 50 --subject-cover 50
    ```
  - step2: calculate AAI
    ```
    python amino_acid_identity.py --in_faa query all_votu.faa --in_blast blastp.tsv --out_tsv aai.tsv
    ```
    我们修改了amino_acid_identity.py line21:`print "parse"`→`print("parse")`; line38:`print "compute"`→`print("compute")`; line52:`print "write"`→`print("write")`以适应python3的语法
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
    我们修改了genus过滤参数中的--min_aai:`--min_aai 40`→`--min_aai 50`, following the parameters in their paper [paper](https://www.nature.com/articles/s41564-021-00928-6)
## For Figure

