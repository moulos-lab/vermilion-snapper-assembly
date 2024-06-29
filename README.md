# Vermilion snapper assembly

The following template with indicative commands describes the general process of
creating the Vermillion snapper assemblies described in the article by 
Roa-Varon et al., 2024. The commands are indicative and no specific path 
structure or file-naming is defined, however the general workflow is depicted
and the process can be repeated.

## Data and files

To assemble, refine and polish the *Vermilion snapper* genome, we are using
PacBio HiFi long reads, an Illumina HiC library and also Illumina Whole Genome
Sequencing reads.

### Filenames

#### WGS

- WGS_FILE_R1.fastq.gz
- WGS_FILE_R2.fastq.gz

#### HiC

- HIC_FILE_R1.fastq.gz
- HIC_FILE_R2.fastq.gz

#### PacBio HiFi

- PACBIO_HIFI.fastq.gz

## Quality control

### PacBio HiFi reads

We are using HiFiAdapterFilt.
Indicative commands:

```
mkdir hififilt
export PATH=$PATH:/PATH/TO/HiFiAdapterFilt/DB
sh pbadapterfilt.sh -p PACBIO_HIFI -t 8
# Will output PACBIO_HIFI.filt.fastq.gz
```

### Illumina HiC and WGS reads 

We are using FastQC to generate QC reports.
Indicative commands:

```
fastqc -t 2 WGS_FILE_R1.fastq.gz WGS_FILE_R2.fastq.gz
fastqc -t 2 HIC_FILE_R1.fastq.gz HIC_FILE_R2.fastq.gz
```

We are using TrimGalore for trimming.
Indicative commands:

```
trim_galore \
  --length 50 \
  --output_dir . \
  --path_to_cutadapt /usr/bin/cutadapt \
  --cores $CORES \
  --paired \
  --retain_unpaired \
  --fastqc \
  --trim-n WGS_FILE_R1.fastq.gz WGS_FILE_R2.fastq.gz
# Will output WGS_FILE_R1.filt.fastq.gz WGS_FILE_R2.filt.fastq.gz

trim_galore \
  --length 50 \
  --output_dir . \
  --path_to_cutadapt /usr/bin/cutadapt \
  --cores $CORES \
  --paired \
  --retain_unpaired \
  --fastqc \
  --trim-n HIC_FILE_R1.filt.fastq.gz HIC_FILE_R2.filt.fastq.gz
# Will output HIC_FILE_R1.filt.fastq.gz HIC_FILE_R2.filt.fastq.gz
```

**From this point forward, we assume that:**

PACBIO_HIFI.fastq.gz = PACBIO_HIFI.filt.fastq.gz
WGS_FILE_R1.fastq.gz = WGS_FILE_R1.filt.fastq.gz 
WGS_FILE_R2.fastq.gz = WGS_FILE_R2.filt.fastq.gz 
HIC_FILE_R1.fastq.gz = HIC_FILE_R1.filt.fastq.gz
HIC_FILE_R2.fastq.gz = HIC_FILE_R2.filt.fastq.gz

## Initial *de novo* assembly

### hifiasm parameter optimization 

We are using Mabs.
Indicative command:

```
mabs-hifiasm.py \
  --pacbio_hifi_reads PACBIO_HIFI.fastq.gz \
  --short_hi-c_reads_R1 HIC_FILE_R1.fastq.gz \
  --short_hi-c_reads_R2 HIC_FILE_R2.fastq.gz \
  --download_busco_dataset actinopterygii_odb10.2021-02-19.tar.gz \
  --threads 32 \
  --output_folder mabs_output
```

`hifiasm` run without Mabs (results are almost identical).
Indicative command:

```
hifiasm -o rsnapper_v1.asm -t32 \
  --h1 HIC_FILE_R1.fastq.gz \
  --h2 HIC_FILE_R2.fastq.gz 
  PACBIO_HIFI.fastq.gz
```

Get fasta from hifiasm output:

```
awk '/^S/{print ">"$2;print $3}' rsnapper_v1.asm.hic.p_ctg.gfa \
  > rsnapper_v1.fasta
awk '/^S/{print ">"$2;print $3}' rsnapper_v1.asm.hic.hap1.p_ctg.gfa \
  > rsnapper_v1.hap1.fasta
awk '/^S/{print ">"$2;print $3}' rsnapper_v1.asm.hic.hap2.p_ctg.gfa \
  > rsnapper_v1.hap2.fasta
```

### QC on initial assembly

We are using QUAST.
Indicative command:

```
quast.py \
  -o quast rsnapper_v1.fasta \
  -pe1 WGS_FILE_R1.fastq.gz \
  -pe2 WGS_FILE_R2.fastq.gz
```

## HiC integration with Juicer and 3D-DNA pipeline
--------------------------------------------------------------------------------

The Juicer and 3D-DNA pipelines are less automated.
Below indicative commands based on the manuals.
Arima2 refers to new option added for the purposes of the manuscript 
(restriction enzyme sequences are not publicly available).

1. Preparation in the working directory

```
ln -s juicer-1.6/CPU/ scripts
cd scripts/common
ln -s juicer-1.6/juicer_tools.2.20.00.jar juicer_tools.jar
cd ../..

mkdir references
cp vsnapper_v1.fasta ./references/
cd references

bwa index vsnapper_v1.fasta
samtools faidx vsnapper_v1.fasta
cut -f1-2 vsnapper_v1.fasta.fai > vsnapper_v1.genome
cd ..

mkdir -p work/fastq
cp HIC_FILE_R1.fastq.gz ./work/fastq/HIC_FILE_R1.fastq.gz
cp HIC_FILE_R2.fastq.gz ./work/fastq/HIC_FILE_R2.fastq.gz

mkdir restriction_sites
cd restriction_sites

python juicer-1.6/misc/generate_site_positions.py Arima2 vsnapper_v1 \
  ./references/vsnapper_v1.fasta
mv vsnapper_v1_Arima2.txt ./restriction_sites/

cd work
scripts/juicer.sh -g vsnapper_v1 \
  -y ../restriction_sites/vsnapper_v1_Arima2.txt \
  -z ./references/vsnapper_v1.fasta \
  -p ./references/vsnapper_v1.genome \
  -a vsnapper_v1_hic \
  -D ./ \
  -t 32
```

2. Run 3D-DNA AFTER putting `lastz` in `$PATH`

```
export PATH=$PATH:lastz-1.04.22/src
run-asm-pipeline.sh \
  --rounds 3 \
  ./references/rsnapper_v1.fasta \
  ./work/aligned/merged_nodups.txt
```

3. At this point, the assembly should be visualized in Juicebox and manually
   curated.

4. After curation, create the final 3D-DNA assembly

```
run-asm-pipeline-post-review.sh \
  -r rsnapper_v1.rawchrom.review.assembly \
  -c 24 \
  --sort-output \
  rsnapper_v1.rawchrom.fasta ./work/aligned/merged_nodups.txt
```

5. Inspect again in Juicebox, repeat curation if necessary

## HiC integration with HiC-Pro and EndHiC pipeline

1. Preparation with HiC-Pro (hicreads directory created according to the HiC-Pro
   manual) - binsize 150k

```
HiC-Pro_3.1.0/bin/utils/digest_genome.py \
  -r [PROPRIETARY_CUT_SITES_HERE] \
  -o rsnapper_v1_hicpro.fasta \
  ./references/rsnapper_v1.fasta

# in hicpro_config.txt BIN_SIZE = 20000 40000 150000 500000 1000000
HiC-Pro_3.1.0/bin/HiC-Pro \
  -i hicreads \
  -o HiCProOut \
  -c config-hicpro.txt
```

2. Run EndHiC - binsize 150k

```
EndHiC/endhic.pl \
  --binsize 150000 \
  ./references/rsnapper_v1.genome \
  ./HiCProOut/hic_results/matrix/HIC_FILE/raw/150000/HIC_FILE_150000_abs.bed \
  ./HiCProOut/hic_results/matrix/HIC_FILE/raw/150000/HIC_FILE_150000.matrix \
  ./HiCProOut/hic_results/matrix/HIC_FILE/iced/150000/HIC_FILE_150000_iced.matrix

EndHiC/cluster2agp.pl \
  z.EndHiC.A.results.summary.cluster \
  ./references/rsnapper_v1.genome > scaffolds.agp

EndHiC/agp2fasta.pl \
  scaffolds.agp \
  ./references/rsnapper_v1.fasta > scaffolds.fa
```

3. Preparation with HiC-Pro (hicreads directory created according to the HiC-Pro
   manual) - binsize 100k

```
# in hicpro_config.txt BIN_SIZE = 25000 50000 100000 200000 500000 1000000
HiC-Pro_3.1.0/bin/HiC-Pro \
  -i hicreads \
  -o HiCProOut \
  -c config-hicpro.txt
```

4. Run EndHiC - binsize 100k

```
EndHiC/endhic.pl \
  --binsize 100000 \
  ./references/rsnapper_v1.genome \
  ./HiCProOut/hic_results/matrix/HIC_FILE/raw/100000/HIC_FILE_150000_abs.bed \
  ./HiCProOut/hic_results/matrix/HIC_FILE/raw/100000/HIC_FILE_150000.matrix \
  ./HiCProOut/hic_results/matrix/HIC_FILE/iced/100000/HIC_FILE_150000_iced.matrix

EndHiC/cluster2agp.pl \
  z.EndHiC.A.results.summary.cluster \
  ./references/rsnapper_v1.genome > scaffolds.agp
  
EndHiC/agp2fasta.pl \
  scaffolds.agp \
  ./references/rsnapper_v1.fasta > scaffolds.fa
```

Telomere evaluation script and report template on GitHub

## Assembly finalization

Polishing with PacBio HiFi reads and Illumina WGS reads.
Indicative commands for all tested assemblies described in the article:

1. Nextpolish2

1.1 Map long reads to the assemblies

```
minimap2 \
  -ax map-hifi \
  -t 16 rsnapper_v1_hifiasm.fasta \
  PACBIO_HIFI.fastq.gz | samtools sort -o hifi_hifiasm.bam

minimap2 \
  -ax map-hifi \
  -t 16 rsnapper_v1_3ddna.fasta \
  PACBIO_HIFI.fastq.gz | samtools sort -o hifi_3ddna.bam

minimap2 \
  -ax map-hifi \
  -t 16 rsnapper_v1_endhic_150k.fasta \
  PACBIO_HIFI.fastq.gz | samtools sort -o hifi_endhic_150k.bam

minimap2 \
  -ax map-hifi \
  -t 16 rsnapper_v1_endhic_100k.fasta \
  PACBIO_HIFI.fastq.gz | samtools sort -o hifi_endhic_100k.bam
```

1.2 Create k-mers with YAK and WGS reads

```
# *unpaired* refer to potentially lost mates after trim_galore which can be 
# useful in k-mer estimations

yak count \
  -t 16 \
  -o k19.yak \
  -k 19 \ 
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz) \
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz)

yak count \
  -t 16 \
  -o k23.yak \
  -k 23 \ 
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz) \
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz)

yak count \
  -t 16 \
  -o k27.yak \
  -k 27 \ 
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz) \
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz)

yak count \
  -t 16 \
  -o k31.yak \
  -k 31 \ 
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz) \
  <(zcat WGS_FILE_R1.fastq.gz) \
  <(zcat WGS_FILE_R1.unpaired.fastq.gz)
```

1.3 Run NextPolish2

```
NextPolish2/target/release/nextPolish2 \
  -t 32 \
  -o rsnapper_v1_hifiasm_polished.fasta \
  -u hifi_hifiasm.bam \
  rsnapper_v1_hifiasm.fasta k19.yak k23.yak k27.yak k31.yak

NextPolish2/target/release/nextPolish2 \
  -t 32 \
  -o rsnapper_v1_3ddna_polished.fasta \
  -u hifi_3ddna.bam \
  rsnapper_v1_3ddna.fasta k19.yak k23.yak k27.yak k31.yak

NextPolish2/target/release/nextPolish2 \
  -t 32 \
  -o rsnapper_v1_endhic_150k_polished.fasta \
  -u hifi_endhic_150k.bam \
  rsnapper_v1_endhic_150k.fasta k19.yak k23.yak k27.yak k31.yak

NextPolish2/target/release/nextPolish2 \
  -t 32 \
  -o rsnapper_v1_endhic_100k_polished.fasta \
  -u hifi_endhic_100k.bam \
  rsnapper_v1_endhic_100k.fasta k19.yak k23.yak k27.yak k31.yak
```

1.4 Visualize EndHiC output to Juibox

```
EndHiC/cluster_to_juciebox_assembly.pl \
  ./references/rsnapper_v1.genome 
  ../output/Round_A.04.summary_and_merging_results/z.EndHiC.A.results.summary.cluster \
  > vsnapper_v1_endhic_100k.assembly

juicer-1.6/misc/generate_site_positions.py \
  Arima2 rsnapper_v1_endhic_100k \
  rsnapper_v1_endhic_100k.fasta

# Create rsnapper_v1_endhic_100k.genome with samtools faidx
juicer_v1/scripts/juicer.sh \
  -g rsnapper_v1_endhic_100k \
  -y rsnapper_v1_endhic_100k_Arima2.txt \
  -z ./references/rsnapper_v1_endhic_100k.fasta \
  -p ./references/rsnapper_v1_endhic_100k.genome \
  -a rsnapper_v1_endhic_hic \
  -D ./ \
  -t 32
```

Then the assembly can be visualized with Juicebox

