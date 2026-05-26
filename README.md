# GABBI pipeline version 1.2.0

<p align="center">
  <img width=450 src="/image/GABBI_logo.png" alt="Genome Alignment-Based Bait Inference"/>
</p>

**GABBI** is a fully automated pipeline for designing target capture baits (or probes) from **whole-genome alignments**.
Probe sets derived from GABBI are expected to target more variable loci than those produced by conventional probe design methods that rely on a single reference genome for mapping reads, thereby providing more sensitive and phylogenetically informative data.

This pipeline was developed and described in the following publication, which presents the first weevil-specific probe set. If you use GABBI, please cite:
> B. Zelvelder, L. Benoit, A. Loiseau, J. Haran, R. Allio (2026) A new method based on genome alignments provides a highly resolutive target enrichment set for weevils (Coleoptera, Curculionoidea); bioRxiv 2026.05.09.724036; doi: [https://doi.org/10.64898/2026.05.09.724036](https://doi.org/10.64898/2026.05.09.724036)

# Requirements
To run this pipeline, the only requirements are the [singularity image](https://cloud.sylabs.io/library/bjzelvelder/pipeline/gabbi) and a Linux environment with ```apptainer``` version >= 1.0.0 or ```singularity``` version >= 3.8.0 (see [Installation](#installation)). This pipeline was developed on Ubuntu 22.04 and is primarily written in shell.

# Table of contents
1. [Installation](#installation)
    - [RECOMMENDED: Downloading the GABBI singularity image](#downloading-the-gabbi-singularity-image)
    - [Building from source (in a singularity image or a sandbox)](#building-from-source)
2. [Pipeline overview](#pipeline-overview)
3. [How to use the GABBI pipeline](#how-to-use-the-gabbi-pipeline)
    - [Preparing input data](#preparing-input-data)
    - [Running the pipeline](#running-the-pipeline)
    - [Reading GABBI outputs](#reading-gabbi-outputs)
4. [What to do next](#what-to-do-next)
    - [Adjusting GABBI conservation thresholds](#adjusting-gabbi-conservation-thresholds)
    - [Running a species tree using the final set of targeted loci](#running-a-species-tree-using-the-final-set-of-targeted-loci)
    - [Generating the final probe set](#generating-the-final-probe-set)
    - [Annotating targeted loci as coding or non-coding sequences for downstream analyses](#annotating-targeted-loci-as-coding-or-non-coding-sequences-for-downstream-analyses)
    - [Testing the probe set _in silico_ on whole genome sequences (WGS)](#testing-the-probe-set-in-silico-on-whole-genome-sequences)
5. [Detailed options](#detailed-options)
6. [Citation](#citation)
7. [References](#references)

---
# Installation
The GABBI pipeline relies on **singularity** to execute within a pre-built container environment (a _singularity image_), eliminating installation issues and incompatibilities. Consequently, the only requirement is a Linux environment with ```apptainer version >= 1.0.0``` or ```singularity version >= 3.0.0```. To install apptainer, instructions are available at: [Installing apptainer](https://apptainer.org/docs/admin/main/installation.html). For clarity, all commands below use the ```singularity``` prefix, but these are fully interchangeable with ```apptainer```.

GABBI is distributed as a Singularity image and has no external dependencies beyond Singularity/Apptainer itself. On HPC clusters managed by SLURM, Singularity is typically available as an environment module that can be searched with ```module search singularity``` and can be loaded with ```module load <singularity module>``` prior to execution.

## Downloading the GABBI singularity image
Because singularity images are too large to be stored on GitHub, the ```.sif``` file is not included in this repository. Instead, the [GABBI image](https://cloud.sylabs.io/library/bjzelvelder/pipeline/gabbi) is hosted on the **Sylabs** remote repository and can be retrieved with the following command:
```
singularity pull --arch amd64 library://bjzelvelder/pipeline/gabbi:v1.2.0
```
The help section can then be accessed with:
```
singularity run-help gabbi_v1.2.0.sif
singularity run gabbi_v1.2.0.sif --help
```

## Building from 'source'
If the above commands failed, or if local modifications to the GABBI scripts are required, the singularity image can be built from the definition file ```GABBI_v1.2.0.def```. 
First, clone the GABBI repository:
```
git clone https://github.com/bjzelvelder/GABBI.git
cd GABBI
```
At this point, the GABBI singularity image can be built locally (which typically takes approximately 45 minutes). Alternatively, a GABBI sandbox can be built for interactive execution.

### Building the GABBI singularity image
Inside GABBI's repository, simply execute:
```
singularity build --fakeroot gabbi_v1.2.0.sif GABBI_v1.2.0.def
```
This generates a ```.sif``` image that can be run identically with:
```
singularity run gabbi_v1.2.0.sif --help
```

### Building a GABBI sandbox
Alternatively, a writable GABBI sandbox stored in ```GABBI_sandbox/``` can be created for interactive use or script modification:
```
singularity build --fakeroot --sandbox GABBI_sandbox/ GABBI_v1.2.0.def
mkdir -p GABBI_sandbox/$PWD
singularity shell -B $PWD --fakeroot --writable GABBI_sandbox/

```
Once inside the singularity shell, the local path is bound via ```-B $PWD```, making local files accessible from within the image. Note that modifications to bound files remain effective even from within the container. The GABBI pipeline can then be launched by sourcing the ```main.sh``` script:
```
source /opt/gabbi/main.sh --help
```
Note: because ```main.sh``` runs ```set -eou pipefail```, errors may terminate the singularity session. To prevent repeated reconnection to the sandbox, run the following command after sourcing ```main.sh```:
```
set +eou pipefail
```
---
# Pipeline overview
GABBI is a fully automated pipeline for designing target capture baits (or probes) from **whole-genome alignments**. Unlike commonly used probe design methods, GABBI does not rely on a **reference taxon**. The resulting targeted loci tend to be more variable than ultraconserved elements (UCEs) and are referred to as shared homologous regions (SHRs). Further details on the rationale and significance of this probe design approach are available in the associated publication ([Zelvelder et al. 2026 pre-print](#citation)).

The GABBI pipeline comprises six phases:
  - **1)** Whole-genome alignment of chromosome-level assemblies using [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) ([Armstrong et al. 2020](#references))
  - **2)** Identification of conserved regions using [MafFilter](https://github.com/jydu/maffilter) and [PhastCons](http://compgen.cshl.edu/phast/phastCons-tutorial.php) ([Hubisz et al. 2011; Dutheil et al. 2014](#references))
  - **3)** Extraction of shared conserved regions across chromosome-level genomes using BLAST (NCBI; [Camacho et al. 2009](#references))
  - **4)** Generation of a temporary probe set from these shared conserved regions, enriched with ancestral sequences inferred by IQ-TREE 3 ([Wong et al. 2025](#references))
  - **5)** Evaluation of the temporary probe set via _in silico_ target capture on additional genomes using [PHYLUCE](https://phyluce.readthedocs.io/en/latest/purpose.html) ([Faircloth 2016](#references))
  - **6)** Extraction of the final SHR set (i.e. targeted loci) enriched with ancestral sequences

<p align="left">
  <img src="/image/GABBI_pipeline.png" alt="GABBI_pipeline"/>
</p>

The primary goal of this pipeline is to make taxon-specific target-capture bait design as straightforward as possible, while remaining flexible enough to accommodate user-specific requirements. For this reason, users are required to provide only **chromosome-level genome assemblies**, a corresponding **guide tree** for the whole-genome alignment (see [preparing input data](#preparing-input-data)), and an additionnal set of **whole-genome assemblies** for probe validation. Alternatively, an existing genome alignment file (HAL format) can be provided as input. Each step of the pipeline is checkpointed to minimise redundant computation upon re-runs, and individual steps can be restarted to fine-tune the probe set with dataset-specific thresholds (see [detailed options](#detailed-options)).

> [!IMPORTANT]
> When provided 48 cpu cores, 4 chromosome-level genomes, and 7 additional weevil genomes, GABBI completed in approximately 13 hours and 30 minutes and produced 20 GB of intermediate files. Given its high degree of parallelisation, allocating as many CPU cores as possible is strongly recommended to reduce computation time, particularly for larger datasets. Several steps also require a lot of memory, so providing sufficient memory is equally recommended (especially on HPC clusters with ressource limitations).

---
# How to use the GABBI pipeline
Once the GABBI singularity image has been obtained and the help command successfully ran, the pipeline is ready for use. This section provides a detailed guide for preparing input data and interpreting GABBI outputs, uing the example dataset available in ```example_data```. This dataset comprises 4 chromosome-level genomes and 7 additional genomes representing a subset of weevils belonging to the Curculioninae subfamily (sensu 'CCCMS'; [Haran et al. 2023](#references)).

## Preparing input data
The first step of the GABBI pipeline is the whole-genome alignment with [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). If a HAL genome alignment file is already available, this step can be bypassed using the ```--hal``` option. Otherwise, GABBI requires two arguments: ```--chr-genomes``` and ```--guide-tree``` to run the genome alignment. In this case, the genomes to be used for probe design must be [downloaded](#downloading-ncbi-genomes) and [organised](#organizing-directory-architecture), and a phylogenetic tree featuring the chromosome-level genomes must be [provided as a guide to the alignment](#getting-a-guide-tree). But first, several considerations regarding the [selection of representative taxa](#choosing-representative-taxa) are detailed below.

### Choosing representative taxa
The primary objective of a probe set is to efficiently hybridise with the fragmented DNA of any taxon within the target taxonomic group. Accordingly, probe sequences should capture as much of the actual variation at each targeted locus as possible. Since representing the full diversity of a taxonomic group in probe design is impractical, a carefully selected subset of representative taxa must be used. Historically, very few genomes were employed to represent large taxonomic scales (e.g., only 6 taxa to represent the ~400,000 species of the Coleoptera UCE probe set; [Faircloth 2017](#references)). The growing availability of complete genome assemblies now allows a substantial improvement in taxonomic representation. The following considerations should be kept in mind:

- Avoid taxonomic redundancy among representative genomes. Certain genera or tribes are often over-represented among available assemblies, while entire families may be missing. Reducing taxonomic redundancy should prevent selected markers from appearing more broadly shared than expected.
- Be aware that including too many genomes can rapidly become computationnaly intensive. Further research is needed to determine an optimum, but providing one or two taxa per tribe is generally sufficient to capture the relevant genetic variation.
- Whole-genome alignments should preferably be conducted on chromosome-level (or very high quality) genome assemblies and validated against a separate set of additional genomes (with ```--add-genomes``` option).
- Testing the final probe set against a third independent set of assembled or unassembled data is also important to validate a probe set (see [What to do next?](#what-to-do-next)).

With these considerations in mind, available assemblies can be browsed at [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/) using filters and column selection to refine the choice of **chromosome-level** and **additional genomes**.

### Downloading NCBI genomes
Once the target genomes have been selected, they can be downloaded manually or automatically using the ```download_genomes_ncbi.sh``` script implemented in the singularity image. Download the NCBI table corresponding to your active selection, ensuring that the **"Assembly Level"** column is included (the script will automatically separate chromosome-level and other assemblies accordingly):

<p align="center">
  <img src="/image/screenshot_ncbi.png" alt="NCBI_screenshot" width="600"/>
</p>

Then, provide the NCBI table as an argument to the script, that will download, rename, and organize the assemblies according to organism names, assembly accessions, and assembly levels:

```
singularity exec gabbi_v1.2.0.sif download_genomes_ncbi.sh [path_to_ncbi_table]
```

### Organizing directory architecture
If genomes were downloaded using the ```download_genomes_ncbi.sh``` script as described above, they should already be appropriately organised. Chromosome-level assemblies that should not be included in the Cactus whole-genome alignment can be moved to the ```additional_genomes/``` folder (see [Choosing representative taxa](#choosing-representative-taxa).

As illustrated by the ```example_data``` directory structure, each genome must reside in its own subdirectory, typically named after its corresponding taxon. The ```download_genomes_ncbi.sh``` script organises assemblies accordingly, but custom assemblies can also be added following the same structure (ensuring the ```.fasta``` or ```.fna``` file extension is used).

High quality, chromosome-level genomes are strongly recommended for Cactus whole-genome alignments and must be stored in a separate directory to be parsed via the ```--chr-level-genomes``` option. The expected directory structure is as follows:

```
+--- chr_level_genomes/
|  +--- Ceutorhynchus_assimilis/
|  |  +--- GCA_917834065.1_PGI_CEUTPL_v4_genomic.fna
|  +--- Ips_nitidus/
|  |  +--- ...
|  +--- Orchestes_rusci/
|  +--- Polydrusus_cervinus/
|
+--- additional_genomes/
|  +--- Anthonomus_rubi/
|  +--- Kuschelorhynchus_macadamiae/
|  +--- Magdalis_cerasi/
|  +--- Oxystoma_pomonae/
|  +--- Pseudeuparius_sepicola/
|  +--- Sitona_lineatus/
|  +--- Xyleborus_glabratus/
```

For clarity, and given that taxonomic redundancy at the species level was avoided in the example dataset, the GenBank ID suffix appended by the ```download_genomes_ncbi.sh``` script was removed using:
```
for i in chr_level_genomes/* additional_genomes/*; do mv $i ${i%_GC*} ;done
```

> [!NOTE]
> Additional genomes provided via the ```--add-genomes``` option can be of any quality; however, low-quality assemblies can reduce SHR recovery and decrease the total number of targeted loci. In such cases, lowering the final [SHR threshold](#detailed-options) may be appropriate.

### Getting a guide tree
The whole-genome alignment step requires a guide tree to run. This phylogenetic tree must be provided as a **NEWICK** file with **taxon names matching exactly** those of chromosome-level genome directories in ```chr_level_genomes/```. The ```chr_level_genomes.nw``` file can be used as a format reference.

For the example dataset, an existing weevil phylogeny was pruned to retain only the relevant taxa, which were then renamed to match the subdirectory names in ```chr_level_genomes/``` exactly. This can be accomplished with tools such as [FigTree](https://github.com/rambaut/figtree) or [Mesquite](https://www.mesquiteproject.org/).

If no phylogenetic inference is available, the [ROADIES](https://github.com/TurakhiaLab/ROADIES) pipeline may be used to generate a rapid species tree from the genomic dataset ([Gupta et al. 2025](#references)).

> _This feature might be incorporated into the GABBI pipeilne in a future release._

## Running the pipeline
Once input data have been prepared, the GABBI pipeline can be executed. Here are a few examples of commands you can run based on the provided ```example_data```:

  - Run the pipeline with minimal options and default parameters, printing logs to STDOUT:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw
```


  - Customise basic pipeline outputs and redirect logs to a file (RECOMMENDED):
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --threads 64 --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output 2>&1 | tee -a GABBI_Curculioninae_v1.log
```
> GABBI-specific log entries can be accessed with ```grep "\[GABBI\]" GABBI_Curculioninae_v1.log```


  - Run GABBI with increased ressources for Cactus (e.g. on a HPC cluster):
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --threads 64 --verbose --cactus-maxDisk 500G --cactus-maxMemory 500G \
    --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output
```
> When running on a HPC cluster, it is advisable to request a large memory allocation or a high-memory node, as default ressources limits may not be sufficient for Cactus, which may end up interrupted.


  - Run (or resume) GABBI on a HAL file:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes \
    --hal GABBI_Curculioninae_v1.hal --threads 64 --verbose \
    --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output  2>&1 | tee -a GABBI_Curculioninae_v1.log
```



  - Run GABBI on a large HAL file, restrincting the number of MAF files to generate based on a subset of taxa (in development):
```
singularity run gabbi_v1.2.0.sif --hal GABBI_Curculioninae_v1.hal \
    --maf-references example_data/maf_references.txt --threads 64 --verbose \
    --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output  2>&1 | tee -a GABBI_Curculioninae_v1.log
```
> Without the ```--add-genomes``` option, GABBI will stop before phase ```05_adding_genomes```. Without the ```--chr-genomes``` option, chromosome-level genomes will be extracted from the HAL file. Only genomes listed in the ```--maf-references``` file will be extracted.


  - Run the GABBI pipeline on a reduced dataset with increased stringency thresholds:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --temp-tax-threshold 100 --temp-allow-dupes 0 \
    --shr-threshold 80 \
    --prefix GABBI_Curculioninae_v2 --out-dir GABBI_Curculioninae_v2_output 2>&1 |tee -a GABBI_Curculioninae_v2.log 
```


  - Adjust the final SHR threshold based on the multifasta table to increase the final number of targeted loci, without rerunning the full pipeline:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --temp-tax-threshold 100 --temp-allow-dupes 0 \
    --shr-threshold 70 --restart 5.5 \
    --prefix GABBI_Curculioninae_v2 --out-dir GABBI_Curculioninae_v2_output 2>&1 |tee -a GABBI_Curculioninae_v2.log 
```

Each step is checkpointed upon successful completion, so re-running the same command will resume the pipeline where it stopped.
If one step fails, the pipeline stops and will restart from the failed step upon re-run. In some cases, issues may not be correctly caught by the checkpointing system, requiring a manual restart from an earlier step once the issue has been corrected. Please report any issues in the corresponding GitHub Issues section.

>[!IMPORTANT]
> If you don't manage to get Cactus running, it may be run independently to fine-tune its parameters or execute it step by step by referring to the [Cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md).
> Once a ```.hal``` file has been produced, GABBI can be resumed using the ```--hal``` option. Ensure that the HAL file is named ```<prefix>.hal``` consistently with the ```--prefix <prefix>``` option.

## Reading GABBI outputs

All GABBI outputs and intermediate files are all stored in the directory specified by ```--output-dir```, or in ```GABBI_output/``` by default. They are organised such that each pipeline phase has its own subdirectory named after the corresponding phase (e.g. 01_cactus_alignment, 02_conserved_loci, etc.).

Files and summary statistics of interest are reported throughout the pipeline execution and can be accessed from the log file using:

```
grep "\[GABBI\]" [log_file]
```

Here is a more detailed list of GABBI outputs generated on the example dataset with the default ```cactus_alignment``` prefix.
- ```01_cactus_alignment/```
  - ```cactus_input.txt```: Input file created using ```--chr-level-genomes``` content and the provided ```--guide-tree```.
  - ```cactus_alignment.hal```: Cactus whole-genome alignment.
- ```02_conserved_loci/```
  - ```maf/```: Whole-genome alignment converted into a Multiple Alignment Format (MAF), using each genome as a reference sequence. Only log files are retained to save disk space. A subset of reference genomes can be provided with ```--maf-references``` option.
  - ```maffilter/```: Each genome has its own subdirectory containing filtered blocks of alignments from each chromosome or scaffold (based on ```--block-size``` and ```--block-length``` options). Only log file are retained to save disk space.
  - ```conserved_loci/*.merge.fasta```: Conserved loci found by PhastCons in each genome.
- ```03_cross_blast/```
  - ```shr_clustering/cactus_alignment.shr_from_blastn.mintax3.dupes0.list```: Filtered results of the cross-BLASTn between PhastCons conserved loci (based on ```--temp-tax-threshold``` and ```--temp-allow-dupes``` options).
- ```04_shr_extraction/```
  - ```shr/```: Temporary SHR sequences obtained from filtered cross-BLAST results.
  - ```cactus_alignment.temp.loci.fasta```: FASTA file containing all temporary targeted loci.
  - ```cactus_alignment.temp.anc.loci.fasta```: FASTA file containing all temporary targeted loci and their ancestral sequences.
- ```05_adding_genomes```
  - ```temp_probes/cactus_alignment.temp.anc.probes.fasta```: Temporary probes generated from temporary targeted loci.
  - ```mapping/```: Results from LASTZ mapping of temporary probes on each chromosome-level and additional genomes.
  - ```final_phyluce_probes/cactus_alignment.70.probe_list-DUPE-SCREENED.fasta```: Final probe set (without ancestral sequences) generated with PHYLUCE.
  - ```consensus_loci/cactus_alignment.70.phyluce.loci.cons.fasta```: Consensus sequences of the PHYLUCE final probe set.
  - ```multifasta_table/cactus_alignment.table```: Table reporting the number of loci conserved by X% of taxa (from 0 to 100%; useful for adjusting the ```--shr-threshold``` option).
- ```06_final_targeted_loci```
  - ```final_alignments/```: _Phylomera_ alignments of the final set of targeted loci. These can be used to easily reconstruct a species tree of all genomes in the dataset (see the following section [Phylogeny of the final set of targeted loci](#run-a-species-tree-using-the-final-set-of-targeted-loci)).
  - ```cactus_alignment.final.anc.loci.fasta```: Final set of targeted SHR and ancestral sequences.
  - ```cactus_alignment.final.anc.loci.cons.fasta```: Consensus sequences of the final set of targeted loci (to be used as references during phylogenomic alignments)

> The asterisk (*) means that multiple files share the same naming pattern.

# What to do next

Once GABBI has completed, a file containing all **targeted loci** from all genomic references, including genomes listed in the ```--chr-genomes``` directory or ```--maf-references``` file, genomes listed in the ```--add-genomes``` directory, and their ancestral sequences (labelled "Node#"), is available in ```GABBI_output/06_final_targeted_loci/prefix.final.anc.loci.fasta```. Note that this file **does not** contain probes (or baits): users are free to design probes independently or to request their bait synthesis provider to generate them from the targeted loci file. _Although this feature might be added in the GABBI pipeilne in a future release._

From this point, the following steps may be considered:
- [Adjusting GABBI conservation thresholds](#adjusting-gabbi-conservation-thresholds) to target more or less loci;
- Verifying the probe set performance by running a [phylogeny of the final set of targeted loci](#running-a-species-tree-using-the-final-set-of-targeted-loci);
- [Generating the final probe set](#generating-your-final-probe-set);
- [Annotating targeted loci as coding or non-coding sequences for downstream analyses](#annotating-targeted-loci-as-coding-or-non-coding-sequences-for-downstream-analyses);
- [Testing the probe set _in silico_ on simulated data from whole genome sequences (WGS)](#testing-the-probe-set-in-silico-on-whole-genome-sequences-WGS).

Recommendations for each of these steps are detailed below.

## Adjusting GABBI conservation thresholds

The default conservation threshold for retaining a locus in the final targeted set is **90%**, meaning that all loci that cannot be found in at least 90% of taxa will be excluded from the probe set. Although this threshold is very stringent, it yielded 4,255 loci for weevils in the [GABBI paper](#citation). If this threshold reduces the number of targeted loci too drastically in a given dataset, it can be lowered by restarting the relevant step using ```--restart 5.5``` option while modifying the ```--shr-threshold``` value in the same command (checkpoints will be detected automatically; see [example commands above](#running-the-pipeline)).

>[!IMPORTANT]
> Previous results will be overwritten. To preserve multiple sets of targeted loci derived from different thresholds (i.e. to compare [multiple species trees](#running-a-species-tree-using-the-final-set-of-targeted-loci)), save the ```06_final_targeted_loci``` directory by moving it to a separate location before re-running.

## Running a species tree using the final set of targeted loci

All targeted loci provided as independent multi-fasta files can be found in ```GABBI_output/06_final_targeted_loci/targeted_loci```. These can be aligned and cleaned to build a supermatrix and run a concatenated phylogenomic analysis of the data that produced our probe set. Note that this analysis relies exclusively on the targeted loci and does not include flanking sequence data. For _in silico_ tests simulating real capture data with reads and potential flanking regions, refer to [Testing the probe set _in silico_ on whole genome sequences](#testing-the-probe-set-in-silico-on-whole-genome-sequences) section.

GABBI includes a built-in tool called **_Phylomera_** that streamlines such phylogenomic analyses. It can be executed through the singularity image, giving all targeted loci as input and inferring a species tree with the substitution model of choice. The example below uses the MFP+MERGE algorythm to minimise biaises that may arise from over-partitionning.

```
singularity exec gabbi_v1.2.0.sif phylomera_v0.8.3.sh \
    --config /opt/gabbi/config/phylomera.conf \
    --input GABBI_output/06_final_targeted_loci/targeted_loci/  \
    --output GABBI_output/06_final_targeted_loci/final_tree  \
    --prefix GABBI_Chryso_10IV2026  \
    --threads 64  \
    --sptree MFP+MERGE  \
    --perc 70  \
    --continue
```
> Please refer to _Phylomera_ ```--help``` section to for parameter details.

Upon completion, IQ-TREE output files will be available in ```GABBI_output/06_final_targeted_loci/final_tree/TREES```.

## Generating the final probe set

Once the final set of targeted loci is satisfactory, probes must be designed for downstream _in silico_ and experimental analyses. As probe design and quality filtering were performed by an external provider in the [GABBI paper](#citation), this step is not currently supported by the GABBI pipeline. Users may design probes independently or submit the targeted loci file to their bait synthesis provider.

Given that GABBI outputs a large number of reference sequences per marker, including highly redundant ancestral sequences, it is highly recommended to remove redundant probes by clustering sequences based on percent identity and overlap. For instance, in the [GABBI paper](#citation), probes sharing more than 75% overlap and 80% identity were removed, reducing the probe set from 721k to 130k probes.

## Annotating targeted loci as coding or non-coding sequences for downstream analyses

To improve downstream phylogenomic inferences, targeted loci can be annotated as coding or non-coding sequences (for codon-partitioned analyses). This step requires at least one well-annotated genome in the dataset. The approach consists of retrieving the coordinates of each targeted locus on an annotated genome and cross-referencing them with the corresponding GFF file to identify overlaps with coding sequences. The detailed procedure is described in the [GABBI paper supplementary data](https://doi.org/10.5281/zenodo.20327231) available on Zenodo.

## Testing the probe set _in silico_ on whole genome sequences

Testing the probe set on simulated target capture data constitutes an important validation step. Any available whole-genome sequences (WGS) can be used to extend the dataset beyond the genomes used for probe design. Several approaches have been proposed in the literature, each suited to different analytical/taxonomical contexts (e.g. [PHYLUCE](#references), [HybPiper](#references), [IBA](#references)). In this section, I will detail the approach used in the [GABBI paper](#citation) using aTRAM [Allen et al. 2015](#references).

First, raw read data must be obtained, either from real WGS data or by simulating them from genome assemblies. In the latter case, ART Illumina ([Huang et al. 2012](#references)) can be used with default parameters, and the resulting reads trimmed with fastp ([Chen et al. 2018](#references)):

```
parallel "
    art_illumina --paired --in {}/{/}.fasta --out {}/{/}.pe150-reads --len 150 --fcov 15 --mflen 500 --sdev 100 -na
    gzip {}/{/}.pe150-reads*
    fastp -i {}/{/}.pe150-reads1.fq.gz -I {}/{/}.pe150-reads2.fq.gz -o {}/{/}.R1.trim.fastq.gz -O {}/{/}.R2.trim.fastq.gz -h {}/fastp.html -j {}/fastp.json
    rm {}/{/}.pe150-reads*
" ::: chr_level_genomes/* additional_genomes/*
```
> Note that GNU parallel is used here to parallelise the process, but these commands can equivalently be executed in a sequential loop.


Performing target capture with aTRAM on large simulated datasets can be computationally intensive, especially with large probe sets. To substancially reduce the number of reads to process (incidentally simulating the hybridization step of target capture), we will be using BBMap ([Bushnell 2014](#references); implemented in the GABBI image) to pre-filter reads against the probe set, with a soft identity threshold of 50%.

```
singularity exec gabbi_v1.2.0.sif bbmap.sh build=1 ref=baits-Moderate-RM25pc-0MT-720061count.fas.clust-75-80
# Usage: singularity exec gabbi_v1.2.0.sif run_BBMap.sh <R1> <R2> <out_R1> <out_R2> <build> <threads>
parallel -j 4 "
    singularity exec gabbi_v1.2.0.sif run_BBMap.sh {}/{/}.R1.trim.fastq.gz {}/{/}.R2.trim.fastq.gz {}/{/}.R1.bbmap.fastq.gz {}/{/}.R2.bbmap.fastq.gz 1 8
" ::: chr_level_genomes/* additional_genomes/*

```
> Here, ```baits-Moderate-RM25pc-0MT-720061count.fas.clust-75-80``` refers to **probe** file, not the targeted loci file (see [Generating the final probe set](#generating-the-final-probe-set) section).
> Note that BBMap has high memory requirements, so avoid launching too many processes in parallel.


Filtered reads are now ready to be processed by aTRAM. The aTRAM script implemented in the GABBI image have been [slightly adjusted](https://github.com/juliema/aTRAM/issues/321) to improve performance on target capture data more with Spades 4.2.0 ([Prjibelski et al. 2020](#references)). First, the aTRAM database must be generated using the ```atram_preprocessor.py``` script:

```
mkdir -p aTRAM_db
parallel "
    singularity exec gabbi_v1.2.0.sif atram_preprocessor.py --blast-db=aTRAM_db/{/} --end-1={}/{/}.R1.bbmap.fastq.gz --end-2={}/{/}.R2.bbmap.fastq.gz --gzip --cpus 1
" ::: chr_level_genomes/* additional_genomes/*
```

The ```aTRAM_db``` directory contains all necessary files for running aTRAM. To parallelise and optimise computation time, the probe file can be split into multiple subsets, each used as a reference for iterative read assembly. While computationally demanding, this approach yields multiple assemblies per marker, providing a notion of coverage useful for identifying duplicate sequences and contaminations. The following commands apply to a HPC cluster with 500 CPUs available per user, but this can be adapted to any ressource limitations. In this example, the probe file is split into 500 subsets, with 25 jobs launching 20 aTRAM instances in parallel, each using a single CPU, across all species.

```
mkdir -p split_probes_500 refseq_files
singularity exec gabbi_v1.2.0.sif split_fasta.py baits-Moderate-RM25pc-0MT-720061count.fas.clust-75-80 split_probes_500 500
for s in {0..480..20};do
    awk -v s="$s" 'FNR>s && FNR<=(s+20)' <(ls split_probes_500/*) > refseq_files/refseq_file$s.txt
done
```

aTRAM can then be run as follows:

- ON A SLURM CLUSTER: retrieve run_aTRAM_v2.sh in [GABBI paper Supplementary files](https://doi.org/10.5281/zenodo.20327231), adjust slurm options, and replace ```singularity run aTRAM.sif``` by ```singularity exec gabbi_v1.2.0.sif atram.py```.
```
mkdir -p slurm-logs capture
for sp in chr_level_genomes/* additional_genomes/*; do
    for ref in refseq_files/*;do
        sbatch -o slurm-logs/slurm-aTRAM_GABBI_${sp##*/}_${ref##*/}.log run_aTRAM_v2.sh capture/${sp##*/} aTRAM_db $ref
    done
done
```

- ON A LOCAL COMPUTER: adjust the number of parallel jobs and spades memory according to available ressources.
```
mkdir -p aTRAM_logs capture
for sp in chr_level_genomes/* additional_genomes/*; do
    mkdir -p capture/${sp##*/}/tmp capture/${sp##*/}/aTRAM
done
parallel -j 80 "
    singularity exec gabbi_v1.2.0.sif atram.py \
        -i 3 -Q {2} -b aTRAM_db/{1/} -o capture/{1/}/aTRAM/ \
        --evalue 1e-3 --word-size 11 --blast-max-target-seqs 100 \
        -a spades --spades-careful --spades-threads 1 --spades-memory 16 \
        --cpus 1 -t capture/{1/}/tmp > aTRAM_logs/aTRAM_GABBI_{1/}_{2/}.log
" ::: chr_level_genomes/* additional_genomes/* ::: split_fasta_500/*
```

aTRAM results are stored in each taxon's subdirectory (e.g. ```capture/Polydrusus_cervinus/aTRAM/```), which contains ```*all_contigs.fasta``` and ```*filtered.fasta``` files corresponding to all spades assemblies and assemblies filtered by back-blasting against the reference, respectively. Filtered assemblies from each marker are then merged per taxon into another dedicated subdirectory ```capture/Polydrusus_cervinus/cons```:

```
regroup_uce() {
    i=$1;
    mkdir -p $i/cons;
    find $i/aTRAM/ -type f -name "*filt*"|egrep -o "uce_[0-9]+_"|sort -u|while read uce;do sed -E "/>/s/>.*(contig_id=.*)/>${i##*/}__${uce}_\1/g" $i/aTRAM/*${uce}*filt* > $i/cons/${i##*/}_${uce}filtered.fasta;done;
    find $i/cons/ -size 0 -print -delete;
}
export -f regroup_uce
parallel -j 80 regroup_uce {} ::: capture/*
```

Pairwise distances between all assemblies of each marker are then computed to detect duplicated sequences and potential contaminations. For this purpose, assemblies are first aligned:

```
find capture -type f | grep "cons" > files_to_align.txt
parallel -j 80 "mafft --auto --adjustdirectionaccurately {} > {.}.mafft.fasta" :::: files_to_align.txt
```
> At this stage, each alignment file contains assemblies obtained from all probes targeting a given locus and should therefore be represented by highly similar sequences.


A consensus sequence is then derived from the largest cluster of nearly identical sequences (allowing up to 5% pairwise distance). If multiple clusters of sequences are detected and the second largest cluster represents more than 10% of all sequences in the alignment, the alignment is discarded. It is important to store the log output from this command as it will be used to remove flagged markers.

```
cons() {
    i=$1
    cd $i/cons
    # Usage: make_consensus_from_mafft_v3.R <all|aligned fasta> <iupac|majority|majseq|clusterize> <distance> [max cluster size (with majseq)]
    singularity exec gabbi_v1.2.0.sif make_consensus_from_mafft_v3.R all majseq 0.05 0.1
}
export -f cons
parallel -j 80 cons {} ::: capture/* > cons.GABBI.log
```

Resulting assemblies can then be merged across all taxa to produce one file per marker:

```
transfer() {
    i=$1
    out=$2
    echo "Transferring $i loci into $out"
    find $i/cons/ -type f -name "*.cons"|egrep -o "uce_[0-9]+_"|while read uce;do
        sed -E "s/>.*/>"${i##*/}"__"${uce}"/g" $i/cons/*${uce}*cons >> $out/${uce%_*}.fasta
    done
}
export -f transfer
parallel -j 80 transfer {} aTRAM_results ::: capture/*
```

During these post-processing steps, alignments containing sequences from multiple origins (either due to duplication or contamination) were removed. Using the consensus command log, markers flagged across multiple taxa can be identified and excluded based on a frequency threshold. During the weevil GABBI probe design, a stringent threshold was used to remove markers duplicated in more than 2 individuals. However, this value can be considerably relaxed when working with real target capture data, expected to be much more prone to contamination than simulated and pre-filtered WGS data.

```
grep "discarded" cons.GABBI.log |egrep -o "uce_[0-9]+_"|sort |uniq -c|sort -n|awk -v "threshold=2" '$1>threshold { print $2 }' > blacklisted_markers.2.GABBI.txt
mkdir aTRAM_results_bl_2
find aTRAM_results/ -type f|egrep -o "uce_[0-9]+"|sort -u|grep -v -f <(sed -E "s/_$/\$/g" blacklisted_markers.2.GABBI.txt )|while read uce;do
    cat aTRAM_results/$uce.fasta > aTRAM_results_bl_2/$uce.fasta
done
```

The resulting directory is finally ready to follow phylogenomic analyses. _Phylomera_ can be run using the consensus reference sequences computed by GABBI in ```GABBI_output/06_final_targeted_loci/cactus_alignment.final.anc.loci.cons.fasta``` and, if available, the [annotated loci file](#annotating-targeted-loci-as-coding-or-non-coding-sequences-for-downstream-analyses).

```
singularity exec gabbi_v1.2.0.sif phylomera_v0.8.3.sh \
    -i aTRAM_results_bl_2 \
    -o phylomera_tax11.no_cds \
    -pre Curculioninae_GABBI.tax11.no_cds \
    -r GABBI_output/06_final_targeted_loci/cactus_alignment.final.anc.loci.cons.fasta \
    -n 11 \
    -g MFP -p 0 \
    -t 80
```
> This command performs alignment cleaning, splits core and flanking regions, and infers gene trees with IQ-TREE using ModelFinder for substitution model selection. To run a concatenated analysis, replace ```-g MFP -p 0``` with, for example, ```-s MFP+MERGE``` (to select the minimum percentage of taxa required to keep a marker interactively) or ```-s MFP+MERGE -p 70``` (to apply the commonly used threshold of 70% taxon completeness). Running this command after the completion of gene trees will **not** overwrite or delete them.

---
# Detailed options
```

    GABBI v1.2.0 — Genome Alignment-Based Bait Inference pipeline
    ======================================================

    Main arguments
    ------------------
      --chr-genomes   DIR   Directory containing chromosome-level genome assemblies to be aligned with Cactus.
                            Each genome should reside in its own sub-directory, named <Taxon_name>.
      --guide-tree    FILE  Newick-format species tree used as a guide topology by Cactus for genome alignment.
    OR
      --hal           FILE  HAL Genome alignment file. If --chr-genomes is provided, sub-directory names must
                            match with genomes names in the HAL file. Otherwise, genomes will be extracted from
                            the HAL file.
      --add-genomes   DIR   Directory containing additional genome assemblies to validate SHR loci for final
                            probe set. Each genome should reside in its own sub-directory, named with
                            a unique <Taxon_name>. Without this option, the pipeline will stop before phase 5.

    Optionnal arguments
    -------------------
      Default values are given in square brackets [].

      Pipeline control
      ----------------
      --prefix        STR   Prefix used throughout the pipeline specific to your dataset [cactus_alignment]
      --out-dir, -o   DIR   Directory that will store GABBI outputs [GABBI_out]
      --restart       FLT   Specify a step to restart when running GABBI again. Valid values are:
                                    1_cactus_alignment, 2.1_hal2maf, 2.2_split_maf, 2.3_maffilter,
                                    2.4_phastcons, 2.5_conserved_loci, 3.1_cross_blast, 3.2_shr_clustering,
                                    4.1_shr_extraction, 4.2_alignments, 4.3_ancestral_seqs, 4.4_temp_loci,
                                    5.1_temp_probes, 5.2_add_genomes_prep, 5.3_lastz, 5.4_multifasta_table,
                                    5.5_final_phyluce_probes, 5.6_consensus_loci, 6.1_extract_targeted_loci,
                                    6.2_align_targeted_loci, 6.3_final_ancestral_seqs, 6.4_final_targeted_loci
                                    [none]
      --stop-before   STR   Interrupt the pipeline before the given phase. Valid phases are numbered from 1 to 6:
                                    01_cactus_alignment, 02_conserved_loci, 03_cross_blast, 04_shr_extraction,
                                    05_adding_genomes, 06_final_targeted_loci [none]
      --verbose             Print additional diagnostic information to stderr during execution.
                            Useful for development and troubleshooting [off]
      --threads, -t   INT   Number of CPU threads to allocate [all available by nproc command]
      --help, -h            Shows this message and exits

      Genome alignment options
      ------------------------
      --cactus-maxDisk   INT   Increase the maximum amount of disk used by cactus, with K, M or G suffix
                               to specify Kilo, Mega or Gigabytes. [300G]
      --cactus-maxCores  INT   Maximum number of cpu used by cactus (too many might cause a cactus to crash) [32]
      --cactus-maxMemory INT   Maximum amount of memory used by cactus (to avoid any issue, the higher is always
                               the better) [128G]
      --maf-references   FILE  File containing a list of reference genomes to generate MAF files from HAL genome
                               alignment. [all chromosome-level genomes]
      --block-size       INT   Minimum number of taxa (including ancestral genomes) required to retain
                               an alignment block [70 % of extent and ancestral genomes]
      --block-length     INT   Minimum length (nt) required to retain an alignment block [200]

      Temporary probe options
      -----------------------
      --cross-blast-ev     FLT    E-value threshold for cross-BLASTn [1E-6]
      --cross-blast-ws     INT    Word size for cross-BLASTn [11]
      --cross-blast-qc     INT    Minimum query coverage (%) for cross-BLASTn [50]
      --temp-tax-threshold INT    Minimum percentage of taxa sharing a locus after phastCons filtering [50]
      --temp-allow-dupes   INT    Maximum number of taxa allowed to carry duplicated loci
                                  after cross-BLASTn filtering [1]

      Final probe design options
      --------------------------
      --shr-threshold        INT  Minimum percentage of taxa sharing a final SHR locus to retain it in
                                  the phyluce probe set [90]
      --final_probes_tiling  INT  Tiling density used by phyluce to make the final phyluce probe set [3]
      --final_probes_masking FLT  Maximum fraction of masked sites authorized in the final phyluce probe set [0.25]
```

# Citation
If you used the **GABBI** pipeline or its dependencies, please cite:
> B. Zelvelder, L. Benoit, A. Loiseau, J. Haran, R. Allio (2026) A new method based on genome alignments provides a highly resolutive target enrichment set for weevils (Coleoptera, Curculionoidea); bioRxiv 2026.05.09.724036; doi: [https://doi.org/10.64898/2026.05.09.724036](https://doi.org/10.64898/2026.05.09.724036)

# References
References cited on this page:

> Allen JM, Huang DI, Cronk QC, Johnson KP. 2015 aTRAM - automated target restricted assembly method: a fast method for assembling loci across divergent taxa from next-generation sequencing data. BMC Bioinformatics 16, 98. [doi:10.1186/s12859-015-0515-2](https://doi.org/10.1186/s12859-015-0515-2)
> 
> Breinholt JW, Earl C, Lemmon AR, Lemmon EM, Xiao L, Kawahara AY. 2018 Resolving relationships among the megadiverse butterflies and moths with a novel pipeline for anchored phylogenomics. Systematic Biology 67, 78–93. [doi:10.1093/sysbio/syx048](https://doi.org/10.1093/sysbio/syx048)
> 
> Bushnell B. 2014 BBMap: a fast, accurate, splice-aware aligner. United States: Ernest Orlando Lawrence Berkeley National Laboratory, Berkeley, CA (US).
>
> Chen S, Zhou Y, Chen Y, Gu J. 2018 fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884–i890. [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
>
> Huang W, Li L, Myers JR, Marth GT. 2012 ART: a next-generation sequencing read simulator. Bioinformatics 28, 593–594. [doi:10.1093/bioinformatics/btr708](https://doi.org/10.1093/bioinformatics/btr708)
>
> Faircloth BC. 2016 PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32, 786–788. [doi:10.1093/bioinformatics/btv646](https://doi.org/10.1093/bioinformatics/btv646)
> 
> Faircloth BC. 2017 Identifying conserved genomic elements and designing universal bait sets to enrich them. Methods Ecol Evol 8, 1103–1112. [doi:10.1111/2041-210X.12754](https://doi.org/10.1111/2041-210X.12754)
> 
> Gupta A, Mirarab S, Turakhia Y. 2025 Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES. Proc. Natl. Acad. Sci. U.S.A. 122, e2500553122. [doi:10.1073/pnas.2500553122](https://doi.org/10.1073/pnas.2500553122)
> 
> Haran J et al. 2023 Phylogenomics illuminates the phylogeny of flower weevils (Curculioninae) and reveals ten independent origins of brood-site pollination mutualism in true weevils. Proc. R. Soc. B. 290, 20230889. [doi:10.1098/rspb.2023.0889](https://doi.org/10.1098/rspb.2023.0889)
>
> Prjibelski A, Antipov D, Meleshko D, Lapidus A, Korobeynikov A. 2020 Using SPAdes De Novo Assembler. CP in Bioinformatics 70, e102. [doi:10.1002/cpbi.102](https://doi.org/10.1002/cpbi.102)


References of tools used by the **GABBI** pipeline. Please consider citing these articles when running GABBI:
> **AMAS:** Borowiec ML. 2016 AMAS: a fast tool for alignment manipulation and computing of summary statistics. PeerJ 4, e1660. [doi:10.7717/peerj.1660](https://doi.org/10.7717/peerj.1660)
> 
> **BLAST+ suite:** Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009 BLAST+: architecture and applications. BMC Bioinformatics 10, 421. [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)
> 
> **Cactus:** Armstrong J et al. 2020 Progressive Cactus is a multiple-genome aligner for the thousand-genome era. Nature 587, 246–251. [doi:10.1038/s41586-020-2871-y](https://doi.org/10.1038/s41586-020-2871-y)
>
> **GNU Parallel:** Tange O. 2011 GNU Parallel - The Command-Line Power Tool. login: The USENIX Magazine, 42-47.
> 
> **HmmCleaner:** Di Franco A, Poujol R, Baurain D, Philippe H. 2019 Evaluating the usefulness of alignment filtering methods to reduce the impact of errors on evolutionary inferences. BMC Evol Biol 19, 21. [doi:10.1186/s12862-019-1350-2](https://doi.org/10.1186/s12862-019-1350-2)
> 
> **IQTREE3:** Wong T et al. 2025 IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models. [doi:10.32942/X2P62N](https://doi.org/10.32942/X2P62N)
> 
> **MACSE v2:** Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. 2018 MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Molecular Biology and Evolution 35, 2582–2584. [doi:10.1093/molbev/msy159](https://doi.org/10.1093/molbev/msy159)
> 
> **MafFilter:** Dutheil JY, Gaillard S, Stukenbrock EH. 2014 MafFilter: a highly flexible and extensible multiple genome alignment files processor. BMC Genomics 15, 53. [doi:10.1186/1471-2164-15-53](https://doi.org/10.1186/1471-2164-15-53)
> 
> **MAFFT:** Katoh K, Standley DM. 2013 MAFFT Multiple Sequence Alignment Software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30, 772–780. [doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)
>
> **OMM_MACSE pipeline:** Scornavacca C, Belkhir K, Lopez J, Dernat R, Delsuc F, Douzery EJP, Ranwez V. 2019 OrthoMaM v10: Scaling-Up Orthologous Coding Sequence and Exon Alignments with More than One Hundred Mammalian Genomes. Molecular Biology and Evolution 36, 861–862. [doi:10.1093/molbev/msz015](https://doi.org/10.1093/molbev/msz015)
> 
> **PHAST:** Hubisz MJ, Pollard KS, Siepel A. 2011 PHAST and RPHAST: phylogenetic analysis with space/time models. Briefings in Bioinformatics 12, 41–51. [doi:10.1093/bib/bbq072](https://doi.org/10.1093/bib/bbq072)
> 
> **PHYLUCE:** Faircloth BC. 2016 PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32, 786–788. [doi:10.1093/bioinformatics/btv646](https://doi.org/10.1093/bioinformatics/btv646)
> 
> **SeqKit:** Shen W, Le S, Li Y, Hu F. 2016 SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11, e0163962. [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)
