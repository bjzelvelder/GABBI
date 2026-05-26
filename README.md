# GABBI pipeline version 1.2.0

<p align="center">
  <img width=450 src="/image/GABBI_logo.png" alt="Genome Alignment-Based Bait Inference"/>
</p>

**GABBI** is a fully automated pipeline to design target capture baits (or probes) using a **whole-genome alignment**.
GABBI-derived probe sets are expected to target more variable loci than usual probe design methods that rely on a base genome to map reads, providing more sensitive and phylogenetically resolutive data.

This pipeline was designed in this publication to produce the first set of weevil-specific probes. If you use GABBI, please cite:
> B. Zelvelder, L. Benoit, A. Loiseau, J. Haran, R. Allio (2026) A new method based on genome alignments provides a highly resolutive target enrichment set for weevils (Coleoptera, Curculionoidea); bioRxiv 2026.05.09.724036; doi: [https://doi.org/10.64898/2026.05.09.724036](https://doi.org/10.64898/2026.05.09.724036)

# Requirements
To run this pipeline, all you need is the [singularity image](https://cloud.sylabs.io/library/bjzelvelder/pipeline/gabbi) and a linux environment with ```apptainer``` version >= 1.0.0 or ```singularity``` version >= 3.8.0 installed (see [Installation](#installation)). This pipeline was designed on Ubuntu 22.04 and mainly written in shell.

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
    - [Tweak GABBI conservation thresholds](#tweak-gabbi-conservation-thresholds)
    - [Run a species tree using the final set of targeted loci](#run-a-species-tree-using-the-final-set-of-targeted-loci)
    - [Generate your final probe set](#generate-your-final-probe-set)
    - [Annotate targeted loci as coding or non-coding sequences for downstream analyses](#annotate-targeted-loci-as-coding-or-non-coding-sequences-for-downstream-analyses)
    - [Test the probe set _in silico_ on whole genome sequences (WGS)](#test-the-probe-set-in-silico-on-whole-genome-sequences)
5. [Detailed options](#detailed-options)
6. [Citation](#citation)
7. [References](#references)

---
# Installation
The GABBI pipeline uses **singularity** to run in a pre-built environment (a so-called _singularity image_) to avoid installation issues and incompatibilities. Thus, all you need is a linux environment with ```apptainer``` version >= 1.0.0 or ```singularity``` version >= 3.0.0 installed. This software can be installed following this tutorial: [Installing apptainer](https://apptainer.org/docs/admin/main/installation.html). For clarity, I will always refer to ```singularity``` rather than ```apptainer``` but all commands that start with ```singularity ...``` can be spelled ```apptainer ...```. 

GABBI is distributed as a Singularity image and has no external dependencies beyond Singularity/Apptainer itself. On HPC clusters managed by SLURM, Singularity is typically available as an environment module that can be searched throught ```module search singularity``` and can be loaded with ```module load <singularity module>``` prior to execution.

## Downloading the GABBI singularity image
Because a singularity image is too heavy to be stored and pulled from github, you will not find the ```.sif``` file in this repository. However, the [GABBI image](https://cloud.sylabs.io/library/bjzelvelder/pipeline/gabbi) is stored on the **Sylabs** remote repository and can simply be pulled with this command:
```
singularity pull --arch amd64 library://bjzelvelder/pipeline/gabbi:v1.2.0
```
You should then be able to see the help section with:
```
singularity run-help gabbi_v1.2.0.sif
singularity run gabbi_v1.2.0.sif --help
```

## Building from 'source'
If previous commands didn't work or if you want to make local modifications to the GABBI scripts, you can build the singularity image yourself using the definition file ```GABBI_v1.2.0.def```. 
First, get GABBI files by cloning the GABBI repository with:
```
git clone https://github.com/bjzelvelder/GABBI.git
cd GABBI
```
At that point, you can build the GABBI singularity image yourself (which should take about 40 minutes). You can also build a GABBI sandbox to run the pipeline interactively (for verboseging purposes).

### Building the GABBI singularity image
Inside GABBI's repository, simply execute:
```
singularity build --fakeroot gabbi_v1.2.0.sif GABBI_v1.2.0.def
```
This will create a ```.sif``` image that can be run the same way with:
```
singularity run gabbi_v1.2.0.sif --help
```

### Building a GABBI sandbox
Alternatively, if you want to make modifications to the GABBI scripts or simply run the pipeline more interactively from within the singularity image (environment), you can make a GABBI sandbox with a writable fakeroot that will be stored in ```GABBI_sandbox/```:
```
singularity build GABBI_sandbox/ GABBI_v1.2.0.def
mkdir -p GABBI_sandbox/$PWD
singularity shell -B $PWD --fakeroot --writable GABBI_sandbox/
```
Once you run singularity shell, you are inside the singularity image. Please be aware that binding your local path with ```-B $PWD``` allows GABBI to access your local files, but you can still make changes to them **even if you are inside the singularity image**.
You can now run the GABBI pipeline by launching the ```main.sh``` script:
```
source /opt/gabbi/main.sh --help
```
Note: because ```main.sh``` runs ```set -eou pipefail```, errors might eject you from the singularity image. To avoid having to connect back to the sandbox repeatidily, run this command after sourcing ```main.sh```:
```
set +eou pipefail
```
---
# Pipeline overview
The GABBI pipeline is a fully automated pipeline to design target capture baits (or probes) using a **whole-genome alignment**. As opposed to other commonly used methods of probe design, GABBI doesn't rely on a **base taxon**. The final set of targeted loci tend to be more variable than ultraconserved elements (UCE) and are referred to as shared homologous regions (SHR). Further details on the context and signficance of this new probe design method can be accessed through the publication that features the GABBI pipeline ([Zelvelder et al. 2026 pre-print](#citation)).

The GABBI pipeline can be segmented into 6 phases:
  - **1)** Aligning chromosome-level genomes using [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) ([Armstrong et al. 2020](#references))
  - **2)** Identifying conserved regions using [MafFilter](https://github.com/jydu/maffilter) and [PhastCons](http://compgen.cshl.edu/phast/phastCons-tutorial.php) ([Hubisz et al. 2011; Dutheil et al. 2014](#references))
  - **3)** Extracting shared conserved regions between chromosome-level genomes using BLAST (NCBI; [Camacho et al. 2009](#references))
  - **4)** Making a temporary probe set out of these shared conserved regions enriched with ancestral sequences using IQ-TREE 3 ([Wong et al. 2025](#references))
  - **5)** Testing the temporary probe set with _in silico_ target capture on additional genomes using [PHYLUCE](https://phyluce.readthedocs.io/en/latest/purpose.html) ([Faircloth 2016](#references))
  - **6)** Extracting the final set of SHR (or targeted loci) enriched with ancestral sequences

<p align="left">
  <img src="/image/GABBI_pipeline.png" alt="GABBI_pipeline"/>
</p>

The goal of this pipeline is to make the design of specific target-capture bait sets as straightforward as possible for any user, while remaining flexible enough to accommodate user-specific requirements. For this reason, you only need to provide **chromosome-level genomes** and a corresponding **guide tree** (see [preparing input data](#preparing-input-data) section) for the genome alignment (although high quality, scaffold-level genomes might work fine for our purpose), and an additionnal set of **whole-genome assemblies**. Alternatively, you can also run this pipeline with a genome alignment file (HAL). Each step of the pipeline is checkpointed to save time on re-run and can be restarted to fine-tune the probe set with specific thresholds that certainly depend on each dataset's specifics (see [detailed options](#detailed-options)).

> [!IMPORTANT]
> When provided 48 cpu cores, 4 chromosome-level genomes and 7 additional genomes of weevils, the GABBI pipeline took roughly 13 hours and 30 minutes to run in total and produced 20 Gb of intermediate files. As it is strongly paralellized, we strongly recommend giving as much cpu cores as possible to GABBI to reduce computational time even further for bigger datasets. Some steps also require a lot of memory, so giving as much memory as possible to GABBI is also recommended (especially on HPC clusters with ressource limitations).

---
# How to use the GABBI pipeline
Once you have fetched the GABBI singularity image and succesfully ran the help command, GABBI is ready to run. To help you prepare input data and reading GABBI output, this section will guide you through the pipeline in greater details using a small examplar dataset available in ```example_data```. This examplar dataset contains 4 chromosome-level genomes and 7 additional genomes of a small sample of weevils that belong to the Curculioninae subfamily (sensu 'CCCMS'; [Haran et al. 2023](#references)).

## Preparing input data
The first step of the GABBI pipeline is the genome alignment with [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus). If you already have a HAL genome alignment file, you can short-circuit this step using the ```--hal``` option. Otherwise, GABBI requires two arguments: ```--chr-genomes``` and ```--guide-tree``` to run the genome alignment. In this case, we need to [download](#downloading-ncbi-genomes) and [organize](#organizing-directory-architecture) the genomes that will be used for the probe design and provide a phylogenetic tree of the chromosome-level genomes to [guide the genomic alignment](#getting-a-guide-tree). But first, there are a few things to note on how to [choose those genomes](#choosing-representative-taxa).

### Choosing representative taxa
Ultimately, the goal of a probe set is to efficiently hybridize with the fragmented DNA of any taxon belonging to the targeted taxonomic group. Thus, we want as much as possible that our probe sequences map the actual variation of each targeted locus. As we cannot represent the entire diversity of our taxonomic target in the probe design (otherwise why would we even bother with probe design?), we have to rely on a drastic subsample of its diveristy by a handful of representative taxa. Historically, very few genomes were used to represent large taxonomic scales (e.g. only 6 taxa represent the 400k species of the Coleoptera UCE probe set; [Faircloth 2017](#references)), but thanks to the growing number of available complete genomes, we can greatly improve this number. Nonetheless, here are a few things you should have in mind:
- Avoid the taxonomic redundancy in representative genomes. Some genus/tribe are often overrepresented among available genomes whereas entire families can be missing. Reducing the taxonomic redundancy in the genome set should prevent selected markers to appear more shared than expected.
- Be aware that providing too many genomes can rapidly become computationnaly intensive. Further research is needed to quantify an optimum, but providing one or two taxa per taxonomic tribe should be enough to cover their genetic variation.
- Genomic alignments should be conducted on chromosome-level (or very high quality) genomes and validated with a separate set of additional genomes (with ```--add-genomes``` option). But testing the final probe set with a third set of assembled or unassembled data can also be useful (see [What to do next?](#what-to-do-next)).

With that in mind, you can check available genomes at [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/) and play with filters and column selection to fine-tune your selection of **chromosome-level** and **additional genomes**.

### Downloading NCBI genomes
Once you have chosen the genomes you want to represent in your probe design, you can download them manually but I recommend you to use the ```download_genomes_ncbi.sh``` script to download and organize them automatically. To do so, you must install the ```ncbi_datasets``` command-line program [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/); if you have conda installed, simply run:
```
conda install -c conda-forge ncbi-datasets-cli
```
Then, download the NCBI table of your active selection, making sure that you have **selected the "Assembly Level" column** on NCBI (the script will separate chromosome-level and other genomes accordingly):

<p align="center">
  <img src="/image/screenshot_ncbi.png" alt="NCBI_screenshot" width="600"/>
</p>

Then, give your NCBI table as an argument to the script that will download, rename and organize genomes depending on organism names, assembly accessions and assembly levels:

```
./download_genomes_ncbi.sh [path_to_ncbi_table]
```

### Organizing directory architecture
If you downloaded your genomes from NCBI using the ```download_genomes_ncbi.sh``` script as above, genomes should already be well organized, but you can still move some chromosome level genomes to the ```additional_genomes/``` folder if you don't want them to be part of cactus whole-genome alignment (see [Choosing representative taxa](#choosing-representative-taxa).

As you can tell by the ```example_data``` directory architecture, each genome must be in its own subfolder, typically named after the taxon name you want to keep during the analysis. The ```download_genomes_ncbi.sh``` script should have organized your genome selection in this way, but you can obvisouly add your own genomic alignments by following the same logic (simply make sure that they have the ```.fasta``` or ```.fna``` file extension).

High quality, chromosome-level genomes are strongly recommended for cactus whole-genome alignments, so they need to be stored in a separate folder as other genomes to be parsed through the ```--chr-level-genomes``` option. Here is how your directory architecture should look like:

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

For clarity, and because I avoided any taxonomic redundancy at the species level in my example dataset, I got rid of the GenBank ID suffix provided by the ```download_genomes_ncbi.sh``` script using this command:
```
for i in chr_level_genomes/* additional_genomes/*; do mv $i ${i%_GC*} ;done
```

> [!NOTE]
> Additional genomes parsed through the ```--add-genomes``` option can be of any quality, but low-quality genomes can impact SHR recovery and reduce the total number of targeted loci. If you have some in your dataset, you might consider lowering the final [SHR threshold](#detailed-options).

### Getting a guide tree
The whole-genome alignment step requires a guide tree to run. This phylogenetic tree must be provided as a **NEWICK** file with **matching taxon names** with the chromosome level genomes provided in the ```chr_level_genomes/``` folder. You can check the ```chr_level_genomes.nw``` format for reference.

Because a weevil phylogeny was already available for the example dataset, I simply pruned an existing tree to only contain the taxa I provided in my dataset and renamed them to match exactly with the ```chr_level_genomes/``` subfolders names. To do so, many tools exist but you can use [FigTree](https://github.com/rambaut/figtree) or [Mesquite](https://www.mesquiteproject.org/).

However, If you are unsure what your phylogenetic tree should look like, you may try to run the [ROADIES](https://github.com/TurakhiaLab/ROADIES) pipeline to compute a fast phylogenetic tree from your genomic dataset ([Gupta et al. 2025](#references)).

> _This feature might be added in the GABBI pipeilne in the future._

## Running the pipeline
Once you have prepared your input data, you can run the GABBI pipeline. Here are a few useful commands you can run based on the ```example_data``` provided:

  - Run the pipeline with minimal options and default arguments, printing logs to STDOUT:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw
```


  - Personnalize basic pipeline outputs and add a log file (RECOMMENDED):
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --threads 64 --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output 2>&1 | tee -a GABBI_Curculioninae_v1.log
```
> Informations outputed by GABBI only can be accessed with ```grep "\[GABBI\]" GABBI_Curculioninae_v1.log```


  - Run GABBI with increased ressources for cactus (e.g. to run on a HPC cluster):
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --threads 64 --verbose --cactus-maxDisk 500G --cactus-maxMemory 500G \
    --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output
```
> If you have access to a HPC cluster, don't hesitate to ask for a lot of memory or a specific node with big memory, as default ressources may not be enough for cactus to run without getting interrupted for memory limits.


  - Run (or resume) GABBI on a HAL file:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes \
    --hal GABBI_Curculioninae_v1.hal --threads 64 --verbose \
    --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output  2>&1 | tee -a GABBI_Curculioninae_v1.log
```



  - Run GABBI on a large HAL file, reducing the number of MAF references to generate:
```
singularity run gabbi_v1.2.0.sif --hal GABBI_Curculioninae_v1.hal \
    --maf-references example_data/maf_references.txt --threads 64 --verbose \
    --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output  2>&1 | tee -a GABBI_Curculioninae_v1.log
```
> Without the ```--add-genomes``` option, GABBI will stop before phase 05_adding_genomes and without the ```--chr-genomes``` option, chromosome-level genomes will be extracted from the HAL file. Only genomes listed in ```--maf-references``` file will be extracted. 


  - Run the GABBI pipeline on a reduced dataset, increasing default stringency thresholds:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --temp-tax-threshold 100 --temp-allow-dupes 0 \
    --shr-threshold 80 \
    --prefix GABBI_Curculioninae_v2 --out-dir GABBI_Curculioninae_v2_output 2>&1 |tee -a GABBI_Curculioninae_v2.log 
```


  - Change the final SHR threshold based on the multifasta table to increase the final number of targeted loci, without having to rerun the entire pipeline:
```
singularity run gabbi_v1.2.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --temp-tax-threshold 100 --temp-allow-dupes 0 \
    --shr-threshold 70 --restart 5.5 \
    --prefix GABBI_Curculioninae_v2 --out-dir GABBI_Curculioninae_v2_output 2>&1 |tee -a GABBI_Curculioninae_v2.log 
```

Each step is checkpointed if it succesfully ran, so running the same command again will resume the pipeline where it stopped.
If one step fails, the pipeline stops and will restart where it failed running the same command. Note that sometimes, issues are not correctly caught by the checkpointing system, so you might have to restart from an anterior step once you found the issue. Please report issues in the corresponding github sections to help me improve the pipeline!

>[!IMPORTANT]
> If you don't manage to get cactus running, you may consider running it separately to fine tune cactus options or run it step by step by referring to [cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md).
> In this case, once you have your ```.hal``` file, you can resume GABBI with the ```--hal``` option. Make sur that it is named ```<prefix>.hal``` as in the ```--prefix <prefix>``` option.

## Reading GABBI outputs

GABBI outputs and temporary files are all stored in the ```--output-dir``` provided or in ```GABBI_output/``` by default. They are organized so that each phase of the pipeline has its own subdirectory named after the phase name (e.g. 01_cactus_alignment, 02_conserved_loci etc.).

Files and statistics of interest are given through the pipeline. They can be accessed in the log file using:

```
grep "\[GABBI\]" [log_file]
```

Here is a more detailed list of some GABBI outputs (generated on the example dataset with the default ```cactus_alignment``` prefix).
- ```01_cactus_alignment/```
  - ```cactus_input.txt```: Input file created using ```--chr-level-genomes``` content and the provided ```--guide-tree```.
  - ```cactus_alignment.hal```: Cactus whole-genome alignment.
- ```02_conserved_loci/```
  - ```maf/```: Whole-genome alignment converted into a "Multiple Alignment Format" (MAF), using each genome as a base genome. Only logs remain to save disk space.
  - ```maffilter/```: Each genome has its own subdirectory containing filtered blocks of alignments of each of their chromosome/scaffolds (based on ```--block-size``` and ```--block-length``` options). Only logs remain to save disk space.
  - ```conserved_loci/*.merge.fasta```: Conserved loci found by phastcons in each genome.
- ```03_cross_blast/```
  - ```shr_clustering/cactus_alignment.shr_from_blastn.mintax3.dupes0.list```: Filtered results of the cross-BLASTn between phastcons conserved loci (based on ```--temp-tax-threshold``` and ```--temp-allow-dupes``` options).
- ```04_shr_extraction/```
  - ```shr/```: List of temporary SHR sequences obtained from filtered cross-BLAST results.
  - ```cactus_alignment.temp.loci.fasta```: FASTA file containing all temporary targeted loci.
  - ```cactus_alignment.temp.anc.loci.fasta```: FASTA file containing all temporary targeted loci and their ancestral sequences.
- ```05_adding_genomes```
  - ```temp_probes/cactus_alignment.temp.anc.probes.fasta```: Temporary probes generated from temporary targeted loci.
  - ```mapping/```: Results from LASTZ mapping of temporary probes on each chromosome-level and additional genomes.
  - ```final_phyluce_probes/cactus_alignment.70.probe_list-DUPE-SCREENED.fasta```: Final probe set (without ancestral sequences) generated with PHYLUCE.
  - ```consensus_loci/cactus_alignment.70.phyluce.loci.cons.fasta```: Consensus sequences of the PHYLUCE final probe set.
  - ```multifasta_table/cactus_alignment.table```: Table listing the number of loci conserved by X% of taxa (from 0 to 100%; useful to modify ```--shr-threshold``` option).
- ```06_final_targeted_loci```
  - ```final_alignments/```: Phylomera alignments of the final set of targeted loci. Can be used to easily build a species tree of all genomes in the dataset (see [future directions](#future-directions)).
  - ```cactus_alignment.final.anc.loci.fasta```: Final set of targeted SHR and ancestral sequences.
  - ```cactus_alignment.final.anc.loci.cons.fasta```: Consensus sequences of the final set of targeted loci (to be used as references during alignments)

> The asterisk (*) means that multiple files have the same pattern.

# What to do next

Once GABBI finished running, you should find a file containing all **targeted loci** from all genomic references, i.e. genomes listed in the ```--chr-genomes``` directory or ```--maf-references``` file supplied, genomes listed in the ```--add-genomes``` directory, and there ancestral sequences (named "Node#") in ```your_GABBI_output/06_final_targeted_loci/prefix.final.anc.loci.fasta```. Note that this file **does not** contain probes (or baits): you are free to design them on your own or ask your bait synthesizer company to generate them from your targeted loci file. _Although this feature might be added in the GABBI pipeilne in the future._

From now on, you may want to: [tweak GABBI conservation thresholds](#tweak-gabbi-conservation-thresholds) to target more or less loci, check if your probe set actually works by running a [phylogeny of the final set of targeted loci](#run-a-species-tree-using-the-final-set-of-targeted-loci), [generate your final probe set](#generate-your-final-probe-set), [annotate targeted loci as coding or non-coding sequences for downstream analyses](#annotate-targeted-loci-as-coding-or-non-coding-sequences-for-downstream-analyses) and/or [test the probe set _in silico_ on simulated data from whole genome sequences (WGS)](#testing-the-probe-set-in-silico-on-whole-genome-sequences-WGS).

Suggestions on how you may proceed for each of these steps are detailed below.

## Tweak GABBI conservation thresholds

As you may have noticed during the execution of the GABBI pipeline, the default conservation score to keep a locus in the final set of targeted loci is **90%**, meaning that all loci that cannot be found in at least 90% of taxa won't be targeted by your probe set. Even if this threshold is very stringent, it allowed us to keep 4,255 loci for weevils in the [GABBI paper](#citation). But if this threshold reduces the number of targeted loci too sharply in your dataset and you want to target more loci, at the cost of being less conservative, you can restart this step by adding the ```--restart 5.5``` option and changing the ```--shr-threshold``` option using the same command-line (checkpoints will be detected automatically, see [example commands above](#running-the-pipeline)).

>[!IMPORTANT]
> Previous results will be deleted, so if you want to save multiple sets of targeted loci based on different thresholds and [run multiple trees](#run-a-species-tree-using-the-final-set-of-targeted-loci) to compare them, save the ```06_final_targeted_loci``` directory by moving it out of GABBI's output.

## Run a species tree using the final set of targeted loci

All targeted loci reduced to one multi-fasta file per marker can be found in ```GABBI_output/06_final_targeted_loci/targeted_loci```. We can align and clean those markers to build a supermatrix and run a concatenated phylogenomic analysis of the data that produced our probe set. Keep in mind that this analysis is done on targeted loci only and doesn't have flanking data. To run _in silico_ tests of the probe set simulating real capture data with reads and potential flanking regions, please refer to [Test the probe set _in silico_ on whole genome sequences](#test-the-probe-set-in-silico-on-whole-genome-sequences) section.

GABBI has a built-in tool that I also designed to make such phylogenomic analysis straightforward, called **_Phylomera_**. We simply run _Phylomera_ through the singularity image, giving all targeted loci as input and asking for a species tree with the model of your choice. Here I use MFP+MERGE to minimize biaises that can result from over-partitionning our dataset.

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
> Please refer to _Phylomera_ ```--help``` section to fine-tune your parameters.

Once it finished running, you will find the IQ-TREE output files in ```GABBI_output/06_final_targeted_loci/final_tree/TREES```.

## Generate your final probe set

If you are satisfied with your final set of targeted loci, you will need to generate the probe set required for _in silico_ analyses. As we relied on an external company to design and clean ours, this step is not currently supported by the GABBI pipeline. You are free to design them on your own or ask your bait synthesizer company to generate them given your targeted loci file.

Because GABBI offers a lot of reference sequences for each marker with ancestral sequences that can be very redundant, we highly encourage you to remove redundant probes by clusterizing probes in your set checking their % of identity and % of overlap. For example, in the [GABBI paper](#citation), we removed redundant probes with over 75% overlap and 80% identity to reduce our probe set to 130k probes (instead of 721k).

## Annotate targeted loci as coding or non-coding sequences for downstream analyses

In order to improve downstream phylogenomic inferences, we also annotated our set of targeted loci as coding or non-coding (for codon-wide partitioning). This step requires you to have at least one finely annotated genome in your dataset. The idea is to fetch the coordinates of each targeted locus on each annotated genome and check with the GFF file if they overlap with coding sequences. We detailed our approach in the [GABBI paper supplementary data](https://doi.org/10.5281/zenodo.20327231) available on Zenodo, if you want to give it a try.

## Test the probe set _in silico_ on whole genome sequences

Lastly, we want to test our probe set on simulated target capture data. To do so, any whole genome sequences (WGS) can be used to get a bigger dataset than the one used to design the probes. Multiple approaches have been proposed in the literature, allowing users to choose the method best suited to their needs (e.g. [PHYLUCE](#references), [HybPiper](#references), [IBA](#references)). In this section, I will detail the approach I used in the [GABBI paper](#citation) using aTRAM [Allen et al. 2015](#references).

First of all, we need to get raw read data, either from real WGS data or simulating them from a genome assembly. In the latter case, we can use ART Illumina ([Huang et al. 2012](#references)) with default options and trim the resulting file with fastp ([Chen et al. 2018](#references)):

```
parallel "
    art_illumina --paired --in {}/{/}.fasta --out {}/{/}.pe150-reads --len 150 --fcov 15 --mflen 500 --sdev 100 -na
    gzip {}/{/}.pe150-reads*
    fastp -i {}/{/}.pe150-reads1.fq.gz -I {}/{/}.pe150-reads2.fq.gz -o {}/{/}.R1.trim.fastq.gz -O {}/{/}.R2.trim.fastq.gz -h {}/fastp.html -j {}/fastp.json
    rm {}/{/}.pe150-reads*
" ::: chr_level_genomes/* additional_genomes/*
```
> Note that we are using GNU parallel here to speed up the process, but this can be translated into a simple loop.


Performing target capture with aTRAM on large simulated datasets can be computationally intensive, especially with a large amount of probes. To drastically reduce the amount of reads to process (incidentally simulating target capture), we will be using BBMap ([Bushnell 2014](#references); implemented in the GABBI image) to filter reads based on our probe set, with a soft threshold of 50% minimum identity.

```
singularity exec gabbi_v1.2.0.sif bbmap.sh build=1 ref=baits-Moderate-RM25pc-0MT-720061count.fas.clust-75-80
# Usage: singularity exec gabbi_v1.2.0.sif run_BBMap.sh <R1> <R2> <out_R1> <out_R2> <build> <threads>
parallel -j 4 "
    singularity exec gabbi_v1.2.0.sif run_BBMap.sh {}/{/}.R1.trim.fastq.gz {}/{/}.R2.trim.fastq.gz {}/{/}.R1.bbmap.fastq.gz {}/{/}.R2.bbmap.fastq.gz 1 8
" ::: chr_level_genomes/* additional_genomes/*

```
> Here, ```baits-Moderate-RM25pc-0MT-720061count.fas.clust-75-80``` is the name of the file containing our **probes**, not targeted loci (see [Generate your final probe set](#generate-your-final-probe-set) section).
> You should also be aware that BBMap requires a lot of memory, so don't launch too many processes in parallel.


Filtered reads are now ready to be processed by aTRAM. I [slightly adjusted](https://github.com/juliema/aTRAM/issues/321) aTRAM scripts to process target capture data more efficiently on Spades 4.2.0 ([Prjibelski et al. 2020](#references)). Those changes are included in the GABBI singularity image as well. We first need to generate the databases for aTRAM to work with using the ```atram_preprocessor.py``` script:

```
mkdir -p aTRAM_db
parallel "
    singularity exec gabbi_v1.2.0.sif atram_preprocessor.py --blast-db=aTRAM_db/{/} --end-1={}/{/}.R1.bbmap.fastq.gz --end-2={}/{/}.R2.bbmap.fastq.gz --gzip --cpus 1
" ::: chr_level_genomes/* additional_genomes/*
```

The ```aTRAM_db``` folder contains all we need to run aTRAM, but we can prepare the probe set to parallelize the process and optimise computation time. The goal is to assemble reads iteratively using each probe sequence as a reference. Although this approach can be kind of overkill, it allows us to get multiple assemblies per marker, which offers some notion of coverage to evaluate duplicate sequences, contaminations and postprocess assemblies accordingly. The following commands apply to a HPC cluster with 500 CPUs available per user, but you can adjust them to your own ressources. My goal is to split the probe file into 500 files that will be processed each on one CPU. 25 jobs will launch 20 instances of aTRAM in parallel, each running on one CPU, repeating that for each species of my dataset.

```
mkdir -p split_probes_500 refseq_files
singularity exec gabbi_v1.2.0.sif split_fasta.py baits-Moderate-RM25pc-0MT-720061count.fas.clust-75-80 split_probes_500 500
for s in {0..480..20};do
    awk -v s="$s" 'FNR>s && FNR<=(s+20)' <(ls split_probes_500/*) > refseq_files/refseq_file$s.txt
done
```

We can now run aTRAM:
```
# ON A SLURM CLUSTER: get run_aTRAM_v2.sh in [GABBI paper Supplementary files](https://doi.org/10.5281/zenodo.20327231), change slurm options and replace "singularity run aTRAM.sif" by "singularity exec gabbi_v1.2.0.sif atram.py"
mkdir -p slurm-logs capture
for sp in chr_level_genomes/* additional_genomes/*; do
    for ref in refseq_files/*;do
        sbatch -o slurm-logs/slurm-aTRAM_GABBI_${sp##*/}_${ref##*/}.log run_aTRAM_v2.sh capture/${sp##*/} aTRAM_db $ref
    done
done

# ON A LOCAL COMPUTER: adjust CPU and spades memory depending on your ressources
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

aTRAM results are stored in each individual's subfolder (i.e. ```capture/Polydrusus_cervinus/aTRAM/```), which includes ```*all_contigs.fasta``` and ```*filtered.fasta``` files that respectively correspond to all spades assemblies, and assemblies blasted back to the reference to remove outliers. For each individual, we will be merging filtered assemblies from each marker in another subfolder ```capture/Polydrusus_cervinus/cons```:

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

Then, we want to compute pairwise distances between all assemblies of each marker to check duplicated sequences and contaminations. To do so, we first need to align them:

```
find capture -type f | grep "cons" > files_to_align.txt
parallel -j 80 "mafft --auto --adjustdirectionaccurately {} > {.}.mafft.fasta" :::: files_to_align.txt
```
> Note that at this point, each alignment file contains assemblies obtained from all probes targeting one locus. Thus, they are supposed to be identical.


Then, get the consensus sequence resulting from the biggest cluster of nearly identical sequences (allowing 95% identity and a second biggest cluster size of 10% of the number of sequences in the alignment). It is important to store the log of this command as it will be useful to remove flagged markers later.

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

Resulting assemblies can now be merged between all individuals to get one file per marker.

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

During these postprocessing steps, we removed alignments that contained sequences from multiple origins (either due to duplication or contamination). Using the consensus command log, we can flag markers that were removed in multiple individuals and decide to remove markers that are found duplicated/contaminated in too many individuals. During the weevil GABBI probe design, I used a very stringent threshold to remove markers duplicated in more than 2 individuals, but this value can be greatly increased with real target capture data, expected to be much more contaminated than simulated and pre-filtered WGS data.

```
grep "discarded" cons.GABBI.log |egrep -o "uce_[0-9]+_"|sort |uniq -c|sort -n|awk -v "threshold=2" '$1>threshold { print $2 }' > blacklisted_markers.2.GABBI.txt
mkdir aTRAM_results_bl_2
find aTRAM_results/ -type f|egrep -o "uce_[0-9]+"|sort -u|grep -v -f <(sed -E "s/_$/\$/g" blacklisted_markers.2.GABBI.txt )|while read uce;do
    cat aTRAM_results/$uce.fasta > aTRAM_results_bl_2/$uce.fasta
done
```

The resulting folder is now ready to follow phylogenomic analyses. You can run _Phylomera_ using the reference sequences computed by GABBI in ```GABBI_output/06_final_targeted_loci/cactus_alignment.final.anc.loci.cons.fasta``` and your [annotated file](#annotate-targeted-loci-as-coding-or-non-coding-sequences-for-downstream-analyses) if you have one.

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
> This command will clean alignments, split core and flanking regions and run gene trees with IQ-TREE using ModelFinder to select the best fitted evolutionary model. To run a concatenated analysis, run the same command, replacing ```-g MFP -p 0``` by, for example, ```-s MFP+MERGE``` (if you want to chose the minimum percentage of taxa required to keep a marker interactively) or ```-s MFP+MERGE -p 70``` (for the typically utilised threshold of 70% spp.)

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
