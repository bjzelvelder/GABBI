# GABBI pipeline version 1.0.0
> _Genome Alignment-Based Bait Inference_


**GABBI** is a fully automated pipeline to design target capture baits (or probes) using a **whole-genome alignment**.
GABBI-derived probe sets are expected to target more variable loci than usual probe design methods that rely on a base genome to map reads, providing more sensitive and phylogenetically resolutive data.

This pipeline was designed in this publication to produce the first set of weevil-specific probes. If you use GABBI, please cite:
> Zelvelder B, ... GABBI: A new method based on genome alignments provides a highly resolutive target enrichment set for weevils (Coleoptera, Curculionoidea), ... doi:X

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
4. [Detailed options](#detailed-options)
5. [Citation](#citation)
6. [References](#references)

---
# Installation
The GABBI pipeline uses **singularity** to run in a pre-built environment (a so-called _singularity image_) to avoid installation issues and incompatibilities. Thus, all you need is a linux environment with ```apptainer``` version >= 1.0.0 or ```singularity``` version >= 3.0.0 installed. This software can be installed following this tutorial: [Installing apptainer](https://apptainer.org/docs/admin/main/installation.html). For clarity, I will always refer to ```singularity``` rather than ```apptainer``` but all commands that start with ```singularity ...``` can be spelled ```apptainer ...```. 

GABBI is distributed as a Singularity image and has no external dependencies beyond Singularity/Apptainer itself. On HPC clusters managed by SLURM, Singularity is typically available as an environment module that can be searched throught ```module search singularity``` and can be loaded with ```module load <singularity module>``` prior to execution.

## Downloading the GABBI singularity image
Because a singularity image is too heavy to be stored and pulled from github, you will not find the ```.sif``` file in this repository. However, the [GABBI image](https://cloud.sylabs.io/library/bjzelvelder/pipeline/gabbi) is stored on the **Sylabs** remote repository and can simply be pulled with this command:
```
singularity pull --arch amd64 library://bjzelvelder/pipeline/gabbi:v1.0.0
```
You should then be able to see the help section with:
```
singularity run-help gabbi_v1.0.0.sif
singularity run gabbi_v1.0.0.sif --help
```

## Building from 'source'
If previous commands doesn't work or that you want make local modifications to the GABBI scripts, you can build the singularity image yourself using the definition file ```.def```. 
First, get GABBI files by cloning the GABBI repository with:
```
git clone https://github.com/bjzelvelder/GABBI.git
cd GABBI
```
At that point, you can build the GABBI singularity image (which should take about 40 minutes) if you only want to build the image and run GABBI. You can also build a GABBI sandbox to run the pipeline interactively (for debugging purposes).

### Building the GABBI singularity image
Inside GABBI's repository, simply execute:
```
singularity build --fakeroot gabbi_v1.0.0.sif GABBI_v1.0.0.def
```
This will create a ```.sif``` image that can be run the same way with:
```
singularity run GABBI_v1.0.0.sif --help
```

### Building a GABBI sandbox
Alternatively, if you want to make modifications to the GABBI scripts or simply run the pipeline more interactively from within the singularity image (environment), you can make a GABBI sandbox with a writable fakeroot that will be stored in ```GABBI_sandbox/```:
```
singularity build GABBI_sandbox/ GABBI_v1.0.0.def
mkdir -p GABBI_sandbox/$PWD
singularity shell -B $PWD --fakeroot --writable GABBI_sandbox/
```
Once you run singularity shell, you are inside the singularity image. Please be aware that binding your local path with ```-B $PWD``` allows GABBI to access your local files, but you can still make changes to them **even if you are inside the singularity image**.
You can now run the GABBI pipeline by launching the ```main.sh``` script:
```
source /opt/gabbi/main.sh --help
```
Note: because of the ```set -eou pipefail``` command inside main.sh, errors might eject you from the singularity image. To avoid having to connect back to the sandbox repeatidily, run this command after sourcing ```main.sh```:
```
set +eou pipefail
```
---
# Pipeline overview
The GABBI pipeline is a fully automated pipeline to design target capture baits (or probes) using a **whole-genome alignment**. As opposed to other commonly used methods of probe design, GABBI doesn't rely on a **base taxon**. The final set of targeted loci tend to be more variable than ultraconserved elements (UCE) and are referred to as shared homologous regions (SHR). Further details on the context and signficance of this new probe design method can be accessed through the publication that features the GABBI pipeline ([Zelvelder et al. 202X](#citation)).

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

The goal of this pipeline is to allow anyone wishing to make a set of specific target capture baits as straightforward as possible, while still being versatile to user specifics. For this reason, you only need to provide **chromosome-level genomes** and a corresponding **guide tree** (see [preparing input data](#preparing-input-data) section) for the genome alignment (although high quality, scaffold-level genomes might work fine for our purpose), and an additionnal set of **whole-genome assemblies**. Each step of the pipeline is checkpointed to save time and can be restarted to fine-tune the probe set with specific thresholds that certainly depend on each dataset specifics (see [detailed options](#detailed-options)).

> [!IMPORTANT]
> When provided 48 cpu cores, 4 chromosome-level genomes and 7 additional genomes of weevils, the GABBI pipeline took roughly 13 hours and 30 minutes to run in total and produced 20 Gb of data. As it is strongly paralellized, we strongly recommend giving as much cpu cores as possible to GABBI to reduce computational time even further for bigger datasets.

---
# How to use the GABBI pipeline
Once you have fetched the GABBI singularity image and succesfully ran the help command, GABBI is ready to run. To help you prepare input data and reading GABBI output, this section will guide you through each option in greater details using a small examplar dataset available in ```example_data```. This examplar dataset contains 4 chromosome-level genomes and 7 additional genomes of a small sample of weevils that belong to the Curculioninae subfamily (sensu 'CCCMS'; [Haran et al. 2023](#references)).

## Preparing input data
GABBI only requires three arguments to run: ```--chr-genomes```, ```--guide-tree```, and ```--add-genomes```. Thus, we need to [download](#downloading-ncbi-genomes) and [organize](#organizing-directory-architecture) the genomes that will be used for the probe design and provide a phylogenetic tree of the chromosome-level genomes to [guide the genomic alignment](#getting-a-guide-tree). But first, there are a few things to note on how to [choose those genomes](#choosing-representative-taxa).

### Choosing representative taxa
Ultimately, the goal of a probe set is to efficiently hybridize with the fragmented DNA of any taxon belonging to the targeted taxonomic group. Thus, we want as much as possible that our probe sequences map the actual variation of each targeted locus. As we cannot represent the entire diversity of our taxonomic target in the probe design (otherwise why would we even bother with probe design?), we have to rely on a drastic subsample of its diveristy by a handful of representative taxa. Historically, very few genomes were used to represent large taxonomic scales (e.g. only 6 taxa represent the 400k species of the Coleoptera UCE probe set; [Faircloth 2017](#references)), but thanks to the growing number of whole-genomes available, we can greatly improve this number. Nonetheless, here are a few things you should have in mind:
- Avoid the taxonomic redundancy in representative genomes. Some genus/tribe are often overrepresented among available genomes whereas entire families can be missing. Reducing the taxonomic redundancy in the genome set should prevent selected markers to appear more shared than expected.
- Be aware that providing too many genomes can rapidly become computationnaly intensive. Further research is needed to quantify an optimum, but providing one or two taxa per taxonomic tribe should be enough to cover their genetic variation.
- Genomic alignment should be conducted on chromosome-level (or very high quality) genomes and validated with a separate set of additional genomes. But testing the final probe set with a third set of genomes/raw data can also be useful!

With that in mind, you can check available genomes at [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/) and play with filters and column selection to fine-tune your selection.

### Downloading NCBI genomes
Once you have chosen the genomes you want to represent in your probe design, you can either download them manually or use the ```download_genomes_ncbi.sh``` script to download them. In this case, you must install the ```ncbi_datasets``` command-line program [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/). If you have conda installed, simply run:
```
conda install -c conda-forge ncbi-datasets-cli
```
Then, download the NCBI table of your active selection, making sure that you have **selected the Assembly Level column** (to let the script organize genomes accordingly):
<p align="center">
  <img src="/image/screenshot_ncbi.png" alt="NCBI_screenshot" width="600"/>
</p>

To download genomes of your downloaded table, give it as an argument to the script:

```
./download_genomes_ncbi.sh [path_to_ncbi_table]
```

### Organizing directory architecture
As you can tell by the ```example_data``` directory architecture, each genome must be in its own subfolder named after the taxon name you want to keep during the analysis. The ```download_genomes_ncbi.sh``` script should have organized your genome selection in this way, but you can obvisouly add your own genomic alignments by following the same logic (simply make sure that they have the ```.fasta``` or ```.fna``` file extension). You can also transfer chromosome level genomes to the ```additional_genomes/``` folder if you don't want them to be part of the genome alignment step (see [Choosing representative taxa](#choosing-representative-taxa).

High quality, chromosome-level genomes are strongly recommended for cactus whole-genome alignments, so they need to be stored in a seperate folder as other genomes to be parsed through the ```--chr-level-genomes``` option. Here is how your directory architecture should look like:

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
for i in chr_level_genomes/* additional_genomes/*; do mv $i ${i%_*} ;done
```

> [!NOTE]
> Additional genomes parsed through the ```--add-genomes``` option can be of any quality, but low-quality genomes can impact SHR recovery and reduce the total number of targeted loci. If you have some in your dataset, you might consider lowering the final [SHR threshold](#detailed-options).

### Getting a guide tree
The whole-genome alignment step requires a guide tree to run. This phylogenetic tree must be provided as a **NEWICK** file with **matching taxon names** with the chromosome level genomes provided in the ```chr_level_genomes/``` folder. You can check the ```chr_level_genomes.nw``` format for reference.

Because a weevil phylogeny was already available for the example dataset, I simply pruned and renamed an existing tree to only contain the taxa I provided in my dataset using [FigTree](https://github.com/rambaut/figtree).

However, if you have no clue on how your phylogenetic tree should look like, you may try to run the [ROADIES](https://github.com/TurakhiaLab/ROADIES) pipeline to compute a fast phylogenetic tree of your dataset ([Gupta et al. 2025](#references)).

> _This feature might be added in the GABBI pipeilne in the future._

## Running the pipeline
Once you have prepared your input data, you can run the GABBI pipeline. Here are a few useful command-lines you can run based on the ```example_data``` provided:

  - Run the pipeline with minimal options and default arguments, printing logs to STDOUT:
```
singularity run gabbi_v1.0.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw
```


  - Personnalize basic pipeline outputs and logs:
```
singularity run gabbi_v1.0.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --threads 64 --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output > GABBI_Curculioninae_v1.log 2>&1
```


  - Run GABBI on a Slurm cluster, increasing default cactus ressources:
```
singularity run gabbi_v1.0.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --threads 128 --debug --cactus-slurm --cactus-maxCores 32G --cactus-maxDisk 300G --cactus-maxMemory 1000G \
    --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output
```

  
  - Run the GABBI pipeline and stop before _in silico_ capture on additional genomes to check temporary probes:
```
singularity run gabbi_v1.0.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --threads 64 --stop-before 05_adding_genomes --prefix GABBI_Curculioninae_v1 --out-dir GABBI_Curculioninae_v1_output > GABBI_Curculioninae_v1.log 2> /dev/null
```
> With this command, you can provide an empty ```additional_genomes/``` folder to start the GABBI pipeline and add your additional genomes later by simply running the pipeline again without the ```--stop-before``` option


  - Run the GABBI pipeline on a reduced dataset, increasing default stringency thresholds:
```
singularity run gabbi_v1.0.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --temp-tax-threshold 100 --temp-allow-dupes 0 --shr-threshold 80 \
    --prefix GABBI_Curculioninae_v2 --out-dir GABBI_Curculioninae_v2_output > GABBI_Curculioninae_v2.log
```


  - Change the final SHR threshold based on the multifasta table to increase the final number of targeted loci, without having to rerun the entire pipeline:
```
singularity run gabbi_v1.0.0.sif --chr-genomes chr_level_genomes --add-genomes additional_genomes --guide-tree chr_level_genomes.nw \
    --temp-tax-threshold 100 --temp-allow-dupes 0 --shr-threshold 70 --restart 5.5_final_phyluce_probes \
    --prefix GABBI_Curculioninae_v2 --out-dir GABBI_Curculioninae_v2_output > GABBI_Curculioninae_v2.log
```

Each step is checkpointed if it succesfully ran, so running the same command again will resume the pipeline where it stopped.
If one step fails, the pipeline stops and will restart where it failed running the same command. Note that sometimes, issues are not correctly caught by the checkpointing system, so you might have to restart from an anterior step once you found the issue. Please report issues in the corresponding github sections to help me improve the pipeline!

>[!IMPORTANT]
> If you don't manage to get cactus running, you may consider running it separately to fine tune cactus options or run it step by step by referring to [cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md).
> In this case, once you have your ```.hal``` file, you can copy it back in the ```GABBI_output/01_cactus_alignment/``` folder, make sur that it is named ```<prefix>.hal``` as in the ```--prefix <prefix>``` option, then manually add a checkpoint for cactus with ```touch GABBI_output/.gabbi_checkpoints/step1.1_cactus_alignment```

## Reading GABBI outputs

GABBI outputs and temporary files are all stored in the ```--output-dir``` provided or in ```GABBI_output/``` by default. They are organized so that each phase of the pipeline has its own subdirectory named after the phase name (e.g. 01_cactus_alignment, 02_conserved_loci etc.).

Files and statistics of interest are given through the pipeline. They can be accessed in your log file using:

```
grep "GABBI" [log_file]
```

Here is a more detailed list of some GABBI outputs (generated on the example dataset with the default ```cactus_alignment``` prefix).
- ```01_cactus_alignment/```
  - ```cactus_input.txt```: Input file created using ```--chr-level-genomes``` content and the provided ```--guide-tree```.
  - ```cactus_alignment.hal```: Cactus whole-genome alignment.
- ```02_conserved_loci/```
  - ```maf/*.maf.gz```: Whole-genome alignment converted into a "Multiple Alignment Format" (MAF), using each genome as a base genome.
  - ```maffilter/```: Each genome has its own subdirectory containing filtered blocks of alignments of each of their chromosome/scaffolds (based on ```--block-size``` and ```--block-length``` options).
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

# Future directions
- Run a species tree /gene trees easily using phylomera_v0.8.3.sh on targeted loci alignments to check the resulting phylogeny
- Run in silico capture on simulated target capture data using art illumina, bbmap and aTRAM or PHYLUCE
- Annotate your probe set using annotated genomes
- Synthesize your probe set (remove redundant sequences)?

---
# Detailed options
```
    GABBI — Genome Alignment-Based Bait Inference pipeline
    ======================================================

    Required arguments
    ------------------
      --chr-genomes   DIR   Directory containing chromosome-level genome assemblies to be aligned with Cactus.
                            Each genome should reside in its own sub-directory, named <Taxon_name>
      --guide-tree    FILE  Newick-format species tree used as a guide topology by Cactus
      --add-genomes   DIR   Directory containing additional genome assemblies to validate SHR loci for final
                            probe set. Each genome should reside in its own sub-directory, named with
                            a unique <Taxon_name> 

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
      --debug               Print additional diagnostic information to stderr during execution.
                            Useful for development and troubleshooting [off]
      --threads, -t   INT   Number of CPU threads to allocate [all available by nproc command]
      --help, -h            Shows this message and exits

      Genome alignment options
      ------------------------
      --cactus-maxDisk  INT   Increase the maximum amount of disk used by cactus, with K, M or G suffix
                              to specify Kilo, Mega or Gigabytes. [50G]
      --cactus-maxCores INT   Maximum number of cpu used by cactus (too many might cause a cactus to crash) [32]
      --block-size      INT   Minimum number of taxa (including ancestral genomes) required to retain
                              an alignment block [70 % of extent and ancestral genomes]
      --block-length    INT   Minimum length (nt) required to retain an alignment block [200]

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
If you used the **GABBI** pipeline, please cite:
> Zelvelder B, ... GABBI: A new method based on genome alignments provides a highly resolutive target enrichment set for weevils (Coleoptera, Curculionoidea), ... doi:X

# References
References cited on this page:
> Faircloth BC. 2017 Identifying conserved genomic elements and designing universal bait sets to enrich them. Methods Ecol Evol 8, 1103–1112. [doi:10.1111/2041-210X.12754](https://doi.org/10.1111/2041-210X.12754)
> 
> Gupta A, Mirarab S, Turakhia Y. 2025 Accurate, scalable, and fully automated inference of species trees from raw genome assemblies using ROADIES. Proc. Natl. Acad. Sci. U.S.A. 122, e2500553122. [doi:10.1073/pnas.2500553122](https://doi.org/10.1073/pnas.2500553122)
> 
> Haran J et al. 2023 Phylogenomics illuminates the phylogeny of flower weevils (Curculioninae) and reveals ten independent origins of brood-site pollination mutualism in true weevils. Proc. R. Soc. B. 290, 20230889. [doi:10.1098/rspb.2023.0889](https://doi.org/10.1098/rspb.2023.0889)

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
