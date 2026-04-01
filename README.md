# GABBI pipeline version 1.0.0
Genome Alignment Based Bait Inference (GABBI)


**GABBI** is a fully automated pipeline to design target capture **probes** (or baits) based on **whole-genome alignments**

# Requirements
To run this pipeline, all you need is a linux environment with ```apptainer``` version >= 1.0.0 or ```singularity``` version >= 3.8.0 installed (see [Installation](#installation)). This pipeline was designed on Ubuntu 22.04 and mainly written in shell.

# Table of contents

# Installation
The GABBI pipeline uses **singularity** to run in a pre-built environment (a so-called _singularity image_) to avoid installation issues and incompatibilities. Thus, all you need is a linux environment with ```apptainer``` version >= 1.0.0 or ```singularity``` version >= 3.0.0 installed. This software can be installed following this tutorial: [Installing apptainer](https://apptainer.org/docs/admin/main/installation.html). For clarity, I will always refer to ```singularity``` rather than ```apptainer``` but all commands that start with ```singularity ...``` can be spelled ```apptainer ...```. 

## Downloading the GABBI singularity image
Because a singularity image is too heavy to be stored on and pulled from github, you will not find the ```.sif``` file in this repository. However, the GABBI image can simply be pulled from the **Sylabs remote repository** with this command:

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
singularity build --fakeroot GABBI_v1.0.0.sif GABBI_v1.0.0.def
```
This will create a ```.sif``` image that can be run the same way with:
```
singularity run-help GABBI_v1.0.0.sif
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

# How to use the GABBI pipeline




# Detailed options
```shell
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
                                    05_final_phyluce_probes, 06_final_targeted_loci [none]
      --debug               Print additional diagnostic information to stderr during execution.
                            Useful for development and troubleshooting [off]
      --threads, -t   INT   Number of CPU threads to allocate [all available by nproc command]
      --help, -h            Shows this message and exits

      Genome alignment options
      ------------------------
      --cactus-maxDisk  INT   Increase the maximum amount of disk used by cactus, with K, M or G suffix
                              to specify Kilo, Mega or Gigabytes. [300G]
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

- 

# References
If you used the **GABBI** pipeline, please consider citing these articles as well.

- 
- 
