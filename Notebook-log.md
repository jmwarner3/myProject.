#This is my Project notebook:2/10
This commit should be to my personal repo and includes the same info as the previous commit, just a different location. 
Data will be given to me from an utside source/lab in order to create a phylogeny for my final project

#ClustalW Notes:2/20
- CLustalW is a Multiple Sequence ALignment (MSA) software that uses 4 main steps in order to determine phylogenetic relationships and distances
- The ClustalW alignment method can be broken down into three main stages:

~Pairwise Alignment and Distance Matrix Calculation: Initially, all pairs of sequences are aligned separately to calculate a distance matrix, which reflects the divergence between each pair of sequences.

~Guide Tree Construction: A guide tree is constructed from the distance matrix. This tree serves as a guide for the final multiple alignment process, with branch lengths proportional to the estimated divergence along each branch.

~Progressive Alignment: Sequences are progressively aligned according to the branching order in the guide tree. Larger groups of sequences are aligned following the branching order. The alignment process involves using dynamic programming algorithms with residue weight matrices and penalties for gap opening and extension.

- ClustalW uses .fasta files as data(input) and provides -tree(NJ Tree) -PIM(output percent identity matrix) and -ALIGN(full multiple alignment) as the output

-#Clustal Command Notes:
To do things:
       -options
           List the command line parameters.

       -help or -check
           Outline the command line params.

       -fullhelp
           Output full help content.

       -align
           Do full multiple alignment.

       -tree
           Calculate NJ tree.

       -pim
           Output percent identity matrix (while calculating the tree).

        -bootstrap=n
           Bootstrap a NJ tree (n= number of bootstraps; def. = 1000).

       -convert
           Output the input sequences in a different file format.

#CLustalW Steps for Future Use and Alignment HW:2/20
1.) Obtain data in a .FASTA file format (you can use clustal omega if not provided in a fasta format)
2.) Provide full path to CLustalW software using COMMAND: $ /usr/local/bin/clustalw2 - infile=input.fasta -tree -pim -type=DNA -case=upper

DATA (sequences)

-INFILE=file.fasta                             :input sequences.
-PROFILE1=file.fasta  and  -PROFILE2=file.fasta  :profiles (old alignment).


                VERBS (do things)

-OPTIONS            :list the command line parameters
-HELP  or -CHECK    :outline the command line params.
-FULLHELP           :output full help content.
-ALIGN              :do full multiple alignment.
-TREE               :calculate NJ tree.
-PIM                :output percent identity matrix (while calculating the tree)
-BOOTSTRAP(=n)      :bootstrap a NJ tree (n= number of bootstraps; def. = 1000).
-CONVERT            :output the input sequences in a different file format.

                PARAMETERS (set things)

***General settings:****
-INTERACTIVE :read command line, then enter normal interactive menus
-QUICKTREE   :use FAST algorithm for the alignment guide tree
-TYPE=       :PROTEIN or DNA sequences
-NEGATIVE    :protein alignment with negative values in matrix
-OUTFILE=    :sequence alignment file name
-OUTPUT=     :GCG, GDE, PHYLIP, PIR or NEXUS
-OUTORDER=   :INPUT or ALIGNED
-CASE        :LOWER or UPPER (for GDE output only)
-SEQNOS=     :OFF or ON (for Clustal output only)
-SEQNO_RANGE=:OFF or ON (NEW: for all output formats)
-RANGE=m,n   :sequence range to write starting m to m+n
-MAXSEQLEN=n :maximum allowed input sequence length
-QUIET       :Reduce console output to minimum
-STATS=      :Log some alignents statistics to file

***Multiple Alignments:***
-NEWTREE=      :file for new guide tree
-USETREE=      :file for old guide tree
-MATRIX=       :Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename
-DNAMATRIX=    :DNA weight matrix=IUB, CLUSTALW or filename
-GAPOPEN=f     :gap opening penalty        
-GAPEXT=f      :gap extension penalty
-ENDGAPS       :no end gap separation pen. 
-GAPDIST=n     :gap separation pen. range
-NOPGAP        :residue-specific gaps off  
-NOHGAP        :hydrophilic gaps off
-HGAPRESIDUES= :list hydrophilic res.    
-MAXDIV=n      :% ident. for delay
-TYPE=         :PROTEIN or DNA
-TRANSWEIGHT=f :transitions weighting
-ITERATION=    :NONE or TREE or ALIGNMENT
-NUMITER=n     :maximum number of iterations to perform
-NOWEIGHTS     :disable sequence weighting

- These steps would provide a NJ tree using -tree command, the least costly alignment of sequences -ALIGN and the PIM using -pim

#Distance & Parsimony Codes:
note: run these commands in R studio 


install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

install.packages("igraph", type = "binary")

library(stats)
library(ade4)
library(ape)

library(adegenet)
library(phangorn)

#DistanceBased
dna <-dna <- fasta2DNAbin(file=Users/jamiewarner/Documents/PCancerVariantData.fasta)


D <- dist.dna(dna, model="TN93")

tre <- nj(D)

tre <- ladderize(tre)

plot(tre, cex=.6)
title("A simple NJ tree")

#ParsimonyBased
dna <- fasta2DNAbin(file=Users/jamiewarner/Documents/PCancerVariantData.fasta)
dna2 <- as.phyDat(dna)

tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)


plot(tre.pars, cex=0.6)


# Package installation

## R packages
This analysis requires `tools`, `parallel`, `data.table`, `devtools`, `phyloRNA` and `beter` packages.
Packages `tools` and `parallel` are part of standard install of R do not need to be installed manually.
Packages `data.table` and `devtools` are available from CRAN. First, start `R` console by typing `R` into terminal:

```{bash}
R
```

and now type:

```{r}
install.packages("data.table")
install.packages("devtools")
```

Packages might have additional system requirements that will not be automatically installed by R. You would need to find them in the error log and install them manually. On `ubuntu`, start terminal and run:
```{bash}
sudo apt install name-of-library
```
Alternatively, you could install these packages with the `conda`, start terminal and type:
```{bash}
conda install r-devtools r-data.table
```
This should automatically install all required system libraries.

### phyloRNA and beter
These packages are not available on CRAN. You can however easily install them directly from github using the `devtools` package. Start the `R` console and type:

```{r}
devtools::install_github("biods/phyloRNA")
devtools::install_github("biods/beter")
```

this will automatically install these packages and required dependencies.

## Python packages
To install the `pysam` package, simply type into the terminal:
```{bash}
pip install pysam
```
# Required files

```

Now type:
```{bash}
cd phyloRNAanalysis
```
to enter the `phyloRNAanalysis` folder.

## Reference genome
First, create `reference` folder in the `phyloRNAanalysis` directory:

```{bash}
mkdir reference
cd reference
```
Now, download and unzip the [GRCh38v15](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) reference genome together with the associated index file:
```{bash}
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```
download and unzip the [annotation](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz):

```{bash}
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
```
as well as set of common variants (no need to unzip these):
```{bash}
wget ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
```

### Filter annotation file:
There is just one more thing. The newest version of `Cellranger` (version 5) does not accept an annotation file with extra scaffolds. This can be fixed by filtering the GTF file according to the fasta index file using the `phyloRNA` package:
```{bash}
Rscript -e 'phyloRNA::filter_gtf("GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf", "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai")'
```

## Data
To download the [GSE163210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163210) data, navigate back to the `phyloRNAanalysis`, create a `data` directory and download the BAM files:
```{bash}
cd ..
mkdir data
```
# Software installation
The analysis requires [R](https://cran.r-project.org/), [python3](https://www.python.org/), [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger), [bamtofastq](https://github.com/10XGenomics/bamtofastq), [GATK](https://gatk.broadinstitute.org/hc/en-us), [VCFtools](https://vcftools.github.io/), [IQtree](http://www.iqtree.org/) and [BEAST2](https://www.beast2.org/)

Required software can be installed independently, but we suggest the [conda](https://docs.conda.io/en/latest/) package manager. `Cellranger` and `bamtofastq` are not available through the `conda` package manager and need to be installed separately.

## Conda
To install the `conda` package manager, either follow the instructions bellow or, for more information, instructions on the conda [website](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).

To download `miniconda` (minimal installation of `conda`), open your terminal and type the following:

```{bash}
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

This will download the `miniconda` into your current directory. Now type:
```{bash}
bash Miniconda3-latest-Linux-x86_64.sh
```

and follow the installation instructions in your terminal. This will install `conda` into your computer. To verify that `conda` is installed, open new terminal window (or type `exec bash`) and type:

```{bash}
conda --help
```

You can now install required software by typing:

```{bash}
conda install r-base gatk4 vcftools iqtree beast2
```

There is no need to install `python3` as it is a part of the `conda` installation.

## Cellranger
Cellranger is a proprietary software and its not available through the `conda` manager.

Visit the `Cellranger` [homepage](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and fill in required information. You will be provided with a download link and [detailed installation instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).

Download provided link with `curl` or `wget`. You will download a compressed archive `cellranger-5.0.1.tar.gz` (or a newer version). To unpack this archive, type:
```{bash}
tar -xzvf cellranger-5.0.1.tar.gz
```

This will unpack the archive in your current directory. Personally, You can now remove `cellranger-5.0.1.tar.gz` by typing:

```{bash}
rm cellranger-5.0.1.tar.gz
```

Now open `.bashrc` or `.profile` in your favourite editor and add `export PATH=~/cellranger-5.0.1:$PATH` or type:

```
echo 'export PATH=~/cellranger-5.0.1:$PATH' >> .bashrc
```

Restart your terminal by typing `exec bash` and type:

```
cellranger --help
```
to verify the installation.

## bamtofastq
To install `bamtofastq` visit the project [github](https://github.com/10XGenomics/bamtofastq) page or follow these instructions:

To download `bamtofastq`, type:

```{bash}
wget -O bamtofastq https://github.com/10XGenomics/bamtofastq/releases/download/v1.3.5/bamtofastq_linux
```

Now you need to add execute (x) permission to `bamtofastq`:

```{bash}
chmod +x bamtofastq
```

Finally, add `bamtofastq` to your PATH like you did with Cellranger:

```{bash}
echo 'export PATH=~/cellranger-5.0.1:$PATH' >> .bashrc
```

Restart your terminal (`exec bash`) and type:

```{bash}
bamtofastq --help
```
to verify the installation.

# Running the analysis
Make sure you have installed required [software](required_software.md), [packages](packages.md) and downloaded necessary [files](required_files.md).

If you want to replicate the analysis as published, navigate to the `phyloRNAanalysis` folder and type:
```{bash}
Rscript run.r
```

This will remap BAM files to reference genome, detect expression levels and SNV, prepare FASTAs and perform phylogenetic reconstruction.

## Running on your own data
If you want to use this for your own data, put your data into the `data` directory and delete the old data.
Then modify the `chemistry`, `densities`, `hdi` and `selection` to suit your needs.

* `chemistry` -- While the automatic detection of chemistry is preferred, `cellranger` might fail to detect chemistry for low-quality data, for this reason, the chemistry was fixed in the analysis. Either change it to `auto` or to chemistry of your data.
* `densities` -- Filter data to a particular data density, set to a single value if you do not care about comparing phylogenies from different data densities.
* `hdi` -- Highest Density Interval method for discretization of expression values. **This method will only work for homogeneous population of cells**. If your population has heterogeneous expression levels, this method of discretization will not work.
* `selection` -- Named vector to select the best performing cells from each sample for the alternative filtering method used in the study.

#' beast.r
#'
#' Functions for phylogenetic reconstruction with beast
import::here("parallel", "mcMap")
import::here("phyloRNA", "mkdir", "corename")
import::here("beter", "process_template")

#' run.r
#'
#' Run the analysis. This includes:
#'
#' Preparation
#' * remapping, demultiplexing, barcode correction and expression counts using Cellranger
#' * cleaning BAM files according to the GATK best practices
#' * adding a sample-specific postfix to cell barcodes
#'
#' Pre-processing:
#' * Expression
#'   -- standardization of genes into mu=0 and sd=1
#'   -- categorization according to empirical 60% and 90% HDI
#' * SNV:
#'   -- as bulk SNV identification and filtering with Mutect2
#'   -- sc SNV identification with vcm.py
#' * stepwise filtration into 20%, 50% and 90% density
#' * alternative filtration into 58 best cells and full dataset, 50% and 90% density
#'
#' Filtering:
#' * stepwise filtration into 20%, 50% and 90% density
#' * subset into 58 best cells and filtering into 50% and 90% density
#'
#' Phylogenetic analysis:
#' * ML with stepwise filtering
#' * ML and BI with alternative filtration
#' * BEAST templates created with the `beter` package
#' * Expression:
#'   -- IQtree: ORDINAL+ASC, ultrafast bootstrap -B 1000
#'   -- BEAST: ordinal from MM, exponential pop growth, coalescent prior, strict clock, two runs
#' * SNV:
#'   -- IQtree: GTR+gamma, ultrafast bootstrap -B 1000
#'   -- BEAST: GTR, exponential pop growth, coalescent prior, strict clock, two runs
#'

# use import::from instead?
import::from("src/prepare.r", "prepare_samples")
import::from("src/snv.r", "detect_snv", "filter_snv")
import::from("src/utils.r", "table2fasta", "fasta2stats")
import::from("src/expr.r", "preprocess_expression", "filter_expression")
import::from("src/beast.r", "beasts")
import::from("src/iqtree.r", "iqtrees")
import::from("phyloRNA", "corename")

main = function(){
    # datasets:
    bam = dir("data", full.names=TRUE)
    outdir = file.path("moravec2021")

    # required reference files:
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    refdir = "reference/ref" # shared cellranger ref
    vcf = "reference/00-common_all.vcf.gz"
    
    # PON and normal samples
    pon = "pon/MDA-MB-231/pon.vcf" # see make_panel_of_normals.r
    normal = file.path( # see make_normal_samples.r
        outdir, "normal",
        c("MDAMB231-MPmt-rep1.prepared.bam", "MDAMB231-MPmt-rep2.prepared.bam")
        )

    # Other settings:
    nthreads = 16
    chemistry = "SC3Pv2"
    densities = c(0.2, 0.5, 0.9)
    hdi = c(0.6, 0.9)
    selection = c("T1" = 20, "T3" = 20, "T2" = 6, "CTC1" = 6, "CTC2" = 6)

    # Preparation:
    prepared = prepare_samples(
        bam, reference, annotation, vcf,
        chemistry = chemistry, nthreads = nthreads,
        outdir = file.path(outdir, "prepare"),
        refdir = refdir
        )

    # SNV part
    vcm = detect_snv(
        bam = prepared$bam,
        barcodes = prepared$barcodes,
        reference = reference,
        normal = normal,
        pon = pon,
        outdir = file.path(outdir, "snv")
        )

    tabdir = file.path(outdir, "snv", "filtered")
    fastadir = file.path(outdir, "snv", "fasta")
    treedir = file.path(outdir, "snv", "tree")

    snv = filter_snv(vcm=vcm, density=densities, prefix="snv", outdir=tabdir)
    snv_fasta = table2fasta(snv, outdir=fastadir)
    fasta2stats(snv_fasta, unknown="N")

    snv_subset = filter_snv(vcm=vcm, selection=selection, prefix="snv_subset",
                            outdir=tabdir)
    snv_subset_fasta = table2fasta(snv_subset, outdir=fastadir)
    fasta2stats(snv_subset_fasta, unknown="N")

    iqtrees(
        c(snv_fasta, snv_subset_fasta),
        model = "TEST",
        bootstrap = 100, parallel = TRUE, nthreads = 16,
        outdir = file.path(treedir, "ML")
        )

    beasts(
        snv_subset_fasta,
        template = file.path("templates", "BDStrictGtr.xml"),
        outdir = file.path(treedir, "BI")
        )

    # expression part
    expr_preprocessed = preprocess_expression(
        h5 = prepared$h5,
        hdi = hdi,
        minGene=0,
        minUMI=0,
        outdir = file.path(outdir, "expr", "prepare"),
        prefix = "all"
        )

    filterdir = file.path(outdir, "expr", "filtered")
    fastadir = file.path(outdir, "expr", "fasta")
    treedir = file.path(outdir, "expr", "tree")
    
    expr = filter_expression(
        expr_preprocessed$discretized, prefix = "expr",
        outdir = filterdir, density = densities
        )
    expr_subset = filter_expression(
        expr_preprocessed$discretized, prefix = "expr_subset",
        outdir = filterdir, selection = selection
        )

    expr_fasta = table2fasta(expr, outdir=fastadir)
    fasta2stats(expr_fasta, unknown="-")

    expr_subset_fasta = table2fasta(expr_subset, outdir=fastadir)
    fasta2stats(expr_subset_fasta, unknown="-")
    expr_zero_fasta = table2fasta(
        expr_subset,
        file.path(fastadir, paste0(corename(expr_subset), "_zero.fasta")),
        outdir = fastadir, zero = "-"
        )

    iqtrees(
        expr_fasta,
        model = "ORDERED+ASC",
        outdir = file.path(treedir, "ML"),
        mc.cores = length(expr_fasta)
        )

    iqtrees(
        c(expr_subset_fasta, expr_zero_fasta),
        model = "ORDERED+ASC",
        bootstrap = 100, parallel = TRUE, nthreads = 16,
        outdir = file.path(treedir, "ML")
        )

    beasts(expr_subset_fasta, outdir = file.path(treedir, "BI"),
           template = file.path("templates", "BDStrictOrdinal.xml"))
    beasts(expr_zero_fasta, outdir = file.path(treedir, "BI"),
           template = file.path("templates", "BDStrictOrdinalZero.xml"))
    }


if(sys.nframe() == 0){
    main()
    }


#' Run BEAST analysis
#'
#' @param fasta a fasta file
#' @param template beter template to generate BEAST XML
#' @param outdir 
#' @param nthreads a number of threads to run the BEAST on
#' @param burnin burnin percentages
#' @param params an additional list of parameters
beast = function(fasta, template, outdir=NULL, nthreads=2, burnin=20, params=NULL){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    if(is.null(params))
        params = list()

    prefix = corename(fasta)

    beastxml = file.path(outdir, paste0(prefix, ".xml"))
    trace = file.path(outdir, paste0(prefix, ".trace"))
    trees = file.path(outdir, paste0(prefix, ".trees"))
    ess = file.path(outdir, paste0(prefix, ".txt"))
    tree = file.path(outdir, paste0(prefix, ".tree"))

    if(!file.exists(beastxml))
        process_template(
            template,
            beastxml,
            alignment = fasta,
            parameters = params
            )

    beast_args = c(
        "-threads", nthreads,
        basename(beastxml)
        )
    if(!file.exists(trace))
        phyloRNA:::systemE("beast", beast_args, dir=outdir)

    log_args = c(
        "-b", burnin,
        basename(trace),
        ">",
        basename(ess)
        )
    if(!file.exists(ess))
        phyloRNA:::systemE("loganalyser", log_args, dir=outdir)

    tree_args = c(
        "-b", burnin,
        "-lowMem",
        basename(trees),
        basename(tree)
        )

    if(!file.exists(tree))
        phyloRNA:::systemE("treeannotator", tree_args, dir=outdir)
    }


beasts = function(fasta, template, outdir=NULL, nthreads=2, burnin=20, params=NULL, mc.cores=1){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    if(is.null(params))
        params = list(list())

    mcMap(
        beast,
        fasta = fasta,
        template = template,
        outdir = file.path(outdir, corename(fasta)),
        nthreads = nthreads,
        burnin = burnin,
        params = params,
        mc.cores = mc.cores
        )
    }
#' expr.r
#'
#' Functions for processing 10x expression data
import::here("utils.r", "filename", "num2char", "mdensity", "write_table", "read_table")
import::here("filter.r", "density_filtering", "subset_filtering")
import::here("phyloRNA",
    "mkdir", "all_files_exist", "corename",
    "expr_merge", "expr_read10xh5",
    "expr_quality_filter", "expr_zero_to_na",
    "expr_normalize", "expr_scale", "expr_discretize",
    "remove_constant", "tab2seq", "write_fasta"
    )


#' Process 10X expression data
#'
#' This function will filter, scale, discretize and transform scRNAseq expression data
#' into a fasta sequence.
#'
#' The 10X expression count matrix can be provided either
#'
#' The function performs following steps:
#' * `phyloRNA::expr_quality_filtering()` -- to remove low-quality cells
#' * `phyloRNA::expr_zero_to_na()` -- transform 0 to NA
#' * `phyloRNA::expr_normalize()` -- normalize expression values accross genes
#' * `phyloRNA::expr_scale()` -- center and scale expression values accross cells
#' * `phyloRNA::expr_discretize()` -- discretize data according to HDI intervals
#' * `phyloRNA::densest_subset()` and `phyloRNA::remove_constant()` -- to filter the data matrix
#' to a chosen density
#' * `phyloRNA::fasta()` -- transform to FASTA sequence
#'
#' @param h5 a path to a expression count matrix stored in a .h5 file
#' or list of such paths
#' @param density **optional** required density or densities of the final matrix
#' @param hdi **optional** a highest density intervals for discretization
#' @param minGene **optional** minimum number of genes per cell
#' @param minUMI **optional** minimum number of UMI per cell
#' @param outdir **optional** an output directory
#' @param prefix **optional** a prefix for file names
#' @param normalize **optional** log-normalize the expression data
#'
#' @return a list of paths of all outputs
preprocess_expression = function(
    h5, hdi=c(0.6,0.9),
    minGene=250, minUMI=500,
    outdir=NULL, prefix=NULL,
    normalize=FALSE
    ){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "data"
    mkdir(outdir)

    result = list(
        intervals = file.path(outdir, paste(prefix, "intervals", "txt", sep=".")),
        discretized = file.path(outdir, paste(prefix, "discretized", "txt", sep="."))
        )

    if(all_files_exist(result))
        return(invisible(result))

    if(length(h5) > 1){
        names = corename(h5)
        data = lapply(h5, expr_read10xh5)
        data = expr_merge(data, names)
        } else {
        data = expr_read10xh5(h5)
        }

    data = process_expression(
        data, hdi = hdi,
        minGene = minGene, minUMI = minUMI,
        trim = FALSE, normalize = normalize,
        intervals = result$intervals
        )

    write_table(data, result$discretized)

    return(invisible(result))
    }


#' Process expression data
#'
#' This function simplifies standard expression data processing, such as quality filtering, scaling,
#' normalization and discretization.
#'
#' @param data an expression matrix
#' @param hdi **optional** a highest density intervals for discretization
#' @param minGene **optional** a minimum amount of represented genes per cell
#' @param minUMI **optional** a minimum amount of total UMI (or count) per cell
#' @param trim **optional** trim empty genes after filtering
#' @param normalize **optional** perform normalization after rescaling
#' @param intervals **optional** a file path to save discretization intervals into file
#' @param unknown **optional** a symbol representing unknown data
#' @return scaled, filtered and discretized count matrix
process_expression = function(
    data, hdi = c(0.6, 0.9),
    minGene = 0, minUMI = 0,
    trim = FALSE, normalize = FALSE,
    intervals = FALSE, unknown = "-"
    ){
    if(minGene > 0 || minUMI > 0 || trim)
        data = expr_quality_filter(data, minGene, minUMI, trim)

    data = expr_zero_to_na(data)

    if(normalize)
        data = expr_normalize(data)

    data = expr_scale(data)
    intervals = calculate_intervals(data, density=hdi, save=intervals)
    data = expr_discretize(data, intervals=intervals, unknown=unknown)
    data
    }


expr2fasta = function(x, fasta, unknown="-", summary=FALSE, process=TRUE, hdi=c(0.6, 0.9)){
    data = x
    if(process)
        data = process_expression(x, hdi, trim=TRUE, unknown=unknown)

    data = remove_constant(data, margin=1, unknown=unknown)
    seq = tab2seq(data, margin=2)

    write_fasta(seq, fasta)

    if(isTRUE(summary))
        summary = file.path(dirname(fasta), paste0(corename(fasta), "_summary.txt"))
    if(is.character(summary))
        count_matrix_summary(data, name=corename(fasta), file=summary)
    }


count_matrix_summary = function(data, name=NULL, file=NULL){
    text = paste0(
        "Sequences: ", ncol(data), "\n",
        "Sites: ", nrow(data), "\n",
        "Unique patterns: ", nrow(unique.matrix(data, MARGIN=1)), "\n",
        "Data density: ", mdensity(data, empty="-")
        )

    if(!is.null(name))
        text = paste0("Name: ", name[1], "\n", text)

    if(!is.null(file))
        writeLines(text, file)

    return(invisible(text))
    }


#' Calculate discretization intervals from data
#'
#' Calculate the discretization intervals empirically using the HDI method.
#'
#' @param data an input data
#' @param density a single or multiple values representing the total density covered by signle
#' or multiple the HDIs
#' @param save **optional** a TRUE or path to save the calculated interval values into a file
#' @return intervals corresponding to the HDI of the density vector
calculate_intervals = function(data, density=c(0.6,0.9), save=FALSE){
    dens = stats::density(data, na.rm=TRUE)

    intervals = lapply(density, function(x) phyloRNA::hdi(dens, 1-x))
    names(intervals) = as.character(density)

    if(isTRUE(save))
        save = "intervals.txt"
    if(is.character(save))
        dput(intervals, save)

    sort(unlist(intervals))
    }


#' Filter Expression dataset
#'
#' Filter expression dataset using two types of filtering approaches
#' see `src/filter.r` for more information
#'
#' @param expr discretized expression data
#' @param selection named list of vector specifying how many cells of each type should be selected
#' @param density desired data density
#' @param outdir an output directory
#' @return a list of filtered files
filter_expression = function(
    expr, prefix, selection=NULL, density=NULL, outdir=NULL
    ){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    if(is.null(selection) && is.null(density))
        stop("Either the selection or the density parameter must be specified.")

    if(!is.null(selection) && is.null(density)){
        file = filename(prefix, outdir=outdir)
        if(file.exists(file))
            return(invisible(file))

        data = read_table(expr)
        file = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            empty = "-",
            outdir = outdir
            )
        return(invisible(file))
        }

    if(!is.null(selection) && !is.null(density)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))
        data = read_table(expr)
        files = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            density = density,
            empty = "-",
            outdir = outdir,
            )
        return(invisible(files))
        }

    if(is.null(selection)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))

        data = read_table(expr)
        files = density_filtering(
            data,
            prefix = prefix,
            density = density,
            empty = "-",
            outdir = outdir
            )
        return(invisible(files))
        }
    }
#' filter.r
#'
#' Functions for the alternative filtering approach
import::here("phyloRNA", "remove_constant", "replace_ordinal", "densest_subset")
import::here("utils.r", "is_empty", "num2char", "write_table")
import::here("phyloRNA",
    "all_files_exist", "mkdir",
    "remove_constant", "replace_ordinal", "densest_subset"
    )

column_density = function(x, empty, sort=TRUE){
    cs = colSums(is_empty(x, empty))
    if(sort)
        cs = sort(cs, decreasing=TRUE)
    cs
    }


#' Divide vector into categories
#'
#' Divide vector into list of categories using the pattern and replace substitution.
#'
#' Pattern can be customized according to specific needs using the `[base::sub]` command.
#' The remainder after character substitution is then used as category.
#'
#' @param x a vector
#' @param pattern a pattern argument for the `[base::sub]` command
#' @param replace a replace argument for the `[base::sub]` command
#' @return a list of vectors for each category
divide_vector = function(x, pattern=".*-", replace=""){
    categories = sub(pattern, replace, x)
    result = list()
    for(category in unique(categories)){
        result[[category]] = x[categories == category]
        }
    result
    }

#' Select number of elements from list
#'
#' Take number of elements, specified by selection, from each category in the input list.
#'
#' @param x a list of vectors from the `divide` function
#' @param selection a list or named vector specifying how many elements should be selected
#' from each category
#' @return a list of the same structure as `x` with number of elements specified by `selection`
select_from_list = function(x, selection){
    result = list()
    categories = names(x)
    for(category in categories){
        result[[category]] = x[[category]][seq_len(selection[[category]])]
        }
    result
    }


subset_rows = function(x, k, empty){
    rs = rowSums( is_empty(x, empty) )
    y = x[rs > k, ]
    y = remove_constant(y, margin=1, empty)
    y
    }


#' Select columns
#'
#' Divide columns according to the pattern. Then select a number of columns according
#' to the selection vector with the least amount of missing data.
#'
#' Pattern a
#'
#' @param x a table
#' @param selection a named vector of selected columns for each subset
#' @param pattern pattern for the sub command
#' @param replace replace for the sub command
#' @param empty missing value specification
#' @return subsetted dataset
select = function(x, selection, pattern=".*-", replace="", empty=NA){
    # These operation work with sorted names
    # rather than whole input data
    print(x)
    print(selection)
    columns = names(column_density(x, empty, sort=TRUE))
    columns = divide_vector(columns, pattern, replace)
    columns = select_from_list(columns, selection)
    columns = unlist(columns)
    x = x[,columns]
    # Filter columns with only single value 
    subset_rows(x, 1, empty)
    }




#' Remove rows with no data
#'
#' Filter dataset by removing rows, down to a particular data density.
#' @param x dataset
#' @param density required density
#' @param empty missing value specification
densest_rows = function(x, density=0.5, empty=NA){
    k = 1
    while(TRUE){
        y = subset_rows(x, k, empty)

        # Filtering would result in a empty dataset
        if(nrow(y) == 0)
            break

        # Stop just bellow desired density
        if(mdensity(y, empty) > density)
            break
 
        x = y
        k = k + 1
        }
    x
    }


density_filenames = function(outdir, prefix, density){
    file.path(outdir, paste0(prefix, num2char(density), ".txt"))
    }


#' Filtering by finding the densest subset
#'
#' Filter data by finding the densest rows and colums.
#' See `[phyloRNA::densest_subset]` for more information.
#'
#' @param x data in tabular format
#' @param density desired data density
#' @param empty missing value specification
#' @param outdir **optional** an output directory, by defautl current working
#' directory is used
#' @param prefix **optional** a prefix for output files, by `filtered` is used
#' @param replace **optional** replace missing value with this character
#' @param rescale **optional** rescale the ordinal scale so that currently present
#' categories form a sequence, this this might be required by some computational
#' software. Note that this redefines the meaning behind categories.
#' If `TRUE`, a numeric sequence starting from `1` is used.
#' Alternatively, user can provide their own ordinal scale to which the object
#' will be rescaled.
#' @return vector of paths for filtered datasets
density_filtering = function(
    x, density=0.5, empty="N",
    outdir=NULL, prefix=NULL,
    replace=NULL, rescale=FALSE
    ){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "filtered"
    mkdir(outdir)

    outfile = density_filenames(outdir, prefix, density)

    if(all_files_exist(outfile))
        return(invisible(outfile))

    for(i in seq_along(density)){
        filtered = densest_subset(x, empty=empty, density=density[i])$result
        filtered = remove_constant(filtered, margin=1, unknown=empty)
        if(!is.null(replace))
            filtered = replace_missing(filtered, empty, replace)
        if(!is.null(rescale)){
            if(length(rescale) == 1 && rescale) # rescale = TRUE
                filtered = replace_ordinal(filtered)
            if(length(rescale) > 1) # e.g.: rescale = 0:10; rescale = letters
                filtered = replace_ordinal(filtered, rescale)
            }
        write_table(filtered, outfile[i])
        }

    invisible(outfile)
    }


replace_missing = function(data, missing, replace){
    replace = as.character(replace)
    if(length(replace) != 1)
        stop("ERROR: Provide exactly one character as a replacement")
    if(nchar(replace) != 1)
        stop("ERROR: Provide exactly one character as a replacement")
    if(is.null(missing) || is.null(replace))
        stop("ERROR: Missing and replace characters cannot be NULL")

    if(is.na(missing))
        data[is.null(data)] = replace
        
    if(!is.null(missing))
        data[data == missing] = replace

    data
    }


#' Filtering by selecting the best cells
#'
#' Filter data by first finding the best cells with the least amount of missing data.
#' Then the rows are filtered to reach desired density.
#'
#' @param x data in tabular format
#' @param selection named list or vector specifying how many cells should be selected of each type
#' or sample
#' @param density desired data density
#' @param empty missing value specification
#' @param outdir an output directory
#' @param prefix a prefix for output files
#' @return vector of paths for filtered datasets
subset_filtering = function(x, selection, density=NULL, empty="N", outdir=NULL, prefix=NULL){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "filtered"
    mkdir(outdir)


    if(is.null(density)){
        file = filename(prefix, outdir=outdir)
        if(file.exists(file))
            return(invisible(file))

        subset = select(x, selection, empty=empty)
        write_table(subset, file)
        return(invisible(file))
        }

    if(!is.null(density)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))
        subset = select(x, selection, empty=empty)
        Map(
            function(d,f) write_table(densest_row(subset, d, empty=empty), f),
            density, files
            )
        return(invisible(files))
        }
    }
#' iqtree.r
#'
#' Functions for phylogenetic reconstruction with iqtree
import::here("parallel", "mcMap", "mclapply")
import::here("phyloRNA", "corename", "mkdir")


#' Run IQtree analysis
#'
#' @param fasta a fasta file
#' @param outdir an output directory, fasta will be copied there
#' @param model a phylogenetic model for IQtree
#' @param ufboot a number of replicates for the ultrafast bootstrap
#' @param bootstrap a number of replicates for the standard bootstrap
#' @param nthreads a number of threads to run the IQtree on
iqtree = function(
    fasta,
    model = NULL, outdir = NULL,
    ufboot = FALSE, bootstrap = FALSE,
    nthreads = "AUTO", parallel = FALSE,
    remake = FALSE
    ){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    cfasta = file.path(outdir, basename(fasta))

    # Early exit if tree exists and not running parallel bootstrap:
    # (there the trees are more complicated)
    tree = paste0(cfasta, ".treefile")
    if(!remake && !parallel && file.exists(tree))
        return(invisible(tree))

    file.copy(fasta, cfasta)

    command = "iqtree"
    args = c(
        "-s", cfasta,
        "-nt", nthreads
        )
    if(!is.null(model))
        args = c(args, "--model", model)
    if(remake)
        args = c(args, "--redo")

    # Run iqtree with the ultrafast bootstrap
    if(ufboot){
        if(is.logical(ufboot)) ufboot = 1000
        args = c(args, "-B", ufboot)
        phyloRNA:::systemE(command, args)
        return(invisible(tree))
        }

    # Run iqtree with the ultrafast bootstrap, but don't parallelize
    if(bootstrap && !parallel){
        if(is.logical(bootstrap)) bootstrap = 100
        args = c(args, "-b", bootstrap)
        phyloRNA:::systemE(command, args)
        return(invisible(tree))
        }

    # run iqtree with parallelized standard bootstrap
    if(bootstrap && parallel){
        if(nthreads == "AUTO")
            nthreads = 4
        trees = iqtree_par(fasta, model, outdir, bootstrap, nthreads, remake)
        return(invisible(trees))
        }

    # run iqtree without boostrap:
    phyloRNA:::systemE(command, args)
    return(invisible(tree))
    }


#' Performs a single bootstrap run with iqtree
#'
#' @param i run id
#' @param fasta a path to fasta file
#' @param a model to be run with
iqboot = function(i, fasta, model = NULL, remake = FALSE){
    command = "iqtree"
    prefix = file.path(dirname(fasta), i)
    tree = paste0(prefix, ".boottrees")

    if(!remake && file.exists(tree))
        return(tree)

    args = c(
        "-s", fasta,
        "-nt", 1,
        "-bo 1",
        "-pre", prefix,
        "-quiet"
        )

    if(remake)
        args = c(args, "--redo")

    if(!is.null(model))
        args = c(args, "--model", model)
    system2(command, args)
    return(tree)
    }


#' Run iqtree with parallelized bootstrap
#'
#' @param fasta an input fasta file
#' @param model an iqtree model, such as HKY, GTR and others
#' @param outdir an output directory
#' @param bootstrap a number of bootstrap replicates
#' @param nthreads a number of threads to run in parallel
#' @return a vector of paths to best tree, bootstrap trees and conensus bootstrap tree
iqtree_par = function(
    fasta,
    model = NULL, outdir = NULL,
    bootstrap = 100, nthreads = 8,
    remake = FALSE
    ){
    command = "iqtree"

    if(is.null(outdir))
        outdir = "."

    bootdir = file.path(outdir, "bootstrap")
    mkdir(outdir)
    mkdir(bootdir)

    core = phyloRNA:::corename(fasta)
    cfasta = file.path(outdir, basename(fasta))
    cbfasta = file.path(bootdir, basename(fasta))

    file.copy(fasta, cfasta)
    file.copy(fasta, cbfasta)

    # Performs a single standard run single run
    tree = iqtree(fasta, model, outdir)

    # perform parallelized bootstrap
    trees = mclapply(
        seq_len(bootstrap),
        iqboot,
        model = model, fasta = cbfasta, remake = remake,
        mc.cores = nthreads
        )

    # merge bootstrap files
    boottrees = file.path(outdir, paste0(core, ".boottrees"))
    
    # prevent appending trees to an existing file
    if(file.exists(boottrees))
        file.remove(boottrees)

    file.append(boottrees, unlist(trees))

    # calculate consensus tree from bootstrap trees
    contree = iqtree_consensus(boottrees)

    # calculate bootstrap support for the best and consensus trees
    bootsup = iqtree_support(contree, boottrees, remake=remake)
    treesup = iqtree_support(tree, boottrees, remake=remake)

    return(c("tree"=treesup, "bootstrap"=boottrees, "consensus"=bootsup))
    }


#' Annotate tree with a boostrap support values
#'
#' Annotate tree with a bootstrap support values using a set of bootstrap trees.
#'
#' @param tree a tree that will be annotated
#' @param bootstrap a set of bootstrap trees
#' @param output **optional** an output file
#' @return a file with splits annotated with a bootstrap support
iqtree_support = function(tree, bootstrap, output=paste0(tree, ".sup"), remake=FALSE){
    if(!remake && file.exists(output))
        return(output)

    command = "iqtree"
    phyloRNA:::systemE(command, c("-sup", tree, "-t", bootstrap))

    annotated = paste0(bootstrap, ".suptree")
    file.rename(annotated, output)
    return(output)
    }


#' Calculate consensus tree from bootstrap trees
#'
#' @param bootstrap a file with bootstrap trees
#' @param output **optional** an output file
#' @return a file with conensus tree
iqtree_consensus = function(bootstrap, output=paste0(bootstrap, ".con"), remake=FALSE){
    if(!remake && file.exists(output))
        return(output)

    command = "iqtree"
    phyloRNA:::systemE(command, c("-con -t", bootstrap))

    consensus = paste0(bootstrap, ".contree")
    file.rename(consensus, output)
    return(output)
    }


iqtrees = function(
    fasta,
    model = NULL, outdir = NULL,
    ufboot = FALSE, bootstrap = FALSE,
    nthreads = "AUTO", parallel = FALSE,
    remake = FALSE, mc.cores=1
    ){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    n = length(fasta)

    if(is.null(model) || length(model) == 1)
        model = rep(list(model), n)

    if(length(outdir) == 1)
        outdir = file.path(outdir, corename(fasta))

    trees = mcMap(
        f=iqtree, fasta, model, outdir, ufboot, bootstrap, nthreads, parallel, remake,
        mc.cores=mc.cores
        )

    invisible(trees)
    }
#' prepare.r
#'
#' Function for preparation of sequences.
import::here("utils.r", "file_sub", "merge_files")
import::here("phyloRNA", "corename")


#' Prepare multiple samples into a single bam file
#'
#' Remap, clean, rename and merge individual sample bam files for further processing, suhc as
#' SNV detection.
#'
#' This function chains together basic processing steps to prepare a previously mapped bam file
#' for further processing. This include:
#'
#' * remapping the bam file to a new reference genome
#' * filtering bam so that only the most well-represented cells are included
#' * cleaning the bam file according to GATK best practices
#' * changing names of the barcodes (cell identifiers) to ensure their uniqueness
#' * merging sample files into a single file
#'
#' It is suggested to specify the `refdir` as it can be shared between analysis using
#' the same reference genome and annotation.
#'
#' @param bams a previously aligned bam files from the same organism that will be remapped
#' to a new reference and merged together
#' @param reference a new genome reference file (".fai")
#' @param annotation an annotation file associated with the reference genome
#' @param vcf a known variants associated with the reference genome that will be masked
#' @param outdir **optional** a general output directory
#' @param refdir **optional** an output directory for the cellranger reference genome files
#' @param chemistry **optiona** 10X chemistry, use only when the automatic detection is failing
#' @param nthreads **optional** a number of threads to use
#' @return a list with a paths to merged bam and barcode files.
prepare_samples = function(
    bams, reference, annotation, vcf,
    obam=NULL, obar=NULL, oh5=NULL,
    outdir=NULL, refdir=NULL,
    chemistry="auto", nthreads=16
    ){

    if(is.null(outdir))
        outdir = "prepare"
    if(is.null(refdir))
        refdir = file.path(outdir, "ref")
    if(is.null(obam))
        obam = file.path(outdir, "all.bam")
    if(is.null(obar))
        obar = file.path(outdir, "all.txt")
    if(is.null(oh5))
        oh5 = file.path(outdir, corename(bams), paste0(corename(bams), ".h5"))

    result = list(
        samples = corename(bams),
        bam = obam,
        barcodes = obar,
        h5 = oh5
        )

    if(file.exists(obam) && file.exists(obar) && all(file.exists(oh5)))
        return(result)

    results = list()
    for(bam in bams)
        results[[bam]] = prepare_sample(
            bam = bam,
            reference = reference,
            annotation = annotation,
            vcf = vcf,
            outdir = file.path(outdir, corename(bam)),
            refdir = refdir,
            chemistry = chemistry,
            nthreads = nthreads
            )

    pcores = unlist(lapply(results, getElement, "prefix"))
    pbams = unlist(lapply(results, getElement, "bam"))
    pbars = unlist(lapply(results, getElement, "barcodes"))
    ph5s = unlist(lapply(results, getElement, "h5"))

    # Defensive programming: Sanity check that the corenames are exactly the same.
    if(all(pcores != result$samples))
        stop("Corenames from samples differ from those in results. This shouldn't happen.")

    # merge prepared bams from all datasets:
    phyloRNA::gatk_MergeSamFiles(pbams, obam)

    # merge prepared barcodes from all datasets:
    merge_files(pbars, obar, overwrite=TRUE)

    return(result)
    }



#' Prepare sample bam file
#'
#' Remap, clean and rename bam file for further processing, such as SNV detection.
#'
#' This function chains together basic processing steps to prepare a previously mapped bam file
#' for further processing. This include:
#'
#' * remapping the bam file to a new reference genome
#' * filtering bam so that only the most well-represented cells are included
#' * cleaning the bam file according to GATK best practices
#' * changing names of the barcodes (cell identifiers) to ensure their uniqueness
#'
#' Function include a number of optional parameters to specify an output folder for each step.
#' `outdir` is a general output directory. If not specified, `prepare` is used. If other optional
#' directories are not specified, they default to a directory sitting in the outdir.
#' * `mapdir` defaults to `map`
#' * `refdir` defaults to `ref`
#' * `cleandir` defaults to `clean`
#'
#' It is suggested to specify the `refdir` as it can be shared between analysis using
#' the same reference genome and annotation.
#'
#' @param bam a previously aligned bam file that will be remapped to a new reference
#' @param reference a new genome reference file (".fai")
#' @param annotation an annotation file associated with the reference genome
#' @param vcf a known variants associated with the reference genome that will be masked
#' @param outdir **optional** a general output directory
#' @param mapdir **optional** an output directory for the remapping step
#' @param refdir **optional** an output directory for the cellranger reference genome files
#' @param cleandir **optional** an output directory for GATK cleaning/preparation steps
#' @param chemistry **optiona** 10X chemistry, use only when the automatic detection is failing
#' @param nthreads **optional** a number of threads to use
#' @return a list with a paths to prepared bam and barcode files.
prepare_sample = function(
    bam, reference, annotation, vcf,
    outdir=NULL, mapdir=NULL, refdir=NULL, cleandir=NULL,
    chemistry = "auto", nthreads=16
    ){
    if(is.null(outdir))
        outdir = "prepare"
    if(is.null(mapdir))
        mapdir = file.path(outdir, "map")
    if(is.null(refdir))
        refdir = file.path(outdir, "ref")
    if(is.null(cleandir))
        cleandir = file.path(outdir, "clean")

    core = phyloRNA::corename(bam)

    bam_cleaned = filename(outdir, core, ".cleaned.bam")
    bam_prepared = filename(outdir, core, ".prepared.bam")

    barcodes_prepared = filename(outdir, core, ".prepared.txt")

    h5_prepared = filename(outdir, core, ".h5")

    result = list(
        prefix = core,
        bam = bam_prepared,
        barcodes = barcodes_prepared,
        h5 = h5_prepared
        )

    # Skip if the final output files already exist
    if(file.exists(result$bam) && file.exists(result$barcodes) && file.exists(result$h5))
        return(result)

    mapped = phyloRNA::remap(
        input = bam,
        reference = reference,
        annotation = annotation,
        outdir = mapdir,
        refdir = refdir,
        chemistry = chemistry,
        nthreads = nthreads,
        copy_h5 = h5_prepared # only thing we want to copy
        )

    phyloRNA::gatk_prepare(
        input = mapped$bam,
        output = bam_cleaned,
        reference = reference,
        vcf = vcf,
        barcodes = mapped$barcodes,
        outdir = cleandir
        )

    pattern = "-1$"
    replace = paste0("-", core)

    phyloRNA::bamtagregex(
        input = bam_cleaned,
        output = bam_prepared,
        tag = "CB",
        pattern = pattern,
        replace = replace
        )

    file_sub(mapped$barcodes, barcodes_prepared, pattern=pattern, replace=replace)

    return(result)
    }


filename = function(dir, core, ext){
    file.path(dir, paste0(core, ext))
    }


merge_h5 = function(inputs, output){
    # well, I can read, merge them, but then not save in the same h5 format.
    # I would need to write a function for that.
    data = lapply(inputs, phyloRNA::expr_read10xh5)
    names = phyloRNA::corenames(inputs)
    data = phyloRNA::expr_merge(data, names)
    }
#!/usr/bin/env python3
"""Add RG tag as a CB tag."""
import argparse
import os
import pysam

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Add RG tag as a CB tag to a SAM/BAM file."""
        )
    parser.add_argument("input", type=str, help="an input SAM/BAM file")
    parser.add_argument("output", type=str, help="an output SAM/BAM file")
    args = parser.parse_args()
    return args

def main():
    """Add RG tag as a CB tag to a SAM/BAM file."""
    args = parse_args()

    imod = mode(args.input)
    omod = mode(args.output)

    with pysam.AlignmentFile(args.input, "r"+imod) as inbam, \
         pysam.AlignmentFile(args.output, "w"+omod, template=inbam) as outbam:
        for read in inbam.fetch(until_eof=True):
            rg2cb(read)
            outbam.write(read)

def rg2cb(read):
    """Add RG tag as a CB tag to a read"""
    text = read.get_tag("RG")
    read.set_tag("CB", text, value_type="Z")

def mode(file):
    """Detect pysam SAM/BAM file mode"""
    ext = os.path.splitext(file)[1]
    if ext == ".bam":
        return "b"
    if ext == ".sam":
        return ""
    raise ValueError(f"File {os.path.basename(file)} must be either SAM or BAM")

if __name__ == "__main__":
    main()
#' snv.r
#'
#' Functions for snv identification and filtering
import::here("utils.r", "filename", "num2char", "read_vcm")
import::here("filter.r", "density_filtering", "subset_filtering")
import::here("phyloRNA", "all_files_exist")

#' Detect SNV for scRNAseq
#'
#' This function will perform SNV detection for scRNAseq data.
#'
#' To detect SNVs on a scRNAseq, this function first detect a high-quality SNVs by treating
#' the sample as a bulk sample, then the most common base pertaining to each cell for every
#' variant is reported and saved to a vcm (variant call matrix) file.
#'
#' @param bam an input sam/bam file
#' @param barcodes a file with cell barcodes
#' @param reference a reference file to which the bam file was mapped
#' @param outdir **optional** a general output directory
#' @param vcfdir **optional** a directory with vcf files
#' @param vcm **optional** an output vcm file
#'
#' @return a path to an alignment table with the most frequent base for every cell at every SNV
#' position
detect_snv = function(
    bam, barcodes, reference,
    normal=NULL, pon=NULL, germline=NULL,
    outdir=NULL, vcfdir=NULL, vcm=NULL,
    nthreads=16
    ){
    core = phyloRNA::corename(bam)

    # Set default parameters:
    if(is.null(outdir))
        outdir = "snv"
    if(is.null(vcfdir))
        vcfdir = file.path(outdir, "vcf")
    if(is.null(vcm))
        vcm = file.path(outdir, paste0(core, ".vcm"))

    phyloRNA::mkdir(outdir)
    phyloRNA::mkdir(vcfdir)

    vcf = file.path(outdir, paste0(core, ".vcf"))

    phyloRNA::gatk_snv(
        bam, reference, vcf,
        normal=normal, pon=pon, germline=germline,
        outdir=vcfdir)
    phyloRNA::vcm(bam, vcf, barcodes, output=vcm, nthreads=nthreads)

    invisible(vcm)
    }




#' Preprocess SNV data
#'
#' Preprocessing SNVs. This include as bulk SNV detection, single-cell SNV detection
#' using bulk SNVs and dataset filtering.
#'
#' @param bam a bam file prepared according to the GATK best practices
#' @param barcodes file with barcodes
#' @param reference a genome reference to which the bam file was mapped
#' @param outdir **optional** an output directory
#' @param nthreads **optional** a number of threads to run on
preprocess_snv = function(
    bam, barcodes, reference, outdir=NULL, nthreads=16
    ){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    vcm = detect_snv(bam, barcodes, reference, outdir=outdir, nthreads=nthreads)

    invisible(vcm)
    }


#' Filter SNV
#'
#' Filter SNV dataset using two types of filtering approaches
#' see `src/filter.r` for more information
#'
#' @param vcm a vcm file
#' @param selection named list or vector specifying how many cells of each type should be selected
#' @param density desired data density
#' @param outdir an output directory
#' @return a list of filtered files
filter_snv = function(vcm, prefix, selection=NULL, density=NULL, outdir=NULL){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    if(is.null(selection) && is.null(density))
        stop("Either the selection or the density parameter must be specified.")


    if(!is.null(selection) && is.null(density)){
        file = filename(prefix, outdir=outdir)
        if(file.exists(file))
            return(invisible(file))

        data = as.data.frame(read_vcm(vcm))
        file = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            empty = "N",
            outdir = outdir
            )
        return(invisible(file))
        }


    if(!is.null(selection) && !is.null(density)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))

        data = as.data.frame(read_vcm(vcm))
        files = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            density = density,
            empty = "N",
            outdir = outdir
            )
        return(invisible(files))
        }

    if(is.null(selection)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))

        data = as.data.frame(read_vcm(vcm))
        files = density_filtering(data, prefix=prefix, density=density, empty="N", outdir=outdir)
        return(invisible(files))
        }
    }
#' sra.r
#'
#' Function to programmatically access to a SRA and GEO databases
import::here("magrittr", "%>%")
import::here("xml2", "read_xml", "as_list")
import::here("phyloRNA", "mkdir", "all_files_exist")
import::here("rentrez", "entrez_search", "entrez_link", "entrez_summary", "extract_from_esummary")

#' Download SRR samples
#'
#' Download 10X fastq files from the SRA database
#'
#' @param srr the SRR (short-read record, presumably) id that identifies reads in the SRA database.
#' SRR ids can be obtained by `get_srr_samples`.
#' @param prefix a sample prefix that identifies downloaded fastq files
#' @param outdir an output directory where the fastq files are downloaded
#' @param cellranger **optional** the files are renamed to the cellranger type
srr_download_sample = function(srr, prefix, outdir=NULL, cellranger=TRUE){
    # For Cellranger, file names need to have this structure:
    # [name]_S1_L001_[read type]_001.fastq.gz
    # S1 -- sample 1
    # L001 -- lane 001
    # read type:
    #     -- I1 -- Index file
    #     -- R1 -- barcodes or actual reads
    #     -- R2 -- barcodes or actual reads
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    if(length(prefix) != 1 || length(srr) != 1 || length(outdir) != 1)
        stop("This function is not vectorized. Provide a single srr, prefix and outdir.")

    if(cellranger){
        read_types = c("I1", "R1", "R2")
        fastqs = paste0(prefix, "_S1_L001_", read_types , "_001.fastq.gz")
        } else {
        fastqs = paste0(prefix, ".fastq.gz")
        }

    fastqs = file.path(outdir, fastqs)
    srr_files = file.path(outdir, paste0(srr, "_", seq_along(fastqs), ".fastq.gz"))

    if(all_files_exist(fastqs))
        return(invisible())

    # Download srr files and check if they exist/were downloaded correctly
    sra_download(srr, outdir)
    if(!all_files_exist(srr_files))
        stop("ERROR: not all files exist.\\n", "Files: ", files)

    file.rename(srr_files, fastqs)

    return(invisible())
    }


sra_download = function(srr, outdir){
    # Check if the srr fastq file or files already exists
    files = dir(outdir, pattern=paste0(srr, ".*\\.fastq.*"))
    if(length(files) != 0)
        return(invisible())

    sra_prefetch(srr, outdir)
    sra_dump(srr, outdir)
    }


sra_prefetch = function(srr, outdir, progress=TRUE){
    command = "prefetch"
    args = c(
        srr,
        "--output-directory", outdir
        )
    if(progress)
        args = c(args, "--progress")

    phyloRNA:::systemE(command, args)
    }


sra_dump = function(srr, outdir){
    command = "fastq-dump"
    args = c(
        srr,
        "--split-files",
        #"--skip-technical", required for 10X
        "--origfmt",
        "--gzip"
        )
    phyloRNA:::systemE(command, args, dir=outdir)
    }


get_srr_samples = function(gse, save=NULL){
    if(!is.null(save) && file.exists(save))
        return(readRDS(save))

    samples = get_gsm_samples(gse)
    samples$srr = sapply(samples$gsm, get_srr)

    if(!is.null(save))
        saveRDS(samples, save)

    samples
    }


get_gsm_samples = function(gse){
    esearch = entrez_search(db="gds", term=paste0(gse, "[ACCN] & gse[ETYP]"))
    esummary = entrez_summary(db="gds", id=esearch$ids)

    samples = esummary$samples
    colnames(samples) = c("gsm", "name")
    samples
    }


get_srr = function(gsm){
    esearch = entrez_search(db="gds", term=paste0(gsm, "[ACCN] gsm[ETYP]"))
    elinks = entrez_link(dbfrom="gds", id=esearch$ids, db="sra")
    esummary = entrez_summary(db="sra", elinks$links$gds_sra)

    run = extract_from_esummary(esummary, "runs")
    attrs = extract_attributes(run)
    srr = attrs[["acc"]]

    srr
    }


extract_attributes = function(x){
    x %>% read_xml %>% as_list() %>% getElement("Run") %>% attributes
    }
#' star.r
#'
#' mapping to the STAR RNA-seq mapping software
import::here("utils.r", "rg2cb")
import::here("phyloRNA", "corename", "mkdir")


#' Prepare sequences using the STAR software
#'
#' Map scRNA-seq sequences using the STAR software.
#'
#' This functions creates a STAR genome index and then maps the fastqs to the reference.
#' Reads are then retaged so that the read-group information is written to the cell barcode tag.
#'
#' @param prefix a prefix for the output bam file
#' @param fastqs one or more fastqs files
#' @param sample a sample information for each fastq file, will be written into readgroup
#' @param flowcell a flowcell information for the whole analysis, will be written into read-group
#' @param reference a fasta reference sequences
#' @param annotatin a gtf annotation file
#' @param overhang a size of allowed overhang, should be equal to the read lengths
#' @param outdir **optional** an output directory
#' @param mapdir **optional** a directory with mapped output
#' @param refdir **optional** a directory with STAR reference index
#' @param manifest **optional** a file path to pre-existing manifest
#' @param gzip **optional** whether the fastq files are compressed
#' @param nthreads **optional** the number of threads to run the analysis on
#' @return a mapped bam file with CB tag
STAR = function(
    prefix, fastqs, sample, flowcell,
    reference, annotation, overhang,
    outdir = NULL, refdir = NULL, mapdir = NULL,
    manifest = NULL, gzip = TRUE,
    nthreads = 16
    ){
    if(is.null(outdir))
        outdir = "."
    if(is.null(refdir))
        refdir = file.path(outdir, "ref")
    if(is.null(mapdir))
        mapdir = file.path(mapdir, "map")
    if(is.null(manifest))
        manifest = file.path(mapdir, "manifest.txt")

    mkdir(outdir)
    mkdir(refdir)
    mkdir(mapdir)

    mapped = file.path(mapdir, "allAligned.sortedByCoord.out.bam")
    tagged = file.path(mapdir, "all.bam")

    if(!file.exists(mapped)){    
        STAR_genome_index(reference, annotation, overhang, refdir, nthreads=nthreads)
        read_groups = make_read_group(sample, flowcell)
        STAR_manifest(fastqs, manifest, read_groups)
        STAR_map(prefix, refdir, manifest=manifest, outdir=mapdir, gzip=gzip)
        }

    if(!file.exists(tagged))
        rg2cb(mapped, tagged)

    return(tagged)
    }


#' Map sequences wih STAR
#'
#' Map scRNA-seq sequences using the STAR software
#'
#' The input fastq files can be provided in two ways, either directly using
#' the `fastq` parameter, or using the `manifest` parameter, which in addition
#' enable specification of the read-groups (see the [STAR_manifest]).
#'
#' @param prefix a file prefix for output files
#' @param reference a path to the STAR genome index, see [STAR_genome_index]
#' @param fastq **optional** fastq files that will be mapped, see details.
#' @param manifest **optional** STAR manifest file, see details
#' @param outdir **optional** an output directory
#' @param gzip **optional** whether the input fastq files are compressed
#' @param nthreads **optional** the number of threads to use
STAR_map = function(
    prefix, reference,
    fastq=NULL, manifest=NULL,
    outdir=NULL, gzip=FALSE, nthreads=16
    ){
    if(is.null(fastq) && is.null(manifest))
        stop("Either fastq or manifest must be specified.")
    if(is.null(fastq) && is.null(manifest))
        stop("Please, specify either fastq or manifest, not both.")

    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    statusfile = file.path(outdir, paste0(prefix, ".finished"))
    if(file.exists(statusfile))
        return(invisible())

    command = "STAR"
    args = c(
        "--genomeDir", reference,
        "--outFileNamePrefix", file.path(outdir, prefix),
        "--outSAMattributes NH HI AS nM NM MD jM jI MC ch RG",
        "--outSAMtype BAM SortedByCoordinate"
        )

    if(!is.null(fastq)){
        id = paste0("ID:", corename(fastq), collapse=" , ")
        args = c(args, "--readFilesIn", paste0(fastq, collapse=","))
        args = c(args, "--outSAMattrRGline", id)
        }

    if(!is.null(manifest))
        args = c(args, "--readFilesManifest", manifest)

    if(gzip)
        args = c(args, "--readFilesCommand gunzip -c")

    phyloRNA:::systemE(command, args)
    file.create(statusfile)
    }


#' Create STAR manifest
#'
#' Create STAR manifest.
#'
#' STAR manifest is an alternative way to specify input `fastq` files for the STAR
#' RNA-seq mapping software. For the specification of read-group tag,
#' see [make_read_group]
#'
#' @param fastq one or more fastq files
#' @param file a character string naming a file
#' @param rg **optional** a read-group tag, see [make_read_group]
STAR_manifest = function(fastq, file, rg=NULL){
    if(file.exists(file))
        return(invisible())

    manifest = fastq

    if(!is.null(rg))
        manifest = paste(fastq, "-", rg, sep="\t")
        
    writeLines(manifest, file)
    }


#' Construct a read-group tag
#'
#' Construct a read-group (RG) tag required by the GATK. See [details] for the
#' description of individual parts of the RG tag.
#'
#' SAM/BAM read-groups tags required by the GATK are:
#' RG tags required by GATK are:
#'    * ID identifies read and is copied to read
#'    * SM identifies sample
#'    * LB identifies library, can have multiple libraries per sample
#'    * PU consist of flowcell barcode, lane and sample barcode
#'    * PL identifies platform, Illumina in this case
#'
#' The the values for the RG tags are:
#'  * ID -- `sample`
#'  * SM -- `sample`
#'  * LB -- `sample`
#'  * PU -- `flowcell:1:sample`
#'  * PL -- "ILLUMINA"
make_read_group = function(sample, flowcell, platform="ILLUMINA"){
    paste(
        paste0("ID:", sample),
        paste0("SM:", sample),
        paste0("LB:", sample),
        paste0("PU:", paste(flowcell, "1", sample, sep=":")),
        paste0("PL:", "ILLUMINA"),
        sep="\t"
        )
    }


#' Generate a STAR genome index
#'
#' Generate a reference STAR genome index.
#'
#' @param reference a path to the reference fasta file
#' @param gtf a path to the annotation gtf file
#' @param overhang **optional** a size of allowed overhang, should be equal
#' to read length
#' @param outdir **optional** an output directory
#' @param nthreads **optional** the number of threads to use
STAR_genome_index = function(reference, gtf, overhang=50, outdir=NULL, nthreads=16){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    statusfile = file.path(outdir, "genome.finished")
    if(file.exists(statusfile))
        return(invisible())
    
    command = "STAR"
    args = c(
        "--runMode genomeGenerate",
        "--genomeDir", outdir,
        "--genomeFastaFiles", reference,
        "--sjdbGTFfile", gtf,
        "--sjdbOverhang", len-1,
        "--runThreadN", nthreads
        )
    phyloRNA:::systemE(command, args)
    file.create(statusfile)
    }
#' stats.r
#'
#' Calculate various summary statistics

# What do I need?
# -- Dimension and the data density of each sample
# -- Dimension and the data density of each filtered dataset
library("Matrix")
library("data.table")
import::here("utils.r", "mdensity", "read_vcm", "write_table", "is_empty")

expr_h5_stats = function(h5, names=NULL, file=NULL){
    if(is.null(names))
        names = phyloRNA::corename(h5)

    # Calculate stats for each dataset
    data = lapply(h5, phyloRNA::expr_read10xh5)
    stats = lapply(data, expr_h5_stats_dataset)
    names(stats) = names

    # calculate stats for the merged dataset
    data = phyloRNA::expr_merge(data, names)
    stats = c(stats, "total"=list(expr_h5_stats_dataset(data)))
    table = do.call(rbind, stats)

    if(!is.null(file))
        write_table(table, file)

    return(invisible(table))
    }

genes_per_cell = function(x){
    mean( colSums(x != 0) )
    }

umi_per_cell = function(x){
    mean( colSums(x) )
    }

umi_per_gene = function(x){
    mean( rowSums(x) )
    }

expr_h5_stats_dataset = function(data){
    stats = list(
        "ncells" = ncol(data),
        "ngenes" = nrow(data),
        "ngenes_detected" = sum(rowSums(data) > 0),
        "genes_per_cell" = genes_per_cell(data),
        "umi" = sum(data),
        "umi_per_gene" = umi_per_gene(data),
        "umi_per_cell" = umi_per_cell(data),
        "density" = mdensity(data, 0)
        )
    return(stats)
    }

vcm_stats = function(vcm, file=NULL){
    data = read_vcm(vcm)

    empty_cells = colnames(data)[colSums(is_empty(data, "N")) == 0]
    stats = list(
        "ncells" = ncol(data),
        "SNVs" = nrow(data),
        "diversity" = diversity(colnames(data)),
        "density" = mdensity(data, "N"),
        "nEmpty" = length(empty_cells),
        "divEmpty" = diversity(empty_cells)
        )
    stats = do.call(rbind, stats)

    if(!is.null(file))
        write_table(stats, file)

    invisible(stats)
    }

filtered_stats = function(files, empty, output){
    stats = list()
    for(file in files){
        name = corename(file)
        data = read_table(file)
        stats[[name]] = filtered_stats_dataset(data, empty)
        }
    stats = do.call(rbind, stats)
    write_table(stats, output) 
    }

fasta_stats = function(fasta, stats=NULL, name=TRUE, unknown="N"){
    if(is.null(stats))
        stats = paste0(tools::file_path_sans_ext(fasta), ".txt")
    if(length(fasta) != length(stats))
        stop("fasta and stats vectors must have the same length")

    if(isTRUE(name))
        name = phyloRNA::corename(fasta)

    n = length(fasta)
    name = rep_len(name, n)
    unknown = rep_len(unknown, n)

    stats = list()
    for(i in seq_along(fasta)){
        seq = phyloRNA::read_fasta(fasta[i])
        tab = phyloRNA::seq2tab(seq)

        text = paste0(
            "Sequences: ", nrow(tab), "\n",
            "Sites: ", ncol(tab), "\n",
            "Unique patterns: ", ncol(unique.matrix(tab, MARGIN=2)), "\n",
            "Data density: ", mdensity(tab, empty=unknown[i]), "\n",
            "Diversity: ", diversity(names(seq))
            )
        if(is.character(name))
            text = paste0("Name: ", name[i], "\n", text)
        stats[[i]] = text
        }

    stats
    }

diversity_stats = function(files, output){
    stats = list()
    for(file in files){
        name = corename(file)
        data = read_table(file)
        stats[[name]] = diversity(colnames(data))
        }
    stats = do.call(rbind, stats)
    write_table(stats, output)
    }

filtered_stats_dataset = function(data, empty){
    stats = list(
        "ncells" = ncol(data),
        "nrow" = nrow(data),
        "density" = mdensity(data, empty)
        )
    return(stats)
    }

diversity = function(names){
    names = sub(".*-", "", names)
    tab = table(names)
    div = paste(names(tab), tab, sep=":", collapse=", ")
    div
    }

expr_stats = function(h5, filtered, outdir=NULL){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    h5_stats = file.path(outdir, "h5_stats.txt")
    if(!file.exists(h5_stats))
        expr_h5_stats(h5, h5_stats)

    expr_stats = file.path(outdir, "expr_stats.txt")
    if(!file.exists(expr_stats))
        filtered_stats(filtered, empty="-", expr_stats)

    expr_diversity = file.path(outdir, "expr_diversity.txt")
    if(!file.exists(expr_diversity))
        diversity_stats(filtered, expr_diversity)
    }

snv_stats = function(vcm, filtered, outdir=NULL){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    vcm_stats = file.path(outdir, "vcm_stats.txt")
    if(!file.exists(vcm_stats))
        vcm_stats(vcm, vcm_stats)

    snv_stats = file.path(outdir, "snv_stats.txt")
    if(!file.exists(snv_stats))
        filtered_stats(filtered, empty="N", snv_stats)

    snv_diversity = file.path(outdir, "snv_diversity.txt")
    if(!file.exists(snv_diversity))
        diversity_stats(filtered, snv_diversity)
    }
#' utils.r
#'
#' shared utility functions
import::here("phyloRNA",
    "remove_constant",
    "write_fasta", "read_fasta",
    "tab2seq", "seq2tab",
    "all_files_exist", "mkdir", "corename"
    )
import::here("data.table", "fread")


#' Write a table
#'
#' Write a table in a particular format. This is a simple wrapper around write.table
#' with a few specified parameters.
#' @param x a matrix or a data frame
#' @param file an output path
write_table = function(x, file){
    write.table(x, file=file, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    }


#' Read a table
#'
#' Read a table in a particular format. This is a simple wrapper around read.table
#' with a few specified parameters.
#' @param file a file in tabular format
#' @return data.frame
read_table = function(file){
    read.table(file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, sep="\t")
    }


#' Convert a numeric value to a character string
#'
#' Converts a numeric value in the format `0.X` into a character string `0X`
#' @param x numeric vector
#' @return character vector
num2char = function(x){
    sub(".", "", as.character(x), fixed=TRUE)
    }


#' Convert a tabular file into a fasta format
#'
#' @param file one or more files in tabular format
#' @param fasta **optional** output path for fasta files
#' @param outdir **optional** an outut directory, if fasta is not specified
#' @param margin whether rows (1) or columns (2) should be concatenated
#' @return a vector of fasta files
table2fasta = function(file, fasta=NULL, outdir=NULL, margin=2, zero=NULL){
    if(!is.null(fasta) && length(fasta) != length(file))
        stop("The file and fasta vectors must have the same length!")
    if(is.null(fasta))
        fasta = paste0(tools::file_path_sans_ext(file), ".fasta")
    if(!is.null(outdir))
        fasta = file.path(outdir, basename(fasta))
    mkdir(outdir)

    if(all_files_exist(fasta))
        return(invisible(fasta))

    for(i in seq_along(file)){
        data = read_table(file[i])
        if(!is.null(zero))
            data[data == zero] = 0
        seq = tab2seq(data, margin=margin)
        write_fasta(seq, file=fasta[i])
        }

    invisible(fasta)
    }


fasta2stats = function(fasta, stats=NULL, name=TRUE, unknown="N"){
    if(is.null(stats))
        stats = paste0(tools::file_path_sans_ext(fasta), ".txt")
    if(length(fasta) != length(stats))
        stop("fasta and stats vector must have the same length")
    if(all_files_exist(stats))
        return(invisible())

    if(isTRUE(name))
        name = corename(fasta)

    n = length(fasta)
    name = rep_len(name, n)
    unknown = rep_len(unknown, n)

    for(i in seq_along(fasta)){
        seq = read_fasta(fasta[i])
        tab = seq2tab(seq)

        text = paste0(
            "Sequences: ", nrow(tab), "\n",
            "Sites: ", ncol(tab), "\n",
            "Unique patterns: ", ncol(unique.matrix(tab, MARGIN=2)), "\n",
            "Data density: ", mdensity(tab, empty=unknown[i])
            )

        if(is.character(name))
            text = paste0("Name: ", name[i], "\n", text)

        writeLines(text, stats[i])
        }
    }

#' Calculate data density of matrix
#'
#' Calculate the data density, that is the proportion of known elements in the matrix.
#' @param x a matrix
#' @param empty an unknown element
#' @return a data density of matrix
mdensity = function(x, empty){
    sum(is_empty(x, empty)) / prod(dim(x))
    }


#' Test if an element is not unknown
#'
#' @param x vector or matrix
#' @param empty an unknown element
#' @return vector or matrix 
is_empty = function(x, empty){
    if(is.na(empty))
        return(!is.na(x))
    x != empty
    }


#' Download file
#'
#' @param url an url from which file is downloaded
#' @param file a character string where downloaded file will be saved
download_file = function(url, file, rewrite=FALSE){
    if(!file.exists(file) || rewrite)
        download.file(url, file)
    }


#' Construct a filename
#'
#' A simple shorthand for construction a file name
filename = function(prefix, suffix="", ext=".txt", outdir="."){
    file.path(outdir, paste0(prefix, suffix, ext))
    }


#' Convert readgroup to a cell barcode
#'
#' Some single-cell detection methodology, such as the `phyloRNA::vcm()` tool, expect that every
#' read is barcoded with a cell-specific barcode. This functions transform non-barcoded single-cell
#' bam file into a barcoded bam file by adding the read-group (RG) to the cell barcode (CB) tag
#' @param input a bam file with reads encoded with RG tag
#' @param output an output bam file with read-group written into the CB tag
rg2cb = function(input, output){
    if(file.exists(output))
        return(invisible())

    command = "python3"
    args = c(
        "src/rg2cb.py",
        input, output
        )
    phyloRNA:::systemE(command, args)
    }


#' Read the vcm file and memoise it
#'
#' This function is memoised (possible reuse) of the vcm file.
#' `data.table::fread()` is used here due to a huge file size.
#' @param vcm a variant call matrix file
#' @return a variant call matrix as a data.table
read_vcm = local({
    memory = list()

    function(vcm){
        if(!is.null(memory[[vcm]]))
            return(memory[[vcm]])

        # using data.tale due to a huge size of the dataset
        data = fread(vcm, header=TRUE)
        # first three columns are not cells (chromosome, position and reference)
        data = data[, -c(1:3)]
        memory[[vcm]] <<- data

        data
        }
    })


#' Convert vcm file to fasta file
#'
#' Convert a vcm fle to fasta file.
#'
#' @param vcm a variant call matrix file
#' @param fasta a fasta file
#' @param selection selected cell barcodes that will be retained
vcm2fasta = function(vcm, fasta, selection=NULL){
    if(file.exists(fasta))
        return(invisible())

    data = read_vcm(vcm)

    if(!is.null(selection)){
        match = selection %in% colnames(data)
        if(!all(match)){
            warning("WARNING: ", sum(!match), " out of ", length(match),
                " requested cells are not present:\n",
                paste0(selection[!match], collapse="\n")
                )
            selection = selection[match]
            }
        data = data[, ..selection] # data.table' subsetting
        }

    data = as.matrix(data)
    data = remove_constant(data)
    seq = tab2seq(data, 2)
    write_fasta(seq, fasta)
    }


#' Substitute a pattern in a file
#'
#' Substitute a pattern over all lines in a file.
#'
#' This is a simple combination of the `sub` replacement function with `readLines` and `writeLines`.
#'
#' @param input an input file
#' @parma output an output file where pattern will be replaced
#' @param pattern a pattern to be replaced
#' @param replace a replacement for pattern
#' @param fixed pattern is a simple character string, not a regular expression
file_sub = function(input, output, pattern, replace, fixed=FALSE){
    lines = readLines(input)
    lines = sub(pattern, replace, lines, fixed=fixed)
    writeLines(lines, output)
    }


#' Merge files
#'
#' Merge multiple files into a single file.
#'
#' @param inputs one or multiple files to merge
#' @param a merged file
#' @param overwrite **optional** if an existing output should be overwritten
merge_files = function(inputs, output, overwrite=FALSE){
    if(file.exists(output) && overwrite)
        file.remove(output)
    if(!file.exists(output))
        file.append(output, inputs)
    }
#' run.r
#'
#' Run the analysis. This includes:
#'
#' Preparation
#' * remapping, demultiplexing, barcode correction and expression counts using Cellranger
#' * cleaning BAM files according to the GATK best practices
#' * adding a sample-specific postfix to cell barcodes
#'
#' Pre-processing:
#' * Expression
#'   -- standardization of genes into mu=0 and sd=1
#'   -- categorization according to empirical 60% and 90% HDI
#' * SNV:
#'   -- as bulk SNV identification and filtering with Mutect2
#'   -- sc SNV identification with vcm.py
#' * stepwise filtration into 20%, 50% and 90% density
#' * alternative filtration into 58 best cells and full dataset, 50% and 90% density
#'
#' Filtering:
#' * stepwise filtration into 20%, 50% and 90% density
#' * subset into 58 best cells and filtering into 50% and 90% density
#'
#' Phylogenetic analysis:
#' * ML with stepwise filtering
#' * ML and BI with alternative filtration
#' * BEAST templates created with the `beter` package
#' * Expression:
#'   -- IQtree: ORDINAL+ASC, ultrafast bootstrap -B 1000
#'   -- BEAST: ordinal from MM, exponential pop growth, coalescent prior, strict clock, two runs
#' * SNV:
#'   -- IQtree: GTR+gamma, ultrafast bootstrap -B 1000
#'   -- BEAST: GTR, exponential pop growth, coalescent prior, strict clock, two runs
#'

# use import::from instead?
import::from("src/prepare.r", "prepare_samples")
import::from("src/snv.r", "detect_snv", "filter_snv")
import::from("src/utils.r", "table2fasta", "fasta2stats")
import::from("src/expr.r", "preprocess_expression", "filter_expression")
import::from("src/beast.r", "beasts")
import::from("src/iqtree.r", "iqtrees")
import::from("phyloRNA", "corename")

main = function(){
    # datasets:
    bam = dir("data", full.names=TRUE)
    outdir = file.path("moravec2021")

    # required reference files:
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    refdir = "reference/ref" # shared cellranger ref
    vcf = "reference/00-common_all.vcf.gz"
    
    # PON and normal samples
    pon = "pon/MDA-MB-231/pon.vcf" # see make_panel_of_normals.r
    normal = file.path( # see make_normal_samples.r
        outdir, "normal",
        c("MDAMB231-MPmt-rep1.prepared.bam", "MDAMB231-MPmt-rep2.prepared.bam")
        )

    # Other settings:
    nthreads = 16
    chemistry = "SC3Pv2"
    densities = c(0.2, 0.5, 0.9)
    hdi = c(0.6, 0.9)
    selection = c("T1" = 20, "T3" = 20, "T2" = 6, "CTC1" = 6, "CTC2" = 6)

    # Preparation:
    prepared = prepare_samples(
        bam, reference, annotation, vcf,
        chemistry = chemistry, nthreads = nthreads,
        outdir = file.path(outdir, "prepare"),
        refdir = refdir
        )

    # SNV part
    vcm = detect_snv(
        bam = prepared$bam,
        barcodes = prepared$barcodes,
        reference = reference,
        normal = normal,
        pon = pon,
        outdir = file.path(outdir, "snv")
        )

    tabdir = file.path(outdir, "snv", "filtered")
    fastadir = file.path(outdir, "snv", "fasta")
    treedir = file.path(outdir, "snv", "tree")

    snv = filter_snv(vcm=vcm, density=densities, prefix="snv", outdir=tabdir)
    snv_fasta = table2fasta(snv, outdir=fastadir)
    fasta2stats(snv_fasta, unknown="N")

    snv_subset = filter_snv(vcm=vcm, selection=selection, prefix="snv_subset",
                            outdir=tabdir)
    snv_subset_fasta = table2fasta(snv_subset, outdir=fastadir)
    fasta2stats(snv_subset_fasta, unknown="N")

    iqtrees(
        c(snv_fasta, snv_subset_fasta),
        model = "TEST",
        bootstrap = 100, parallel = TRUE, nthreads = 16,
        outdir = file.path(treedir, "ML")
        )

    beasts(
        snv_subset_fasta,
        template = file.path("templates", "BDStrictGtr.xml"),
        outdir = file.path(treedir, "BI")
        )

    # expression part
    expr_preprocessed = preprocess_expression(
        h5 = prepared$h5,
        hdi = hdi,
        minGene=0,
        minUMI=0,
        outdir = file.path(outdir, "expr", "prepare"),
        prefix = "all"
        )

    filterdir = file.path(outdir, "expr", "filtered")
    fastadir = file.path(outdir, "expr", "fasta")
    treedir = file.path(outdir, "expr", "tree")
    
    expr = filter_expression(
        expr_preprocessed$discretized, prefix = "expr",
        outdir = filterdir, density = densities
        )
    expr_subset = filter_expression(
        expr_preprocessed$discretized, prefix = "expr_subset",
        outdir = filterdir, selection = selection
        )

    expr_fasta = table2fasta(expr, outdir=fastadir)
    fasta2stats(expr_fasta, unknown="-")

    expr_subset_fasta = table2fasta(expr_subset, outdir=fastadir)
    fasta2stats(expr_subset_fasta, unknown="-")
    expr_zero_fasta = table2fasta(
        expr_subset,
        file.path(fastadir, paste0(corename(expr_subset), "_zero.fasta")),
        outdir = fastadir, zero = "-"
        )

    iqtrees(
        expr_fasta,
        model = "ORDERED+ASC",
        outdir = file.path(treedir, "ML"),
        mc.cores = length(expr_fasta)
        )

    iqtrees(
        c(expr_subset_fasta, expr_zero_fasta),
        model = "ORDERED+ASC",
        bootstrap = 100, parallel = TRUE, nthreads = 16,
        outdir = file.path(treedir, "ML")
        )

    beasts(expr_subset_fasta, outdir = file.path(treedir, "BI"),
           template = file.path("templates", "BDStrictOrdinal.xml"))
    beasts(expr_zero_fasta, outdir = file.path(treedir, "BI"),
           template = file.path("templates", "BDStrictOrdinalZero.xml"))
    }


if(sys.nframe() == 0){
    main()

