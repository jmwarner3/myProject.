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


