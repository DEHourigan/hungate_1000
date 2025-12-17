# Hungate alignments

# Set working directory
setwd("/data/users/david/hungate_410_asmash/analysis")

library(msa)
library(ape)
library(ggtree)
library(ggplot2)
library(tinytex)
library(data.table)

<<<<<<< HEAD
# Function to align sequences and generate LaTeX and PDF files
=======
# Function to align sequences 
>>>>>>> Initial commit: add Hungate_alignments.R and core_peptide_alignments
align_and_generate <- function(input_file, output_tex, output_fasta, paper_width, paper_height, consensus_threshold = 50) {
    seq <- readAAStringSet(input_file, format = "fasta")
    aln <- msa(seq, method = "ClustalO")
    msaPrettyPrint(aln, file = output_tex, output = "tex",
                                 showNames = "left", showNumbering = "none", showLogo = "top",
                                 showConsensus = "bottom", consensusThreshold = consensus_threshold,
                                 shadingMode = "identical", logoColors = "rasmol",
                                 verbose = FALSE, askForOverwrite = FALSE, paperWidth = paper_width, paperHeight = paper_height,
                                 alFile = output_fasta)
    tinytex::pdflatex(output_tex)
}

<<<<<<< HEAD
# Alignments and LaTeX generation
align_and_generate("alignment_folder/pediocin_like_no_leader.faa", 
                                     "alignment_folder/pediocin_like_no_leader.tex", 
                                     "alignment_folder/pediocin_like_no_leader.fasta", 
                                     10, 3.5)

align_and_generate("core_peptide_alignments/nisin_like.faa", 
                                     "core_peptide_alignments/nisin_like.tex", 
                                     "core_peptide_alignments/nisin_like.fasta", 
                                     10, 3.5)

align_and_generate("core_peptide_alignments/non-nisin-like-class-i-corepeptide.faa", 
                                     "core_peptide_alignments/lantibiotic_non_nisin.tex", 
                                     "core_peptide_alignments/lantibiotic_non_nisin.fasta", 
                                     16, 3.5)

align_and_generate("core_peptide_alignments/circular_bacteriocin.faa", 
                                     "core_peptide_alignments/circular_bacteriocin.tex", 
                                     "core_peptide_alignments/circular_bacteriocin.fasta", 
                                     14, 3.5)

align_and_generate("core_peptide_alignments/IIB/IIB_groupA.faa", 
                                     "core_peptide_alignments/IIB_bacteriocin.tex", 
                                     "core_peptide_alignments/IIB_bacteriocin.fasta", 
                                     20, 9, consensus_threshold = 10)

# Lassopeptide alignments
setwd("core_peptide_alignments/lassopeptide")
align_and_generate("lassopeptide_cores.faa", "lassopeptide.tex", "lassopeptide.fasta", 20, 9, consensus_threshold = 10)
align_and_generate("lasso_clas_I.faa", "lasso_clas_I.tex", "lasso_clas_I.fasta", 8, 3)
align_and_generate("lasso_clas_II.faa", "lasso_clas_II.tex", "lasso_clas_II.fasta", 10, 6)
align_and_generate("GAD_motif.faa", "GAD_motif.tex", "GAD_motif.fasta", 10, 6)

# Lanthipeptide II alignments
setwd("../lanthi_II")
align_and_generate("filtered_35.faa", "lanthi_class_II.tex", "filtered_35.fasta", 20, 9, consensus_threshold = 10)
=======
# Alignments 
align_and_generate("core_peptide_alignments/pediocin/pediocin_like_no_leader.faa", 
                                     "core_peptide_alignments/pediocin/pediocin_like_no_leader.tex", 
                                     "core_peptide_alignments/pediocin/pediocin_like_no_leader.fasta", 
                                     7, 3)

align_and_generate("core_peptide_alignments/nisin/nisin_like.faa", 
                                     "core_peptide_alignments/nisin/nisin_like.tex", 
                                     "core_peptide_alignments/nisin/nisin_like.fasta", 
                                     10, 3.5)

align_and_generate("core_peptide_alignments/non-nisin/non-nisin-like-class-i-corepeptide.faa", 
                                     "core_peptide_alignments/non-nisin/lantibiotic_non_nisin.tex", 
                                     "core_peptide_alignments/non-nisin/lantibiotic_non_nisin.fasta", 
                                     16, 3.5)

align_and_generate("core_peptide_alignments/IIC/circular_bacteriocin.faa", 
                                     "core_peptide_alignments/IIC/circular_bacteriocin.tex", 
                                     "core_peptide_alignments/IIC/circular_bacteriocin.fasta", 
                                     14, 3.5)

align_and_generate("core_peptide_alignments/IIB/IIB_groupA.faa", 
                                     "core_peptide_alignments/IIB/IIB_bacteriocin.tex", 
                                     "core_peptide_alignments/IIB/IIB_bacteriocin.fasta", 
                                     20, 9, consensus_threshold = 10)

align_and_generate("core_peptide_alignments/lassopeptide/lassopeptide_cores.faa",
				   "core_peptide_alignments/lassopeptide/lassopeptide.tex",
				   "core_peptide_alignments/lassopeptide/lassopeptide.fasta",
				   20, 9, consensus_threshold = 10)

align_and_generate("core_peptide_alignments/lassopeptide/lasso_clas_I.faa",
				   "core_peptide_alignments/lassopeptide/lasso_clas_I.tex",
				   "core_peptide_alignments/lassopeptide/lasso_clas_I.fasta",
				   8, 3)

align_and_generate("core_peptide_alignments/lassopeptide/lasso_clas_II.faa",
				   "core_peptide_alignments/lassopeptide/lasso_clas_II.tex",
				   "core_peptide_alignments/lassopeptide/lasso_clas_II.fasta",
				   10, 6)

align_and_generate("core_peptide_alignments/lassopeptide/GAD_motif.faa",
				   "core_peptide_alignments/lassopeptide/GAD_motif.tex",
				   "core_peptide_alignments/lassopeptide/GAD_motif.fasta",
				   10, 6)


align_and_generate("core_peptide_alignments/lanthi_II/filtered_35.faa", 
				   "core_peptide_alignments/lanthi_II/lanthi_class_II.tex", 
				   "core_peptide_alignments/lanthi_II/filtered_35.fasta",
				   20, 9, consensus_threshold = 10)

align_and_generate("core_peptide_alignments/ranthipeptide/ranthipeptide_90.faa", 
	"core_peptide_alignments/ranthipeptide/ranthipeptide_90.tex", 
	"core_peptide_alignments/ranthipeptide/ranthipeptide_90.fasta",
	20, 9, consensus_threshold = 100)

align_and_generate("core_peptide_alignments/cyclic_lactone_ai/cai_90.faa", 
	"core_peptide_alignments/cyclic_lactone_ai/cai_90.tex", 
	"core_peptide_alignments/cyclic_lactone_ai/cai_90.fasta", 20, 12)


align_and_generate("core_peptide_alignments/IIB/IIB_groupA_core_only.faa", 
	"core_peptide_alignments/IIB/IIB_groupA.tex", 
	"core_peptide_alignments/IIB/IIB_groupA.fasta", 20, 12)


align_and_generate("core_peptide_alignments/angicin/angicin.faa",
	"core_peptide_alignments/angicin/angicin.tex", 
	"core_peptide_alignments/angicin/angicin.fasta", 20, 12)
>>>>>>> Initial commit: add Hungate_alignments.R and core_peptide_alignments

# Tree plotting
lan2_tree <- read.tree("raxml_out.raxml.supportFBP")
lan2_tree_plot <- ggtree(lan2_tree, size = 0.3) +
    geom_tiplab(size = 1) +
    theme_tree() +
    theme(legend.position = "none") +
    expand_limits(x = c(0, 10))

ggsave(lan2_tree_plot, filename = "lan2.png", width = 5, height = 5, units = "cm", dpi = 300)

<<<<<<< HEAD
# Additional alignments
setwd("../ranthipeptide")
align_and_generate("ranthipeptide_90.faa", "ranthipeptide_90.tex", "ranthipeptide_90.fasta", 20, 9, consensus_threshold = 100)

setwd("../cyclic_lactone_ai")
align_and_generate("cai_90.faa", "cai_90.tex", "cai_90.fasta", 20, 12)

setwd("../IIB")
align_and_generate("IIB_groupA_core_only.faa", "IIB_groupA.tex", "IIB_groupA.fasta", 20, 12)

setwd("../angicin")
align_and_generate("angicin.faa", "angicin.tex", "angicin.fasta", 20, 12)
=======
>>>>>>> Initial commit: add Hungate_alignments.R and core_peptide_alignments
