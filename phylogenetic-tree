library(ape)
library(phytools)
library(adephylo)



# read aln files
system("rm *.aln concat*")
fas_lst <- list.files(pattern = ".fa")
aln <- list()
all_seqnames <- NULL
for(i in 1:length(fas_lst)){
	aln[[i]] <- read.dna(fas_lst[[i]], format="fasta")

	# clean seq names
	attributes(aln[[i]])$dimnames[[1]] <- substr(attributes(aln[[i]])$dimnames[[1]], 9, nchar(attributes(aln[[i]])$dimnames[[1]]))
	all_seqnames <- c(all_seqnames, attributes(aln[[i]])$dimnames[[1]])

}

all_seqnames <- unique(all_seqnames)

# check for dup seq names -- and rm shorter sequence from alignment
#for(i in 1:2){
for(i in 1:length(fas_lst)){
	print(paste0("Now doing aln ", i, " -- ", fas_lst[i]))
	cur_aln <- aln[[i]]
	
	# DB seq
	bin_seq_pos <- which(grepl("Bin_", attributes(cur_aln)$dimnames[[1]]))
	db_seq_pos <- which(grepl("GCA_", attributes(cur_aln)$dimnames[[1]]))

	# find dups
	cur_seq_names <- attributes(cur_aln)$dimnames[[1]]
	dup_pos <- which(duplicated(cur_seq_names))
	print(paste0("- found ", length(dup_pos), " dup pairs"))
	
	# ID shorter seq
	seq_dups_2_rm <- NULL
	for(j in 1:length(dup_pos)){
		print(paste0("- removing pair ", j))
		cur_pos <- dup_pos[j]
		cur_dup_pos <- which(cur_seq_names == cur_seq_names[cur_pos])
		
		if(length(cur_dup_pos) > 2){
			print("!!!!!")
		}

		#seq1 <- cur_aln[dup_pos[1],]
		#n_indels_seq1 <- which(as.character(seq1) == "-")
		#seq2 <- cur_aln[dup_pos[2],]
		#n_indels_seq2 <- which(as.character(seq2) == "-")
		#paste(as.character(seq1), collapse="")
		
		# tmp aln with only DB seqs and the 2 dups
		#tmp_aln <- cur_aln[c(db_seq_pos, cur_dup_pos), ] ## too slow...
		subset_aln <- 500
		tmp_aln <- cur_aln[c(sample(db_seq_pos, subset_aln, replace=F), cur_dup_pos), ]
		attributes(tmp_aln)$dimnames[[1]][subset_aln+1] <- paste0(attributes(tmp_aln)$dimnames[[1]][subset_aln+1], "dup1")
		attributes(tmp_aln)$dimnames[[1]][subset_aln+2] <- paste0(attributes(tmp_aln)$dimnames[[1]][subset_aln+2], "dup2")
		tmp_dist <- dist.dna(tmp_aln, model="JC69")
		# fix NA distances
		tmp_dist[which(is.na(tmp_dist))] <- 1e-6
		tmp_dist[which(tmp_dist < 1e-6)] <- 1e-6
		tmp_dist[which(tmp_dist > 999)] <- 999
		#tmp_tr <- bionj(tmp_dist)
		tmp_tr <- nj(tmp_dist)
		tmp_tr_MP_rooted <- midpoint.root(tmp_tr)
		# fix <0 brlen
		tmp_tr_MP_rooted$edge.length[tmp_tr_MP_rooted$edge.length < 0] <- 0.0
		dist_pair <- distRoot(tmp_tr_MP_rooted, tips = c(subset_aln+1, subset_aln+2))
		max_dist_pair <- which.max(dist_pair)
		seq_dups_2_rm <- c(seq_dups_2_rm, cur_dup_pos[max_dist_pair])
	}
	# update aln[[i]]
	aln[[i]] <- cur_aln[-seq_dups_2_rm, ]

	# ID missing sequences in each alignment
	miss_seq <- setdiff(all_seqnames, attributes(aln[[i]])$dimnames[[1]])
	
	# fill in missing sequences
	write.dna(aln[[i]], paste0(fas_lst[i], ".aln"), format = "fasta")
	for(j in 1:length(miss_seq)){
		cat(paste0("
		> ", miss_seq[j], "
		", paste(rep("N", dim(aln[[i]])[2]), collapse="")) , file = paste0(fas_lst[i], ".aln"), sep = "", append = T)
	}
}

clean_aln <- list()
for(i in 1:length(fas_lst)){
	clean_aln[[i]] <- read.dna(paste0(fas_lst[i], ".aln"), format="fasta")
}

cat_aln <- do.call(cbind, clean_aln)
write.dna(cat_aln, "concat.aln")


            ##############
            # Trimming Alignment Files
            ##############
           
#remove NCBI genomes with more than 10% gaps over the whole alignment, and metagenomic bins with more than 50% gaps
alignment_length <- nchar(outfile[2])
too_short_seq_ncbi <- which(nchar(as.character(outfile)) -nchar( gsub("-", "", outfile, fixed=T)) > alignment_length*0.1)
too_short_seq_ncbi <- too_short_seq_ncbi[too_short_seq_ncbi > max(grep("Bin", outfile)) + 1]
too_short_seq_bins <- which(nchar(as.character(outfile)) -nchar( gsub("-", "", outfile, fixed=T)) > alignment_length*0.5)
too_short_seq_bins <- too_short_seq_bins[too_short_seq_bins <= max(grep("Bin", outfile)) + 1]
too_short_seq <- c(too_short_seq_bins, too_short_seq_ncbi)
outfile <- outfile[-sort(c(too_short_seq - 1, too_short_seq))]
write.table(outfile, file="C:/Users/matti/VBShare/concatenated_RP.fasta", quote=F, row.names=F, col.names=F)

#make trees with FastTree: 
#FastTreeMP -gtr -gamma -nt concatenated_RP.fasta > concatenated_RP.tre

#read in the tree and genomes files
concatenated_tree <- read.tree(file="C:/Users/matti/VBShare/concatenated_RP2.tre")
ncbi_taxonomy <- read.table("C:/Users/matti/VBShare/genomes_taxonomy.txt")
colnames(ncbi_taxonomy) <- c("genome", "ncbi_id")
ncbi_taxonomy <- ncbi_taxonomy[which(ncbi_taxonomy$genome %in% concatenated_tree$tip.label),]

#extract taxonomy for the NCBI genomes and save the file (so it can later be read in)
taxize_classification <- classification(ncbi_taxonomy$ncbi_id, db="ncbi")
taxize_classification <- do.call(rbind, taxize_classification)
taxize_classification$genome <- gsub(".[0-9]+$", "", rownames(taxize_classification))
taxize_classification <- taxize_classification[-3]
taxize_classification <- taxize_classification[-which(taxize_classification$name %in% "cellular organisms"),]
taxize_classification <- reshape(taxize_classification, timevar = "rank", idvar = "genome", direction = "wide")[1:8]
taxize_classification <- taxize_classification[-1]
colnames(taxize_classification) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxize_classification) <- ncbi_taxonomy$genome
write.table(taxize_classification, file="C:/Users/matti/VBShare/NCBI_genomes_tax_table.txt", quote=F, sep="\t")

#construct a tax table and otu table for the tree of both ncbi and metagenome genomes
ncbi_genomes_taxonomy <- read.csv("C:/Users/matti/VBShare/NCBI_genomes_tax_table.txt", quote="", sep="\t")
ncbi_genomes_taxonomy$label <- rownames(ncbi_genomes_taxonomy)
bin_tip_names <- concatenated_tree$tip.label[which(grepl("Bin", concatenated_tree$tip.label))]
ncbi_genomes_taxonomy <- rbind.fill(ncbi_genomes_taxonomy, data.frame(label=bin_tip_names))
ncbi_genomes_taxonomy$source <- ifelse(grepl("GCA", ncbi_genomes_taxonomy$label),"NCBI","metagenome") 

#plot a tree with ggtree
new_ggtree <- as_data_frame(phy_tree(concatenated_tree))
all_taxonomy <- ncbi_genomes_taxonomy
all_taxonomy$fulltax <- apply(all_taxonomy[,2:7], 1, paste, collapse = "_")
all_taxonomy[all_taxonomy$source == "metagenome",]$fulltax <-  all_taxonomy[all_taxonomy$source == "metagenome",]$label
all_taxonomy$fulltax <- gsub("NA_", "", all_taxonomy$fulltax)
new_ggtree <- full_join(new_ggtree, all_taxonomy, by="label")

#change the phylum level to class level for Proteobacteria
levels(new_ggtree$Phylum) <- c(levels(new_ggtree$Phylum), levels(new_ggtree$Class[new_ggtree$Phylum %in% "Proteobacteria"]))
new_ggtree$Phylum[new_ggtree$Phylum %in% "Proteobacteria"] <- new_ggtree$Class[new_ggtree$Phylum %in% "Proteobacteria"]

new_ggtree2 <- as.treedata(new_ggtree)

image <- ggtree(new_ggtree2, ladderize = T) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab(aes(label=new_ggtree$Phylum)) + geom_tippoint(aes(alpha=new_ggtree$source), color="red", shape=20, size=5) + scale_alpha_manual(values=c("1", "0"))
ggsave(file = "C:/Users/matti/VBShare/CPR_phyla_image4.pdf", plot=image, units="mm", width=2000, height=4000, limitsize = F)
