#Created by N.R.Ratib on March 22, 2020
##This script cleans up the complied mutations output (.tsv) from the Breseq compare function to generate a presence-absence matrix of mutations in clones and plot that information as heatmap panels separated by time points. 

library(prodlim)
library(plotrix)

df <- read.delim(file='/Volumes/GoogleDrive/My Drive/Job Search/sample code/T1_compare.tsv', sep='\t', header=TRUE, as.is=TRUE)

#data columns from .tsv output we want to keep
keep <- c("aa_new_seq", "aa_position", "aa_ref_seq", "codon_new_seq", "codon_position", "codon_ref_seq", "del_end", "del_start", "duplication_size", "gene_name", "gene_position", "gene_product", "gene_strand", "genes_inactivated", "genes_overlapping", "genes_promoter", "insert_position", "mediated", "mutation_category", "new_seq", "position", "repeat_length", "repeat_name", "repeat_new_copies", "repeat_ref_copies", "repeat_seq", "size", "snp_type", "title")

#make sure all 29 data columns are present
match(keep, colnames(df))

#subset dataframe to keep only columns we want
df <- df[,match(keep, colnames(df))]

#multiple clones have the same mutation so this identifies all unique mutations among clones
urows <- unique(df[,1:ncol(df)-1])

#matching the unique mutations to clones
rmatch <- sapply(1:nrow(df), function(x){
	row.match(df[x,1:ncol(df)-1], urows)
})

#clone names for mutation presence-absence matrix
ind_names <- unique(df$title)

#empty matrix where columns are clones and rows are unique mutations
df_temp <- matrix(0, ncol=length(ind_names), nrow=nrow(urows))

#add data to matrix where 0 = mutation absent in clone and 1 = present
df_snps <- sapply(1:length(ind_names), function(x){
	ind <- which(df$title == ind_names[x])
	c <- rep(0, nrow(df_temp))
	c[rmatch[ind]] <- 1
	df_temp[,x] <- c
})

#generate mutation names from columns in df
gene_names <- sapply(1:nrow(urows), function(x){
	paste(urows$gene_name[x], "_", urows[x,3], urows[x,2], urows[x,1], sep="")
})

#name columns and rows in matrix
colnames(df_snps) <- ind_names
rownames(df_snps) <- gene_names

#makes a list of the time point each clone belongs to 
temp_names <- sapply(colnames(df_snps), function(i){
  name <- unlist(strsplit(i, '_'))[1]
})

#makes a list of unique time points
time_points <- unique(temp_names)


df_zo <- df_snps
rownames(df_zo) <- 1:nrow(df_zo)

freq_counts <- sapply(time_points, function(j){
  ixs <- which(temp_names == j) #columns from time point
  df_temp <- df_zo[,ixs] #subset df_zo into timepoint
  df_temp <- df_temp[rowSums(df_temp == 1) > 0,] #pull out rows that have mutations
  snp_freqs <- rowSums(df_temp == 1) #count frequency of 
  sort_vals <- sort(snp_freqs, index.return=T)
  rev(as.numeric(row.names(df_temp)[sort_vals$ix]))
})

temp_ixs <- sapply(1:length(freq_counts), function(i){
  prevs_snps <- freq_counts[1:(i-1)]
  prevs_snps <- unlist(prevs_snps)
  ixs <- which(!(freq_counts[[i]] %in% prevs_snps))
  unlist(freq_counts[[i]])[ixs]
})
snp_order <- c(unlist(freq_counts[[1]]), unlist(temp_ixs))
snp_names_reorder <- gene_names[snp_order]

str_comparisons <- expand.grid(ind_names, ind_names, stringsAsFactors=F)

distances <- apply(str_comparisons, 1, function(i){
    vec1 <- df_snps[,i[1]]
    vec2 <- df_snps[,i[2]]
    ixs <- which(!(vec1 == 'N' | vec2 == 'N'))
    sum(vec1[ixs] != vec2[ixs])
})

# calculating distance measures
dist_matrix <- matrix(distances, nrow = ncol(df_zo)) / nrow(df_zo)
rownames(dist_matrix) <- colnames(df_zo)
colnames(dist_matrix) <- colnames(df_zo)


# making cluster object here
hc_obj <- hclust(dist(dist_matrix))

# extracting strain order from the cluster object
strain_order_all <- sapply(hc_obj$labels[hc_obj$order], function(j){
	unlist(strsplit(j, split='_'))[1]
	})
	
# breaking up strain order by time points 
strain_order_names <- sapply(time_points, function(x){
	names(strain_order_all[which(strain_order_all == x)])
	})

inds <- sapply(time_points, function(x){
	sum(temp_names == x)
	})

strain_order <- match(unlist(strain_order_names), colnames(df_zo))

df_image <- df_zo[snp_order,strain_order]
colnames(df_image) <- colnames(df_zo)[strain_order]
rownames(df_image) <- snp_names_reorder

write.table(df_image, file='df_image.txt', sep='\t')
write.table(urows[snp_order,], file='mutations.txt', sep='\t')

multi_panel_layout <- matrix(c(rep(unlist(
  sapply(1:length(inds), function(i){
    rep(i, ceiling(inds/10)[i])
  })), each=10)), nrow=10)

r_inds <- cumsum(c(1, inds))

#######################
image_plot <- function(k, colors, col_order) {
  genotypes <- cutree(hc_obj, k)
  match_genotypes <- match(genotypes, col_order)
  names(match_genotypes) <- names(genotypes)
  col_pal <- c(rep('gray95', length(colors)), colors, rep('white', length(colors)))
  col_pal <- matrix(col_pal, nrow=length(colors))
  file_name <- paste('image_strains_', k, '.pdf', sep='')
  pdf(file_name, height=20, width=40)
  layout(multi_panel_layout)
  par(mar=c(2,1,4,1), oma=c(2,9,2,2))
  jnk <- sapply(1:length(time_points), function(i){
    ixs <- r_inds[i]:(r_inds[(i+1)]-1) # -1 because it's start1 to start2-1
    df_temp <- df_image[,ixs]+1 # +1 to change 0 to 1 for color matrix index
    geno_temp <- match_genotypes[match(colnames(df_temp), names(match_genotypes))]
    color_matrix <- sapply(1:ncol(df_temp), function(x){
      sapply(1:nrow(df_temp), function(y){
        col_pal[geno_temp[x],df_temp[y,x]]
        })
      })
    snp_names_numbered <- sapply(1:length(snp_names_reorder), function(i){paste(snp_names_reorder[i], i, sep=' - ')})
    color2D.matplot(df_temp, cellcolors=color_matrix, 
      xaxt='n', yaxt='n', main=time_points[i], axes=F, border=NA, cex.main=2)
    if (i == 1){
      axis(2, at=1:nrow(df_image)-0.5, labels = rev(snp_names_numbered), cex.axis=0.5,las=2)
      }
     else {
     	axis(2, at=1:nrow(df_image)-0.5, labels = rev(1:length(snp_names_reorder)), cex.axis=0.5,las=2)
     }
  })
  dev.off()
}



