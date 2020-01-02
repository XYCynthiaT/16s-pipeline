# setwd(dirname(parent.frame(2)$ofile))
# working dir: ~/fasting
# input dir: ~/s5Dada2
# output dir: ~/s5Dada2

pkgs <- c("Biostrings", "dplyr")
for (pkg in pkgs) {
    library(pkg, character.only = T)
}

loger = function(msg){
    t = strftime(Sys.time(), "%Y-%m-%d %H-%M-%S")
    cat("[ ", t, " ] ", msg, "\\n")
}

load(snakemake@input[[1]])
load(snakemake@input[[2]])

# Arrange rows of seqtab.nochim -------------------------------------------

loger("Arrange rows of seqtab.nochim")
orderIndex <- rownames(seqtab.nochim) %>%
    as.integer() %>%
    order
 seqtab.nochim <- seqtab.nochim[orderIndex,]   

# Combine reverse complement ASV counts -----------------------------------

loger("Combine reverse complement ASV counts")
duplicate <- intersect(dna, reverseComplement(dna))
index <- vector("integer")
rc_index <- vector("integer")
for (i in seq_along(duplicate)) {
    if(sum(length(index), length(rc_index)) != length(duplicate)) {
        if(!(i %in% rc_index)){
            index = append(index, i)
            j = which(reverseComplement(duplicate[i]) == duplicate)
            rc_index = append(rc_index, j)
        }
    }
}
indexInRaw <- which(colnames(seqtab.nochim) %in% duplicate[index])
rc_indexInRaw <- which(colnames(seqtab.nochim) %in% duplicate[rc_index])
uniqSeqtab <- seqtab.nochim[,indexInRaw]+seqtab.nochim[,rc_indexInRaw]
uniqSeqtab <- cbind(uniqSeqtab, seqtab.nochim[,-append(indexInRaw, rc_indexInRaw)]) %>%
    t()
seqs <- rownames(uniqSeqtab) %>%
    DNAStringSet()
rownames(uniqSeqtab) <- paste0("ASV", 1:nrow(uniqSeqtab))


# Feature data (fdata) ----------------------------------------------------

loger("Clean up taxaid for feature data")
orderedTaxid <- taxid[as.character(seqs),]
# <Inspecting steps>
seqs_taxid <- rownames(orderedTaxid) %>%
    DNAStringSet()
print("Does DNA seqs(rowname) of taxid match DNA seqs(rowname) of seqtab.nochim?")
print(identical(seqs, seqs_taxid)) # check the rows of fdata matches the rows of edata
# <inspecting steps end>
rownames(orderedTaxid) <- paste0("ASV", 1:nrow(orderedTaxid))

# Export to .csv ----------------------------------------------------------

loger("Export to .csv files")
write.csv(uniqSeqtab, file = "s6TidyData/ASV_samplCounts.csv")
write.csv(orderedTaxid, file = "s6TidyData/ASV_taxanomy.csv")
