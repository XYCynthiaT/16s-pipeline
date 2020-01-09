library(dada2)
library(dplyr)
library(DECIPHER)

loger = function(msg){
    t = strftime(Sys.time(), "%Y-%m-%d %H-%M-%S")
    cat("[ ", t, " ] ", msg, "\n")
}

load(snakemake@input[[1]])
load(snakemake@input[[2]])

# --> Step10. Assign taxonomy
loger("Assign taxonomy, take a long time...")
dna = seqtab.nochim %>%
    getSequences %>%
    DNAStringSet() # Create a DNAStringSet from ASVs
ids = IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE) #use all processors
ranks = c("domain", "phylum", "class", "order", "family", "genus", "species")
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid = t(sapply(ids, function(x){
    m = match(ranks, x$rank)
    taxa = x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA 
    taxa 
}))
colnames(taxid) = ranks
rownames(taxid) = getSequences(seqtab.nochim)

save(dna, taxid, file = "output/s5DADA2/DADA2_taxaid.rda")
