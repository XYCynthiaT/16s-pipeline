library(dada2)
library(dplyr)
library(ggplot2)
library(DECIPHER)
packageVersion("dada2")

loger = function(msg){
    t = strftime(Sys.time(), "%Y-%m-%d %H-%M-%S")
    cat("[ ", t, " ] ", msg, "\\n")
}

# --> Step1. Getting ready
loger("Getting ready")
path <- snakemake@params[[1]]
truncLen_r1 <- snakemake@params[[2]]
truncLen_r2 <- snakemake@params[[3]]
# cutLens <- read.csv(snakemake@input[[1]])

# Read in fastq file names by r1 and r2
loger("Read in fastq file names")

fnfs = list.files(
    path, pattern=".*.r1.fastq$", full.names = TRUE
) %>%
    sort()
fnrs = list.files(
    path, pattern = ".*.r2.fastq$", full.names = TRUE
) %>%
    sort()
# Extract sample id"
# fastq file name format: {sample}.r1.fq
sampleID = fnfs %>%
    basename() %>% 
    strsplit("\\.") %>% 
    sapply(`[`, 1)

# --> Step2. Inspect read quality profiles
# Done by fastq

#  --> Step3. Filter and trim
# Place filtered files in filtered/ subdirectory
loger("Filter and trim")
filfs = file.path(path, "filtered", paste0("fasting_", sampleID, "_r1_filt.fq.gz"))
filrs = file.path(path, "filtered", paste0("fasting_", sampleID, "_r2_filt.fq.gz"))
names(filfs) = sampleID
names(filrs) = sampleID
out = filterAndTrim(
    fnfs, filfs, fnrs, filrs, 
    truncLen=c(truncLen_r1, truncLen_r2),
    # truncLen=c(cutLens$r1, cutLens$r2), # Ref to FastQC
    maxN = 0, # DADA2 requires no Ns
    maxEE=c(2, 2), # the maximum number of “expected errors” allowed in a read
    truncQ=2, 
    compress=TRUE, 
    multithread=TRUE # On Windows set F
)

# removing samples that is empty after filtering
loger("Removing samples that is empty after filtering")
filfs_reduced = c()
filrs_reduced = c()
for(i in seq_along(filfs)){
    if(!file.exists(filfs[i]) | !file.exists(filrs[i])){
        next
    }
    filfs_reduced = c(filfs_reduced, filfs[i])
    filrs_reduced = c(filrs_reduced, filrs[i])
}
filfs = filfs_reduced
filrs = filrs_reduced

# --> Step4. Learn the Error Rates
loger("Learn error rates")
errF = learnErrors(filfs, multithread = TRUE)
errR = learnErrors(filrs, multithread = TRUE)
loger("Generating plots of error rates")
ploterrF = plotErrors(errF, nominalQ = TRUE)
ploterrR = plotErrors(errR, nominalQ = TRUE)
ggsave("output/s5DADA2/img/ploterrF.png", plot = ploterrF, width = 8, height = 8, units = "in")
ggsave("output/s5DADA2/img/ploterrR.png", plot = ploterrR, width = 8, height = 8, units = "in")

# --> Step5. Sample Inference
loger("Sample inference, take some time...")
dadaFs = dada(filfs, err=errF, multithread=TRUE)
dadaRs = dada(filrs, err=errR, multithread=TRUE)

# --> Step6. Merge paired reads
loger("Merge paired reads")
mergers = mergePairs(dadaFs, filfs, dadaRs, filrs, verbose=TRUE)
# Inspect the merger df from the first sample
#head(mergers[[1]])

# --> Step7. Construct sequence table
loger("Construct seq table")
seqtab = makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
loger("Inspect distribuion of sequence lengths")
seqtab %>% 
    getSequences() %>%
    nchar() %>%
    table() 
# Sequences that are much longer or shorter than expected may be the result of non-specific priming.

# --> Step8. Remove chimeras
loger("Remove chimera")
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
print("passing rate:")
print(sum(seqtab.nochim)/sum(seqtab))

save(seqtab.nochim, file = "output/s5DADA2/DADA2_seqtab_nochim.rda")