#!/bin/sh
# properties = {"type": "single", "rule": "DADA2", "local": false, "input": ["s3Combine"], "output": ["s5DADA2/DADA2_seqtab_nochim.rda", "s5DADA2/img/ploterrF.png", "s5DADA2/img/ploterrR.png"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 1, "cluster": {"name": "16s-pipeline.DADA2", "time": "1-0", "nodes": 1, "ntasks": 1, "mem": "2G", "partition": "high", "output": "slurmout/DADA2/%A.out", "error": "slurmout/DADA2/%A.err"}}
cd /home/xinyut/workflows/16s-pipeline && \
/home/xinyut/miniconda3/envs/16s-pipeline/bin/python3.7 \
-m snakemake s5DADA2/DADA2_seqtab_nochim.rda --snakefile /home/xinyut/workflows/16s-pipeline/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /home/xinyut/workflows/16s-pipeline/.snakemake/tmp.k1_5wywn s3Combine --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules DADA2 --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/home/xinyut/workflows/16s-pipeline/.snakemake/tmp.k1_5wywn/1.jobfinished" || (touch "/home/xinyut/workflows/16s-pipeline/.snakemake/tmp.k1_5wywn/1.jobfailed"; exit 1)

