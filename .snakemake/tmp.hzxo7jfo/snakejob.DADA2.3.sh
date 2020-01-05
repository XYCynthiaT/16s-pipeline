#!/bin/sh
# properties = {"type": "single", "rule": "DADA2", "local": false, "input": ["s3Combine"], "output": ["s5DADA2/DADA2_seqtab_nochim.rda", "s5DADA2/img/ploterrF.png", "s5DADA2/img/ploterrR.png"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 3, "cluster": {"name": "16s-pipeline.DADA2.", "time": "1-0", "n": 1, "ntasks": "1", "mem": "2G", "partition": "med", "output": "slurmout/DADA2/.out", "error": "slurmout/DADA2/.err"}}
cd /home/xinyut/workflows/16s-pipeline && \
/home/xinyut/miniconda3/envs/16s-pipeline/bin/python3.7 \
-m snakemake s5DADA2/DADA2_seqtab_nochim.rda --snakefile /home/xinyut/workflows/16s-pipeline/microbiome.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /home/xinyut/workflows/16s-pipeline/.snakemake/tmp.hzxo7jfo s3Combine --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules DADA2 --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/home/xinyut/workflows/16s-pipeline/.snakemake/tmp.hzxo7jfo/3.jobfinished" || (touch "/home/xinyut/workflows/16s-pipeline/.snakemake/tmp.hzxo7jfo/3.jobfailed"; exit 1)

