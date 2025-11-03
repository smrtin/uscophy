#!/usr/bin/env python3
import sys
import subprocess
from snakemake.utils import read_job_properties
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print("Usage: submit.py <jobscript>")
        sys.exit(1)

    jobscript = sys.argv[1]

    # Parse Snakemake job properties from the job script
    job_properties = read_job_properties(jobscript)

    # Extract properties with defaults
    job_name = job_properties.get("rule", "snakemake_job")
    threads = job_properties.get("threads", 1)
    
    Path("Cluster_job_logs").mkdir(parents=True, exist_ok=True)
    # Build qsub command with parsed properties

    qsub_cmd = [
    "qsub",
    "-terse",
    "-o", "Cluster_job_logs" ,
    "-e", "Custer_job_logs" ,
    "-V", 
    "-S", "/bin/bash", 
    "-cwd" ,
    "-j", "n" ,
    "-q" , "medium.q,fast.q,large.q,small.q",
    "-N", job_name,
    "-pe", "smp", str(threads),
    jobscript
    ]

    print("Submitting job with command:", " ".join(qsub_cmd))

    # Submit the job and capture output and errors
    try:
        result = subprocess.run(qsub_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print("Submission output:", result.stdout.strip())
    except subprocess.CalledProcessError as e:
        print("Error submitting job:", e.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

