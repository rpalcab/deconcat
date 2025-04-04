# Deconcat
Pipeline for the identification, resolution and validation of biological plasmid multimers.

## Installation
1. Install and activate conda environment with required programs and Python packages.
   
```
conda env create -n deconcat --file deconcat_env.yml
```

2. Run program with test data. Check execution is successful.

```
python3 bin/deconcat.py --fasta_file /path/to/sample.fasta --fastq_file /path/to/sample.fastq.gz
```

## Usage
```
usage: deconcat.py [-h] --fasta_file FASTA_FILE [--fastq_file FASTQ_FILE]
                   [--out_path OUT_PATH] [--threads THREADS]
                   [--max_len MAX_LEN]

Identification, resolution and biological confirmation of plasmid multimers.

optional arguments:
  -h, --help            show this help message and exit

Input arguments:
  --fasta_file FASTA_FILE
                        Path to plasmid fasta file. REQUIRED
  --fastq_file FASTQ_FILE
                        Path to long read fastq.gz

Output arguments:
  --out_path OUT_PATH   Output path. Default: plasmid_file_monomer

Tunning arguments:
  --threads THREADS     Number of threads (default: 6)
  --max_len MAX_LEN     Maximum plasmid length (default: 200000)
```

## Arguments
### Required arguments:
- `--fasta_file`: Sample assembly in FASTA format.

### Optional arguments:
- `--fastq_file`: Path to a long-read fastq.gz file for the biological multimer validation step.
- `--out_path OUT_PATH`: Output directory path (default: `sample_monomer`).
- `--threads THREADS`: Number of threads to use (default: 6).
- `--max_len MAX_LEN`: Maximum plasmid length (default: 200,000 bp). Plasmids longer than this size will be discarded.

## Output
After running the program, a folder named either as specified or in the format `"<sample_name>_monomer"` will be created. Inside this folder, you will find:

- **Assembly FASTA (`_corr.fasta`)**: The resolved multimer sequence in its monomer form.
- **Monomer FASTA (`_consensus.fasta`)** : The resolved multimer sequence in its monomer form.
- **Filtered Reads (`.filtered.fastq.gz`)**: A compressed FASTQ file containing reads longer than 150% of the monomer length, supporting the biological existence of the monomer.
- **Unfiltered Reads (`.sorted.fastq.gz`)**: A compressed FASTQ file containing all reads mapped to the multimer, regardless of their length.
- **Read Length Plot (.png)**: A plot visualizing the lengths of reads mapping to the multimer.
- **Subfolder `blastn`**: Contains intermediate files from the BLASTN execution and logs of the process.
- **Subfolder `monomers`**:
  - A multi-FASTA file of the multimer split into individual monomer sequences, with coordinates of the unresolved plasmid in the header.
  - An alignment file (`.aln`) of complete monomers used to reconstruct the consensus.
- **Subfolder `multimer_mapping`**: Contains intermediate files from the mapping process and program logs.

## Schematic workflow
![deconcat_stacked](https://github.com/user-attachments/assets/0f0545e1-573e-4b95-9513-98148ae9144c)
