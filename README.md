# deconcat
Pipeline for the identification, resolution and validation of biological plasmid multimers.

## Installation
1. Install and activate conda environment with required programs and Python packages.
   
```
conda env create -n deconcat --file deconcat_env.yml
```

2. Download the replication genes database from the MOB-suite database, [manually](https://zenodo.org/records/10304948/files/data.tar.gz?download=1) or via command-line:
```
# Download DB
wget https://zenodo.org/records/10304948/files/data.tar.gz
# Decompress file and move to folder 
tar -xf data.tar.gz; cd data/
# Remove all databases but replication genes DB
ls | grep -xv "rep.dna.fas" | xargs rm
# Move DB to custom path
mv rep.dna.fas /new/db/path/
```

3. Run program with test data. Check execution is successful.

```
python3 deconcat.py --fasta_file /path/to/sample_plasmid.fasta --assembly_file /path/to/sample.fasta --fastq_file /path/to/sample.fastq.gz --circ_db /path/to/rep.dna.fas
```

## Usage
```
usage: deconcat.py [-h] --fasta_file FASTA_FILE [--assembly_file ASSEMBLY_FILE] [--fastq_file FASTQ_FILE] [--out_path OUT_PATH] --circ_db CIRC_DB
                   [--threads THREADS] [--max_len MAX_LEN] [--plasmid_id PLASMID_ID]

options:
  -h, --help            show this help message and exit
  --fasta_file FASTA_FILE
                        Path to plasmid fasta file
  --assembly_file ASSEMBLY_FILE
                        Path to complete assembly fasta file
  --fastq_file FASTQ_FILE
                        Path to long read fastq.gz
  --out_path OUT_PATH   Output path (default: plasmid_file_monomer)
  --circ_db CIRC_DB     Circlator DB path
  --threads THREADS     Number of threads (default: 6)
  --max_len MAX_LEN     Maximum plasmid length (default: 200000)
  --plasmid_id PLASMID_ID
                        Plasmid id (contig name) in complete assembly fasta (default: same id as in plasmid fasta file)
```

## Input
### Required Inputs:
- `--fasta_file`: A single-FASTA file of the plasmid to be analyzed (multi-FASTA files are not supported).
- `--circ_db`: The path to a database of Rep genes.

### Optional Parameters:
- `--assembly_file ASSEMBLY_FILE`: Path to a complete assembly FASTA file for the biological multimer validation step.
- `--fastq_file FASTQ_FILE`: Path to a long-read fastq.gz file for the biological multimer validation step.
- `--out_path OUT_PATH`: Output directory path (default: `plasmid_file_monomer`).
- `--threads THREADS`: Number of threads to use (default: 6).
- `--max_len MAX_LEN`: Maximum plasmid length (default: 200,000 bp). Plasmids longer than this size will be discarded.
- `--plasmid_id PLASMID_ID`: Plasmid ID (contig name) in the complete assembly FASTA file (default: same ID as in the plasmid FASTA file). Required for step 3. The contig ID in the plasmid FASTA must match the plasmid contig ID in the full genome assembly FASTA file.

## Output
After running the program, a folder named either as specified or in the format `"<sample_name>_monomer"` will be created. Inside this folder, you will find:

- **Consensus FASTA**: The resolved multimer sequence in its monomer form.
- **Report File**: A readable text file containing basic information about the process.
- **Filtered Reads (`fastq.gz`)**: A compressed FASTQ file containing reads longer than 110% of the monomer length, supporting the biological existence of the monomer.
- **Read Length Plot**: A plot visualizing the lengths of reads mapping to the multimer.
- **Subfolder `blastn`**: Contains intermediate files from the BLASTN execution and logs of the process.
- **Subfolder `circlator`**: Contains intermediate files from Circlator execution and logs.
- **Subfolder `monomers`**:
  - A multi-FASTA file of the multimer split into individual monomer sequences, with positions relative to the recircularized plasmid.
  - An alignment file (`.aln`) of complete monomers used to reconstruct the consensus.
- **Subfolder `multimer_mapping`**: Contains intermediate files from the mapping process and program logs.

## Schematic workflow
![deconcat_stacked](https://github.com/user-attachments/assets/0f0545e1-573e-4b95-9513-98148ae9144c)

## Contact
Comments and suggestions are always welcome at: rpalcab@gmail.com

## Next steps
To do:
- [x] Emitir informe resumen
- [x] Generar logs de programas intermedios
- [x] Subir archivos de prueba
- [x] Subir conda env
- [x] Redactar README con info básica (workflow) y manual de uso
- [ ] Organizar en módulos
- [ ] Pasar a nextflow
