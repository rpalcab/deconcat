# %%
import os
import sys
import gzip
from Bio import SeqIO
from Bio import Align
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# %%
def read_fasta(fasta_file, max_len=200000, query_len=1200):
    d_fasta = {}
    if os.path.isfile(fasta_file) and os.path.getsize(fasta_file) > 0:
        n_conts = 0
        pl_len = 0
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            d_fasta[seq_record.id] = seq_record.seq
            n_conts += 1
            query = str(seq_record.seq[0:query_len])
            if len(seq_record) > pl_len:
                pl_len = len(seq_record)
        if n_conts == 1 and pl_len < max_len:
            return True, query, d_fasta
        return False, 'Multifasta or maximum length exceeded', d_fasta
    return False, 'Empty fasta or wrong path', d_fasta

# %%
def sliding_window(pos_list, window=2):   
    positions = []
    for i in range(len(pos_list)- window + 1):
        positions.append(pos_list[i:i+window])
    return positions

# %%
def extract_monomers(df_reps, d_fasta, contig_id):
    start_positions = list(df_reps[9])
    start_positions.append(len(d_fasta[contig_id])+1)

    positions = sliding_window(start_positions, 2)
    monomers = {}
    largest_len = 0
    largest_id = ""
    for start, end in positions:
        mono_id = f'{contig_id}_{start}_{end-1}'
        mono_seq = d_fasta[contig_id][start-1:end-1]
        monomers[mono_id] = mono_seq
        if len(mono_seq) > largest_len:
            largest_len = len(mono_seq)
            largest_id = mono_id

    return monomers, largest_id

# %%
def similarity_check(monomers, largest_id):
    list_check = [i for i in monomers.keys() if i != largest_id]
    l_complete = [largest_id]
    l_partial = []
    status = True

    for i in list_check:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aln = aligner.align(monomers[largest_id], monomers[i])[0]
    
        identity = aln.score / len(monomers[largest_id])
        aligned_query = aln.aligned[0]
        aligned_length = sum([end - start for start, end in aligned_query])
        qcov = aligned_length / len(monomers[largest_id])

        if identity >= 0.9:
            if qcov >= 0.9:
                l_complete.append(i)
            else:
                l_partial.append(i)
        else:
            status = False

    return status, l_complete, l_partial


# %%
def extract_read_lengths_and_filter_reads(fastq_gz_file, monomer_length):
    read_lengths = []
    filtered_reads = []
    threshold_length = monomer_length * 1.1
    with gzip.open(fastq_gz_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_length = len(record.seq)
            read_lengths.append(read_length)
            if read_length > threshold_length:
                filtered_reads.append(record)
    return read_lengths, filtered_reads

# %%
def plot_read_lengths(read_lengths, output_file, monomer_length, max_length):
    plt.figure(figsize=(10, 6))
    plt.hist(read_lengths, bins=100, alpha=0.75, edgecolor='black')
    plt.axvline(x=monomer_length, color='red', linestyle='--', linewidth=2)
    plt.text(monomer_length, max(plt.ylim()) * 0.9, f'X={monomer_length}', color='red', ha='right')
    plt.axvline(x=max_length, color='blue', linestyle='--', linewidth=2)
    plt.text(max_length, max(plt.ylim()) * 0.8, f'X={max_length}', color='blue', ha='right')
    plt.yscale('log')
    plt.title('Distribution of Read Lengths')
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()

# %%
def save_filtered_reads(filtered_reads, output_file):
    with gzip.open(output_file, "wt") as handle:
        SeqIO.write(filtered_reads, handle, "fastq")

# %%
def parse_arguments():

    parser = argparse.ArgumentParser(description="Identification, resolution and biological confirmation of plasmid multimers.")
    
    # Input paths
    parser.add_argument("--fasta_file", type=str, required=True, 
                        help="Path to plasmid fasta file")
    
    parser.add_argument("--assembly_file", type=str, default=None, 
                        help="Path to complete assembly fasta file")
    
    parser.add_argument("--fastq_file", type=str, default=None, 
                        help="Path to long read fastq.gz")
    
    parser.add_argument("--out_path", type=str, default=None, 
                        help="Output path (default: plasmid_file_monomer")

    parser.add_argument("--circ_db", type=str, required=True, 
                        help="Circlator DB path")

    parser.add_argument("--threads", type=int, default=6, 
                        help="Number of threads (default: 6)")
    
    parser.add_argument("--max_len", type=int, default=200000, 
                        help="Maximum plasmid length (default: 200000)")
    
    parser.add_argument("--plasmid_id", type=str, default=None, 
                        help="Plasmid id (contig name) in complete assembly fasta (default: same id as in plasmid fasta file)")
    
    args = parser.parse_args()
    
    if args.out_path is None:
        args.out_path = os.path.splitext(args.fasta_file)[0] + '_monomer'

    return args

# %%
def main():
    args = parse_arguments()

    ## 0. inputs
    # max plasmid length
    max_len=args.max_len
    # plasmid fasta
    fasta_file = args.fasta_file
    fasta_basename = os.path.splitext(os.path.basename(fasta_file))[0]
    # assembly fasta
    assembly_file = args.assembly_file
    # comprobar que el fasta tiene una sola seq menor que max_len
    check_status, query, d_fasta = read_fasta(fasta_file, max_len)
    if check_status != True:
        print(f'Please check fasta file: {check_status}')
        sys.exit(1)
    # fastq
    fastq_file = args.fastq_file
    # output folder
    out_path = args.out_path
    if not os.path.exists(out_path): os.mkdir(out_path)
    # Databases
    # circlator db
    circ_db = args.circ_db
    # threads
    threads = args.threads
    # plasmid id
    plasmid_id = args.plasmid_id

    print('Repetition identification')
    ## 1. identificación de las repeticiones

    # 1.1 recircularizar fasta
    circ_path = os.path.join(out_path, 'circlator')
    circ_fasta = os.path.join(circ_path, f'{fasta_basename}.fasta')
    if not os.path.exists(circ_path): os.mkdir(circ_path)
    circ_cmd = ['circlator', 'fixstart', '--genes_fa', circ_db, fasta_file, os.path.join(circ_path, fasta_basename)]
    subprocess.run(circ_cmd)
    check_status, query, d_fasta = read_fasta(circ_fasta, max_len)

    # 1.2 creación de blastdb a partir de fasta circularizado
    blast_path = os.path.join(out_path, 'blastn')
    if not os.path.exists(blast_path): os.mkdir(blast_path)
    makeblastdb_cmd = ['makeblastdb', '-in', circ_fasta, '-dbtype', 'nucl', '-out', circ_fasta]
    subprocess.run(makeblastdb_cmd)

    # 1.3 guardar query en un archivo
    query_fasta = os.path.join(blast_path, 'query.fasta')
    with open(query_fasta, 'w') as f:
        f.write(f">query_sequence\n{query}\n")

    # 1.4 localizar repeticiones con blastn
    blast_out = os.path.join(blast_path, f'{fasta_basename}_blast.out')
    blastn_cmd = ['blastn', '-query', query_fasta, '-db', circ_fasta, '-strand', 'plus', 
                '-outfmt', '6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send sframe evalue bitscore',
                '-perc_identity', '90', '-qcov_hsp_perc', '95']
    sort_cmd = ['sort', '-n', '-k', '10,11']
    sp = subprocess.run(blastn_cmd, check=True, capture_output=True)
    with open(blast_out, 'w') as of: 
        subprocess.run(sort_cmd, input=sp.stdout, stdout=of)

    ## 2. creación del monómero

    monomer_path = os.path.join(out_path, 'monomers')
    if not os.path.exists(monomer_path): os.mkdir(monomer_path)

    # 2.1 Extracción de monómeros
    df_reps = pd.read_table(blast_out, header=None)
    contig_id = list(d_fasta.keys())[0]

    monomers, largest_id = extract_monomers(df_reps, d_fasta, contig_id)
    if len(monomers.keys()) == 1:
        print('The plasmid introduced is a monomer')
        sys.exit(0)

    # 2.2 Comprobación de similitud
    """
        pairwise de todos los monómeros frente al de mayor tamaño.
        si la identidad es menor de 90%, no son monómeros
        si es mayor y el qcov es menor del 90%, no usamos en alineamiento múltiple
        si es mayor, añadimos al alineamiento múltiple
    """
    status, l_complete, l_partial = similarity_check(monomers, largest_id)
    if status == False:
        print('The plasmid introduced is a monomer')
        sys.exit(0)
    
    # 2.3 Alineamiento múltiple de las secuencias que pasan los criterios "completa"
    monomer_file = os.path.join(monomer_path, 'complete_monomers.fasta')
    with open(monomer_file, 'w') as of: 
        for i in l_complete:
            of.write(f'>{i}\n')
            of.write(f'{monomers[i]}\n')

    aln_file = os.path.join(monomer_path, 'complete_monomers.aln')
    mafft_cmd = ['mafft', '--adjustdirectionaccurately', '--thread', str(threads), monomer_file]
    with open(aln_file, 'w') as of:
        subprocess.run(mafft_cmd, stdout=of)

    # 2.4 Generación del consenso
    cons_file = os.path.join(out_path, f'{fasta_basename}_consensus.fasta')
    cons_cmd = ['cons', '-sequence', aln_file, '-outseq', cons_file, '-name', fasta_basename]
    subprocess.run(cons_cmd)

    print('The plasmid introduced is a multimer: ')
    print(f'Find monomer consensus on: {cons_file}')
    print(f'Find divided multimer on: {monomer_file}')
    # print(f'Find statistical report on: {}')

    if assembly_file != None and fastq_file != None:
        print('Validating multimer according to long reads')
        ## 3. validación del multímero biológico
        # Mapeo
        map_path = os.path.join(out_path, 'multimer_mapping')
        if not os.path.exists(map_path): os.mkdir(map_path)

        if args.plasmid_id is None:
            plasmid_id = list(d_fasta.keys())[0]

        # mapeo original
        sam_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.sam')
        minimap_cmd = ['minimap2', '--secondary=no', '-t', str(threads), '-ax', 'lr:hq', '-o', sam_file, assembly_file, fastq_file]
        subprocess.run(minimap_cmd)

        # paso a bam eliminando lecturas no mapeadas
        bam_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.bam')
        samtools_filter_cmd = ['samtools', 'view', '-h', '-@', str(threads), '-F', '4', '-bS', sam_file]
        with open(bam_file, 'w') as of:
            subprocess.run(samtools_filter_cmd, stdout=of)

        # ordenar
        bam_sorted_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.sorted.bam')
        samtools_sort_cmd = ['samtools', 'sort', '--threads', str(threads), bam_file, '-o', bam_sorted_file]
        subprocess.run(samtools_sort_cmd)

        # generar index
        samtools_index_cmd = ['samtools', 'index', '-@', str(threads), bam_sorted_file]
        subprocess.run(samtools_index_cmd)

        # extraer solo lecturas que mapean con el plásmido en cuestión
        sam_concat_file = os.path.join(map_path, f'{fasta_basename}_concat.sam')
        samtools_extract_cmd = ['samtools', 'view', '-h', bam_sorted_file, plasmid_id]
        with open(sam_concat_file, 'w') as of:
            subprocess.run(samtools_extract_cmd, stdout=of)

        # pasar sam a bam
        bam_concat_file = os.path.join(map_path, f'{fasta_basename}_concat.bam')
        samtools_bam_cmd = ['samtools', 'view', '-h', '-@', str(threads), '-F', '4', '-bS', sam_concat_file]
        with open(bam_concat_file, 'w') as of:
            subprocess.run(samtools_bam_cmd, stdout=of)

        # generar index
        samtools_index_cmd = ['samtools', 'index', '-@', str(threads), bam_concat_file]
        subprocess.run(samtools_index_cmd)

        fastq_sorted_file = os.path.join(map_path, f'{fasta_basename}_concat.sorted.fastq.gz')
        samtools_fastq_cmd = ['samtools', 'fastq',  bam_concat_file]
        gzip_cmd = ['gzip']
        sp = subprocess.run(samtools_fastq_cmd, check=True, capture_output=True)
        with open(fastq_sorted_file, 'w') as of:
            subprocess.run(gzip_cmd, input=sp.stdout, stdout=of)

        # Plot
        monomer_len = len(monomers[largest_id])
        max_len = len(d_fasta[list(d_fasta.keys())[0]])
        read_lengths, filtered_reads = extract_read_lengths_and_filter_reads(fastq_sorted_file, monomer_len)
        output_file_plot = os.path.join(out_path, f'{fasta_basename}_read_plot.png')
        plot_read_lengths(read_lengths, output_file_plot, monomer_len, max_len)

        print(f"Plot saved as {output_file_plot}")
    
    else:
        print('Skipping multimer validation.')
        print('To validate, introduce fastq.gz and complete assembly fasta.')
    
    print('End of pipeline.')
    sys.exit(0)

# %%
if __name__ == "__main__":
    main()


