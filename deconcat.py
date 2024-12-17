#!/usr/bin/env python

import os
import sys
import gzip
from Bio import SeqIO
from Bio import Align
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import csv

# %%
def read_fasta(fasta_file, max_len=200000, query_len=1200):
    d_fasta = {}
    if os.path.isfile(fasta_file) and os.path.getsize(fasta_file) > 0:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if len(seq_record) < max_len and len(seq_record) > query_len:
                query = str(seq_record.seq[0:query_len])
                d_fasta[seq_record.id] = [seq_record.seq, query]
            else:
                d_fasta[seq_record.id] = [seq_record.seq, '']
        return d_fasta
    exit('Empty fasta or wrong path')

def write_fasta(output_file, d_seq):
    with open(output_file, "w") as output_handle:
        for k, v in d_seq.items():
            output_handle.write(k + "\n")
            output_handle.write(v + "\n")
    return None

def create_blastdb(out_path, contig, sequence):
    blast_path = os.path.join(out_path, 'blastn')
    if not os.path.exists(blast_path): os.mkdir(blast_path)
    contig_path = os.path.join(blast_path, f"{contig}.fasta")
    write_fasta(contig_path, {f'>{contig}': sequence})
    mkblast_err = os.path.join(blast_path, f'makeblastdb_{contig}.log')
    makeblastdb_cmd = ['makeblastdb', '-in', contig_path, '-dbtype', 'nucl', '-out', contig_path]
    with open(mkblast_err, 'w') as errf:
        subprocess.run(makeblastdb_cmd, stdout=errf, stderr=errf)
    return blast_path, contig_path

def repetition_search(blast_path, contig_path, query_fasta, fasta_basename, contig):
    blast_out = os.path.join(blast_path, f'{fasta_basename}_{contig}_blast.out')
    blast_err = os.path.join(blast_path, f'blast_{contig}.log')
    blastn_cmd = ['blastn', '-query', query_fasta, '-db', contig_path, '-strand', 'plus', 
                '-outfmt', '6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send sframe evalue bitscore',
                # '-perc_identity', '90', '-qcov_hsp_perc', '95']
                '-perc_identity', '80', '-qcov_hsp_perc', '80']
    sort_cmd = ['sort', '-n', '-k', '10,11']
    sp = subprocess.run(blastn_cmd, check=True, capture_output=True)
    with open(blast_out, 'w') as stdf, open(blast_err, 'w') as errf: 
        subprocess.run(sort_cmd, input=sp.stdout, stdout=stdf, stderr=errf)
    return blast_out

# %%
def sliding_window(pos_list, window=2):   
    positions = []
    for i in range(len(pos_list)- window + 1):
        positions.append(pos_list[i:i+window])
    return positions

# %%
def extract_monomers(df_reps, sequence, contig_id):
    start_positions = list(df_reps[9])
    start_positions.append(len(sequence)+1)

    positions = sliding_window(start_positions, 2)
    monomers = {}
    largest_len = 0
    largest_id = ""
    for start, end in positions:
        mono_id = f'{contig_id}_{start}_{end-1}'
        mono_seq = sequence[start-1:end-1]
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
    
    # parser.add_argument("--assembly_file", type=str, default=None, 
    #                     help="Path to complete assembly fasta file")
    
    parser.add_argument("--fastq_file", type=str, default=None, 
                        help="Path to long read fastq.gz")
    
    parser.add_argument("--out_path", type=str, default=None, 
                        help="Output path (default: plasmid_file_monomer")

    # parser.add_argument("--circ_db", type=str, required=True, 
    #                     help="Circlator DB path")

    parser.add_argument("--threads", type=int, default=6, 
                        help="Number of threads (default: 6)")
    
    parser.add_argument("--max_len", type=int, default=200000, 
                        help="Maximum plasmid length (default: 200000)")
    
    # parser.add_argument("--plasmid_id", type=str, default=None, 
    #                     help="Plasmid id (contig name) in complete assembly fasta (default: same id as in plasmid fasta file)")
    
    args = parser.parse_args()
    
    if args.out_path is None:
        args.out_path = os.path.splitext(args.fasta_file)[0] + '_monomer'

    return args

# #%%
# def write_report(report_path, multimer, monomer_len, max_len, monomers):
#     d_output = {'is_multimer': multimer,
#                 'original_len': max_len,
#                 'monomer_len': monomer_len,
#                 'multiplicity': len(monomers.keys())}
    
#     with open(report_path, mode='w', newline='') as of:
#         csv_writer = csv.DictWriter(of, fieldnames=d_output.keys())
#         csv_writer.writeheader()
#         csv_writer.writerow(d_output)
#     return None

# # %%
def main():
    args = parse_arguments()

    ## 0. inputs
    # max plasmid length
    max_len=args.max_len
    # plasmid fasta
    fasta_file = args.fasta_file
    fasta_basename = os.path.splitext(os.path.basename(fasta_file))[0]
    # procesar fasta y determinar contigs mayores que el query (1.2kb) y menores que longitud máxima (200kb default)
    d_fasta = read_fasta(fasta_file, max_len)
    # fastq
    fastq_file = args.fastq_file
    # output folder
    out_path = args.out_path
    if not os.path.exists(out_path): os.makedirs(out_path, exist_ok=True)
    report_path = os.path.join(out_path, f'{fasta_basename}_report.txt')
    # threads
    threads = args.threads

    ## 1. identificación de las repeticiones
    print()
    print('Step 1: Repetition identification')
    print()

    # 1.1 Procesar contigs uno a uno
    d_contig_corr = {}
    multimers = []
    for contig in d_fasta.keys():
        print(contig)
        sequence = str(d_fasta[contig][0])
        query = d_fasta[contig][1]
        if d_fasta[contig][1] == '':
            print('Skipping... Too large or too small.')
            d_contig_corr[f'>{contig}'] = sequence
        else:
            print('Checking...')
            # 1.1 creación de blastdb a partir de fasta circularizado
            blast_path, contig_path = create_blastdb(out_path, contig, sequence)

            # 1.2 guardar query en un archivo
            query_fasta = os.path.join(blast_path, f'query_{contig}.fasta')
            write_fasta(query_fasta, {f'>query_sequence': query})

            # 1.3 localizar repeticiones con blastn
            blast_out = repetition_search(blast_path, contig_path, query_fasta, fasta_basename, contig)

            ## 2. creación del monómero
            # 2.1 Extracción de monómeros
            df_reps = pd.read_table(blast_out, header=None)

            monomers, largest_id = extract_monomers(df_reps, sequence, contig)
            monomer_len = len(monomers[largest_id])
            # max_len = len(d_fasta[list(d_fasta.keys())[0]])
            if len(monomers.keys()) == 1:
                print("solo un monómero")
                print('Skipping... Not a multimer')
                d_contig_corr[f'>{contig}'] = sequence
        #         write_report(report_path, 'No', monomer_len, max_len, monomers)
        #         print(f'Find statistical report on: {report_path}')
        #         sys.exit(0)
            else:
            # 2.2 Comprobación de similitud
                """
                    pairwise de todos los monómeros frente al de mayor tamaño.
                    si la identidad es menor de 90%, no son monómeros
                    si es mayor y el qcov es menor del 90%, no usamos en alineamiento múltiple
                    si es mayor, añadimos al alineamiento múltiple
                """
                print("similarity_check")
                status, l_complete, l_partial = similarity_check(monomers, largest_id)
                if status == False:
                    d_contig_corr[f'>{contig}'] = sequence
                    print('Skipping... Not a multimer')
                    # write_report(report_path, 'No', max_len, max_len, monomers)
                    # print(f'Find statistical report on: {report_path}')
                    # sys.exit(0)
                else:
                    # 2.3 Alineamiento múltiple de las secuencias que pasan los criterios "completa"
                    multimers.append(contig)
                    print('Multimer!')
                    print()
                    print('Step 2: Monomer creation')

                    monomer_path = os.path.join(out_path, 'monomers')
                    if not os.path.exists(monomer_path): os.mkdir(monomer_path)
                    monomer_file = os.path.join(monomer_path, 'complete_monomers.fasta')
                    with open(monomer_file, 'w') as of: 
                        for i in l_complete:
                            of.write(f'>{i}\n')
                            of.write(f'{monomers[i]}\n')

                    aln_file = os.path.join(monomer_path, 'complete_monomers.aln')
                    aln_err = os.path.join(monomer_path, 'mafft.log')
                    mafft_cmd = ['mafft', '--adjustdirectionaccurately', '--thread', str(threads), monomer_file]
                    with open(aln_file, 'w') as stdf, open(aln_err, 'w') as errf:
                        subprocess.run(mafft_cmd, stdout=stdf, stderr=errf)

                    # 2.4 Generación del consenso
                    cons_file = os.path.join(out_path, f'{contig}_consensus.fasta')
                    cons_cmd = ['em_cons', '-sequence', aln_file, '-outseq', cons_file, '-name', contig]
                    subprocess.run(cons_cmd)
                    print('Multimer resolved')
                    d_cons = read_fasta(cons_file, query_len=0)
                    d_contig_corr[f'>{contig}'] = str(d_cons[contig][0])

    output_file = os.path.join(out_path, f'{fasta_basename}_corr.fasta')
    write_fasta(output_file, d_contig_corr)
    print(f'Written in {output_file}')

#     # 2.5 Devolver report
#     write_report(report_path, 'Yes', monomer_len, max_len, monomers)

#     print(f'Find monomer consensus at: {cons_file}')
#     print(f'Find divided multimer at: {monomer_file}')
#     print(f'Find statistical report on: {report_path}')
#     print()

    ## 3. validación del multímero biológico
    if fastq_file != None and len(multimers) > 0:
        print('Step 3: Biological multimer validation')
        print()

        # Mapeo
        print('Mapping reads to assembly')
        map_path = os.path.join(out_path, 'multimer_mapping')
        if not os.path.exists(map_path): os.mkdir(map_path)

        # mapeo original
        sam_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.sam')
        sam_err = os.path.join(map_path, f'minimap.log')
        minimap_cmd = ['minimap2', '--secondary=no', '-t', str(threads), '-ax', 'lr:hq', '-o', sam_file, fasta_file, fastq_file]
        with open(sam_err, 'w') as errf:
            subprocess.run(minimap_cmd, stderr=errf)

        # paso a bam eliminando lecturas no mapeadas
        bam_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.bam')
        smt_err = os.path.join(map_path, f'samtools.log')
        samtools_filter_cmd = ['samtools', 'view', '-h', '-@', str(threads), '-F', '4', '-bS', sam_file]
        with open(bam_file, 'w') as stdf, open(smt_err, 'w') as errf:
            subprocess.run(samtools_filter_cmd, stdout=stdf, stderr=errf)

        # ordenar
        bam_sorted_file = os.path.join(map_path, f'{fasta_basename}_whole_genome.sorted.bam')
        samtools_sort_cmd = ['samtools', 'sort', '--threads', str(threads), bam_file, '-o', bam_sorted_file]
        with open(smt_err, 'a') as errf:
            subprocess.run(samtools_sort_cmd, stderr=errf)

        # generar index
        samtools_index_cmd = ['samtools', 'index', '-@', str(threads), bam_sorted_file]
        with open(smt_err, 'a') as errf:
            subprocess.run(samtools_index_cmd, stderr=errf)

        # extraer solo lecturas que mapean con el plásmido en cuestión
        for multimer in multimers:
            print(f'Filtering reads mapping plasmid: {multimer}')
            print()
            sam_concat_file = os.path.join(map_path, f'{fasta_basename}_{multimer}.sam')
            samtools_extract_cmd = ['samtools', 'view', '-h', bam_sorted_file, multimer]
            with open(sam_concat_file, 'w') as stdf, open(smt_err, 'a') as errf:
                subprocess.run(samtools_extract_cmd, stdout=stdf, stderr=errf)

            # pasar sam a bam
            bam_concat_file = os.path.join(map_path, f'{fasta_basename}_{multimer}.bam')
            samtools_bam_cmd = ['samtools', 'view', '-h', '-@', str(threads), '-F', '4', '-bS', sam_concat_file]
            with open(bam_concat_file, 'w') as stdf, open(smt_err, 'a') as errf:
                subprocess.run(samtools_bam_cmd, stdout=stdf, stderr=errf)

            # generar index
            samtools_index_cmd = ['samtools', 'index', '-@', str(threads), bam_concat_file]
            with open(smt_err, 'a') as errf:
                subprocess.run(samtools_index_cmd, stderr=errf)

            fastq_sorted_file = os.path.join(out_path, f'{fasta_basename}_{multimer}.sorted.fastq.gz')
            samtools_fastq_cmd = ['samtools', 'fastq',  bam_concat_file]
            gzip_cmd = ['gzip']
            sp = subprocess.run(samtools_fastq_cmd, check=True, capture_output=True)
            with open(fastq_sorted_file, 'w') as stdf, open(smt_err, 'a') as errf:
                subprocess.run(gzip_cmd, input=sp.stdout, stdout=stdf, stderr=errf)

            # Plot
            print('Plotting reads length')
            read_lengths, filtered_reads = extract_read_lengths_and_filter_reads(fastq_sorted_file, monomer_len)
            output_file_plot = os.path.join(out_path, f'{fasta_basename}_{multimer}_read_plot.png')
            multimer_len = len(d_fasta[multimer][0])
            monomer_len = len(d_contig_corr[f'>{multimer}'])
            plot_read_lengths(read_lengths, output_file_plot, monomer_len, multimer_len)

            print(f"Plot saved as {output_file_plot}")
            print()
        
    else:
        print('Skipping multimer validation.')
        print('To validate, introduce fastq.gz and complete assembly fasta.')
        print()
        
    print('End of pipeline.')
    sys.exit(0)

# %%
if __name__ == "__main__":
    main()


