from pandas import read_csv, DataFrame, ExcelWriter
import re
from collections import defaultdict
from Bio import SeqIO

fasta_file = snakemake.input['fasta']
aln_file = snakemake.input['aln']
data_file = snakemake.input['data']
clstr_file = snakemake.input['clstr']

output_file = str(snakemake.output)

def read_clstr(clstr_file):
    clusters = defaultdict(list)
    reps = {}
    with open(clstr_file) as file:
        for line in file:
            if line.startswith('>Cluster'):
                cluster = line
            else:
                match = re.match('(\\d+)\t(\\d+)aa, >(.+)[.][.][.] (.+)', line)
                num, seq_len, seq_name, identity = match.groups()
                if identity == '*':
                    reps[cluster] = seq_name
                clusters[cluster].append(seq_name)
    rep_data = []
    clust_data = []
    for cluster, labels in clusters.items():
        rep = reps[cluster]
        rep_data.append({ 'label': rep, 'cluster_members': ','.join(labels) })
        for label in labels:
            clust_data.append({ 'label': label, 'representative': rep })
    return DataFrame(clust_data).set_index('label'), DataFrame(rep_data).set_index('label')

def read_fasta(fasta_file):
    data = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        data.append({ 'label': record.id, 'sequence': record.seq })
    return DataFrame(data).set_index('label')

clstr, reps = read_clstr(clstr_file)
fasta = read_fasta(fasta_file)
aln = read_fasta(aln_file).rename(columns = { 'sequence': 'stripped_sequence' })

data = read_csv(data_file, index_col = 'label').join(fasta)
seq_data = data.join(clstr)
rep_data = aln.join(data)

writer = ExcelWriter(output_file, engine = 'xlsxwriter')
seq_data.to_excel(writer, sheet_name = 'Sequences')
rep_data.to_excel(writer, sheet_name = 'Representatives')
writer.close()
