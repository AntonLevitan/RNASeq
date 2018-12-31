from Bio import SeqIO
try:
    from urllib.request import urlretrieve, urlopen
except ImportError:
    from urllib2 import urlopen
    from urllib import urlretrieve
from io import BytesIO
import gzip
import os
import subprocess
import pandas as pd
import pysam
import argparse
from normalization import rpkm, tpm


def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(species_flag, dest='sp', required=True)
    parser.add_argument(fastq_flag, dest='fastq', required=True)
    parser.add_argument(matched_norm_flag, dest='mn')

    return parser.parse_args()


fastq_flag = '--fastq'
species_flag = '--sp'
DATA_DIRECTORY = 'Data' + os.sep
matched_norm_flag = '--mn'

file_name_prefix = str(arguments().fastq)[:-6]

fastq_file = file_name_prefix + '.fastq'
bam_file = file_name_prefix + '.bam'
count_file = file_name_prefix + '.txt'
csv_file = file_name_prefix + '.csv'


if arguments().sp == 'ca':
    CHR_URL = 'http://candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_chromosomes.fasta.gz'
    GFF_URL = 'http://candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff'
    DIPLOID_FILENAME = 'A22_diploid_current_chromosomes.fasta'
    CHR_FILENAME = 'haplotype_A_A22_current_chromosomes.fasta'
    CURRENT_FEATURES_GFF = 'A22_current_features.gff'
elif arguments().sp == 'cg':
    CHR_URL = 'http://www.candidagenome.org/download/sequence/C_glabrata_CBS138/current/C_glabrata_CBS138_current_chromosomes.fasta.gz'
    GFF_URL = 'http://www.candidagenome.org/download/gff/C_glabrata_CBS138/C_glabrata_CBS138_current_features.gff'
    CHR_FILENAME = 'C_glabrata_CBS138_current_chromosomes.fasta'
    CURRENT_FEATURES_GFF = 'C_glabrata_CBS138_current_features.gff'
else:
    raise ValueError('Unknown species flag specified')


if not os.path.exists(DATA_DIRECTORY):
        os.makedirs(DATA_DIRECTORY)

if not os.path.exists(DATA_DIRECTORY + CURRENT_FEATURES_GFF):
    urlretrieve(GFF_URL, DATA_DIRECTORY + CURRENT_FEATURES_GFF)
    print('downloaded {} file'.format(CURRENT_FEATURES_GFF))


def download_gzip(url, filename):

    download = urlopen(url)
    zipped = BytesIO(download.read())
    unzipped = gzip.GzipFile(fileobj=zipped)

    with open(DATA_DIRECTORY + filename, 'wb') as file:
        file.write(unzipped.read())

    unzipped.close()
    zipped.close()
    download.close()


if arguments().sp == 'ca':
    if not os.path.exists(DATA_DIRECTORY + DIPLOID_FILENAME):
        download_gzip(CHR_URL, DIPLOID_FILENAME)
        print('downloaded {} file'.format(DIPLOID_FILENAME))

    if not os.path.exists(DATA_DIRECTORY + CHR_FILENAME):
        chro = list(SeqIO.parse(DATA_DIRECTORY + DIPLOID_FILENAME, "fasta"))
        haplotype_A = []

        for i in range(len(chro)):
            if 'A' in chro[i].id:
                haplotype_A.append(chro[i])

        SeqIO.write(haplotype_A, DATA_DIRECTORY + CHR_FILENAME, 'fasta')
        print('created {} file'.format(CHR_FILENAME))
elif arguments().sp == 'cg':
    if not os.path.exists(DATA_DIRECTORY + CHR_FILENAME):
        download_gzip(CHR_URL, CHR_FILENAME)
        print('downloaded {} file'.format(CHR_FILENAME))
else:
    raise ValueError('Unknown species flag specified')

if not os.path.exists(DATA_DIRECTORY + bam_file):
    subprocess.call(['STAR', '--runThreadN 8 ', '--runMode genomeGenerate', '--genomeDir ' + DATA_DIRECTORY,
                     '--genomeFastaFiles ' + DATA_DIRECTORY + CHR_FILENAME, '--sjdbGTFtagExonParentTranscript ID'])
    subprocess.call(['STAR', '--runThreadN 8 ', '--outSAMstrandField intronMotif', '--genomeDir ' + DATA_DIRECTORY,
                     '--readFilesIn ' + fastq_file, '--outSAMtype BAM SortedByCoordinate',
                     '--outFileNamePrefix ' + DATA_DIRECTORY + file_name_prefix, '--alignIntronMin 30', '--alignIntronMax 1000'])
    subprocess.call(['mv', DATA_DIRECTORY + file_name_prefix + 'Aligned.sortedByCoord.out.bam', DATA_DIRECTORY + bam_file])
    pysam.index(DATA_DIRECTORY + bam_file)

# TODO: change this to a direct python implementation using HTSeq library
# https://media.readthedocs.org/pdf/htseq/release_0.10.0/htseq.pdf
if not os.path.exists(count_file):
    count_logfile = open(count_file, 'w')
    proc = subprocess.Popen(['htseq-count', '--stranded=no', '-t', 'gene', '-i', 'ID', '-f', 'bam',
                             DATA_DIRECTORY + bam_file, DATA_DIRECTORY + CURRENT_FEATURES_GFF], stdout=count_logfile)
    proc.wait()
    # proc.kill()


def mapping_reads(data, output_data):

    data = pd.read_csv(data, delimiter="\t", header=None)
    reads = data.iloc[:-5, :]
    reads.columns = ['id', 'read_count']
    reads['read_count'] = reads['read_count'] + 1

    if arguments().sp == 'ca':
        reads['chromosome'] = reads['id'].str.split('_', 1).str[0]
        reads = reads[reads['id'].str.split('_', 2).str[2] == 'A']
        location = pd.read_csv(DATA_DIRECTORY + CURRENT_FEATURES_GFF, sep='\t', skiprows=25)
        location = location[location.iloc[:, 2] == 'gene']
        location['name'] = location.iloc[:, 8].str.split(';', 1).str[0].str.split('=', 1).str[1]
        location = location[location['name'].str.split('_',  2).str[2] == 'A']
        reads['start'] = list(location.iloc[:, 3])
        reads['end'] = list(location.iloc[:, 4])
    elif arguments().sp == 'cg':
        location = pd.read_csv(DATA_DIRECTORY + 'cg_features_coords.csv')
        reads = reads.merge(location)
        reads = reads. drop('Unnamed: 0', axis=1).rename(columns={'chrom': 'chromosome'})
        reads['chromosome'] = reads['chromosome'].str.split('_', 1).str[0]
    reads['length'] = list(abs(reads['start'] - reads['end']))
    reads['rpkm'] = list(rpkm(reads['read_count'], reads['length']))
    reads['tpm'] = list(tpm(reads['read_count'], reads['length']))
    print('')
    print('mapped to features: ' + str(sum(reads['read_count'])))
    print('total reads: ' + str(sum(data.iloc[:, 1])))
    print('ratio: ' + str(sum(reads['read_count']) / sum(data.iloc[:, 1])))
    reads.to_csv(output_data)
    return reads


if not os.path.exists(csv_file):
    mapping_reads(count_file, csv_file)


def matched_normalization():

    reference_file_prefix = arguments().mn
    reference_txt = reference_file_prefix + '.txt'
    reference_csv = reference_file_prefix + '.csv'
    excel_file = file_name_prefix + '_' + reference_file_prefix + '.xlsx'
    
    reference = mapping_reads(reference_txt, reference_csv)
    examined = mapping_reads(count_file, csv_file)
    
    examined['ref_reads'] = reference['read_count']
    examined['ref_rpkm'] = reference['rpkm']
    examined['ref_tpm'] = reference['tpm']
    examined = examined.replace(0, 1)
    examined['reads_ratio'] = examined['read_count'] / examined['ref_reads']
    examined['rpkm_ratio'] = examined['rpkm'] / examined['ref_rpkm']
    examined['tpm_ratio'] = examined['tpm'] / examined['ref_tpm']
    examined = examined[examined['reads_ratio'] < examined['reads_ratio'].quantile(.999)]
    avg = examined.groupby(['chromosome']).mean()
    print(avg)
    examined.to_excel(excel_file)
    return examined


if arguments().mn:
    matched_normalization()
