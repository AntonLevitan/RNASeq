from Bio import SeqIO
from urllib.request import urlretrieve, urlopen
from io import BytesIO
import gzip
import os
import subprocess

A22_CHR_URL = 'http://candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_chromosomes.fasta.gz'
A22_GFF_URL = 'http://candidagenome.org/download/gff/C_albicans_SC5314/Assembly22/C_albicans_SC5314_A22_current_features.gff'
DATA_DIRECTORY = 'genome' + os.sep
DIPLOID_FILENAME = 'A22_diploid_current_chromosomes.fasta'
HAPLOTYPE_A_FILENAME = 'haplotype_A_A22_current_chromosomes.fasta'
A22_CURRENT_FEATURES_GFF = 'A22_current_features.gff'

file_name_prefix = input('Input raw fastq filename here: ')

fastq_file = file_name_prefix + '.fastq'
bam_file = file_name_prefix + '.bam'
count_file = file_name_prefix + '.txt'

if not os.path.exists(A22_CURRENT_FEATURES_GFF):
    urlretrieve(A22_GFF_URL, A22_CURRENT_FEATURES_GFF)
    print('downloaded {} file'.format(A22_CURRENT_FEATURES_GFF))


def download_gzip(url, filename):

    if not os.path.exists(DATA_DIRECTORY):
        os.makedirs(DATA_DIRECTORY)

    download = urlopen(url)
    zipped = BytesIO(download.read())
    unzipped = gzip.GzipFile(fileobj=zipped)

    with open(DATA_DIRECTORY + filename, 'wb') as file:
        file.write(unzipped.read())

    unzipped.close()
    zipped.close()
    download.close()


if not os.path.exists(DATA_DIRECTORY + DIPLOID_FILENAME):
    download_gzip(A22_CHR_URL, DIPLOID_FILENAME)
    print('downloaded {} file'.format(DIPLOID_FILENAME))

if not os.path.exists(DATA_DIRECTORY + HAPLOTYPE_A_FILENAME):
    chro = list(SeqIO.parse(DATA_DIRECTORY + DIPLOID_FILENAME, "fasta"))
    haplotype_A = []

    for i in range(len(chro)):
        if 'A' in chro[i].id:
            haplotype_A.append(chro[i])

    SeqIO.write(haplotype_A, DATA_DIRECTORY + HAPLOTYPE_A_FILENAME, 'fasta')
    print('created {} file'.format(HAPLOTYPE_A_FILENAME))

if not os.path.exists(bam_file):
    subprocess.call(['STAR', '--runThreadN 12', '--runMode genomeGenerate', '--genomeDir ' + DATA_DIRECTORY, '--genomeFastaFiles ' + DATA_DIRECTORY + HAPLOTYPE_A_FILENAME, '--sjdbGTFtagExonParentTranscript ID'])
    subprocess.call(['STAR', '--runThreadN 12', '--outSAMstrandField intronMotif', '--genomeDir ' + DATA_DIRECTORY, '--readFilesIn ' + fastq_file, '--outSAMtype BAM SortedByCoordinate', '--alignIntronMin 30', '--alignIntronMax 1000'])
    subprocess.call(['mv', 'Aligned.sortedByCoord.out.bam', bam_file])
    subprocess.call(['samtools', 'index', bam_file])

if not os.path.exists(count_file):
    count_logfile = open(count_file, 'w')
    proc = subprocess.Popen(['htseq-count', '-t', 'gene', '-i', 'ID', '-f', 'bam', bam_file, A22_CURRENT_FEATURES_GFF], stdout=count_logfile)
    proc.wait()
    proc.kill()
