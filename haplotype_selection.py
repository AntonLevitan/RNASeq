from Bio import SeqIO

records = list(SeqIO.parse("C_albicans_SC5314_A22_current_chromosomes.fasta", "fasta"))

haplotype_A = [records[0], records[2], records[4], records[6], records[8], records[10], records[12], records[15]] 

SeqIO.write(haplotype_A, 'haplotype_A_A22_current chromosomes.fasta', 'fasta')