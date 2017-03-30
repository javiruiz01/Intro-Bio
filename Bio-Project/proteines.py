from codeGenetique import codeGenetique


def traduction(gene):
    result = ''
    i = 0
    while (i < len(gene)):
        codon = gene[i:i + 3]
        if codon.endswith("\n"):
            break
        if (len(codon) % 3 != 0):
            break
        result += codeGenetique[codon]
        i += 3
    return result


def newOrder(fichier):
    f = open(fichier, 'r')
    rawData = ''
    for line in f:
        if line.startswith(">"):
            continue
        rawData += line
    data = rawData.replace('\n', '')
    f.close()
    return data


def main():
    print "ANNOTATION PAR HOMOLOGIE AVEC BLAST --> Coding genes"
    print "\tFirst genome of Vibrio Cholerae"
    gene_1_codants = open("./fasta_files_1_genome/genes_codants_Vibrio_cholerae.fasta", "r")
    proteines_1 = open("./blast_files/proteines_codants_Vibrio_cholerae_1.fasta", "a+")
    for line in gene_1_codants:
        if line.startswith(">"):
            proteines_1.write(line)
            continue
        proteines_1.write(traduction(line) + "\n")
    gene_1_codants.close()
    proteines_1.close()
    print "\t\tNew file created with coding genes in proteines for first genome"

    print "\tSecond genome of Vibrio Cholerae"
    gene_2_codants = open("./fasta_files_2_genome/genes_codants_Vibrio_cholerae_2.fasta")
    proteines_2 = open("./blast_files/proteines_codants_Vibrio_cholerae_2.fasta", "a+")
    for line in gene_2_codants:
        if line.startswith(">"):
            proteines_2.write(line)
            continue
        proteines_2.write(traduction(line) + "\n")
    gene_2_codants.close()
    proteines_2.close()
    print "\t\tNew file created with coding genes in proteines for second genome"

    print "\tCreating new file with proteines from E. coli"
    gene_codants = open("./fasta_files_1_genome/genes_codants_Escherichia.coli.fasta", "r")
    proteines_coli = open("./blast_files/proteines_codants_Coli.fasta", "a+")
    for line in gene_codants:
        if line.startswith(">"):
            proteines_coli.write(line)
            continue
        proteines_coli.write(traduction(line) + "\n")
    gene_codants.close()
    proteines_coli.close()
    print "\t\tNew file created with coding genes in proteines for Coli"

    print "ANNOTATION PAR HOMOLOGIE AVEC BLAST --> Non coding genes"
    print "\tFirst genome of Vibrio Cholerae"
    gene_1_non_codants = open("./fasta_files_1_genome/genes_non_codants_Vibrio_cholerae.fasta", "r")
    proteines_1_non_codants = open("./blast_files/proteines_non_codants_cholerae_1.fasta", "a+")
    for line in gene_1_non_codants:
        if line.startswith(">"):
            proteines_1_non_codants.write(line)
            continue
        proteines_1_non_codants.write(traduction(line) + "\n")
    gene_1_non_codants.close()
    proteines_1_non_codants.close()
    print "\t\tNew file created with non coding genes in proteines for first genome"

    print "\tSecond genome of Vibrio Cholerae"
    gene_2_non_codants = open("./fasta_files_2_genome/genes_non_codants_Vibrio_cholerae_2.fasta", "r")
    proteines_2_non_codants = open("./blast_files/proteines_non_codants_cholerae_2.fasta", "a+")
    for line in gene_2_non_codants:
        if line.startswith(">"):
            proteines_2_non_codants.write(line)
            continue
        proteines_2_non_codants.write(traduction(line) + "\n")
    gene_2_non_codants.close()
    proteines_2_non_codants.close()
    print "\t\tNew file created with non coding genes in proteines for second genome"

    print "\tSecond genome of Vibrio Cholerae"
    coli_non_codants = open("./fasta_files_1_genome/genes_non_codants_Escherichia.coli.fasta", "r")
    proteines_coli_non_codants = open("./blast_files/proteines_non_codants_Coli.fasta", "a+")
    for line in coli_non_codants:
        if line.startswith(">"):
            proteines_coli_non_codants.write(line)
            continue
        proteines_coli_non_codants.write(traduction(line) + "\n")
    coli_non_codants.close()
    proteines_coli_non_codants.close()
    print "\t\tNew file created with non coding genes in proteines for E. coli"


main()
