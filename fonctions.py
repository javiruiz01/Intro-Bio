from codeGenetique import codeGenetique
import sys
import os.path
import subprocess


def pourcentage(seq):  # 6.1 --> DONE
    percent = 100.0 / len(seq)
    nA = float(seq.count("A")) * percent
    print ('A = ' + str(nA) + ' %')
    nC = float(seq.count("C")) * percent
    print ('C = ' + str(nC) + ' %')
    nG = float(seq.count("G")) * percent
    print ('G = ' + str(nG) + ' %')
    nT = float(seq.count("T")) * percent
    print ('T = ' + str(nT) + ' %\n')


def nombreCG(seq):  # 6.2 --> DONE
    nombre = float(seq.count("CG"))
    nCG = (nombre / len(seq))
    return [nombre, round(nCG, 1)]


def listeCodons(seq):  # 6.3 --> DONE
    result = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if (len(codon) % 3 != 0):
            break
        result.append(seq[i:i + 3])
    return result


def orf_aux(seq):  # 6.4 AUX
    result = []
    start = 'ATG'
    stop = ['TAA', 'TAG', 'TGA']
    begin = False
    # seq = seq[offset: len(seq)]
    # print ("From AUX, the sequence is: " + seq)
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if (len(codon) % 3 != 0):
            break
        if codon == start and begin is False:
            begin = True
        if not begin:
            continue
        result.append(codon)
        if codon in stop:
            break
    return result


def orfs(seq):  # 6.4 --> DONE
    cds = []
    seqs = [seq, complementaire(seq)]
    for s in seqs:
        # print("This is the sequence: " + s)
        for i in range(0, 3):
            # print ("i = " + str(i))
            seq_aux = seq[i:len(seq)]
            result = orf_aux(seq_aux)
            # print ("Result = ", result[0:-1])
            cds.append(result[:])
    return cds


def complementaire(seq):  # 6.5 --> DONE
    dic = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    result = ''
    for i in seq:
        result += dic[i]
    return result


def traduction(gene):  # 6.6 --> DONE
    result = []
    i = 0
    while (i < len(gene)):
        codon = gene[i:i + 3]
        if (len(codon) % 3 != 0):
            break
        result.append(codeGenetique[codon])
        i += 3
    return result


def pourcentageFichier(fichier):  # 7.1 --> DONE
    data = newOrder(fichier)
    pourcentage(data)


def newOrder(fichier):  # 7.1 AUX
    f = open(fichier, 'r')
    rawData = ''
    for line in f:
        if line.startswith(">"):
            continue
        rawData += line
    data = rawData.replace('\n', '')
    return data


def genes(fichier1, fichier2):  # 7.2 --> DONE
    f1 = open(fichier1, 'r')
    genome = newOrder(fichier2)
    next(f1)
    nouveauFichier = 'fasta_' + os.path.basename(fichier1)[0:-4] + '.txt'
    fasta = open(nouveauFichier, 'a+')
    for line in f1:
        replicon = line.split("\t")[2]
        accesion = line.split("\t")[3]
        orientation = line.split("\t")[4]
        identifiant = line.split("\t")[5]
        data = genome[int(replicon) - 1:int(accesion) - 3]
        if orientation == '-':
            data = data[::-1]
        # Faut chercher ATG?
        data = orf_aux(data)
        data = ''.join(data)
        # Maintenant on fait le fichier avec format ASTA
        if len(data) != 0:
            fasta.write(">" + identifiant + "\n")
            fasta.write(data + "\n")
    fasta.close()
    return nouveauFichier

def pourcentageFASTA(fichier):  # 7.3
    f = open(fichier, 'r')
    id = ''
    for line in f:
        if line.startswith(">"):
            id += line
            continue
        print(id)
        pourcentage(line)
        id = ''


def writeToFile(fichier, trad, id):  # Aux
    nouveauFichier = fichier.replace('fasta', 'fastaProteines')
    f = open(nouveauFichier, 'a+')
    f.write(id + "\n")
    for item in trad:
        f.write(item)
    f.write("\n\n")
    return nouveauFichier


def traductionGenes(fichier):  # 7.4
    f = open(fichier, 'r')
    id = ''
    for line in f:
        if line.startswith(">"):
            id = line
            continue
        trad = traduction(line[0:len(line) - 1])
        nouveauFichier = writeToFile(fichier, trad, id)
        id = ''
    return nouveauFichier


def tailleMoyenne(fichier):  # 7.5
    f = open(fichier, 'r')
    proteinCounter = 0
    acc = 0
    for line in f:
        if line.startswith(">"):
            proteinCounter += 1
            continue
        # line = line[0:len(line) - 1]
        acc += len(line[0:len(line) - 1])
    taille = acc / proteinCounter
    print("Taille moyenne des proteines: " + str(taille))


def pourcentageCodante(fichier, total): # 7.6
    # PRE: fichier des genes + nombre total des nucleotides du genome
    data = newOrder(fichier)
    percentage = len(data) / total
    print(str(percentage))


def main():
    # pourcentage("ACCTGGACAT") # 6.1
    # for item in nombreCG("AACGTGGCACG"): # 6.2
    #     print item
    # for item in listeCodons("AACGTGGCA"): # 6.3
    #     print item
    # Seq= "CCTCAGCGAGGACAGCAAGGGACTAGCCAGGAGGGAGAACAGAAACTCCAGAACATCTTGGAAATAGCTCCCAGAAAAGCAAGCAGCCAACCAGGCAGGTTCTGTCCCTTTCACTCACTGGCCCAAGGCGCCACATCTCCCTCCAGAAAAGACACCATGAGCACAGAAAGCATGATCCGCGACGTGGAACTGGCAGAAGAGGCACTCCCCCAAAAGATGGGGGGCTTCCAGAACTCCAGGCGGTGCCTATGTCTCAGCCTCTTCTCATTCCTGCTTGTGGCAGGGGCCACCACGCTCTTCTGTCTACTGAACTTCGGGGTGATCGGTCCCCAAAGGGATGAGAAGTTCCCAAATGGCCTCCCTCTCATCAGTTCTATGGCCCAGACCCTCACACTCAGATCATCTTCTCAAAATTCGAGTGACAAGCCTGTAGCCCACGTCGTAGCAAACCACCAAGTGGAGGAGCAGCTGGAGTGGCTGAGCCAGCGCGCCAACGCCCTCCTGGCCAACGGCATGGATCTCAAAGACAACCAACTAGTGGTGCCAGCCGATGGGTTGTACCTTGTCTACTCCCAGGTTCTCTTCAAGGGACAAGGCTGCCCCGACTACGTGCTCCTCACCCACACCGTCAGCCGATTTGCTATCTCATACCAGGAGAAAGTCAACCTCCTCTCTGCCGTCAAGAGCCCCTGCCCCAAGGACACCCCTGAGGGGGCTGAGCTCAAACCCTGGTATGAGCCCATATACCTGGGAGGAGTCTTCCAGCTGGAGAAGGGGGACCAACTCAGCGCTGAGGTCAATCTGCCCAAGTACTTAGACTTTGCGGAGTCCGGGCAGGTCTACTTTGGAGTCATTGCTCTGTGAAGGGAATGGGTGTTCATCCATTTCTACCCAGCCCCCACTCTGACCCCTTTACTCTGACCCCTTTATTGTCTACTCCTCAGAGCCCCCAGTCTGTATCCTTCTAACTTAGAAAGGGGATTATGGCTCAGGGTCCAACTCTGTGCTCAGAGCTTTCAACAACTACTCAGAAACACAAGATGCTGGGACAGTGACCTGGACTGTGGGCCTCTCATGCACCACCATCAAGGACTCAAAGGGCTTTCCGAATTCACTGGAGCCTCGAATGTCCATTCCTGAGTTCTGCAAAGGGAGAGTGGTCAGGTTGCCTCTGTCTCAGAATGAGGCTGGATAAGATCTCAGGCCTTCCTACCTTCAGACCTTTCCAGATTCTTCCCTGAGGTGCAATGCACAGCCTTCCTCACAGAGCCAGCCCCCCTCTATTTATATTTGCACTTATTATTTATTATTTATTTATTATTTATTTATTTGCTTATGAATGTATTTATTTGGAAGGCCGGGGTGTCCTGGAGGACCCAGTGTGGGAAGCTGTCTTCAGACAGACATGTTTTCTGTGAAAACGGAGCTGAGCTGTCCCCACCTGGCCTCTCTACCTTGTTGCCTCCTCTTTTGCTTATGTTTAAAACAAAATATTTATCTAACCCAATTGTCTTAATAACGCTGATTTGGTGACCAGGCTGTCGCTACATCACTGAACCTCTGCTCCCCACGGGAGCCGTGACTGTAATCGCCCTACGGGTCATTGAGAGAAATAA"
    # for item in orfs(Seq): # 6.4
    #     print item
    # print complementaire('AACGTGGCA') # 6.5
    # for item in traduction('') # 6.6
    #     print item
    # pourcentageFichier(sys.argv[1])
    # genes('../data/Escherichia.coli/Escherichia.coli.tab', '../data/Escherichia.coli/Escherichia.coli.genome')
    # pourcentageFASTA('fasta.txt')
    # traductionGenes('fasta.txt')
    # En parametre passer les fichiers genome et tout
    # tailleMoyenne('fastaProteines.txt')
    # final('fasta.txt', 58021)

    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print ('Usage: python ' + sys.argv[0] + ' PATH_TO_TAB_FILE PATH_TO_GENOME')
    else:
        # Creation d'un fichier au format fasta avec les genes d'un genome
        print('Extraction des genes a partir d\'un genome')
        fastaGenes = genes(sys.argv[1], sys.argv[2])
        print('Nouveau fichier: ' + fastaGenes)

        # Creation d'un fichier avec les proteines des diferents fichiers
        print('Traduction des genes a proteines')
        fastaProteines = traductionGenes(fastaGenes)
        print('Nouveau fichier: ' + fastaProteines)
        tailleMoyenne(sys.argv[1])

        # Pourcentage de la region codante
        # print('Pourcentage de la region codante de une spece')
        # total = subprocess.Popen("grep -v \">\" " + sys.argv[2] + " | wc - c", shell=True)
        # pourcentageCodante(fastaGenes, total)

main()
