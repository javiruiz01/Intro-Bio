import os.path

# Trouver pour l'organisme chosi (Vibrio cholerae):
#   Nombre de chromosomes
#   Nombre de plasmides
#   Longueurs
#   Longueur totale du genome
#   Pourcentage en GC globale
#   Composition de nucleotides

hello = "holita"


def newOrder(fichier):
    f = open(fichier, 'r')
    rawData = ''
    for line in f:
        if line.startswith(">"):
            continue
        rawData += line
    data = rawData.replace('\n', '')
    return data


def nChromosomes(genome):
    total = 0
    for line in genome:
        if line.startswith(">"):
            total += 1
    return total


def nNucleotidesGenome(sequence):
    return len(sequence)


# def nNucleotidesChromosome(genome):
#     lengths = []
#     total = 0
#     for line in genome:
#         if line.startswith(">"):
#             lengths.append(total)
#             total = 1
#         else:
#             total += len(line) - 1
#     return lengths


def pourcentage(sequence):
    length = len(sequence)
    D = {
        'A': (float(sequence.count('A')) / length),
        'C': (float(sequence.count('C')) / length),
        'G': (float(sequence.count('G')) / length),
        'T': (float(sequence.count('T')) / length)
    }
    return D


def orfs(seq):
    cds = []
    pos = []
    for i in range(0, 3):
        seq_aux = seq[i:len(seq)]
        result, positions = orf_aux(seq_aux)
        cds.append(result)
        pos.append(positions)
    # Maintenant, il faut les garder dans un fichier
    newFile = printToFile(cds, pos)
    return newFile


def printToFile(cds, pos):
    newFile = open('orfFile.txt', 'a+')
    counter = 1
    for i, j in zip(cds, pos):
        for codon, position in zip(i, j):
            # newFile.write("> ORF numero " + str(counter) + "\n")
            # newFile.write("> Position = " + str(position) + "\n")
            newFile.write(str(position) + "\t" + str(len(codon) + position) + "\t" + codon + "\n")
            counter += 1
            # newFile.write(codon + "\n")
    return newFile.name


def orf_aux(seq):  # 6.4 AUX
    orf = ''
    result = []
    start = 'ATG'
    stop = ['TAA', 'TAG', 'TGA']
    begin = False
    positions = []
    position = 0
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) % 3 != 0:
            break
        if codon == start and begin is False:
            position = i
            begin = True
        if not begin:
            continue
        orf += codon
        if codon in stop:
            begin = False
            result.append(orf)
            positions.append(position)
            orf = ''
    return result, positions


def genes(fichier1, fichier2):
    f1 = open(fichier1, 'r')
    genome = newOrder(fichier2)
    length = len(genome)
    next(f1)
    nouveauFichier = 'genes' + os.path.basename(fichier1)[0:-4] + '.txt'
    fasta = open(nouveauFichier, 'a+')
    for line in f1:
        replicon = line.split("\t")[2]
        accesion = line.split("\t")[3]
        orientation = line.split("\t")[4]
        identifiant = line.split("\t")[5]
        data = genome[int(replicon) - 1:int(accesion) - 3]
        if orientation == '-':
            # data = data[::-1]
            continue
        # Faut chercher ATG?
        # data = orf_aux(data)
        data = ''.join(data)
        # Maintenant on fait le fichier avec format FASTA
        if len(data) != 0:
            # fasta.write(">" + identifiant + "\n")
            fasta.write(
                str(int(replicon) - 1) + "\t" + str(int(accesion) - 3) + "\t" + data + "\n"
            )
    fasta.close()
    return nouveauFichier, length


def genesToList(fichierGenes, length):
    f1 = open(fichierGenes, 'r')
    result = [0] * length
    for line in f1:
        gene = line.split("\t")[2]
        if not gene.startswith("ATG"):
            continue
        replicon = line.split("\t")[0]
        accesion = line.split("\t")[1]
        for j in range(int(replicon), int(accesion)):
            result[j] = 1
    return result


def evaluation(fichierGenes, fichierSeuil, length):
    genes = genesToList(fichierGenes, length)
    seuil = seuilToList(fichierSeuil, length)
    return compare_intervalle(genes, seuil)


def compare_intervalle(genome, orfs):
    vp = vn = fn = fp = 0
    for i in range(0, len(genome)):
        gen = genome[i]
        orf = orfs[i]
        if gen + orf == 2:
            vp += 1
        elif gen + orf == 0:
            vn += 1
        elif gen > orf:
            fn += 1
        else:
            fp += 1
    result = [[vn, fp], [fn, vp]]
    return result


def seuilToList(fichierSeuil, length):
    f1 = open(fichierSeuil, 'r')
    result = [0] * length
    for line in f1:
        if line.startswith("START"):
            continue
        start = line.split("\t")[0]
        end = line.split("\t")[1]
        for j in range(int(start), int(end)):
            result[j] = 1
    return result


def senSpe(matrice):
    vraip = matrice[1][1]
    sen = vraip / float(vraip + matrice[1][0])
    vn = matrice[0][0]
    spe = vn / float(vn + matrice[0][1])
    vp = vraip / float(vraip + matrice[0][1])
    return sen, spe, vp


def main():
    cholerae = open("./data/Vibrio_cholerae_genome.fasta", 'r')
    coli = open("./data/Escherichia.coli.fasta", 'r')

    # Nombre de chromosomes:
    # cholerae_txt = newOrder(cholerae)
    # coli_txt = newOrder(coli)
    # total1 = nChromosomes(cholerae)
    # total2 = nChromosomes(coli)
    # print "Nombre de chromosomes: \n\tCholerae = " + str(total1) + "\n\tColi = " + str(total2)

    # Nombre de nucleotides dans le genome(longueur):
    # print nNucleotides(cholerae_txt)
    # print nNucleotides(coli_txt)

    # Nombre de nucleotides dans les chromosomes:
    # lengths = nNucleotidesChromosome(cholerae)
    # for i in lengths:
    #     print "Hello" + str(i)

    # Pourcentage de nucleotides:
    # cholerae_txt = newOrder(cholerae)
    # coli_txt = newOrder(coli)
    # D = pourcentage(cholerae_txt)
    # for i in D:
    #     print str(i) + " = " + str(D[i])

    # trouver les orfs
    data = newOrder("./data/Vibrio_cholerae_genome.fasta")
    if not os.path.exists("orfFile.txt"):
        orfs(data)  # Ours

    if not os.path.exists("genesVibrio_cholerae.txt"):
        genes("./data/Vibrio_cholerae.tab", "./data/Vibrio_cholerae_genome.fasta")

    result = evaluation("genesVibrio_cholerae.txt", "orfFile.txt", len(data))
    print result
    sen, spe, vp = senSpe(result)
    print ("SEN = " + str(sen) + "\nSPE = " + str(spe) + "\nVP = " + str(vp))


main()
