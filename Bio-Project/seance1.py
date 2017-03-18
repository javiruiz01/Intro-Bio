# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os.path
import shutil

# On travaille avec les genomes Vibrio Cholerae et E. Coli
coli = "./data/Escherichia.coli.fasta"
cholerae = "./data/Vibrio_cholerae_genome.fasta"
tabColi = "./data/Escherichia.coli.tab"
tabCholerae = "./data/Vibrio_cholerae.tab"

spaces = "\s\s\s\s"


###################
#   <Seance 1>    #
###################

# Change un fichier fasta en un String sans saut de ligne
def newOrder(fichier):  # 7.1 AUX
    f = open(fichier, 'r')
    rawData = ''
    for line in f:
        if line.startswith(">"):
            continue
        rawData += line
    data = rawData.replace('\n', '')
    f.close()
    return data


# Nombre de chromosomes
def nChromosomes(fastaFile):
    f = open(fastaFile, 'r')
    total = 0
    for line in f:
        if line.startswith(">"):
            total += 1
    f.close()
    return total


# Longueur de chromosomes
def lengthChromosome(fastaFile):
    lengths = []
    if nChromosomes(fastaFile) == 1:
        lengths.append(len(newOrder(fastaFile)))
        return lengths
    else:
        f = open(fastaFile, 'r')
        length = 0
        for line in f:
            if line.startswith(">"):
                if length == 0:
                    continue
                lengths.append(length)
                length = 0
                continue
            length += len(line) - 1
        lengths.append(length)
        f.close()
        return lengths


# Longueur totale du genome
def totalLength(fastaFile):
    return len(newOrder(fastaFile))


# Composition en nucleotides
def pourcentage(fastaFile):
    sequence = newOrder(fastaFile)
    length = len(sequence)
    D = {
        'A': (float(sequence.count('A')) / length),
        'C': (float(sequence.count('C')) / length),
        'G': (float(sequence.count('G')) / length),
        'T': (float(sequence.count('T')) / length)
    }
    return D


# Découpez le génome en blocs 1kbp non chevauchants
def decoupage(fastaFile, tabFile):
    fasta = newOrder(fastaFile)
    tab = open(tabFile, 'r')
    blocs = []
    bloc = fasta[0:1000]
    total = 1000  # 1kbp
    for line in tab:
        if line.startswith("#"):
            continue
        end = line.split()[3]
        if int(total) < int(end):
            bloc += fasta[int(total):int(end)]
            blocs.append(bloc)
            total = int(end) + 1000
            bloc = fasta[int(end):total]
    tab.close()
    return blocs


# Pourcentage de la chaine GC dans une sequence
def gcPercent(sequence):
    percent = (float(sequence.count('GC')) / len(sequence))
    return percent


# Pourcentage de G et de C dans une sequence et pas dans un fichier
def gcBlocPercentage(sequence):
    D = {
        'G': (float(sequence.count('G')) / len(sequence)),
        'C': (float(sequence.count('C')) / len(sequence))
    }
    return D


# Histogramme de la distribution de GCs
def histogrammeGCBlocs(blocs, name):
    percents = []
    for bloc in blocs:
        d = gcBlocPercentage(bloc)
        percents.append(d['C'] + d['G'])
    plt.hist(percents)
    plt.title(name)
    plt.show()


# ORFs du genome choisi
def orfs(sequence):
    cds = []
    pos = []
    seqs = [sequence, complementaire(sequence)[
                      ::-1]]  # TODO: Attention, je crois qu'on doit inverser lorsqu'on fait le complementaire
    for seq in seqs:
        for i in range(0, 3):
            seq_aux = seq[i:len(seq)]
            result, positions = orf_aux(seq_aux)
            cds.append(result)
            pos.append(positions)
    # Creation fichier FASTA
    return printToFile(cds)


# Fonction qui calcule le complementaire d'une sequence donnee
def complementaire(seq):
    dic = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    result = ''
    for i in seq:
        # TODO: Faire attention a ça, il ya des "Y" dans le gene
        if i == "Y" or i == "S" or i == "M":
            i = "C"
        if i == "W" or i == "K" or i == "N":
            i = "T"
        if i == "R":
            i = "A"
        result += dic[i]
    return result


# Fonction auxiliare pour calculer les ORFs
def orf_aux(sequence):
    orf = ''
    result = []
    start = 'ATG'
    stop = ['TAA', 'TAG', 'TGA']
    begin = False
    positions = []
    position = 0
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
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


# Fonction qui cree un fichier FASTA a partir des ORFs trouves dans une sequence
# On n'utilise pas le position pour le moment
def printToFile(cds):
    newFile = open('orfFasta.fsa', 'a+')
    counter = 1
    for i in cds:
        for codon in i:
            newFile.write('>Sequence_' + str(counter) + "[organism=Vibrio Cholerae]\n")
            counter += 1
            newFile.write(codon + "\n")
    newFile.close()
    return newFile.name


# Extraire sequences codantes des genes a partir du fichier tab et du genome
def genes(genomeFile, tabFile):
    tab = open(tabFile, 'r')
    genome = newOrder(genomeFile)
    next(tab)  # On saute la premiere ligne
    fastaFile = open("genes_codants_" + os.path.basename(tabFile)[0:-4] + ".fasta", 'a+')
    fastaNonCodants = open("genes_non_codants_" + os.path.basename(tabFile)[0:-4] + ".fasta", 'a+')
    counter = 0
    counter_non_codants = 0
    codonStart = False
    for line in tab:
        start = line.split("\t")[2]  # Start
        end = line.split("\t")[3]  # End
        orientation = line.split("\t")[4]
        # identifiant = line.split("\t")[5]
        if orientation == '-':
            data = genome[int(start) - 1: int(end)][::-1]
            data = complementaire(data)
        else:
            data = genome[int(start) - 1:int(end) - 3]
        data = ''.join(data)
        nStopTAA = data.count("TAA")
        nStopTAG = data.count("TAG")
        nStopTGA = data.count("TGA")
        if (not data.startswith("ATG")) or not ((int(nStopTAA) + int(nStopTAG) + int(nStopTGA)) != 0):
            counter_non_codants += 1
            fastaNonCodants.write('>Sequence_non_codante' + str(counter_non_codants) + "\n")
            fastaNonCodants.write(data + "\n")
            continue
        # if len(data) != 0:
        counter += 1
        fastaFile.write('>Sequence_' + str(counter) + "\n")
        fastaFile.write(data + "\n")
    fastaFile.close()
    tab.close()
    return fastaFile


# Creation d'un fichier avec les genes non codants, cad,
# les parties du gene qui ne sont pas dans le fichier tab
def nonCodants(genomeFile, tabFile):
    tab = open(tabFile, 'r')
    genome = newOrder(genomeFile)
    next(tab)
    fastaFile = open("genes_non_codants_" + os.path.basename(tabFile)[0:-4] + ".fasta", 'a+')
    counter = 0
    position = 0
    for line in tab:
        start = line.split("\t")[2]
        data = genome[int(position): int(start)]
        position = line.split("\t")[3]  # Position est la position de la fin du gene codant
        if len(data) != 0:
            counter += 1
            fastaFile.write('Sequence_non_codante_' + str(counter) + "\n")
            fastaFile.write(data + "\n")
    fastaFile.close()
    tab.close()
    return fastaFile


# Calculer pour chaque gene son pourcentage en GC et mettre a jour le fichier tab
def updateTabFile(genes_codants):
    f = open(genes_codants, 'rw+')
    newFile = open("tmp.fasta", 'a+')
    for line in f:
        if line.startswith(">"):
            newFile.write(line)
            continue
        d = gcBlocPercentage(line[0:len(line) - 1])
        percent = float(d['C'] + d['G']) * 100
        newFile.write(str(percent) + "\t" + line)
    f.close()
    newFile.close()
    os.rename(f.name, f.name + "_old")
    os.rename("tmp.fasta", f.name)
    return


# Histogramme des pourcentages en GC des gènes
def histogramGCgenes(fastaFile, name):
    f = open(fastaFile, 'r')
    percents = []
    for line in f:
        if line.startswith(">"):
            continue
        percents.append(float(line.split("\t")[0]))
    plt.hist(percents)
    plt.title(name)
    plt.show()


# Fonction pour ranger tous les noveaus fichier crees
def cleanUp():
    if os.path.exists("fasta_files"):
        shutil.rmtree('./fasta_files/')
    os.makedirs("./fasta_files")
    for file in os.listdir("./"):
        if file.endswith("fasta") or file.endswith("fsa") or file.endswith("fasta_old"):
            os.rename(file, "./fasta_files/" + file)


def seance1():
    # Nombre de chromosomes
    totalColi = nChromosomes(coli)
    totalCholerae = nChromosomes(cholerae)
    print "Total number of chromosomes:\n\tE.Coli = " + str(totalColi) + "\n\tVibrio Cholerae = " + str(totalCholerae)

    # Longueur de chromosomes
    i = 0
    lengthColi = lengthChromosome(coli)
    for length in lengthColi:
        i += 1
        print "Length of the " + str(i) + " chromosome from E. Coli = " + str(length)
    i = 0
    lengthCholerae = lengthChromosome(cholerae)
    for length in lengthCholerae:
        i += 1
        print "Length of the " + str(i) + " chromosome from Vibrio Cholerae = " + str(length)

    # Longueur totale du genome
    print "Total length of E. Coli = " + str(totalLength(coli))
    print "Total length of Vibirio Cholerae = " + str(totalLength(cholerae))

    # Composition en nucleotides
    dColi = pourcentage(coli)
    for item in dColi:
        print "Percentage of " + item + " in E. Coli = " + str(dColi[item])
    dCholerae = pourcentage(cholerae)
    for item in dColi:
        print "Percentage of " + item + " in Vibrio Cholerae = " + str(dCholerae[item])

    # Pourcentage en GC global
    # Cas ou GC veut dire la chaine GC
    # print "Pourcentage en GC global de E. Coli = " + str(gcPercent(newOrder(coli)))
    # print "Pourcentage en GC global de Vibrio Cholerae = " + str(gcPercent(newOrder(cholerae)))
    # Cas ou on veut la somme des G et des C
    print "Percentage of GC in E. Coli = " + str(dColi['G'] + dColi['C'])
    print "Percentage of GC in Vibrio Cholerae = " + str(dCholerae['G'] + dCholerae['C'])

    # Découpez le génome en blocs 1kbp non chevauchants
    blocsColi = decoupage(coli, tabColi)
    print "Number of blocs that we will be analysing from E. Coli = " + str(len(blocsColi))
    blocsCholerae = decoupage(cholerae, tabCholerae)
    print "Number of blocs that we will be analysing from Vibrio Cholerae = " + str(len(blocsCholerae))

    # Histogramme de la distribution de GC
    print "Histogram with the distribution of GC in each bloc, uncomment to see graphics."
    # histogrammeGCBlocs(blocsColi, "E. Coli")
    # histogrammeGCBlocs(blocsCholerae, "Vibrio Cholerae")

    # ORFs du genome choisi et creation du fichier FASTA
    if os.path.exists("orfFasta.fsa"):
        os.remove("orfFasta.fsa")
    print "Creation of FASTA file for submission at NCBI's glimmer tool for Vibrio Cholerae"
    print str(len(orfs(newOrder(cholerae)))) + " = Length of the ORFs"

    # Extraire sequences codantes des genes a partir du fichier tab
    if os.path.exists("genes_codants_Vibrio_cholerae.fasta"):
        os.remove("genes_codants_Vibrio_cholerae.fasta")
    print "Creation of FASTA file with coding genes for Vibrio Cholerae"
    genes(cholerae, tabCholerae)

    # print "Attention, on a 74122 genes codants si on calcule les ORFs, mais avec le fichier tab, on en trouve que 3504"
    # print "On va supposer que c'est bon car on trouve 37275 ATG dans le genome du Vibrio Cholerae"
    # print "Et donc, si on prend en compte que on doit calculer aussi son complementaire, ça fait du sens"

    # Creation d'un fichier avec les genes non codants, cad,
    # les parties du gene qui ne sont pas dans le fichier tab
    if os.path.exists("genes_non_codants_Vibrio_cholerae.fasta"):
        os.remove("genes_non_codants_Vibrio_cholerae.fasta")
    print "Creation of FASTA file with non-coding genes for Vibrio Cholerae"
    nonCodants(cholerae, tabCholerae)

    # Calculer pour chaque gene son pourcentage en GC et mettre a jour le fichier tab
    print "Updating FASTA file with coding genes and adding percentage in GC nucleotides"
    updateTabFile("genes_codants_Vibrio_cholerae.fasta")

    # Histogramme des pourcentages en GC des gènes
    print "Histogram with the distribution of GC in each coding gene for Vibrio Cholerae, uncomment to see graphic"
    # histogramGCgenes("genes_codants.fasta", "Pourcentage en GC des genes codants de Vibrio Cholerae")

    print "Now we do everything we need to do to compare it to E. Coli, we will be using the tab file"
    if os.path.exists("genes_codants_Escherichia.coli.fasta"):
        os.remove("genes_codants_Escherichia.coli.fasta")
    print "\tCreation of FASTA file with coding genes for E. Coli"
    genes(coli, tabColi)

    print "\tUpdating FASTA file with coding genes and adding percentage in GC nucleotides"
    updateTabFile("genes_codants_Escherichia.coli.fasta")

    print "\tHistogram with the distribution of GC in each coding gene for E. Coli, uncomment to see graphic"
    # histogramGCgenes("genes_codants_Escherichia.coli.fasta", "Pourcentage en GC des genes codants de E. Coli")

    print "Cleaning up"
    cleanUp()


####################
#   </Seance 1>    #
####################


##################
#   <Seance 2>   #
##################


# Histogramme generique
def histogramme(list, name):
    plt.hist(list)
    plt.title(name)
    plt.show()


# Longueurs des genes predits avec Glimmer
def lengthGlimmer(glimmerFile, name):
    f = open(glimmerFile, 'r')
    lengths = []
    for line in f:
        if not line.startswith("orf"):
            continue
        lengths.append(abs(int(line.split("    ")[1]) - int(line.split("    ")[2])))
    # histogramme(lengths, name)
    return lengths


# Longueur des genes du fichier tab
def lengthTabFile(tabFile, name):
    f = open(tabFile, 'r')
    lengths = []
    for line in f:
        if line.startswith("#"):
            continue
        lengths.append(abs(int(line.split("\t")[2]) - int(line.split("\t")[3])))  # Start - End
    # histogramme(lengths, name)
    return lengths


# Fonction qui cree une liste de 1s et 0s en comparant le fichier de prediction (ORFs ou Glimmer) et le genome
def binaryGlimmer(glimmerFile, genomeFile):  # TODO: On ne se preoccupe pas du frame, on s'interesse qu'aux positions
    glimmer = open(glimmerFile, 'r')
    genome = newOrder(genomeFile)
    binary = [0] * len(genome)
    for line in glimmer:
        if not line.startswith("orf"):
            continue
        start = line.split("    ")[1]
        end = line.split("    ")[2]
        if int(start) > int(end):
            start, end = end, start
        for i in range(int(start), int(end)):
            binary[i] = 1
    return binary


# Fonction qui cree une liste de 1s et 0s en comparant le fichier .tab et le genome
def binaryTab(tabFile, genomeFile):
    tab = open(tabFile, 'r')
    genome = newOrder(genomeFile)
    binary = [0] * len(genome)
    for line in tab:
        if line.startswith("#"):
            continue
        start = line.split("\t")[2]
        end = line.split("\t")[3]
        for i in range(int(start), int(end)):
            binary[i] = 1
    return binary


# Comparation des annotations (binary) issues de glimmer et du fichier tab
def compare_intervalle(glimmerBinary, tabBinary):
    vp = vn = fn = fp = 0
    length = len(glimmerBinary)
    for i in range(0, length):
        glimmer = glimmerBinary[i]
        tab = tabBinary[i]
        if glimmer + tab == 2:
            vp += 1
        elif glimmer + tab == 0:
            vn += 1
        elif glimmer > tab:
            fn += 1
        else:
            fp += 1
    return [[vn, fp], [fn, vp]]


# Fonction qui, a partir du resultat (matrice) de la fonction compare_intervalle
# Donne la specifité, la sensibilite et le taux de vrais poistifs
def senSpe(matrice):
    vraip = matrice[1][1]
    sen = vraip / float(vraip + matrice[1][0])
    vn = matrice[0][0]
    spe = vn / float(vn + matrice[0][1])
    vp = vraip / float(vraip + matrice[0][1])
    return sen, spe, vp


####################
#   </Seance 2>    #
####################

def main():
    seance1()

    ###################
    #   <Seance 2>    #
    ###################

    # Produire un fichier d'annotation des gènes en utilisant Glimmer -> glimmer_data.txt
    # Longueurs des genes predits avec Glimmer
    print "Calculating lengths of ORFs from the Glimmer file for Vibrio Cholerae, uncomment to see graphics"
    lengthGlimmer("./glimmer_data.txt", "Lengths of Glimmer file for Vibrio Cholerae, uncomment to see graphics")

    # Longueur des genes du fichier tab
    print "Calculating lengths from tab file of Vibrio Cholerae"
    lengthTabFile(tabCholerae, "Lengths of tab file for Vibrio Cholerae")

    # Maintenant la meme chose pour E. Coli
    print "We do the same thing for E. Coli, uncomment to see graphics"
    print "\tCalculating lengths of ORFs from the Glimmer file for E. Coli, uncomment to see graphics"
    lengthGlimmer("glimmer_data_coli.txt", "Lengths of tab file for E. Coli")

    print "Creating binary list from the Glimmer file of Vibrio Cholerae and its genome"
    glimmerBinary = binaryGlimmer("./glimmer_data.txt", cholerae)

    print "Creating binary list from the .tab file of Vibrio Cholerae and its genome"
    tabBinary = binaryTab(tabCholerae, cholerae)

    print "Comparing both glimmer and tab binaries: "
    annotation = compare_intervalle(glimmerBinary, tabBinary)
    print "\tTrue negatives: " + str(annotation[0][0])
    print "\tFalse positives: " + str(annotation[0][1])
    print "\tFalse negatives: " + str(annotation[1][0])
    print "\tTrue positives: " + str(annotation[1][1])

    print "And lastly, we look at the sensitivity and the specifity:"
    sen, spe, vp = senSpe(annotation)
    print "\tSensibility = " + str(sen)
    print "\tSpecificity = " + str(spe)
    print "\tRate of true positives = " + str(vp)

    ####################
    #   </Seance 2>    #
    ####################


main()
