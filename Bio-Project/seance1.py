# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os.path
import shutil

# On travaille avec les genomes Vibrio Cholerae et E. Coli
coli = "./data/Escherichia.coli.fasta"
cholerae = "./data/Vibrio_cholerae_genome.fasta"
tabColi = "./data/Escherichia.coli.tab"
tabCholerae = "./data/Vibrio_cholerae.tab"


# Change un fichier fasta en un String sans saut de ligne
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
    print len(blocs)
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


def atPercent(sequence):
    D = {
        'A': (float(sequence.count('A')) / len(sequence)),
        'T': (float(sequence.count('T')) / len(sequence))
    }
    return D


# Histogramme de la distribution de GCs
def histogrammeGCBlocs(blocs, name):
    percents = []
    positions = []
    position = 0
    for bloc in blocs:
        positions.append(position)
        position += len(bloc)
        d = gcBlocPercentage(bloc)
        percents.append(d['C'] + d['G'])
    print len(blocs)
    x_axis_list = range(len(positions))
    plt.plot(x_axis_list, percents)
    plt.title(name)
    plt.savefig("./histogrammes/" + name + ".png")
    plt.show()

def histogrammeGCBlocs_hist(blocs, name):
    percents = []
    for bloc in blocs:
        d = gcBlocPercentage(bloc)
        percents.append(d['C'] + d['G'])
    print len(blocs)
    plt.plot(percents)
    plt.title(name)
    plt.savefig("./histogrammes/" + name + "_hist.png")
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
        if i == "Y" or i == "S" or i == "M" or i == "W" or i == "K" or i == "N" or i == "R":
            continue
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
    fastaGeneric = open("genes_" + os.path.basename(tabFile)[0:-4] + ".fasta", "a+")
    globalCounter = 0
    counter = 0
    counter_non_codants = 0
    codonStart = False
    for line in tab:
        start = line.split("\t")[2]  # Start
        end = line.split("\t")[3]  # End
        orientation = line.split("\t")[4]
        identifiant = line.split("\t")[5]
        if orientation == '-':
            data = genome[int(start) - 1: int(end)][::-1]
            if ("Y" in data) or ("S" in data) or ("M" in data) or ("W" in data) or ("K" in data) or ("N" in data) or (
                        "R" in data):
                continue
            data = complementaire(data)
        else:
            data = genome[int(start) - 1:int(end) - 3]
            if ("Y" in data) or ("S" in data) or ("M" in data) or ("W" in data) or ("K" in data) or ("N" in data) or (
                        "R" in data):
                continue
        data = ''.join(data)
        nStopTAA = int(data.count("TAA"))
        nStopTAG = int(data.count("TAG"))
        nStopTGA = int(data.count("TGA"))
        globalCounter += 1
        fastaGeneric.write(">" + identifiant + "\n")
        fastaGeneric.write(data + "\n")
        if (not data.startswith("ATG")) or not ((nStopTAA + nStopTAG + nStopTGA) != 0):
            counter_non_codants += 1
            fastaNonCodants.write(">" + identifiant + "\n")
            fastaNonCodants.write(data + "\n")
            continue
        counter += 1
        fastaFile.write(">" + identifiant + "\n")
        fastaFile.write(data + "\n")
    fastaFile.close()
    tab.close()
    fastaGeneric.close()
    return fastaGeneric


# Calculer pour chaque gene son pourcentage en GC et mettre a jour le fichier tab
def updateFastaFile(genes_codants, name):
    f = open(genes_codants, 'r')
    newFile = open("genes_updated_" + name + ".fasta", 'a+')
    for line in f:
        if line.startswith(">"):
            newFile.write(line)
            continue
        d = gcBlocPercentage(line[0:len(line) - 1])
        percent = float(d['C'] + d['G']) * 100
        newFile.write(str(percent) + "\t" + line)
    f.close()
    newFile.close()
    # os.rename(genes_codants, genes_codants + "_old")
    # os.rename("tmp.fasta", "genes.fasta")
    return newFile.name


# Update le fichier .tab, en rajoutant le pourcentage de gc dans le genome, prendre fichier updated
def updateTabFile(genes_codants, tabFile):
    f = open(genes_codants, "r")
    tab = open(tabFile, "r")
    newFile = open(tabFile + "_updated", "a+")
    d = dict()
    tmp = ''
    for line in f:
        if line.startswith(">"):
            tmp = line[1:-1]
            continue
        d[tmp] = line.split("\t")[0]
        tmp = ''
    for line in tab:
        if line.startswith("#"):
            newFile.write(line[0:len(line) - 2] + "\tPerc_GC\n")
            continue
        geneId = line.split("\t")[5]
        if geneId not in d.keys():
            continue
        newFile.write(line[0:len(line) - 2] + "\t" + d[geneId] + "\n")
    f.close()
    tab.close()
    newFile.close()
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
    plt.savefig("./histogrammes/" + name + ".png")
    plt.show()


# Fonction pour ranger tous les noveaus fichier crees
def cleanUp():
    if os.path.exists("fasta_files_1_genome"):
        shutil.rmtree('./fasta_files_1_genome/')
    os.makedirs("./fasta_files_1_genome")
    for file in os.listdir("./"):
        if file.endswith("fasta") or file.endswith("fsa") or file.endswith("fasta_old") or file.endswith("updated"):
            os.rename(file, "./fasta_files_1_genome/" + file)


# Histogramme generique
def histogramme(list, name):
    plt.hist(list)
    plt.title(name)
    plt.savefig("./histogrammes/" + name + ".png")
    plt.show()


# Longueurs des genes predits avec Glimmer
def lengthGlimmer(glimmerFile, name):
    f = open(glimmerFile, 'r')
    lengths = []
    for line in f:
        if not line.startswith("orf"):
            continue
        lengths.append(abs(int(line.split("    ")[1]) - int(line.split("    ")[2])))
    histogramme(lengths, name)
    return lengths


# Longueur des genes du fichier tab
def lengthTabFile(tabFile, name):
    f = open(tabFile, 'r')
    lengths = []
    for line in f:
        if line.startswith("#"):
            continue
        lengths.append(abs(int(line.split("\t")[2]) - int(line.split("\t")[3])))  # Start - End
    histogramme(lengths, name)
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


# question 9, tab file est le fichier tab full, name est le nom du fichier quon va creer
def composition_nuc(tab_file, genome_1, genome_2, name):
    tab = open(tab_file, "r")
    gen_1 = newOrder(genome_1)
    gen_2 = newOrder(genome_2)
    newFile = open(name, "a+")
    newFile.write("Chromosome\tStart\tEnd\tFeature\tPercGC\tPercA\tPercT\n")
    for line in tab:
        if line.startswith("#"):
            continue
        start = line.split("\t")[2]  # Start
        end = line.split("\t")[3]  # End
        chromosome = line.split("\t")[0]
        name = line.split("\t")[1]
        if chromosome == "I":
            data = gen_1[int(start) - 1:int(end) - 3]
        else:
            data = gen_2[int(start) - 1: int(end) - 3]
        Dgc = gcBlocPercentage(data)
        Dat = atPercent(data)
        newFile.write(
            name + "\t" + str(start) + "\t" + str(end) + "\tNONE\t" + str(Dgc['C'] + Dgc['G']) + "\t" + str(
                Dat['A']) + "\t"
            + str(Dat['T']) + "\n")
    tab.close()
    newFile.close()
    return


# Decoupage en blocs de 50 bp, ils peuvent se chevaucher
def blocs50(fasta_file, chromosome_name):
    nucleotides = newOrder(fasta_file)
    blocs = []
    bloc = ''
    counter = 0
    for i in range(0, len(nucleotides)):
        if counter == 50:
            blocs.append(bloc)
            bloc = ''
            counter = 0
        bloc += nucleotides[i]
        counter += 1
    blocs.append(bloc)
    newFile = open("50bp.igv", "a+")
    if os.path.getsize("./50bp.igv") == 0:
        newFile.write("Chromosome\tStart\tEnd\tFeature\t%GC\n")
    position = 0
    for bloc in blocs:
        d = gcBlocPercentage(bloc)
        end = str(position + len(bloc))
        newFile.write(
            chromosome_name + "\t" + str(position) + "\t" + end + "\tNONE\t" + str(float(d['G'] + d['C'])) + "\n")
        position = int(end) + 1
    newFile.close()
    return newFile


def blocs50_v0(tab_file, genome1, genome2):
    tab = open(tab_file, "r")
    gen1 = newOrder(genome1)
    gen2 = newOrder(genome2)
    newFile = open("./50bp.igv", "a+")
    if os.path.getsize("./50bp.igv") == 0:
        newFile.write("Chromosome\tStart\tEnd\tFeature\t%GC\n")
    for line in tab:
        if line.startswith("#"):
            continue
        start = line.split("\t")[2]  # Start
        end = line.split("\t")[3]  # End
        chromosome = line.split("\t")[0]
        name = line.split("\t")[1]
        def_end = int(end)
        if chromosome == "I":
            data = gen1[int(start) - 1:int(end) - 3]
        else:
            data = gen2[int(start) - 1: int(end) - 3]
        for i in range(0, len(data), 50):
            if i == 0:
                position = int(start) + i
            else:
                position = int(end) + 1
            bloc = data[i:int(i) + 50]
            bloc = bloc[0: (bloc.count("A") + bloc.count("T") + bloc.count("G") + bloc.count("C"))]
            if len(bloc) != 50:
                end = def_end
            else:
                end = position + len(bloc)
            d = gcBlocPercentage(bloc)
            newFile.write(
                name + "\t" + str(position) + "\t" + str(end) + "\tNONE\t" + str(float(d['G'] + d['C'])) + "\n"
            )
        position = 0
    tab.close()
    newFile.close()
    return newFile


# Deuxieme fichier igv avec le pourcentage de GC des regions annotees
# Faire attention, glimmer file et genome file doivent etre du meme chromosome
def geneAnnotesIGV(glimmer_file, genome_file, chromosome_name):
    glimmer = open(glimmer_file, "r")
    genome = newOrder(genome_file)
    newFile = open("./genes_annotes.igv", "a+")
    if os.path.getsize("./genes_annotes.igv") == 0:
        newFile.write("Chromosome\tStart\tEnd\tFeature\t%GC\n")
    for line in glimmer:
        if not line.startswith("orf"):
            continue
        start = line.split("    ")[1]
        end = line.split("    ")[2]
        if int(start) > int(end):
            start, end = end, start
        d = gcBlocPercentage(complementaire(genome[int(start):int(end)]))
        newFile.write(
            chromosome_name + "\t" + start + "\t" + end + "\tNONE\t" + str(float(d['G'] + d['C'])) + "\n")
    newFile.close()
    glimmer.close()
    return newFile


def main():
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
    blocsCholerae_1 = decoupage(cholerae, tabCholerae)
    blocsCholerae_2 = decoupage("./data/Vibrio_cholerae_genome_2.fasta", "./data/Vibrio_cholerae_2.tab")
    print "Number of blocs that we will be analysing from Vibrio Cholerae = " + str(len(blocsCholerae_1))

    # Histogramme de la distribution de GC
    print "Histogram with the distribution of GC in each bloc, uncomment to see graphics."
    # histogrammeGCBlocs(blocsColi, "Distribution of GC in E. Coli by 1000 kbp")
    # histogrammeGCBlocs(blocsCholerae_1, "Distribution of GC in Vibrio Cholerae first chromosome by 1000 kbp")
    # histogrammeGCBlocs(blocsCholerae_2, "Distribution of GC in Vibrio Cholerae second chromosome by 1000 kbp")
    histogrammeGCBlocs_hist(blocsColi, "Distribution of GC in E. Coli by 1000 kbp (histogram)")
    histogrammeGCBlocs_hist(blocsCholerae_1, "Distribution of GC in Vibrio Cholerae first chromosome by 1000 kbp (histogram)")
    histogrammeGCBlocs_hist(blocsCholerae_2, "Distribution of GC in Vibrio Cholerae second chromosome by 1000 kbp (histogram)")

    # ORFs du genome choisi et creation du fichier FASTA
    if os.path.exists("orfFasta.fsa"):
        os.remove("orfFasta.fsa")
    print "Creation of FASTA file for submission at NCBI's glimmer tool for Vibrio Cholerae"
    print str(len(orfs(newOrder(cholerae)))) + " = Length of the ORFs"

    # Extraire sequences codantes des genes a partir du fichier tab
    if os.path.exists("genes_codants_Vibrio_cholerae.fasta"):
        os.remove("genes_codants_Vibrio_cholerae.fasta")
    print "Creation of FASTA file with coding genes and another FASTA file for non-coding genes for Vibrio Cholerae"
    genesFileGenericCholerae = genes(cholerae, tabCholerae)

    # print "Attention, on a 74122 genes codants si on calcule les ORFs, mais avec le fichier tab, on en trouve que 3504"
    # print "On va supposer que c'est bon car on trouve 37275 ATG dans le genome du Vibrio Cholerae"
    # print "Et donc, si on prend en compte que on doit calculer aussi son complementaire, ça fait du sens"

    # Calculer pour chaque gene son pourcentage en GC et mettre a jour le fichier tab
    print "Updating FASTA file with coding genes and adding percentage in GC nucleotides"
    updatedFile_Cholerae = updateFastaFile(genesFileGenericCholerae.name, "cholerae")
    updateTabFile(updatedFile_Cholerae, "./data/Vibrio_cholerae.tab")

    # Histogramme des pourcentages en GC des gènes
    print "Histogram with the distribution of GC in each coding gene for Vibrio Cholerae, uncomment to see graphic"
    histogramGCgenes("genes_updated_cholerae.fasta", "Pourcentage en GC des genes de Vibrio Cholerae - Chromosome 1")

    print "Now we do everything we need to do to compare it to E. Coli, we will be using the tab file"
    if os.path.exists("genes_codants_Escherichia.coli.fasta"):
        os.remove("genes_codants_Escherichia.coli.fasta")
    print "\tCreation of FASTA file with coding genes and another for non-coding genes for E. Coli"
    genesFileGenericColi = genes(coli, tabColi)

    print "\tUpdating FASTA file with coding genes and adding percentage in GC nucleotides"
    updateFile_Coli = updateFastaFile(genesFileGenericColi.name, "coli")
    updateTabFile(updateFile_Coli, "./data/Escherichia.coli.tab")

    print "\tHistogram with the distribution of GC in each coding gene for E. Coli, uncomment to see graphic"
    histogramGCgenes("genes_updated_coli.fasta", "Pourcentage en GC des genes de E. Coli")

    print "Cleaning up"
    cleanUp()

    # Produire un fichier d'annotation des gènes en utilisant Glimmer -> glimmer_data.txt
    # Longueurs des genes predits avec Glimmer
    print "Calculating lengths of ORFs from the Glimmer file for Vibrio Cholerae, uncomment to see graphics"
    # lengthGlimmer("./glimmer_data.txt", "Lengths of Glimmer file for Vibrio Cholerae Chromosome 1")
    # lengthGlimmer("./glimmer_data_2.txt", "Lengths of Glimmer file for Vibrio Cholerae Chromosome 2")

    # Longueur des genes du fichier tab
    print "Calculating lengths from tab file of Vibrio Cholerae"
    # lengthTabFile(tabCholerae, "Lengths of tab file for Vibrio Cholerae chromosome 1")
    # lengthTabFile("./data/Vibrio_cholerae_2.tab", "Length of tab file for Vibrio Cholerae chromosome 2")

    # Maintenant la meme chose pour E. Coli
    print "We do the same thing for E. Coli, uncomment to see graphics"
    print "\tCalculating lengths of ORFs from the Glimmer file for E. Coli, uncomment to see graphics"
    # lengthGlimmer("glimmer_data_coli.txt", "Lengths of Glimmer file for E. Coli")
    lengthTabFile(tabColi, "Lengths of tab file for E. Coli")

    print "Creating binary list from the Glimmer file of Vibrio Cholerae and its genome - Chromosome 1"
    glimmerBinary = binaryGlimmer("./glimmer_data.txt", cholerae)

    print "Creating binary list from the Glimmer file of Vibrio Cholerae and its genome - Chromosome 2"
    glimmerBinary_2 = binaryGlimmer("./glimmer_data_2.txt", "./data/Vibrio_cholerae_genome_2.fasta")

    print "Creating binary list from the .tab file of Vibrio Cholerae and its genome - Chromosome 1"
    tabBinary = binaryTab(tabCholerae, cholerae)
    print "Creating binary list from the .tab file of Vibrio Cholerae and its genome - Chromosome 2"
    tabBinary_2 = binaryTab("./data/Vibrio_cholerae_2.tab", "./data/Vibrio_cholerae_genome_2.fasta")

    print "Comparing both glimmer and tab binaries: Chromosome 1"
    annotation = compare_intervalle(glimmerBinary, tabBinary)
    print "Comparing both glimmer and tab binaries: Chromosome 2"
    annotation_2 = compare_intervalle(glimmerBinary_2, tabBinary_2)
    print "\t\t\tChromosome 1:"
    print "\tTrue negatives: " + str(annotation[0][0])
    print "\tFalse positives: " + str(annotation[0][1])
    print "\tFalse negatives: " + str(annotation[1][0])
    print "\tTrue positives: " + str(annotation[1][1])

    print "\t\t\tChromosome 2:"
    print "\tTrue negatives: " + str(annotation_2[0][0])
    print "\tFalse positives: " + str(annotation_2[0][1])
    print "\tFalse negatives: " + str(annotation_2[1][0])
    print "\tTrue positives: " + str(annotation_2[1][1])

    print "And lastly, we look at the sensitivity and the specifity - Chromosome 1:"
    sen, spe, vp = senSpe(annotation)
    print "\tSensibility = " + str(sen)
    print "\tSpecificity = " + str(spe)
    print "\tRate of true positives = " + str(vp)

    print "And lastly, we look at the sensitivity and the specifity - Chromosome 2:"
    sen_2, spe_2, vp_2 = senSpe(annotation_2)
    print "\tSensibility = " + str(sen_2)
    print "\tSpecificity = " + str(spe_2)
    print "\tRate of true positives = " + str(vp_2)

    print "Making nucleotides composition --> new File IVG.ivg"
    if os.path.exists("./IVG.ivg"):
        os.remove("./IVG.ivg")
    composition_nuc("./data/Vibrio_cholerae_full.tab", "./data/Vibrio_cholerae_genome.fasta",
                    "./data/Vibrio_cholerae_genome_2.fasta", "IVG.ivg")

    print "Making new file with percentage of GC in blocs of 50 bp of Vibrio Cholerae, first chromosome"
    if os.path.exists("./50bp.igv"):
        os.remove("./50bp.igv")
    blocs50("./data/Vibrio_cholerae_genome.fasta", "NC_002505.1")
    # blocs50_v0("./data/Vibrio_cholerae_full.tab", "./data/Vibrio_cholerae_genome.fasta",
    # "./data/Vibrio_cholerae_genome_2.fasta")

    print "\tCompleting the file with the second chromosome"
    blocs50("./data/Vibrio_cholerae_genome_2.fasta", "NC_002506.1")

    print "Making new file with percentage of GC in annotated genomes for Vibrio Cholerae, first chromosome"
    if os.path.exists("./genes_annotes.igv"):
        os.remove("./genes_annotes.igv")
    geneAnnotesIGV("./glimmer_data.txt", "./data/Vibrio_cholerae_genome.fasta", "NC_002505.1")

    print "\tCompleting the file with the second chromosome"
    geneAnnotesIGV("./glimmer_data_2.txt", "./data/Vibrio_cholerae_genome_2.fasta", "NC_002506.1")


main()
