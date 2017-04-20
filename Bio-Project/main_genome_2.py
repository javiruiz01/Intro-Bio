# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os.path
import shutil
import seance1

# On travaille avec les genomes Vibrio Cholerae et E. Coli
coli = "./data/Escherichia.coli.fasta"
cholerae = "./data/Vibrio_cholerae_genome_2.fasta"
tabColi = "./data/Escherichia.coli.tab"
tabCholerae = "./data/Vibrio_cholerae_2.tab"


# Fonction pour ranger tous les noveaus fichier crees
def cleanUp():
    if os.path.exists("fasta_files_2_genome"):
        shutil.rmtree('./fasta_files_2_genome/')
    os.makedirs("./fasta_files_2_genome")
    for file in os.listdir("./"):
        if file.endswith("fasta") or file.endswith("fsa") or file.endswith("fasta_old"):
            os.rename(file, "./fasta_files_2_genome/" + file)


def main():
    # print "Readying main python script"
    # f = open("seance1.py", "wb")
    # f.seek(-len(os.linesep), )
    print "GENOME 2 DE VIBRIO CHOLERAE"
    # Nombre de chromosomes
    totalColi = seance1.nChromosomes(coli)
    totalCholerae = seance1.nChromosomes(cholerae)
    print "Total number of chromosomes:\n\tE.Coli = " + str(totalColi) + "\n\tVibrio Cholerae = " + str(totalCholerae)

    # Longueur de chromosomes
    i = 0
    lengthColi = seance1.lengthChromosome(coli)
    for length in lengthColi:
        i += 1
        print "Length of the " + str(i) + " chromosome from E. Coli = " + str(length)
    i = 0
    lengthCholerae = seance1.lengthChromosome(cholerae)
    for length in lengthCholerae:
        i += 1
        print "Length of the " + str(i) + " chromosome from Vibrio Cholerae = " + str(length)

    # Longueur totale du genome
    print "Total length of E. Coli = " + str(seance1.totalLength(coli))
    print "Total length of Vibirio Cholerae = " + str(seance1.totalLength(cholerae))

    # Composition en nucleotides
    dColi = seance1.pourcentage(coli)
    for item in dColi:
        print "Percentage of " + item + " in E. Coli = " + str(dColi[item])
    dCholerae = seance1.pourcentage(cholerae)
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
    blocsColi = seance1.decoupage(coli, tabColi)
    print "Number of blocs that we will be analysing from E. Coli = " + str(len(blocsColi))
    blocsCholerae = seance1.decoupage(cholerae, tabCholerae)
    print "Number of blocs that we will be analysing from Vibrio Cholerae = " + str(len(blocsCholerae))

    # Histogramme de la distribution de GC
    print "Histogram with the distribution of GC in each bloc, uncomment to see graphics."
    # histogrammeGCBlocs(blocsColi, "Distribution of GC in E. Coli")
    # histogrammeGCBlocs(blocsCholerae, "Distribution of GC in Vibrio Cholerae")

    # ORFs du genome choisi et creation du fichier FASTA
    if os.path.exists("orfFasta.fsa"):
        os.remove("orfFasta.fsa")
    print "Creation of FASTA file for submission at NCBI's glimmer tool for Vibrio Cholerae"
    print str(len(seance1.orfs(seance1.newOrder(cholerae)))) + " = Length of the ORFs"

    # Extraire sequences codantes des genes a partir du fichier tab
    if os.path.exists("genes_codants_Vibrio_cholerae.fasta"):
        os.remove("genes_codants_Vibrio_cholerae.fasta")
    print "Creation of FASTA file with coding genes and another FASTA file for non-coding genes for Vibrio Cholerae"
    genesFileGenericCholerae = seance1.genes(cholerae, tabCholerae)

    # print "Attention, on a 74122 genes codants si on calcule les ORFs, mais avec le fichier tab, on en trouve que 3504"
    # print "On va supposer que c'est bon car on trouve 37275 ATG dans le genome du Vibrio Cholerae"
    # print "Et donc, si on prend en compte que on doit calculer aussi son complementaire, ça fait du sens"

    # Calculer pour chaque gene son pourcentage en GC et mettre a jour le fichier tab
    print "Updating FASTA file with coding genes and adding percentage in GC nucleotides"
    updateFileCholerae = seance1.updateFastaFile(genesFileGenericCholerae.name, "cholerae_2")
    seance1.updateTabFile(updateFileCholerae, "./data/Vibrio_cholerae_2.tab")

    # Histogramme des pourcentages en GC des gènes
    print "Histogram with the distribution of GC in each coding gene for Vibrio Cholerae, uncomment to see graphic"
    seance1.histogramGCgenes("genes_updated_cholerae_2.fasta", "Pourcentage en GC des genes de Vibrio Cholerae - Chromosome 2")

    print "Now we do everything we need to do to compare it to E. Coli, we will be using the tab file"
    if os.path.exists("genes_codants_Escherichia.coli.fasta"):
        os.remove("genes_codants_Escherichia.coli.fasta")
    print "\tCreation of FASTA file with coding genes and another for non-coding genes for E. Coli"
    genesFileGenericColi = seance1.genes(coli, tabColi)

    print "\tUpdating FASTA file with coding genes and adding percentage in GC nucleotides"
    seance1.updateFastaFile(genesFileGenericColi.name, "coli")

    print "\tHistogram with the distribution of GC in each coding gene for E. Coli, uncomment to see graphic"
    # histogramGCgenes("genes_codants_Escherichia.coli.fasta", "Pourcentage en GC des genes codants de E. Coli")

    cleanUp()

    # Produire un fichier d'annotation des gènes en utilisant Glimmer -> glimmer_data.txt
    # Longueurs des genes predits avec Glimmer
    print "Calculating lengths of ORFs from the Glimmer file for Vibrio Cholerae, uncomment to see graphics"
    seance1.lengthGlimmer("./glimmer_data.txt",
                          "Lengths of Glimmer file for Vibrio Cholerae, uncomment to see graphics")

    # Longueur des genes du fichier tab
    print "Calculating lengths from tab file of Vibrio Cholerae"
    seance1.lengthTabFile(tabCholerae, "Lengths of tab file for Vibrio Cholerae")

    # Maintenant la meme chose pour E. Coli
    print "We do the same thing for E. Coli, uncomment to see graphics"
    print "\tCalculating lengths of ORFs from the Glimmer file for E. Coli, uncomment to see graphics"
    seance1.lengthGlimmer("glimmer_data_coli.txt", "Lengths of tab file for E. Coli")

    print "Creating binary list from the Glimmer file of Vibrio Cholerae and its genome"
    glimmerBinary = seance1.binaryGlimmer("./glimmer_data_2.txt", cholerae)

    print "Creating binary list from the .tab file of Vibrio Cholerae and its genome"
    tabBinary = seance1.binaryTab(tabCholerae, cholerae)

    print "Comparing both glimmer and tab binaries: "
    annotation = seance1.compare_intervalle(glimmerBinary, tabBinary)
    print "\tTrue negatives: " + str(annotation[0][0])
    print "\tFalse positives: " + str(annotation[0][1])
    print "\tFalse negatives: " + str(annotation[1][0])
    print "\tTrue positives: " + str(annotation[1][1])

    print "And lastly, we look at the sensitivity and the specifity:"
    sen, spe, vp = seance1.senSpe(annotation)
    print "\tSensibility = " + str(sen)
    print "\tSpecificity = " + str(spe)
    print "\tRate of true positives = " + str(vp)


main()
