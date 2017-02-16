import math

from codeGenetique import codeGenetique
from freqCodons import freq_par_aa
import matplotlib.pyplot as plt
import os.path


def pourcentage(sequence):
    length = len(sequence)
    D = {
        'A': (float(sequence.count('A')) / length),
        'C': (float(sequence.count('C')) / length),
        'G': (float(sequence.count('G')) / length),
        'T': (float(sequence.count('T')) / length)
    }
    return D


def newOrder(fichier):  # 7.1 AUX
    f = open(fichier, 'r')
    rawData = ''
    for line in f:
        if line.startswith(">"):
            continue
        rawData += line
    data = rawData.replace('\n', '')
    return data


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
            newFile.write("> ORF numero " + str(counter) + "\n")
            newFile.write("> Position = " + str(position) + "\n")
            counter += 1
            newFile.write(codon + "\n")
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


def histogramme(fichier):
    f = open(fichier, "r")
    lengths = []
    for line in f:
        if line.startswith(">"):
            continue
        lengths.append(len(line))
    plt.hist(lengths, range(0, 1000))
    plt.show()


def seuilToTab(fichier):  # On prend comme seuil 201
    f = open(fichier, "r")
    newFile = open("seuil.tab", "a+")
    newFile.write("START\tEND\tORF\n")
    for line in f:
        if line.startswith(">"):
            if line.find("Position"):
                line_aux = line
                pos = line_aux.replace("> Position = ", "")
                pos = pos.replace("\n", "")
            continue
        length = len(line)
        if length > 201:
            fin = int(pos) + len(line)
            newFile.write(str(pos) + "\t" + str(fin) + "\t" + line)


def senSpe(matrice):
    vraip = matrice[1][1]
    sen = vraip / float(vraip + matrice[1][0])
    vn = matrice[0][0]
    spe = vn / float(vn + matrice[0][1])
    vp = vraip / float(vraip + matrice[0][1])
    return sen, spe, vp


def ecrit_intervalle(debut, fin, length):
    result = [0] * length
    for i in range(0, len(debut)):
        start = debut[i] - 1
        end = fin[i]
        for i in range(start, end):
            result[i] = 1
    return result


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
                identifiant + "\t" + str(int(replicon) - 1) + "\t" + str(int(accesion) - 3) + "\t" + data + "\n"
            )
    fasta.close()
    return nouveauFichier, length


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


def genesToList(fichierGenes, length):
    f1 = open(fichierGenes, 'r')
    result = [0] * length
    for line in f1:
        gene = line.split("\t")[3]
        if not gene.startswith("ATG"):
            continue
        replicon = line.split("\t")[1]
        accesion = line.split("\t")[2]
        for j in range(int(replicon), int(accesion)):
            result[j] = 1
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


def evaluation(fichierGenes, fichierSeuil, length):
    genes = genesToList(fichierGenes, length)
    seuil = seuilToList(fichierSeuil, length)
    return compare_intervalle(genes, seuil)


def compteCodons(sequences):
    # On resoit une liste de frequences
    total = 0
    for item in sequences:
        for i in range(0, len(item), 3):
            codon = item[i:i + 3]
            if not (len(codon) % 3 == 0):
                break
            freq_par_aa[codon] += 1
            total += 1
    total += 64
    # D = {key:value for key, value in freq_par_aa.items()}
    D = freq_par_aa
    for i in D:
        D[i] = (float(freq_par_aa[i]) / total) * 100
        print ("KEY: " + i + " VALUE: " + str(D[i]))
    return D


def logproba(liste_lettres, tuple_frequences):  # tuple frequences cest avec la fonction frequences
    total = 0
    for i in liste_lettres:
        total += math.log(float(tuple_frequences[i]), 2)
    return total


def logprobagene(sequence, frequences):  # frequences c'est le dictionaire de compteCodons
    total = 0
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        if not (len(codon) % 3 == 0):
            break
        total += math.log(frequences[codon], 2)
    return total


def isItGene(sequence, dictLettres, dictCodons):
    logLettres = logproba(sequence, dictLettres)
    logCodon = logprobagene(sequence, dictCodons)
    if logCodon - logLettres > 1:
        print "Oui :)"
    else:
        print "Non!!!!!!!!!!!!!!!!!"
    return


def graphique(genome, fenetre):
    print "Holita"

def main():
    # data = newOrder("data/Escherichia.coli/Escherichia.coli.genome")
    # orfFile = orfs(data)
    # histogramme(orfFile)
    # seuilToTab(orfFile)

    # matrice = [[2.0, 3.0], [4.0, 5.0]]
    # sen, spe, vp = senSpe(matrice)
    # print("SEN = " + str(sen) + "\nSPE = " + str(spe) + "\nVP = " + str(vp))

    # intervalle = ecrit_intervalle([3, 11, 16], [8, 13, 21], 22)
    # print intervalle

    # genome = [0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,1,1,1,0]
    # orfs = [0,0,0,1,1,1,1,0,0,0,1,1,1,0,0,0,0,1,1,1,1,0]
    # result = compare_intervalle(genome, orfs)
    # print result

    # genes("data/Escherichia.coli/Escherichia.coli.tab", "data/Escherichia.coli/Escherichia.coli.genome")

    # data = newOrder("data/Escherichia.coli/Escherichia.coli.genome")
    # result = evaluation("genesEscherichia.coli.txt", "seuil.tab", len(data))
    # print "4,641,652 == " + str(len(data))
    # print "[[vn, fp], [fn, vp]]"
    # print result
    # sen, spe, vp = senSpe(result)
    # print ("SEN = " + str(sen) + "\nSPE = " + str(spe) + "\nVP = " + str(vp))
    # # Pour E.coli:
    # # 4, 641, 652 == 4641652
    # # [[vn, fp], [fn, vp]]
    # # [[1906908, 1027207], [6606, 1700931]]
    # # SEN = 0.996131269776
    # # SPE = 0.649909086726
    # # VP = 0.623476891565


    # On va utiliser le fichier seuil.txt, avec les genes commensant par START et finissant par STOP
    f1 = open("seuil.tab", 'r')
    sequences = []
    for line in f1:
        if line.startswith("START"):
            continue
        sequence = line.split("\t")[2]
        sequences.append(sequence[0:len(sequence) - 1])
    # sequences = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGAC"
    DCodons = compteCodons(sequences)
    data = newOrder("data/Escherichia.coli/Escherichia.coli.genome")
    DLettres = pourcentage(data)

    for seq in sequences:
        isItGene(seq, DLettres, DCodons)

    print ('Holita')


main()
