from codeGenetique import codeGenetique
import matplotlib.pyplot as plt
import os.path


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
    # seqs = [seq, complementaire(seq)]
    # for s in seqs:
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
    sen = vraip / (vraip + matrice[1][0])
    vn = matrice[0][0]
    spe = vn / (vn + matrice[0][1])
    vp = vraip / (vraip + matrice[0][1])
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
            data = data[::-1]
        # Faut chercher ATG?
        # data = orf_aux(data)
        data = ''.join(data)
        # Maintenant on fait le fichier avec format FASTA
        if len(data) != 0:
            # fasta.write(">" + identifiant + "\n")
            fasta.write(
                identifiant + "\t" + str(int(replicon) - 1) + "\t" + str(int(accesion) - 3) + "\t" + data + "\n")
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


def evaluation(fichierGenes, fichierSeuil, length):
    result = [0] * length
    f1 = open(fichierGenes, "r")
    f2 = open(fichierSeuil, "r")
    for line in f1:
        print "holita que tal"

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

    print ('Holita')


main()
