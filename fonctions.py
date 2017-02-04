from codeGenetique import codeGenetique

def pourcentage(seq): # 6.1
    longueur = len(seq)
    nA = (float(seq.count("A")) / longueur) * 100
    print ('A = ' + str(nA) + ' %')
    nC = (float(seq.count("C")) / longueur) * 100
    print ('C = ' + str(nC) + ' %')
    nG = (float(seq.count("G")) / longueur) * 100
    print ('G = ' + str(nG) + ' %')
    nT = (float(seq.count("T")) / longueur) * 100
    print ('T = ' + str(nT) + ' %\n')

def nombreCG(seq): # 6.2
    nombre = float(seq.count("CG"))
    nCG = (nombre/len(seq))
    return [nombre, nCG]

def listeCodons(seq): # 6.3
    result = []
    i = 0
    while (i < len(seq)):
        result.append(seq[i:i+3]) 
        i += 3   
    return result

def orfs(seq): # 6.4
    result = []
    start = 'ATG'
    stop = ['TAA', 'TAG', 'TGA']
    i = 0
    while (i < len(seq)):
        codon = seq[i:i+3]
        result.append(codon)
        if codon == start:
            result = []
            result.append(codon)
        if codon in stop:
            break
        i += 3
    return result

def complementaire(seq): # 6.5
    dic = {
        'A' : 'T',
        'T' : 'A',
        'C' : 'G',
        'G' : 'C'
    }
    result = ''
    for i in seq:
        result += dic[i]
    return result

def traduction(gene): # 6.6
    result = []
    i = 0
    while (i < len(gene)):
        codon = gene[i:i+3]
        if (len(codon)%3 == 0):
            result.append(codeGenetique[codon])
        i += 3
    return result

def pourcentageFichier (fichier): # 7.1
    data = newOrder(fichier)
    pourcentage(data)

def newOrder(fichier): # Aux
    f = open(fichier, 'r')
    rawData = ''
    for line in f:
        if line.startswith(">"):
            continue
        rawData += line
    data = rawData.replace('\n', '')
    return data

def genes(fichier1, fichier2): # 7.2
    f1 = open(fichier1, 'r')
    genome = newOrder(fichier2)
    next(f1)
    fasta = open('fasta.txt','a+')
    for line in f1:
        replicon = line.split("\t")[2]
        accesion = line.split("\t")[3]
        orientation = line.split("\t")[4]
        identifiant = line.split("\t")[5]
        data = genome[int(replicon) - 1:int(accesion) - 3]
        if orientation == '-':
            data = data[::-1]
        # Maintenant on fait le fichier avec format ASTA
        fasta.write(">" + identifiant + "\n")
        fasta.write(data + "\n")

def pourcentageFASTA(fichier): # 7.3
    f = open(fichier, 'r')
    id = ''
    for line in f:
        if line.startswith(">"):
            id += line
            continue
        print(id)
        pourcentage(line)
        id = ''   

def writeToFile(trad, id): # Aux
    f = open('fastaProteines.txt', 'a+')
    f.write(id + "\n")
    for item in trad:
        f.write(item)
    f.write("\n\n")

def traductionGenes(fichier): # 7.4
    f = open(fichier, 'r')
    id = ''
    for line in f:
        if line.startswith(">"):
            id += line
            continue
        trad = []
        trad = traduction(line[0:len(line) - 1])
        writeToFile(trad, id)
        id = ''
    return

def tailleMoyenne(fichier): # 7.5
    f = open(fichier, 'r')
    proteinCounter = 0
    acc = 0
    for line in f:
        if line.startswith(">"):
            proteinCounter += 1
            continue
        line = line[0:len(line) - 1]
        acc += len(line)
    taille = acc/proteinCounter
    print("Taille moyenne: " + str(taille))

def final(fichier, total):
    data = newOrder(fichier)
    percentage = len(data)/total
    print(str(percentage))

def main():
    # pourcentage("ACCTGGACAT")
    # for item in nombreCG("AACGTGGCACG"):
    #     print item
    # for item in listeCodons("AACGTGGCA"):
    #     print item
    # for item in orfs("ATGATGGCCCTHTAA"):
    #     print item
    # print complementaire('AACGTGGCA')
    # for item in traduction('')
    #     print item
    # pourcentageFichier('../data/Escherichia.coli/Escherichia.coli.genome')
    # genes('../data/Escherichia.coli/Escherichia.coli.tab', '../data/Escherichia.coli/Escherichia.coli.genome')
    # pourcentageFASTA('fasta.txt')
    # traductionGenes('fasta.txt')
    # En parametre passer les fichiers genome et tout
    # tailleMoyenne('fastaProteines.txt')
    final('fasta.txt', 58021)

main()
