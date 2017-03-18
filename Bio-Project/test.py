# -*- coding: utf-8 -*-
import os.path

tabColi = "./data/Escherichia.coli.tab"
coli = "./data/Escherichia.coli.fasta"
# coli = "./test.fasta"

stop = ['TAA', 'TAG', 'TGA']


# stop = ['AAT', 'GAT', 'AGT']


def newOrder(fichier):  # Works fine
    f = open(fichier, 'r')
    rawData = ''
    for line in f:
        if line.startswith(">"):
            continue
        rawData += line
    data = rawData.replace('\n', '')
    f.close()
    return data


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
        if i not in dic:
            continue
        result += dic[i]
    return result


def main():
    print "Removing previous test file"
    if os.path.exists("testFile.txt"):
        os.remove("testFile.txt")

    print "Commencing test"
    f = open(tabColi)
    genome = newOrder(coli)
    # print "Genome position [0:200] = " + genome[0:200]
    new = open("testFile.txt", 'a+')
    counter = 0
    for line in f:
        if line.startswith("#"):
            continue
        start = line.split("\t")[2]  # Start
        end = line.split("\t")[3]  # End
        orientation = line.split("\t")[4]
        if orientation != '-':
            continue
        counter += 1
        data = genome[int(start) - 1: int(end)]
        data = data[::-1]
        data = complementaire(data)
        data = ''.join(data)
        new.write("Sequence_" + str(counter) + " dans la position: [" + str(start) + ":" + str(end) + "]\n")
        new.write(data + "\n")
        # if counter == 1:
        #     lineTest = genome[0:int(accesion)]
    f.close()
    new.close()

    print "Testing if positions are correct"
    data = genome[35377:36162]
    print data
    # print "Data from 0 to first coding sequence:"
    # for i in range(0, len(data)):
    #     print "Position [" + str(i) + "] = " + data[i]
    print "Raw data, without the stop codon and with the start codon: " + data

    print "Testing how many coding genes start by ATG"
    test = open("testFile.txt", 'r')
    counter = 0
    for line in test:
        if line.endswith("ATG"):
            counter += 1
    print "\tNº of ATG starting lines: " + str(counter)
    test.close()

    test = open("testFile.txt", 'r')
    print "Testing how many coding genes start by GTA"
    counter = 0
    for line in test:
        if line.endswith("GTA"):
            counter += 1
    print "\tNº of GTA starting lines: " + str(counter)
    test.close()

    test = open("testFile.txt", 'r')
    print "Nº of lines that start with a stop codon, looking at it from the other side"
    counter = 0
    for line in test:
        if line.endswith("AAT") or line.startswith("GTA") or line.startswith("AGT"):
            counter += 1
    print  "\tCounter = " + str(counter)
    test.close()


main()
