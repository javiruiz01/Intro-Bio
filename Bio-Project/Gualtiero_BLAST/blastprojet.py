import matplotlib.pyplot as plt


def keep_sup_tenpow3(fileblast,savefile): #garde les aligenements qui on un evalue inferieur a 10e-3
    print'Question : 3.a'
    print'deleting sequences with an E-value superior to 10e-3'
    blastpresult = open(fileblast, 'r')
    resultData = ''
    for line in blastpresult:
        line = line.split()
        if(float(line[10])< 0.001):
            for word in line:
                resultData += str(word) + "\t"
            resultData += "\n"
    savfile = open(savefile, 'w')
    savfile.write(resultData)
    blastpresult.close()
    savfile.close()
    print'done'


def couverture(fileblast,fileprot,resultfile):
    print'Question : 3.b/c'
    print'calculating the coverage of the coding sequences'
    print'the sequences that have a coverage over 80% are kept'
    fblast = open(fileblast, 'r')
    fprot = open(fileprot, 'r')
    fresult = open(resultfile, 'w')
    res = ''
    nb_valid = 0
    blastlines = fblast.readlines()
    protlines = fprot.readlines()
    for blastline in blastlines:
        blastline = blastline.split()
        for protline in protlines:
            protline = protline.split()
            if protline[0].startswith(">"):
                nameseq = list()
                nameprot = protline[0].replace('>', '')
                nameprot = nameprot.replace('\n', '')
                nameseq .append(nameprot)
            else:
                nameseq.append(blastline[0])
                if(blastline[0] == nameseq[0]):
                    lenprot = len(nameseq[1])-1
                    lenblast = int(blastline[7])-int(blastline[6])
                    cover = lenblast / float(lenprot)
                    if (cover > 0.8):
                        res += str(blastline[0]+'\t'+blastline[1])+"\n"
                        nb_valid += 1
    fresult.write(res)
    fblast.close()
    fprot.close()
    fresult.close()
    print "number of valid sequences  = ", nb_valid
    print'done'

def aligements_significatif(filetab,filecorrespondace,resultfile):
    print'Question : 3.d'
    print'adding a column in the .tab file in first position'
    print'the column represents the number of significant alignements for each sequence'
    ftab = open(filetab, 'r')
    fcorr = open(filecorrespondace, 'r')
    tablines = ftab.readlines()
    corrlines = fcorr.readlines()
    res = ''
    for tabline in tablines:
        tablinelist = tabline.split()
        corrcount = 0
        for corrline in corrlines:
            corrline =corrline.split()
            if(tablinelist[5]==corrline[0]):
                corrcount += 1
        res += str(corrcount) + '\t' + tabline
    fresult = open(resultfile, 'w')
    fresult.write(res)
    ftab.close()
    fcorr.close()
    fresult.close()
    print'done'

def merge_alignements(alignement_1,aligenement_2,resultfile):
    print'Question : 3.d followup'
    print'merging the alignements of the two chromosomes of Vibrio cholerae'
    falign_1 = open(alignement_1, 'r')
    falign_2 = open(aligenement_2, 'r')
    fresult = open(resultfile, 'w')
    alignlines_1 = falign_1.readlines()
    alignlines_2 = falign_2.readlines()
    res = ''
    for line in range(1,len(alignlines_1)):
        alignline_1 = alignlines_1[line].split()
        if(int(alignline_1[0]) != 0):
            res += alignlines_1[line]
        else:
            res += alignlines_2[line]
    fresult.write(res)
    falign_1.close()
    falign_2.close()
    fresult.close()
    print'done'

def gcPercentage(sequence):
    G = float(sequence.count('G'))/len(sequence)
    C = float(sequence.count('G'))/len(sequence)
    return G + C

def no_alignement(filetab,filefasta):
    print('Question : 3.e')
    print('searching sequences wtih no alignement')
    no_alig = 0
    ismyprot = False
    perncentagelist = list()
    ftab = open(filetab, 'r')
    ffasta = open(filefasta,'r')
    tablines = ftab.readlines()
    fastalines = ffasta.readlines()
    for tabline in tablines:
        tablinelist = tabline.split()
        if(int(tablinelist[0])==0):
            no_alig += 1
            for fastaline in fastalines:
                fastaline = fastaline.split()
                if fastaline[0].startswith(">"):
                    nameprot = fastaline[0].replace('>', '')
                    nameprot = nameprot.replace('\n', '')
                    if(tablinelist[6] == nameprot):
                        ismyprot = True
                else:
                    if(ismyprot == True):
                        perncentagelist.append(gcPercentage(fastaline[0]))
                        ismyprot = False
    print "number of unaligned sequences  = ",no_alig
    print('generating histogram...')
    plt.hist(perncentagelist)
    plt.title("percentage of GC in non coding sequences of cholerae")
    plt.xlabel("percentage")
    plt.ylabel("number of sequences")
    plt.savefig('Fig/nGCcholerae')
    ffasta.close()
    ftab.close()
    print'done'

def Onilne_aligement():
    print('Question : 4.a/b')
    print('most observed family when aligning non coding sequences against nr database : g-proteobacteria')


def add_COG_id_function(corrfile,COGfile,resultfile):
    print('Question : 5/6')
    print('searching sequences wtih no alignement')
    fcorr = open(corrfile, 'r')
    fcog = open(COGfile,'r')
    res = ''
    corrlines = fcorr.readlines()
    coglines = fcog.readlines()
    for corrline in corrlines:
        corrlinelist = corrline.split('\t')
        for cogline in coglines:
            coglinelist = cogline.split('\t')
            if(corrlinelist[1].replace('\n','') == coglinelist[6]):
                corrline = corrline.replace('\n','')
                res += corrline +'\t'+ coglinelist[0]+'\t'+coglinelist[16]+'\n'
                break
    fresult = open(resultfile,'w')
    fresult.write(res)
    fcog.close()
    fcorr.close()
    fresult.close()
    print('there are 30 genes that do not have a classe fonctionnelle')
    print('and 88 others that are in the S class')
    print('done')

def histoClassCOG(blastfile):
    print('Question : 7')
    print('generating the histogram of COG ')
    fblast = open(blastfile, 'r')
    blastlines = fblast.readlines()
    ClassCOG = list()
    for blastline in blastlines:
        blastline = blastline.split()
        ClassCOG.append(blastline[3])
    dict = {}
    for lettre in ClassCOG:
        if(dict.has_key(lettre)==False):
            dict[lettre] = ClassCOG.count(lettre)
    plt.bar(range(len(dict)), dict.values(), align='center')
    plt.xticks(range(len(dict)), dict.keys())
    plt.title("COG functional Classes of cholerae")
    plt.xlabel("COG functional Classes")
    plt.ylabel("number of COG")
    plt.savefig('Fig/nCOGCholerae')
    print('done')




#keep_sup_tenpow3('Blast/cholerae_blast_2','Blast/cholerae_tenpow_2')
#couverture('Blast/cholerae_tenpow_1','cholerae_proteines/proteines_codants_Vibrio_cholerae_1.fasta','Blast/correspondance_1.blast')
#aligements_significatif('Data/Vibrio_cholerae.tab','Blast/correspondance_2.blast','Blast/correspondance_aligements_2.blast')
#merge_alignements('Blast/correspondance_aligements_1.blast','Blast/correspondance_aligements_2.blast','Blast/correspondance_aligements_full')
#no_alignement('Blast/correspondance_aligements_full','cholerae_nucleotides_codants/genes_codants_Vibrio_cholerae_full.fasta')
#Onilne_aligement()
#add_COG_id_function('Blast/correspondance_full.blast','Data/Escherichia_coli_COG_annotation.tsv','Blast/corresp_COG_full.blast')
histoClassCOG('Blast/corresp_COG_full.blast')
