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