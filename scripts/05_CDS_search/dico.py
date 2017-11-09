import string

def dico(F1):    
    dicoco = {}
    while 1:
        next2 = F1.readline()
        if not next2:
            break
        if next2[0] == ">":
            fasta_name_query = next2[:-1]
            Sn = string.split(fasta_name_query, "||")
            fasta_name_query = Sn[0]
            next3 = F1.readline()
            fasta_seq_query = next3[:-1]
            dicoco[fasta_name_query]=fasta_seq_query    
    return(dicoco)
