import string, itertools

def dico(F1):
    dicoco = {}
    with open(F1, "r") as file:
        for name, query in itertools.zip_longest(*[file]*2):
            if name[0] == ">":
                fasta_name_query = name[:-1]
                Sn = str.split(fasta_name_query, "||")
                fasta_name_query = Sn[0]                
                fasta_seq_query = query[:-1]
                dicoco[fasta_name_query] = fasta_seq_query    
    return(dicoco)
