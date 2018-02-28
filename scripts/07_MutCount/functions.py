import itertools

def simplify_fasta_name(fasta_name,LT):
    for abbreviation in LT:
        if abbreviation in fasta_name:
            new_fasta_name = abbreviation

    return(new_fasta_name)

## Generates bash, with key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)
def dico(fasta_file,LT):
    #count_fastaName = 0
    bash1 = {}
    with open(fasta_file, "r") as file:
        for name, query in itertools.izip_longest(*[file]*2):
            if not name:
                break
            if name[0] == ">":
                #count_fastaName += 1
                fasta_name = name[1:-1]
                sequence = query[:-1]
                if fasta_name not in bash1.keys():
                    fasta_name = simplify_fasta_name(fasta_name, LT)
                    bash1[fasta_name] = sequence
                else :
                    print fasta_name

    kk = bash1.keys()
    key0 = kk[0]
    seq0 = bash1[key0]
    ln_seq = len(seq0)
    
    return(bash1)