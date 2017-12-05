def simplify_fasta_name(fasta_name,LT):
    for abbreviation in LT:
        if abbreviation in fasta_name:
            new_fasta_name = abbreviation

    return(new_fasta_name)

## Generates bash, with key = fasta name; value = sequence (WITH GAP, IF ANY, REMOVED IN THIS FUNCTION)
def dico(fasta_file,LT):

    count_fastaName=0
    F1 = open(fasta_file, "r")
    
    bash1 = {}
    while 1:
        nextline = F1.readline()
        #print nextline
        if not nextline :
            break
        
        if nextline[0] == ">":
            count_fastaName = count_fastaName + 1
            fasta_name = nextline[1:-1]
            nextline = F1.readline()
            sequence = nextline[:-1]
            
            if fasta_name not in bash1.keys():
                fasta_name = simplify_fasta_name(fasta_name,LT)  ### DEF 0 ###
                bash1[fasta_name] = sequence
            else:
                print fasta_name

    # Find alignment length
    kk = bash1.keys()
    key0 = kk[0]
    seq0 = bash1[key0]
    ln_seq = len(seq0)

    F1.close()
    
    return(bash1)
