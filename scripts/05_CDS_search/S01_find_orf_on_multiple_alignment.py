#!/usr/bin/env python
# coding: utf8
# Author: Eric Fontanillas
# Modification: 03/09/14 by Julie BAFFARD
# Last modification : 10/09/21 by Charlotte Berthelier

"""
Description: Predict potential ORF on the basis of 2 criteria + 1 optional criteria

- CRITERIA 1

Get the longest part of the sequence alignemen without codon stop "*",
and test in the 3 potential ORF and check with a Blast the best coding sequence

- CRITERIA 2

This longest part should be > 150nc or 50aa

- CRITERIA 3 [OPTIONNAL]

A codon start "M" should be present in this longuest part, before the last 50 aa

OUTPUTs:
"05_CDS_aa" & "05_CDS_nuc" => NOT INCLUDE THIS CRITERIA
"06_CDS_with_M_aa" & "06_CDS_with_M_nuc" => INCLUDE THIS CRITERIA
"""

import os
import re
import argparse

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from dico import dico


def code_universel(file1):
    """ Creates bash for genetic code (key : codon ; value : amino-acid) """
    bash_code_universel = {}

    with open(file1, "r") as file:
        for line in file.readlines():
            item = str.split(line, " ")
            length1 = len(item)
            if length1 == 3:
                key = item[0]
                value = item[2][:-1]
                bash_code_universel[key] = value
            else:
                key = item[0]
                value = item[2]
                bash_code_universel[key] = value

    return bash_code_universel


def multiple3(seq):
    """
    Tests if the sequence is a multiple of 3, and if not removes extra-bases
    Possible to lost a codon, when I test ORF (as I will decay the ORF)
    """

    multiple = len(seq) % 3
    if multiple != 0:
        return seq[:-multiple], multiple
    return seq, multiple


def detect_methionine(seq_aa, ortho, minimal_cds_length):
    """ Detects if methionin in the aa sequence """

    size = len(seq_aa)
    cutoff_last_50aa = size - minimal_cds_length

    # Find all indices of occurances of "M" in a string of aa
    list_indices = [pos for pos, char in enumerate(seq_aa) if char == "M"]

    # If some "M" are present, find whether the first "M" found is not in the
    # 50 last aa (indice < CUTOFF_Last_50aa) ==> in this case: maybenot a CDS
    if list_indices != []:
        first_m = list_indices[0]
        if first_m < cutoff_last_50aa:
            ortho = 1  # means orthologs found

    return ortho


def reverse_complement2(seq):
    """ Reverse complement DNA sequence """
    seq1 = 'ATCGN-TAGCN-atcgn-tagcn-'
    seq_dict = {seq1[i]: seq1[i + 6]
                for i in range(24) if i < 6 or 12 <= i <= 16}

    return "".join([seq_dict[base] for base in reversed(seq)])


def simply_get_orf(seq_dna, gen_code):
    """ Generate the ORF sequence from DNA sequence """
    seq_by_codons = [seq_dna.upper().replace('T', 'U')[i:i + 3]
                     for i in range(0, len(seq_dna), 3)]
    seq_by_aa = [gen_code[codon] if codon in gen_code.keys()
                 else '?' for codon in seq_by_codons]

    return ''.join(seq_by_aa)


def find_good_orf_criteria3(bash_aligned_nc_seq, bash_code_universel, \
        minimal_cds_length):
    """
    Multiple sequence based : Based on the alignment of several sequences (orthogroup)
    Criteria 1 : Get the segment in the alignment with no codon stop
    """

    # 1 - Get the list of aligned aa seq for the 3 ORF:
    bash_of_aligned_aa_seq_3orf = {}
    bash_of_aligned_nuc_seq_3orf = {}
    best_longuest_subsequence_list_position = []

    for fasta_name in bash_aligned_nc_seq.keys():
        # Get sequence, chek if multiple 3, then get 6 orfs
        sequence_nc = bash_aligned_nc_seq[fasta_name]
        new_sequence_nc = multiple3(sequence_nc)
        new_sequence_rev = reverse_complement2(new_sequence_nc)
        # For each seq of the multialignment => give the 6 ORFs (in nuc)
        bash_of_aligned_nuc_seq_3orf[fasta_name] = [new_sequence_nc,
                                                    new_sequence_nc[1:-2],
                                                    new_sequence_nc[2:-1],
                                                    new_sequence_rev,
                                                    new_sequence_rev[1:-2],
                                                    new_sequence_rev[2:-1]]

        seq_prot_orf1 = simply_get_orf(new_sequence_nc, bash_code_universel)
        seq_prot_orf2 = simply_get_orf(
            new_sequence_nc[1:-2], bash_code_universel)
        seq_prot_orf3 = simply_get_orf(
            new_sequence_nc[2:-1], bash_code_universel)
        seq_prot_orf4 = simply_get_orf(new_sequence_rev, bash_code_universel)
        seq_prot_orf5 = simply_get_orf(
            new_sequence_rev[1:-2], bash_code_universel)
        seq_prot_orf6 = simply_get_orf(
            new_sequence_rev[2:-1], bash_code_universel)

        # For each seq of the multialignment => give the 6 ORFs (in aa)
        bash_of_aligned_aa_seq_3orf[fasta_name] = [
            seq_prot_orf1,
            seq_prot_orf2,
            seq_prot_orf3,
            seq_prot_orf4,
            seq_prot_orf5,
            seq_prot_orf6]

    # 2 - Test for the best ORF (Get the longuest segment in the alignment
    # with no codon stop ... for each ORF ... the longuest should give the
    # ORF)
    best_max = 0

    for i in [0, 1, 2, 3, 4, 5]:  # Test the 6 ORFs
        orf_aligned_aa = []
        orf_aligned_nuc = []

        # 2.1 - Get the alignment of sequence for a given ORF
        # Compare the 1rst ORF between all sequence => list them in
        # orf_aligned_aa // them do the same for the second ORF, and them the
        # 3rd
        for fasta_name in bash_of_aligned_aa_seq_3orf:
            orf_sequence = bash_of_aligned_aa_seq_3orf[fasta_name][i]
            aa_length = len(orf_sequence)
            # List of all sequences in the ORF nb "i" =
            orf_aligned_aa.append(orf_sequence)

        for fasta_name in bash_of_aligned_nuc_seq_3orf:
            orf_sequence = bash_of_aligned_nuc_seq_3orf[fasta_name][i]
            # List of all sequences in the ORF nb "i" =
            orf_aligned_nuc.append(orf_sequence)

        # 2.2 - Get the list of sublist of positions whithout codon stop in the alignment
        # For each ORF, now we have the list of sequences available
        # Next step is to get the longuest subsequence whithout stop
        # We will explore the presence of stop "*" in each column of the
        # alignment, and get the positions of the segments between the
        # positions with "*"
        j = 0  # Start from first position in alignment
        list_of_list_subsequences = []
        list_positions_subsequence = []
        while j < aa_length:
            column = []
            for seq in orf_aligned_aa:
                column.append(seq[j])
            j = j + 1
            if "*" in column:
                # Add previous list of positions
                list_of_list_subsequences.append(list_positions_subsequence)
                # Re-initialyse list of positions
                list_positions_subsequence = []
            else:
                list_positions_subsequence.append(j)

        # 2.3 - Among all the sublists (separated by column with codon stop
        # "*"), get the longuest one (BETTER SEGMENT for a given ORF)
        longuest_subsequence_list_position = []
        maxi = 0
        for sublist in list_of_list_subsequences:
            if maxi and minimal_cds_length < len(sublist):
                maxi = len(sublist)
                longuest_subsequence_list_position = sublist

        # 2.4. - Test if the longuest subsequence start exactly at the
        # beginning of the original sequence (i.e. means the ORF maybe
        # truncated)
        if longuest_subsequence_list_position != []:
            if longuest_subsequence_list_position[0] == 0:
                cds_maybe_truncated = 1
            else:
                cds_maybe_truncated = 0
        else:
            cds_maybe_truncated = 0

        # 2.5 - Test if this BETTER SEGMENT for a given ORF,
        # is the better than the one for the other ORF (GET THE BEST ORF)
        # Test whether it is the better ORF
        if maxi > best_max:
            best_max = maxi
            best_orf = i + 1
            best_longuest_subsequence_list_position = longuest_subsequence_list_position

    # 3 - ONCE we have this better segment (BEST CODING SEGMENT)
    # ==> GET THE STARTING and ENDING POSITIONS (in aa position and in nuc position)
    # And get the INDEX of the best ORF [0, 1, or 2]
    if best_longuest_subsequence_list_position != []:
        pos_min_aa = best_longuest_subsequence_list_position[0]
        pos_min_aa = pos_min_aa - 1
        pos_max_aa = best_longuest_subsequence_list_position[-1]
        best_orf_bash_of_aligned_aa_seq = {}
        best_orf_bash_of_aligned_aa_seq_coding = {}
        for fasta_name in bash_of_aligned_aa_seq_3orf:
            # cause list going from 0 to 2 in LIST_3_ORF, while the ORF nb is
            # indexed from 1 to 3
            index_best_orf = best_orf - 1
            seq = bash_of_aligned_aa_seq_3orf[fasta_name][index_best_orf]
            seq_coding = seq[pos_min_aa:pos_max_aa]
            best_orf_bash_of_aligned_aa_seq[fasta_name] = seq
            best_orf_bash_of_aligned_aa_seq_coding[fasta_name] = seq_coding

        # 4 - Get the corresponding position (START/END of BEST CODING SEGMENT)
        # for nucleotides alignment
        # And Blast the best coding sequence <=> longuest sequence to check if
        # it's the better one
        pos_min_nuc = pos_min_aa * 3
        pos_max_nuc = pos_max_aa * 3

        best_orf_bash_aligned_nc_seq = {}
        best_orf_bash_aligned_nc_seq_coding = {}
        for fasta_name in bash_aligned_nc_seq.keys():
            seq = bash_of_aligned_nuc_seq_3orf[fasta_name][index_best_orf]
            seq_coding = seq[pos_min_nuc:pos_max_nuc]
            best_orf_bash_aligned_nc_seq[fasta_name] = seq
            best_orf_bash_aligned_nc_seq_coding[fasta_name] = seq_coding
            seq_cutted = re.sub(r'^.*?[a-zA-Z]', '', seq)
            sequence_for_blast = (fasta_name + '\n' + seq_cutted + '\n')
            good_orf_found = False
            try:
                result_handle = NCBIWWW.qblast(
                    "blastn", "nt", sequence_for_blast, expect=0.001, hitlist_size=1)
                blast_records = NCBIXML.parse(result_handle)
            except BaseException:
                good_orf_found = False
            else:
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.expect < 0.001:
                                good_orf_found = True

            if good_orf_found is False:
                # Blast other possible CDS and keep best evalue, if no one keep
                # the longuest sequence
                for i in range(0, 5):
                    sequence = bash_of_aligned_nuc_seq_3orf[fasta_name][i]
                    if sequence == seq:
                        pass
                    else:
                        sequence = re.sub(r'^.*?[a-zA-Z]', '', sequence)
                        sequence_to_blast = (
                            fasta_name + '\n' + sequence + '\n')
                        try:
                            result_handle = NCBIWWW.qblast(
                                "blastn", "nt", sequence_to_blast, expect=0.001, hitlist_size=1)
                            blast_records = NCBIXML.parse(result_handle)
                        except BaseException:
                            pass
                        else:
                            for blast_record in blast_records:
                                for alignment in blast_record.alignments:
                                    for hsp in alignment.hsps:
                                        if hsp.expect < 0.001:
                                            # New good CDS found
                                            best_orf_bash_aligned_nc_seq[fasta_name] = sequence
                                            seq_coding = sequence[pos_min_nuc:pos_max_nuc]
                                            best_orf_bash_aligned_nc_seq_coding[fasta_name] = \
                                                seq_coding

            else:
                pass

    else:  # no CDS found
        best_orf_bash_aligned_nc_seq = {}
        best_orf_bash_aligned_nc_seq_coding = {}
        best_orf_bash_of_aligned_aa_seq = {}
        best_orf_bash_of_aligned_aa_seq_coding = {}

    # Check whether their is a "M" or not, and if at least 1 "M" is present,
    # that it is not in the last 50 aa

    best_orf_bash_of_aligned_aa_seq_cds_with_m = {}
    best_orf_bash_of_aligned_nuc_seq_cds_with_m = {}

    ortho = 0
    for fasta_name in best_orf_bash_of_aligned_aa_seq_coding:
        seq_aa = best_orf_bash_of_aligned_aa_seq_coding[fasta_name]
        ortho = detect_methionine(
            seq_aa, ortho, minimal_cds_length)  # DEF6 ###

    # CASE 1: A "M" is present and correctly localized (not in last 50 aa)
    if ortho == 1:
        best_orf_bash_of_aligned_aa_seq_cds_with_m = best_orf_bash_of_aligned_aa_seq_coding
        best_orf_bash_of_aligned_nuc_seq_cds_with_m = best_orf_bash_aligned_nc_seq_coding

    # CASE 2: in case the CDS is truncated, so the "M" is maybe missing:
    if ortho == 0 and cds_maybe_truncated == 1:
        best_orf_bash_of_aligned_aa_seq_cds_with_m = best_orf_bash_of_aligned_aa_seq_coding
        best_orf_bash_of_aligned_nuc_seq_cds_with_m = best_orf_bash_aligned_nc_seq_coding

    # CASE 3: CDS not truncated AND no "M" found in good position (i.e. before the last 50 aa):
        # => the 2 bash "CDS_with_M" are left empty ("{}")

    return best_orf_bash_aligned_nc_seq, best_orf_bash_aligned_nc_seq_coding, \
        best_orf_bash_of_aligned_nuc_seq_cds_with_m, best_orf_bash_of_aligned_aa_seq, \
        best_orf_bash_of_aligned_aa_seq_coding, best_orf_bash_of_aligned_aa_seq_cds_with_m, \
        best_longuest_subsequence_list_position


def write_output_file(results_dict, name_elems, path_out):
    """
    Output the best ORF in proteic and nucleic format
    """
    if results_dict != {}:
        name_elems[3] = str(len(results_dict.keys()))
        new_name = "_".join(name_elems)

        out1 = open("%s/%s" % (path_out, new_name), "w")
        for fasta_name in results_dict.keys():
            seq = results_dict[fasta_name]
            out1.write("%s\n" % fasta_name)
            out1.write("%s\n" % seq)
        out1.close()


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "codeUniversel",
        help="File describing the genetic code (code_universel_modified.txt")
    parser.add_argument(
        "min_cds_len",
        help="Minmal length of a CDS (in amino-acids)",
        type=int)
    parser.add_argument(
        "min_spec", help="Minimal number of species per alignment")
    parser.add_argument("list_files", help="File with all input files names")
    args = parser.parse_args()

    minimal_cds_length = int(args.min_cds_len)  # in aa number
    bash_code_universel = code_universel(args.codeUniversel)
    minimum_species = int(args.min_spec)

    # Inputs from file containing list of species
    list_files = []
    with open(args.list_files, 'r') as file:
        for line in file.readlines():
            list_files.append(line.strip('\n'))

    # Directories for results
    dirs = ["04_best_orf_nuc", "04_best_orf_aa", "05_CDS_nuc",
            "05_CDS_aa", "06_CDS_with_M_nuc", "06_CDS_with_M_aa"]
    for directory in dirs:
        os.mkdir(directory)

    count_file_processed, count_file_with_cds, count_file_without_cds, \
        count_file_with_cds_plus_m = 0, 0, 0, 0
    count_file_with_cds_and_enought_species, count_file_with_cds_m_and_enought_species = 0, 0

    # ! : Currently, files are named "orthogroup_x_y_sequences.fasta, where x is the number
    # of the orthogroup (not important, juste here to make a distinct name),
    # and y is the number of sequences/species in the group. These files are
    # outputs of blastalign, where species can be removed. y is then modified.
    name_elems = ["orthogroup", "0", "with", "0", "species.fasta"]

    # by fixing the counter here, there will be some "holes" in the outputs directories
    # (missing numbers), but the groups between directories will correspond
    # n0 = 0

    for file in list_files:
        # n0 += 1

        count_file_processed = count_file_processed + 1
        nb_gp = file.split('_')[1]  # Keep trace of the orthogroup number
        fasta_file_path = "./%s" % file
        bash_fasta = dico(fasta_file_path)

        best_orf_nuc, best_orf_nuc_coding, best_orf_nuc_cds_with_m, best_orf_aa, \
                best_orf_aa_coding, best_orf_aa_cds_with_m, longuest_subsequence_list_position \
                = find_good_orf_criteria3(bash_fasta, bash_code_universel, minimal_cds_length)
        name_elems[1] = nb_gp

        # Update counts and write group in corresponding output directory
        if best_orf_nuc != {}:
            count_file_with_cds += 1
            if len(best_orf_nuc.keys()) >= minimum_species:
                count_file_with_cds_and_enought_species += 1
                write_output_file(best_orf_nuc, name_elems,
                                  dirs[0])  # OUTPUT best_orf_nuc
                # The most interesting
                write_output_file(best_orf_aa, name_elems, dirs[1])
        else:
            count_file_without_cds += 1

        if best_orf_nuc_coding != {} and len(
                best_orf_nuc_coding.keys()) >= minimum_species:
            write_output_file(best_orf_nuc_coding, name_elems, dirs[2])
            write_output_file(best_orf_aa_coding, name_elems, dirs[3])

        if best_orf_nuc_cds_with_m != {}:
            count_file_with_cds_plus_m += 1
            if len(best_orf_nuc_cds_with_m.keys()) >= minimum_species:
                count_file_with_cds_m_and_enought_species += 1
                write_output_file(best_orf_nuc_cds_with_m, name_elems, dirs[4])
                write_output_file(best_orf_aa_cds_with_m, name_elems, dirs[5])

    print("*************** CDS detection ***************")
    print("\nFiles processed: %d" % count_file_processed)
    print("\tFiles with CDS: %d" % count_file_with_cds)
    print("\tFiles wth CDS and more than %s species: %d" %
          (minimum_species, count_file_with_cds_and_enought_species))
    print("\t\tFiles with CDS plus M (codon start): %d" %
          count_file_with_cds_plus_m)
    print(
        "\t\tFiles with CDS plus M (codon start) and more than %s species: %d" %
        (minimum_species, count_file_with_cds_m_and_enought_species))
    print("\tFiles without CDS: %d \n" % count_file_without_cds)
    print("")


if __name__ == '__main__':
    main()
