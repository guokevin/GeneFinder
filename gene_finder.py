# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

KEVIN GUO
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        Added additional tests for all possible cases
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    else:
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        Tests given are enough since there are no additional cases or extremes to check for
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CGAATGTAGGGTAAAAACT")
    'AGTTTTTACCCTACATTCG'
    """
    complement = ""
    for i in range(-1,-len(dna)-1,-1):
        complement += get_complement(dna[i])

    return complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        Additional tests ensure that the function checks for validity of assumed input and for cases when there is no in frame stop codon
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATAGATAGG")
    DNA sequence does not begin with a start codon
    >>> rest_of_ORF("ATGAGAAGG")
    'ATGAGAAGG'

    Start codon - ATG
    First in Frame Stop Codon - TAG, TAA, TGA
    """

    threes = split_list(dna)
    # print 'a'
    if(threes[0] != 'ATG'):
        print("DNA sequence does not begin with a start codon")
        return
    # print 'b'
    frame = earliest_frame_stop(threes)
    if(frame != -1):
        return dna[:(3*(frame))]
    else:
        return dna;
    # print 'c'


def earliest_frame_stop(threes):
    for i in range(0, len(threes)):
        if(threes[i] == 'TGA' or threes[i] == 'TAA' or threes[i] == 'TAG'):
            return i
    return -1

def split_list(dna):
    threes = []
    for i in range(0,len(dna)/3):
        templist = dna[:3]
        threes.append(templist)
        dna = dna[3:]
    # print dna
    # if(len(dna)%3 != 0):
    #     threes.append(dna[3:])
    return threes

def open_frame(threes):
    for i in range(0, len(threes)):
        if(threes[i] == 'ATG'):
            return i
    return -1

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        Added code to see what would happen if inserted program had no ATG to begin with, had no ATG at all, and increased the vigor of the test
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("CAT")
    []
    >>> find_all_ORFs_oneframe("CATGAATGTAGATAGATGTGCCCC")
    ['ATGTGCCCC']
    >>> find_all_ORFs_oneframe("ATGCATAGATAGCCCATGTGCCCATGACCAATGCCC")
    ['ATGCATAGA', 'ATGTGCCCA', 'ATGCCC']
    """
    threes = split_list(dna)
    orfs = []

    while(len(threes) > 0):
        # print len(threes)
        if(open_frame(threes) == -1):
            return orfs
        start_frame = open_frame(threes)
        threes = threes[start_frame:]
        dna = dna[3*start_frame:]
        orfs.append(rest_of_ORF(dna))
        end_frame = earliest_frame_stop(threes)

        # print'd'

        if (end_frame == -1):
            return orfs

        threes = threes[end_frame:]
        dna = dna[3*end_frame:]

        start_frame = open_frame(threes)
        threes = threes[start_frame:]
        dna = dna[3*start_frame:]

        

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        Tests to ensure that it finds no false orfs
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("GCATAATTAG")
    []
    """
    return find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        Ensures it finds the other strands' orfs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("AGTGTAAAACAT")
    ['ATGTTTTACACT']
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

        Ensure it find longest ORF of both strands

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("AGTGTAAAACAT")
    'ATGTTTTACACT'
    """
    maximum = ''
    for orf in find_all_ORFs_both_strands(dna):
        # print orf
        if len(orf) > len(maximum):
            maximum = orf

    return maximum
    # return max(find_all_ORFs_both_strands(dna), key = len)

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
    """
    maximum = ''
    for i in range(0,num_trials):
        dna = shuffle_string(dna)
        # print "dna:",dna
        # print i
        temp = longest_ORF(dna)
        if len(temp) > len(maximum):
            maximum = temp
        # print "maximum:",maximum
    # print len(maximum)
    return len(maximum)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        Added case for starting with not ATG to see result just to ensure it doesn't check for codons

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("CCCGCTTT")
        'PA'
    """
    dna_list = split_list(dna)

    aa = ""
    for i in range(0, len(dna_list)):
        aa += aa_table[dna_list[i]]

    return aa


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1000)
    print threshold
    #threshold = 100

    orfs = find_all_ORFs_both_strands(dna)

    return_orfs = []

    for i in range(0, len(orfs)):
        if (len(orfs[i]) > threshold):
            return_orfs.append(orfs[i])

    # return return_orfs
    for i in range(0, len(return_orfs)):
        return_orfs[i] = coding_strand_to_AA(return_orfs[i])
    
    return return_orfs

if __name__ == "__main__":
    from load import load_seq
    # dna = load_seq("./data/X73525.fa")
    # print dna
    # print(gene_finder(dna))
    # gene_finder(dna)
    import doctest
    doctest.testmod()
    doctest.run_docstring_examples(find_all_ORFs_both_strands,globals(),verbose=True)

    #print(longest_ORF_noncoding('ATGTAG',30))