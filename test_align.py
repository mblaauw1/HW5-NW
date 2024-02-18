# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw_instance = NeedlemanWunsch("../HW5-NW/substitution_matrices/BLOSUM62.mat", float(-10), float(-1))
    results = nw_instance.align(seq1, seq2)
    print(results)
    if seq1[0]==None:
        return "Alignment sequences cannot be blank"
    if seq2[0]==None:
        return "Alignment sequences cannot be blank"
    test_list = results[4]
 
    # printing original list
    print("The original list : " + str(test_list))
 
    # using any() + list comprehension
    # to Search in Matrix
    res = any(None in sub for sub in test_list)
 
    # printing result
    print("Does Matrix contain None value ? : " + str(res))
    assert res== False




def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw_instance = NeedlemanWunsch("../HW5-NW/substitution_matrices/BLOSUM62.mat", float(-10), float(-1))
    results = nw_instance.align(seq3, seq4)
    if seq1[0]==None:
        return "Alignment sequences cannot be blank"
    if seq2[0]==None:
        return "Alignment sequences cannot be blank"
    pass
    assert np.sum(results[5]) < results[0]



