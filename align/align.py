# Importing Dependencies
import numpy as np
from typing import Tuple
from collections import OrderedDict

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = OrderedDict()  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        seqA=seqA
        seqB=seqB
        working_matrix=np.zeros((len(seqA)+3,len(seqB)+3))
        storage_matrix={}
        working_matrix[0,0]=0
        working_matrix[0,:]=self.gap_open
        working_matrix[:,0]=self.gap_open
        print (working_matrix)
        #may need to correct if open gap error should be propogated to beginning states
        countx=0
        county=0
        Quit=0
        storage_matrix[(0,0)]=["Stop","nogap"]
        for i in range(len(seqB)):
            storage_matrix[(0,i+1)]=[(0,i),"gap"]
            i+=1
        for i in range(len(seqA)):
            storage_matrix[(i+1,0)]=[(i,0),"gap"]
            i+=1
        while Quit!=1:
            up = working_matrix[(countx + 1, county)]
            key_to_check=(countx+1,county)
            if key_to_check in storage_matrix.keys():
                if storage_matrix[(countx+1,county)][1]=="gap":
                    up+=self.gap_open
                elif storage_matrix[(countx+1,county)][1]=="nogap":
                    up+=self.gap_extend
            key_to_check=(countx,county+1)
            if key_to_check in storage_matrix.keys():
                left=working_matrix[(countx,county+1)]
                if storage_matrix[(countx,county+1)][1]=="gap":
                    left+=self.gap_open
                elif storage_matrix[(countx,county+1)][1]=="nogap":
                    left+=self.gap_extend
            if countx<len(seqA):
                valuex=countx
            else:
                valuex=len(seqA)-1
            if county<len(seqB):
                valuey=county
            else:
                valuey=len(seqB)-1
            dictionary_return=self.sub_dict[(seqA[valuex],seqB[valuey])]
            print(dictionary_return)
            if countx<len(seqA) and county<len(seqB):
                diagonal=working_matrix[(countx,county)]+dictionary_return
            else: 
                diagonal = -10000000
            if county<len(seqB):
                up=up
            else: 
                up = -10000000
            if countx<len(seqA):
                left=left
            else: 
                left = -10000000
            direction=np.max([diagonal, up, left])
            print(direction)
            
            if direction == diagonal:
                working_matrix[(countx+1, county+1)]=diagonal
                storage_matrix[(countx+1, county+1)]=[(countx,county),"nogap","diagonal"]
                
            elif direction ==left:
                working_matrix[(countx+1, county+1)]=left
                storage_matrix[(countx+1, county+1)]=[(countx,county+1),"gap","left"]
                
            elif direction ==up:
                working_matrix[(countx+1, county+1)]=up
                storage_matrix[(countx+1, county+1)]=[(countx+1,county),"gap","up"]
                
            if direction == diagonal:
                countx+=1
                county+=1
            elif direction ==left:
                countx+=1
            elif direction ==up:
                county+=1
            print (countx)
            print (county)
            if countx>len(seqA):
                if county>len(seqB):
                    Quit=1 
        self.storage_matrix=storage_matrix
        self.working_matrix=working_matrix
        

           

       
       
       
        """
        TODO
        build M matrix
            define initial first row parameters
            define initial first column parameters
            dict_sub=_read_sub_matrix
            define penalties for alignment_score, gap_open, gap_extend
        at each position, look to up, diagonal, and left
            compare max values for up+0, left+0, diagonal+penalty
                calculate values for up+0, left+0

            
            
            
            add max value in matrix
                if multiple are equal, set equal to random of the equal values
            add position pulled from to a list for traceback once matrix is completed
        Read backtrace
            if neither i nor j is the same, align s1[i] and s2[j]
            if j is the same, align gap in s2[j]
            if i is the same, align gap in s1[i]
            
            
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        pass

        
        # TODO: Implement global alignment here
        pass      		
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        storage_matrix = self.storage_matrix
        working_matrix = self.working_matrix
        end_position = (len(self._seqA), len(self._seqB))
        
        path_taken = []
        seqA_alignment = []
        seqB_alignment = []
        alignment_score = 0
        
        current_position = end_position
        while True:
            while current_position != (0, 0):
                key_to_check=current_position
                if key_to_check not in storage_matrix.keys():
                    current_position=(current_position[0]-1,current_position[1])
                    key_to_check=current_position
                    if key_to_check not in storage_matrix.keys():
                        current_position=(current_position[0]+1,current_position[1]-1)
                        key_to_check=current_position
                        if key_to_check not in storage_matrix.keys():
                            current_position=(current_position[0]-2,current_position[1]+1)
                            key_to_check=current_position
                            if key_to_check not in storage_matrix.keys():
                                current_position=(current_position[0]+2,current_position[1]-2)

                if len(current_position)>2:
                    current_position=current_position
                    print(current_position)
                prev_position, alignment_cue= storage_matrix[current_position]
                
                if alignment_cue == "left":
                    seqB_alignment.append("_")
                elif alignment_cue == "up":
                    seqA_alignment.append("_")
                elif alignment_cue == "diagonal":
                    seqA_alignment.append(self._seqA[current_position[0] - 1])
                    seqB_alignment.append(self._seqB[current_position[1] - 1])

                alignment_score += working_matrix[current_position]
                path_taken.append(current_position)
                current_position = prev_position

            seqA_alignment = ''.join(seqA_alignment[::-1])
            seqB_alignment = ''.join(seqB_alignment[::-1])

            return (alignment_score, seqA_alignment, seqB_alignment, storage_matrix, working_matrix)

    
    
    
    """def _backtrace(self) -> Tuple[float, str, str]:
        storage_matrix=self.storage_matrix
        working_matrix=self.working_matrix
        end_value=working_matrix[len(self._seqA),len(self._seqB)]
        #start_value=storage_matrix[(1,1)][1]
        #path_taken=[start_value]
        path_taken=[]
        seqA_alignment=[]
        seqB_alignment=[]
        alignment_score=0
        #key_to_check=(len(self._seqA),len(self._seqB))
        #if key_to_check in storage_matrix.keys():
        storage_matrix_list = list(storage_matrix.items())
        key_to_check = storage_matrix_list[-1]
        print(key_to_check)
        print(storage_matrix)
        key_value=key_to_check[1]
        next_direction = key_value[0]
        
        alignment_cue = key_value[2]


        if alignment_cue=="left":
            seqB_alignment+="_"
        if alignment_cue=="up":
            seqA_alignment+="_"
        if alignment_cue=="diagonal":
            seqA_alignment=self._seqA[next_direction[0]-1]
            seqB_alignment=self._seqB[next_direction[1]-1]
        if next_direction != type(None):
            path_taken+=next_direction
            alignment_score+=working_matrix[next_direction]

            
        while next_direction != "Stop":
            next_direction=storage_matrix[next_direction][0]
            if alignment_cue=="left":
                seqB_alignment+="_"
            elif alignment_cue=="up":
                seqA_alignment+="_"
            elif alignment_cue=="diagonal":
                seqA_alignment=self._seqA[next_direction[0]-1]
                seqB_alignment=self._seqB[next_direction[1]-1]
            path_taken+=next_direction
            alignment_score+=working_matrix[next_direction]
        seqA_alignment=seqA_alignment[::-1]
        seqB_alignment=seqB_alignment[::-1]"""
        
        
       
        
    """TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB"""
        
 



def read_fasta(fasta_file: str) -> Tuple[str, str]:
    
    """DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header




#check indexing on array function