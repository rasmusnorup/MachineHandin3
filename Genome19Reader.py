import numpy as np


noStates = 19


init_probs_19_state = [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]

#States: 0:N 1-2-3:Start 4-5-6:Code 7-8-9:Stop 10-11-12:Rstart 13-14-15:RCode 16-17-18:Rstop

def translate_observations_to_indices(obs):
    mapping = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    return [mapping[symbol.lower()] for symbol in obs]


def convertAnnToState(ann):
    state = []
    isCoding = False
    isRCoding = False
    i=0
    while i<len(ann):
        if ann[i] == 'N':
            state.append(0)
            i += 1
        elif ann[i] == 'C':
            if isCoding == False:
                isCoding = True
                state.extend([1,2,3])
            elif ann[i+3] != "C":
                isCoding = False
                state.extend([7,8,9])
            else:
                state.extend([4,5,6])
            i += 3
        elif ann[i] == 'R':
            if isRCoding == False:
                isRCoding = True
                state.extend([16, 17, 18])
            elif ann[i + 3] != "R":
                isRCoding = False
                state.extend([10, 11, 12])
            else:
                state.extend([13, 14, 15])
            i += 3
    return state

def translate_indices_to_path(indices):
    mapping = ['C', 'C', 'C', 'N', 'R', 'R', 'R']
    return ''.join([mapping[i] for i in indices])

def read_fasta_file(filename):
    """
    Reads the given FASTA file f and returns a dictionary of sequences.

    Lines starting with ';' in the FASTA file are ignored.
    """
    sequences_lines = {}
    current_sequence_lines = None
    with open(filename) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            if line.startswith('>'):
                sequence_name = line.lstrip('>')
                current_sequence_lines = []
                sequences_lines[sequence_name] = current_sequence_lines
            else:
                if current_sequence_lines is not None:
                    current_sequence_lines.append(line)
    sequences = {}
    for name, lines in sequences_lines.items():
        sequences[name] = ''.join(lines)
    return sequences[list(sequences.keys())[0]]

genome1 = read_fasta_file("genome1.fa")
trueann1 = read_fasta_file("true-ann1.fa")
print(convertAnnToState(trueann1))
