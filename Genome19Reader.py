import numpy as np


noStates = 19


init_probs_19_state = [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]

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

def countEmissionProbs(genome, annotation):
    genome = translate_observations_to_indices(genome)
    annotation = convertAnnToState(annotation)
    probs = make_table(noStates,4)
    for i in range(0,len(genome)):

        probs[annotation[i]][genome[i]] += 1
    for n in range(0,noStates):
        total = sum(probs[n])
        for m in range(0,4):
            probs[n][m] = probs[n][m]/total
    return probs

def countTransisionProbs(annotation):
    annotation = convertAnnToState(annotation)
    probs = make_table(noStates,noStates)
    for i in range(0,len(annotation)-1):
        probs[annotation[i]][annotation[i+1]] += 1
    for n in range(0,noStates):
        total = sum(probs[n])
        for m in range(0,noStates):
            probs[n][m] = probs[n][m]/total
    return probs

def viterbi(initialProb, transProbs, emProbs, genome):
    genome = translate_observations_to_indices(genome)
    T = len(genome)
    T1 = make_table(noStates, T)
    T2 = make_table(noStates, T)

    for i in range(0,noStates):
        if (emProbs[i][genome[0]] != 0 and initialProb[i] != 0):
            T1[i][0] = np.log(emProbs[i][genome[0]])+np.log(initialProb[i])
        T2[i][0] = 0

    p = 0
    for i in range(1,T):
        p+=1
        if p==10000:
            p=0
            print(str(int(i/T*100)) + "%")
        for j in range(0,noStates):
            best = -np.inf
            bestk = 0

            for k in range(0,noStates):
                if transProbs[k][j] != 0:
                    temp = T1[k][i-1]+np.log(transProbs[k][j])
                    if temp>best:
                        best = temp
                        bestk = k
            if emProbs[j][genome[i]] != 0:
                T1[j][i] = best+np.log(emProbs[j][genome[i]])
            else:
                T1[j][i] = best
            T2[j][i] = bestk
    z=[0]*T
    x=[0]*T
    best = 0
    bestk = 0
    for k in range(0, noStates):
        temp = T1[k][T-1]
        if temp>best:
            best = temp
            bestk = k
    x[T-1] = bestk
    for j in range(1,T+1):
        i=T-j
        x[i-1] = T2[x[i]][i]
    return x

def translate_indices_to_path(indices):
    mapping = ['N','C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R']
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
def compute_accuracy(true_ann, pred_ann):
    if len(true_ann) != len(pred_ann):
        return -1.0
    return sum(1 if true_ann[i] == pred_ann[i] else 0
               for i in range(len(true_ann))) / len(true_ann)
def make_table(m, n):
    """Make a table with `m` rows and `n` columns filled with zeros."""
    return [[0] * n for _ in range(m)]
def translate_observations_to_indices(obs):
    mapping = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    return [mapping[symbol.lower()] for symbol in obs]

def saveFasta(filename, annotation):
    file = open(filename,'w')
    file.write( ">" + filename + "\n")
    for i in range(0,len(annotation),60):
        


genome1 = read_fasta_file("genome1.fa")
genome2 = read_fasta_file("genome2.fa")
genome3 = read_fasta_file("genome3.fa")
genome4 = read_fasta_file("genome4.fa")
genome5 = read_fasta_file("genome5.fa")
genomes = genome1 + genome2 + genome3 + genome4 + genome5


trueann1 = read_fasta_file("true-ann1.fa")
trueann2 = read_fasta_file("true-ann2.fa")
trueann3 = read_fasta_file("true-ann3.fa")
trueann4 = read_fasta_file("true-ann4.fa")
trueann5 = read_fasta_file("true-ann5.fa")
trueanns = trueann1 + trueann2 + trueann3 + trueann4 + trueann5

emProbs = countEmissionProbs(genomes,trueanns)
transProbs = countTransisionProbs(trueanns)

result = viterbi(init_probs_19_state, transProbs, emProbs, genome5)
result = translate_indices_to_path(result)
print(result)
print(compute_accuracy(trueann5,result))
