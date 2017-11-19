import numpy as np

class hmm:
    def __init__(self, init_probs, trans_probs, emission_probs):
        self.init_probs = init_probs
        self.trans_probs = trans_probs
        self.emission_probs = emission_probs
noStates = 7


init_probs_7_state = [0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00]

trans_probs_7_state = [
    [0.00, 0.00, 0.90, 0.10, 0.00, 0.00, 0.00],
    [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00],
    [0.00, 0.00, 0.05, 0.90, 0.05, 0.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00],
    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00],
    [0.00, 0.00, 0.00, 0.10, 0.90, 0.00, 0.00],
]

emission_probs_7_state = [
    #   A     C     G     T
    [0.30, 0.25, 0.25, 0.20],
    [0.20, 0.35, 0.15, 0.30],
    [0.40, 0.15, 0.20, 0.25],
    [0.25, 0.25, 0.25, 0.25],
    [0.20, 0.40, 0.30, 0.10],
    [0.30, 0.20, 0.30, 0.20],
    [0.15, 0.30, 0.20, 0.35],
]

hmm_7_state = hmm(init_probs_7_state, trans_probs_7_state, emission_probs_7_state)



def translate_observations_to_indices(obs):
    mapping = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    return [mapping[symbol.lower()] for symbol in obs]


def convertAnnToState(ann):
    state = []
    cCount = 2
    rCount = 4
    for i in range(0,len(ann)):
        if ann[i] == 'N':
            state.append(3)
        if ann[i] == 'C':
            state.append(cCount)
            cCount -= 1
            if cCount==-1:
                cCount = 2
        if ann[i] == 'R':
            state.append(rCount)
            rCount += 1
            if rCount == 7:
                rCount = 4
    return state


def translate_indices_to_path(indices):
    mapping = ['C', 'C', 'C', 'N', 'R', 'R', 'R']
    return ''.join([mapping[i] for i in indices])


def translate_indices_to_observations(indices):
    mapping = ['a', 'c', 'g', 't']
    return ''.join(mapping[idx] for idx in indices)

def translate_path_to_indices(path):
    return list(map(lambda x: int(x), path))


def make_table(m, n):
    """Make a table with `m` rows and `n` columns filled with zeros."""
    return [[0] * n for _ in range(m)]

def compute_accuracy(true_ann, pred_ann):
    if len(true_ann) != len(pred_ann):
        return 0.0
    return sum(1 if true_ann[i] == pred_ann[i] else 0
               for i in range(len(true_ann))) / len(true_ann)

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

def countEmissionProbs(genome, annotation):
    probs = make_table(noStates,4)
    for i in range(0,len(genome)):

        probs[annotation[i]][genome[i]] += 1
    for n in range(0,noStates):
        total = sum(probs[n])
        for m in range(0,4):
            probs[n][m] = probs[n][m]/total

    return probs

def countTransisionProbs(annotation):
    probs = make_table(noStates,noStates)
    for i in range(0,len(annotation)-1):
        probs[annotation[i]][annotation[i+1]] += 1
    for n in range(0,noStates):
        total = sum(probs[n])
        for m in range(0,noStates):
            probs[n][m] = probs[n][m]/total

    return probs

def viterbi(initialProb, transProbs, emProbs, genome):
    T = len(genome)
    T1 = make_table(noStates, T)
    T2 = make_table(noStates, T)

    for i in range(0,noStates):
        T1[i][0] = emProbs[i][genome[0]]*initialProb[i]
        T2[i][0] = 0

    p = 0
    for i in range(1,T):
        p+=1
        if p==10000:
            p=0
            print(str(int(i/T*100)) + "%")
        for j in range(0,noStates):
            best = 0
            bestk = 0
            for k in range(0,noStates):
                temp = T1[k][i-1]*transProbs[k][j]*emProbs[j][genome[i]]

                #print("T1:" + str(T1[k][i-1]) +"  TransP:"+ str(transProbs[k][j]) +"   emProb:"+ str(emProbs[j][genome[i]]))

                if temp>best:
                    best = temp
                    bestk = k
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



genome1 = translate_observations_to_indices(read_fasta_file("genome1.fa")+read_fasta_file("genome2.fa"))
ann1 = convertAnnToState(read_fasta_file("true-ann1.fa")+read_fasta_file("true-ann2.fa"))
emProbs = countEmissionProbs(genome1,ann1)
transProbs = countTransisionProbs(ann1)
result = viterbi(init_probs_7_state,transProbs,emProbs,genome1[0:100000])
print(result)
