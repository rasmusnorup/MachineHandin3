import numpy as np
import itertools

class hmm:
    def __init__(self, init_probs, trans_probs, emission_probs):
        self.init_probs = init_probs
        self.trans_probs = trans_probs
        self.emission_probs = emission_probs
noStates = 7


init_probs_7_state = [0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00]
init_probs_codon_state = [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]

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

                #print("T1:" + str(T1[k][i-1]) +"  TransP:"+ str(transProbs[k][j]) +"   emProb:"+ str(emProbs[j][genome[i]]))

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

def codonCountTransmissionProbs(states):
    probs = make_table(noStates,68)

    for i in range(0,len(states)-1):
        probs[states[i]][states[i+1]] += 1
    for n in range(0,noStates):
        total = sum(probs[n])
        for m in range(0,noStates):
            probs[n][m] = probs[n][m]/total
    return probs

def codonCountEmissionProbs(genomeIndices, states):
    result = make_table(noStates,68)
    for i in range(0, len(states)):
        result[states[i]][genomeIndices[i]] += 1
    for n in range(0,noStates):
        total = sum(result[n])
        for m in range(0,68):
            result[n][m] = result[n][m]/total
    return result



def codonAnotationToStates(annotation):
    result = []
    i = 0
    isCoding = False
    isRCoding = False
    while i<len(annotation):
        if annotation[i] == "N":
            result.append(0)
            i+=1
        elif annotation[i] == "C":
            if isCoding == False:
                result.append(1)
                isCoding = True
            else:
                if annotation[i+3] != "C":
                    result.append(3)
                    isCoding = False
                else:
                    result.append(2)
            i += 3
        elif annotation[i] == "R":
            if isRCoding == False:
                result.append(6)
                isRCoding = True
            else:
                if annotation[i + 3] != "R":
                    result.append(4)
                    isRCoding = False
                else:
                    result.append(5)
            i += 3
    return result

def codonGenomeToIndices(genome,states):
    result = []
    perms = [''.join(i) for i in itertools.product("ACGT", repeat=3)]
    singleMap = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    mapping = dict(zip(perms, range(4, 68)))
    j = 0
    for i in range(0, len(states)):
        if states[i] == 0:
            result.append(singleMap[genome[j]])
            j += 1
        else:
            result.append(mapping[genome[j] + genome[j + 1] + genome[j + 2]])
            j += 3
    return result

# state 0=noncoding, 1=start, 2=code, 3=stop, 4=Rstart, 5=Rcode, 6=Rstop
# Remember that Rstop comes before Rstart



genome1 = read_fasta_file("genome1.fa")

#genome2 = translate_observations_to_indices(read_fasta_file("genome2.fa"))
trueann1 = read_fasta_file("true-ann1.fa")
truestates1 = codonAnotationToStates(trueann1)
genomeIndices1 = codonGenomeToIndices(genome1,truestates1)

#print(codonAnotationToStates(trueann1))
emProbs = codonCountEmissionProbs(genomeIndices1,truestates1)
transProbs = codonCountTransmissionProbs(truestates1)
print(viterbi(init_probs_codon_state,transProbs,emProbs,genomeIndices1))
#trueann2 = read_fasta_file("true-ann2.fa")
#ann1 = convertAnnToState(read_fasta_file("true-ann1.fa"))
#emProbs = countEmissionProbs(genome1,ann1)
#transProbs = countTransisionProbs(ann1)
#result = viterbi(init_probs_7_state,transProbs,emProbs,genome2)
#learnann=translate_indices_to_path(result)
#print(compute_accuracy(trueann2,learnann))
#print(result)
#print(translate_indices_to_path(result))