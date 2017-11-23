def saveFasta(file,annname, annotation):
    file.write( ">" + annname + "\n")
    i= 0
    string = ""
    while i < len(annotation):
        for j in range(0,60):
            if i+j < len(annotation):
                string = string + annotation[i+j]
        file.write(string + "\n")
        i = i+60
        string = ""


def saveConcFasta(filename, filenames, annotations):
    file = open(filename + ".fa", 'w')
    for i in range(0,len(filenames)):
        saveFasta(file, filenames[i], annotations[i])
    file.close()


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

pred = [ read_fasta_file("pred-ann6.fa"), read_fasta_file("pred-ann7.fa"), read_fasta_file("pred-ann8.fa"), read_fasta_file("pred-ann9.fa"), read_fasta_file("pred-ann10.fa")]
filenames = ["pred-ann6", "pred-ann7", "pred-ann8", "pred-ann9", "pred-ann10"]
saveConcFasta("pred-all", filenames, pred)