import re

STOP_M = ['TAA', 'TAG', 'AGA', "AGG"]
START_M = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']

START_S = ["ATG", "TTG", "CTG"]
STOP_S = ["TAA", "TAG", "TGA"]

gencode_s = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}


gencode_m = {'ATA': 'M', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': '*', 'AGG': '*',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': 'W', 'TGG': 'W'}


def find_first_in_register_stop(seq, STOP):
    i = 0
    while i < len(seq) - 2 and seq[i:i+3] not in STOP:
        i += 3
    if i < len(seq) - 2:
        return i + 3
    else:
        return -1


def translate(seq, START, STOP, genome):
    i = 0
    if genome == "standart":
        while i < len(seq) + 2 and seq[i:i + 3] not in START:
            i += 1
        prot = []
        while i < len(seq) - 2 and seq[i:i + 3] not in STOP:
            prot.append(gencode_s[seq[i:i + 3]])
            i += 3
        return ''.join(prot)
    else:
        while i < len(seq) + 2 and seq[i:i + 3] not in START:
            i += 1
        prot = []
        while i < len(seq) - 2 and seq[i:i + 3] not in STOP:
            prot.append(gencode_m[seq[i:i + 3]])
            i += 3
        return ''.join(prot)


def find_all_starts(seq, START):
    starts = []
    for i in range(len(START)):
        i = seq.find(START[i])
        while i >= 0:
            starts.append(i)
            i = seq.find(START, i + 1)
        return tuple(starts)


def orf_finder(seq, starts, stops):
    for i in range(len(seq) - 2):
        if seq[i:i + 3] in starts:
            start = i
            stop = find_first_in_register_stop(seq[start + 3:], stops)
            if stop > 0:
                return seq[start:stop + start], stop + start + 3
    return -1, -1


def find_all_orfs(seq, starts, stops):
    orfs = []
    i = 0
    while i < len(seq):
        orf = orf_finder(seq[i:], starts, stops)
        if orf[0] == -1:
            return orfs
        else:
            orfs.append(orf[0])
            i += orf[1]
    return orfs


def read_data(file):
    return file.readlines()[1].replace("\n", "").replace("\r", "")


def write_data(file, seq_n, seq_a, count, bool):
    if bool == True:
        file.write(">ORF_" + str(count) + "_nucl\n" + seq_n + "\n")
        file.write(">ORF_" + str(count) + "_prot\n" + seq_a + "\n")
    else:
        file.write(None)


def main():
    input_name = "input.fa"
    genome = "mitochondrial"

    input = open(input_name, 'r')
    output = open('output.fa', 'w')
    seq = read_data(input)

    if genome == "standard":
        count = 1
        orfs = find_all_orfs(seq, START_S, STOP_S)
        for i in orfs:
            write_data(output, i, translate(i, START_S, STOP_S, genome), count, True)
            count += 1
    elif genome == "mitochondrial":
        orfs = find_all_orfs(seq, START_M, STOP_M)
        for i in orfs:
            write_data(output, i, translate(i, START_M, STOP_M, genome), 1, True)
    else:
        write_data("None", 0, '', 0, False)


if __name__ == '__main__':
    main()
