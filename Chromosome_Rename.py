import argparse

parser = argparse.ArgumentParser(description=\
        'fixes dumb coverage files')
parser.add_argument('--f1', help='input coverage file')
parser.add_argument('--f2', help='output coverage file')
args = parser.parse_args()

f = open(args.f1, 'r')
f2 = open(args.f2, 'w')
real_ids = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
for line in f:
        if line.split('\t')[0] in real_ids:
                line = line.split('\t')
                line[0] = 'chr'+line[0]
                f2.write('\t'.join(line))
f.close()
f2.close()