import argparse

parser = argparse.ArgumentParser(description=\
        'removes chr from chromosomes files')
parser.add_argument('--f1', help='input coverage file')
parser.add_argument('--f2', help='output coverage file')
args = parser.parse_args()

f = open(args.f1, 'r')
f2 = open(args.f2, 'w')
real_ids = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrY']
for line in f:
        if line.split('\t')[0] in real_ids:
                line = line.split('\t')
                line[0] = line[0].replace('chr', '')
                f2.write('\t'.join(line))
f.close()
f2.close()