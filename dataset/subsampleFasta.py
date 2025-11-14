'''
python subsampleFasta.py [fasta file] [train tab]
'''

from sys import argv

def parseTab(fNm):
	f = open(fNm)

	fids = {}
	for i in f:
		fid = i.strip('\n').split('\t')[0]
		fids[fid] = 0

	f.close()

	return fids

def parseFasta(fNm, fids):
	f = open(fNm)

	currSeq = ''
	currGen = ''
	for i in f:
		if i[0] == '>':
			if currSeq in fids:
				print '>' + currSeq
				print currSeq
			currSeq = ''
			currGen = i.strip().split('>')[1]
			continue
		currSeq += i.strip()

	f.close()

def main():
	fids = parseTab(argv[2])
	parseFasta(argv[1], fids)

if __name__ == '__main__':
	main()
