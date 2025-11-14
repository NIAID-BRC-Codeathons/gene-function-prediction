'''
python parseGenomes.py [files list] [tabular out pref]
'''

from sys import argv,stderr
import os
from glob import glob

uHsh = {
	'hypothetical': 0,
	'uncharacterized gene': 0,
	'uncharacterized protein': 0,
	'unknown gene': 0,
	'unknown protein': 0
}

def err(s):
	stderr.write(s)

def getFLst(fNm):
	f = open(fNm)

	fLst = []
	for i in f:
		i = i.strip('\n')
		fLst.append(i)

	f.close()

	return fLst

def parseHeader(f):
	hHsh = {}

	ln = f.readline().strip('\n').split('\t')
	for i in range(0,len(ln)):
		hHsh[ln[i]] = i

	return hHsh

def parseFeatureTab(fNm):
	f = open(fNm)
	fidLabHsh = {}

	hHsh = parseHeader(f)
	flg = False
	for i in f:
		i = i.strip('\n').split('\t')

		fid = i[hHsh['patric_id']]
		prd = i[hHsh['product']]
		plf = i[hHsh['plfam_id']]
		pgf = i[hHsh['pgfam_id']]

		fidLabHsh[fid] = [plf, pgf, prd]

	f.close()

	return fidLabHsh

def addSeq(fidSeqHsh, fid, seq):
	if fid != "" and seq != "":
		fidSeqHsh[fid] = seq

def parseFasta(fNm):
	f = open(fNm)

	fidSeqHsh = {}

	fid = ''
	seq = ''
	for i in f:
		if i[0] == '>':
			addSeq(fidSeqHsh, fid, seq)
			fid = '|'.join(i.strip().split()[0].split('>')[1].split('|')[:2])
			seq = ''
			continue
		seq += i.strip()

	addSeq(fidSeqHsh, fid, seq)

	f.close()

	return fidSeqHsh

def filtLabs(fidLabHsh):
	dLst = []
	for i in fidLabHsh:
		prd = fidLabHsh[i][-1]
		prd = prd.lower()

		for j in uHsh:
			if j in prd:
				dLst.append(i)
				break

	for i in dLst:
		del fidLabHsh[i]

def mergeFidLabHsh(fxxLabhsh, fidLabHsh, fxxSeqHsh):
	for i in fxxSeqHsh:
		if i not in fidLabHsh:
			continue

		seq = fxxSeqHsh[i]
		lab = ','.join(fidLabHsh[i])

		if seq not in fxxLabhsh:
			fxxLabhsh[seq] = {}
		if lab not in fxxLabhsh[seq]:
			fxxLabhsh[seq][lab] = 0
		fxxLabhsh[seq][lab] += 1

def parseGenomes(fLst):
	ffnLabHsh = {}
	faaLabHsh = {}
	fidFfnHsh = {}
	fidFaaHsh = {}
	cnt = 0
	inc = len(fLst) / 50.
	err("Parsing features and fasta files...\n\t")
	for i in fLst:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		gid = os.path.basename(i[:-1])
		gFLst = glob(i+'*')
		# print gid
		# for j in gFLst:
		# 	print j

		fTbFNm = i + gid + '.PATRIC.features.tab'
		ffnFNm = i + gid + '.PATRIC.ffn'
		faaFNm = i + gid + '.PATRIC.faa'

		fidLabHsh = parseFeatureTab(fTbFNm)
		ffnSeqHsh = parseFasta(ffnFNm)
		faaSeqHsh = parseFasta(faaFNm)

		# filtLabs(fidLabHsh)

		mergeFidLabHsh(ffnLabHsh, fidLabHsh, ffnSeqHsh)
		mergeFidLabHsh(faaLabHsh, fidLabHsh, faaSeqHsh)
		fidFfnHsh.update(ffnSeqHsh)
		fidFaaHsh.update(faaSeqHsh)

		# break
	err('\n')

	return ffnLabHsh, faaLabHsh, fidFfnHsh, fidFaaHsh

def printFxxLabHsh(fxxLabhsh, fxxFidHsh, fNm):
	err("Printing: " + fNm)

	f = open(fNm, 'w')

	for i in fxxLabhsh:
		arr = [i]
		arr.append(fxxFidHsh[i][0].replace('|','_'))
		for j in fxxLabhsh[i]:
			arr.append(j + ',' + str(fxxLabhsh[i][j]))
		aStr = '\t'.join(arr)
		f.write(aStr + '\n')

	f.close()

	err('\n')

def makeRevHsh(hsh):
	rHsh = {}
	for i in hsh:
		val = i
		key = hsh[i]
		if key not in rHsh:
			rHsh[key] = []
		rHsh[key].append(val)

	return rHsh

def printIDEqs(fxxFidHsh, fNm):
	err("Printing: " + fNm)

	f = open(fNm, 'w')

	for i in fxxFidHsh:
		sids = fxxFidHsh[i]
		f.write('\t'.join(sids) + '\n')

	f.close()

	err('\n')

def printFasta(fNm, sid, seq):
	f = open(fNm, 'w')

	f.write('>' + sid + '\n')
	f.write(seq)

	f.close()

def printFastas(fxxFidHsh, dNm):
	if dNm[-1] != '/':
		dNm += '/'

	if not os.path.exists(dNm):
		os.mkdir(dNm)

	err("Printing Fastas: " + dNm + '\n\t')
	cnt = 0
	inc = len(fxxFidHsh) / 50.
	for i in fxxFidHsh:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		seq = i
		sid = fxxFidHsh[i][0].replace('|','_')
		fNm = dNm + sid + '.fasta'
		printFasta(fNm, sid, seq)
		# break
	err('\n')

def main():
	fLst = getFLst(argv[1])
	# fLst = fLst[:100]
	ffnLabHsh, faaLabHsh, fidFfnHsh, fidFaaHsh = parseGenomes(fLst)
	ffnFidHsh = makeRevHsh(fidFfnHsh)
	faaFidHsh = makeRevHsh(fidFaaHsh)
	printFxxLabHsh(ffnLabHsh, ffnFidHsh, argv[2] + '.ffn.tab')
	printFxxLabHsh(faaLabHsh, faaFidHsh, argv[2] + '.faa.tab')
	printIDEqs(ffnFidHsh, argv[2] + '.ffn.ids.equiv.tab')
	printIDEqs(faaFidHsh, argv[2] + '.faa.ids.equiv.tab')

if __name__ == '__main__':
	main()
