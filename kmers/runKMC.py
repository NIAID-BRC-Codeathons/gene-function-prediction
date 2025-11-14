'''
python runKMC.py [fxx.fasta] [k] [out name pref]
fxx.fasta: fasta file where we want to run KMC for each entry in 
           the fasta

KMC main loop
	for each line in tab:
	- write temp fasta file per sequence
	- run kmc on temp fasta
	- read kmc output into hash
	- get master hash of kmers

Index master hash of kmers

Matrix main loop
	for each entry in hash
	- turn hash entry into vector
	- add vector into matrix with sequence ID
'''

from sys import argv,stderr
import os

PROT = True
# PROT = False

def err(s):
	stderr.write(s)

def parseKMC(fNm):
	f = open(fNm)

	kHsh = {}
	for i in f:
		i = i.lower().strip('\n').split('\t')
		kHsh[i[0]] = i[1]

	f.close()

	return kHsh

def getAllKmrs(fNm, k):
	k = int(k)
	if PROT:
		err("Getting all Kmers...\n\t")
		f = open(fNm)

		kHsh = {}
		currSeq = ''
		for i in f:
			if i[0] == '>':
				for i in range(0,len(currSeq)):
					sSeq = currSeq[i:i+k]
					if len(sSeq) < k:
						break
					if sSeq not in kHsh:
						kHsh[sSeq] = 0
				currSeq = ''
				continue
			currSeq += i.strip()

		f.close()
	else:
		err("Getting all Kmers...\n\t")
		err("Running KMC...\n\t")
		cmnd = ['kmc.sh', str(k), fNm, 'temp.fasta', '.', '1', '> /dev/null']
		os.system(' '.join(cmnd))

		err("Reading KMC...\n\t")
		kHsh = parseKMC('temp.fasta.' + str(k) + '.kmrs')

	err("Indexing kmers...")
	cnt = 0
	for i in sorted(kHsh):
		kHsh[i] = cnt
		cnt += 1
	err('\n')

	return kHsh

def rComp(sSeq):
	sSeq = list(sSeq[::-1])
	for i in range(0,len(sSeq)):
		if sSeq[i] == 'a':
			sSeq[i] = 't'
		elif sSeq[i] == 't':
			sSeq[i] = 'a'
		elif sSeq[i] == 'g':
			sSeq[i] = 'c'
		else:
			sSeq[i] = 'g'

	return ''.join(sSeq)

def runKMC(currGen, currSeq, k):
	if currGen == '' or currSeq == '':
		return {}

	# f = open('temp.fasta', 'w')
	# f.write('>' + currGen + '\n' + currSeq)
	# f.close()

	# cmnd = ['kmc.sh', str(k), 'temp.fasta', 'temp.fasta', '.', '1', '> /dev/null']
	# os.system(' '.join(cmnd))

	# return parseKMC('temp.fasta.' + str(k) + '.kmrs')

	currSeq = currSeq.lower()

	k = int(k)
	kHsh = {}
	for i in range(0,len(currSeq)):
		sSeq = currSeq[i:i+k]
		if len(sSeq) < k:
			break

		if not PROT:
			rcSSeq = rComp(sSeq)

			if rcSSeq < sSeq:
				sSeq = rcSSeq

		if sSeq not in kHsh:
			kHsh[sSeq] = 0
		kHsh[sSeq] += 1

	return kHsh

def kHshToArr(kHsh, allKmrs):
	arr = ['0'] * len(allKmrs)
	for i in kHsh:
		val = kHsh[i]
		if i not in allKmrs:
			i = rComp(i)
		if i not in allKmrs:
			continue
		ind = allKmrs[i]
		arr[ind] = str(val)

	return arr

def parseFasta(fNm, allKmrs, k):
	err("Running individual KMCs...")
	f = open(fNm)

	# kmrHsh = {}
	currSeq = ''
	currGen = ''
	for i in f:
		if i[0] == '>':
			kHsh = runKMC(currGen, currSeq, k)
			if len(kHsh) > 0:
				kmrArr = kHshToArr(kHsh, allKmrs)
				# kmrHsh[currGen] = kmrArr
				print '\t'.join([currGen] + kmrArr)#[:10]
			currSeq = ''
			currGen = i.strip().split('>')[1]

			# if len(kmrHsh) > 100:
			# 	break

			continue
		currSeq += i.strip()
	f.close()

	kHsh = runKMC(currGen, currSeq, k)
	if len(kHsh) > 0:
		kmrArr = kHshToArr(kHsh, allKmrs)
		print '\t'.join([currGen] + kmrArr)#[:10]
		# kmrHsh[currGen] = kmrArr

	err('\n')

	# return kmrHsh

def printKmerHsh(kmrHsh, fNm):
	f = open(fNm, 'w')

	for i in kmrHsh:
		arr = [i] + kmrHsh[i]
		f.write('\t'.join(arr) + '\n')

	f.close()

def printAllKmrs(allKmrs, fNm):
	f = open(fNm, 'w')

	for i in sorted(allKmrs):
		arr = [i, str(allKmrs[i])]
		f.write('\t'.join(arr) + '\n')

	f.close()

def main():
	allKmrs = getAllKmrs(argv[1], argv[2])
	printAllKmrs(allKmrs, argv[3] + '.' + str(argv[2]) + '.kmrs.ind.tab')
	kmrHsh = parseFasta(argv[1], allKmrs, argv[2])
	# printKmerHsh(kmrHsh, argv[3] + '.' + str(argv[2]) + '.kmrs.mat.tab')

if __name__ == '__main__':
	main()
