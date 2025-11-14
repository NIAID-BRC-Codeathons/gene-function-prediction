'''
python filterTab.py [feats.fxx.tab] [out pref]
'''

from sys import argv

hypoSyn = {
	'hypothetical':0,
	'unknown':0,
	'uncharacterized':0
}

def getCnts(arr, cnts):
	hshCnts = {}
	for i in range(0,len(arr)):
		if arr[i] not in hshCnts:
			hshCnts[arr[i]] = 0
		hshCnts[arr[i]] += cnts[i]

	return hshCnts

def parseTab(fNm):
	f = open(fNm)

	# allFuns = {}
	plfTab = []
	pgfTab = []
	funTab = []
	fidSeq = []
	for i in f:
		i = i.strip('\n').split('\t')
		# print i

		seq = i[0]
		fig = i[1]
		plfs = []
		pgfs = []
		funs = []
		cnts = []
		for j in range(2,len(i)):
			i[j] = i[j].split(',')
			plfs.append(i[j][0])
			pgfs.append(i[j][1])
			funs.append(','.join(i[j][2:-1]))
			cnts.append(int(i[j][-1]))

		flg = False
		for j in funs:
			for k in hypoSyn:
				if k.lower() in j:
					flg = True
					break
			if flg:
				break
		if flg:
			continue

		plfCnts = getCnts(plfs, cnts)
		pgfCnts = getCnts(pgfs, cnts)
		funCnts = getCnts(funs, cnts)
		totCnts = sum(cnts)

		mPLF = plfCnts[sorted(plfCnts,key = lambda x: plfCnts[x], reverse = True)[0]]
		mPGF = pgfCnts[sorted(pgfCnts,key = lambda x: pgfCnts[x], reverse = True)[0]]
		mFun = funCnts[sorted(funCnts,key = lambda x: funCnts[x], reverse = True)[0]]

		flg = False
		if mPLF / float(totCnts) > 0.9:
			plf = sorted(plfCnts,key = lambda x: plfCnts[x], reverse = True)[0]
			# print plf
			plfTab.append([seq, fig, plf])
			flg = True
		if mPGF / float(totCnts) > 0.9:
			pgf = sorted(pgfCnts,key = lambda x: pgfCnts[x], reverse = True)[0]
			# print pgf
			pgfTab.append([seq, fig, pgf])
			flg = True
		if mFun / float(totCnts) > 0.9:
			fun = sorted(funCnts,key = lambda x: funCnts[x], reverse = True)[0]
			# print fun
			funTab.append([seq, fig, fun])
			flg = True

		if flg:
			fidSeq.append([fig, seq])

	f.close()

	return plfTab, pgfTab, funTab, fidSeq

def printTab(tab, fNmPref):
	f = [
		open(fNmPref + '.seq.fasta', 'w'),
		open(fNmPref + '.lab.tab', 'w')
	]

	for i in tab:
		seq, fig, lab = i
		f[0].write('>' + fig + '\n' + seq + '\n')
		f[1].write('\t'.join([fig, lab.lower()]) + '\n')

	for i in f:
		i.close()

def printFullFasta(tab, fNm):
	f = open(fNm, 'w')

	for i in tab:
		f.write('>' + i[0] + '\n' + i[1] + '\n')

	f.close()

def main():
	plfTab, pgfTab, funTab, fidSeq = parseTab(argv[1])
	printTab(plfTab, argv[2] + '.plf')
	printTab(pgfTab, argv[2] + '.pgf')
	printTab(funTab, argv[2] + '.fun')
	printFullFasta(fidSeq, argv[2] + '.allSeqs.fasta')

if __name__ == '__main__':
	main()
