'''
python oneVAll.py [tab] [out pref]
'''

from sys import argv,stderr
import os

def err(s):
	stderr.write(s)

def parseTab(fNm):
	f = open(fNm)

	allFun = {}
	tab = []
	for i in f:
		i = i.strip('\n').split('\t')
		tab.append(i)
		allFun[i[-1]] = 0

	f.close()

	cnt = 0
	for i in sorted(allFun):
		allFun[i] = cnt
		cnt += 1

	return tab, allFun

def printAllFun(fNm, allFun):
	f = open(fNm, 'w')

	for i in sorted(allFun, key = lambda x: allFun[x]):
		f.write(i + '\t' + str(allFun[i]) + '\n')

	f.close()

def makeOneVAll(tab, fun):
	oTab = []
	for i in tab:
		if i[-1] == fun:
			oTab.append([i[0], '1'])
		else:
			oTab.append([i[0], '0'])

	return oTab

def printTab(fNm, oTab):
	f = open(fNm, 'w')

	for i in oTab:
		f.write('\t'.join(i) + '\n')

	f.close()

def printFunTabs(tab, allFun, dNm):
	if dNm[-1] != '/':
		dNm += '/'

	fNames = []
	inc = len(allFun) / 50.
	cnt = 0
	for i in allFun:
		if cnt >= inc:
			err('=')
			cnt = 0
		cnt += 1

		oTab = makeOneVAll(tab, i)
		fNm = dNm + 'l' + str(allFun[i]) + '.tab'
		printTab(fNm, oTab)
		fNames.append(fNm)
		# break
	err('\n')

	return fNames

def main():
	tab, allFun = parseTab(argv[1])

	printAllFun(argv[2] + '.fun.enc.tab', allFun)

	oDir = argv[2] + 'fun.tabs'
	if not os.path.exists(oDir):
		os.mkdir(oDir)

	fNames = printFunTabs(tab, allFun, oDir)

	for i in fNames:
		print i

if __name__ == '__main__':
	main()





