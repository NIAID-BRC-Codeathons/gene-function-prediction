'''
python subsampleLabsTab.py [lab tab file]
'''

from sys import argv

MAXSAM = 50

def parseTab(fNm):
	f = open(fNm)

	tab = []
	for i in f:
		tab.append(i.strip().split('\t'))

	f.close()

	return tab

def filterTab(tab):
	labCnts = {}
	for i in tab:
		lab = i[1]
		if lab not in labCnts:
			labCnts[lab] = 0
		if labCnts[lab] < MAXSAM:
			print '\t'.join(i)
		labCnts[lab] += 1

def main():
	tab = parseTab(argv[1])
	filterTab(tab)

if __name__ == '__main__':
	main()
