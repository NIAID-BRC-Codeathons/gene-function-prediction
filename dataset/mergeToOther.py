'''
python mergeToOther.py [lab tab file] [filt]
'''

from sys import argv

def parseTab(fNm, filt):
	f = open(fNm)

	tab = []
	labCnts = {}
	for i in f:
		i = i.strip('\n').split('\t')
		tab.append(i)
		if i[-1] not in labCnts:
			labCnts[i[-1]] = 0
		labCnts[i[-1]] += 1

	f.close()

	for i in tab:
		if labCnts[i[-1]] < filt:
			i[-1] = 'OTHER'
		print '\t'.join(i)

def main():
	parseTab(argv[1], int(argv[2]))

if __name__ == '__main__':
	main()
