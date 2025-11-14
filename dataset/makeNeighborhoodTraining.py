'''
python makeNeighborhoodTraining.py [typ.nei.X.lst] [typ.pid.7.lst]
typ = fun, plf, or pgf
X = neighborhood size
'''

from sys import argv

def parseTab(fNm):
	f = open(fNm)

	fHsh = {}
	allF = {}
	for i in f:
		i = i.strip('\n').split('\t')
		mid = int(len(i)/2.0)

		lab = i[mid]
		feat = i[:mid] + i[mid+1:]

		feat = '\t'.join(feat)

		# print lab, feat

		if feat not in fHsh:
			fHsh[feat] = {}
		if lab not in fHsh[feat]:
			fHsh[feat][lab] = 0
		fHsh[feat][lab] += 1

		if i[0] not in allF:
			allF[i[0]] = 0

	f.close()

	cnt = 0
	for i in sorted(allF):
		allF[i] = cnt
		cnt += 1

	return fHsh, allF

def main():
	fHsh, allF = parseTab(argv[1])

	cnt = 0
	for i in fHsh:
		print cnt
		for j in fHsh[i]:
			print '\t',fHsh[i][j]
		cnt += 1

if __name__ == '__main__':
	main()
