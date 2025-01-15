import sys,os,string,glob
PATH = '/home/DATA/YUHL/WES/Batch_macrogen1/CNV/exomedepth/'
filelist=os.listdir(PATH)
fp = [file for file in filelist if file.endswith('_dedup.csv')]
#fp = glob.glob('*_dedup.csv')
os.chdir(PATH)

for fname in fp:
	sample = string.split(fname,'_dedup.csv')[0]
	Sample = string.split(sample,'.')
	Rename = '-'.join(Sample)
	print(Rename)
	os.system('mv '+fname+' '+Rename+'_dedup.csv')

