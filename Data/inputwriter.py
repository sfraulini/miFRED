import sys
import os
if len(sys.argv)<3:
    print('Syntax:',sys.argv[0],'pathtofolder(nofinalslash)','binfile_name(with_.txt)')
    exit()
new=open(sys.argv[2],'w')
for fa in os.listdir(sys.argv[1]):
    ind = fa.index('.f')
    genome = fa[:ind]
    f=open(sys.argv[1]+'/'+fa)
    for line in f:
        if line[0]=='>':
            print(line.strip().split()[0][1:]+'\t'+genome,file=new)
new.close()
