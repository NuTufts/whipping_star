import os,sys

fout = open("rerun_list_dm2_um4.txt",'w')

NJOBS = 5000

for i in range(NJOBS):
    print(i," ",0,file=fout)
    
fout.close()
    
