import os,sys

fit_result_dir = "fit_results"
fout = os.listdir(fit_result_dir)

MAXJOBS=5000

finished = []
incomplete = []

nread = 0
for f in fout:
    f = f.strip()
    jobid = int(f.split(".")[-2].split("_")[-1])
    
    fx = open(fit_result_dir+"/"+f.strip())
    fl = fx.readlines()

    print(jobid,f)
        
    if nread>0 and nread%200==0:
        print(jobid,f)
        
    dexp = {}
    imax = 0
    for l in fl:
        l = l.strip()
        info = l.split()
        try:
            iexp = int(info[0])
        except:
            continue
        if iexp>imax:
            imax = iexp
        dexp[iexp] = True
    nexperiments = len(dexp)
    if nexperiments==1000:
        finished.append(jobid)
    else:
        incomplete.append( (jobid,imax) )
        
finished.sort()
incomplete.sort()
        
print("INCOMPLETE LIST")
covered_reruns = {}
for (jobid,imax) in incomplete:
    #print(jobid," ",imax+1,file=frerun)
    covered_reruns[jobid] = imax+1

frerun = open('rerun_w_starts.txt','w')    
for i in range(MAXJOBS):
    if i not in covered_reruns and i not in finished:
        print(i," ",0,file=frerun)

rerun_keys = covered_reruns.keys()
for i in rerun_keys:
    print(i," ",covered_reruns[i],file=frerun)
frerun.close()

ffinished = open("finished.txt","w")
for i in finished:
    print(i,file=ffinished)

print("NUM FINISHED: ",len(finished))
        

