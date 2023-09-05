import os,sys

CLEAN_LOGS = False
CLEAN_RESULTS = True
CLEAN_CORE_DUMPS = False

ffinished = open("finished.txt",'r')

lfin = ffinished.readlines()

for l in lfin:
    l = l.strip()
    print(l)
    gridid = int(l)
    jobdir="job_fc2dprofile_%d"%(gridid)

    logfile = jobdir+"/log_fc2dprofile_%d.txt"%(gridid)

    if CLEAN_LOGS and os.path.exists(logfile):
        print("removing logfile: ",logfile)
        os.system("rm %s"%(logfile))

    if CLEAN_RESULTS and os.path.exists(jobdir)==True:
        print("removing jobdir: ",jobdir)
        os.system("rm -r %s"%(jobdir))

    if CLEAN_CORE_DUMPS:
        print("removing core dumps, if any")
        os.system("rm %s/core.*"%(jobdir))

    


