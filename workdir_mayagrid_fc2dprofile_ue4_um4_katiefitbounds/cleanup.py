import os,sys

CLEAN_LOGS = True
CLEAN_RESULTS = False

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

