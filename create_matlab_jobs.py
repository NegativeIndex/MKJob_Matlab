#!/usr/bin/env python3
import os
import re
import subprocess
import glob
import numpy as np
import shutil 
import sys
sys.path.insert(0,'/Users/wdai11/function')
import flatten
import itertools

######################### 
# global variables
##########################
class common:
    fname="run_angles.m"     
    otherfiles=("gd_function.m",)
    nslot=4
    memory=None  #None/64/128/192/256/384/512

######################### 
# generate one job
##########################
def find_line_numbers(lines,regexp):
    numbers=[]
    for i,line in enumerate(lines):
        # print(line)
        if re.match(regexp, line):
            numbers.append(i)
    return numbers

######################### 
# generate one job
##########################
def copy_other_files(dname):
    # copy files
    for ff in common.otherfiles:
        print('Copy file '+ff)
        shutil.copy(ff,dname)

def create_job_file(dname):
    nslot=common.nslot
    memory=common.memory
    # generate job file
    os.chdir(dname)
    fname2="dwt-matlab.job"
    shutil.copy2('/Users/wdai11/bin/dwt-matlab-job-file.job',
                 fname2)
    
    # read file
    with open(fname2, 'r') as f:
        jobfile=f.readlines()

    # edit the jobfile
    res=subprocess.check_output(["unique-number"]).decode("utf-8").rstrip()
    bname=re.sub('.m\s*$','',common.fname)
  
    for idx, line in enumerate(jobfile):
        line=re.sub("CURRENTDIRECTORY", os.getcwd(), line)
        line=re.sub("MYJOBNAME", "dwt-matlab-"+res, line)
        line=re.sub("BASENAME",bname,line)
        line=re.sub("#\$ -pe smp \d+","#$ -pe smp {}".format(nslot),line)
        if memory:
            line=re.sub('^#+\$\s+-l\s+\d+G=true','#$ -l {}G=true'.format(
                memory),line)
        else:
            line=re.sub('^#+\$\s+-l\s+(\d+G=true)','##$ -l \\1'.format(
                memory),line)
        jobfile[idx]=line

    print('Create dwt-matlab.job')
    ##
    with open(fname2, 'w') as f:
        for line in jobfile:
            f.write("{}".format(line))
       
    # generate job.begin
    open("job.begin", 'a').close()
    print("Copy job.begin")
    
########################## 
# main function
##########################
# prepare m file
path=os.path.abspath("./")
fname=common.fname

mlines=[]
with open(fname, "rt") as fin:
    mlines=fin.readlines()

line1=find_line_numbers(mlines,'%% parameters')
line2=find_line_numbers(mlines,r'% end of parameters')

lmin=line1[0]
i=0
while line2[i]<lmin: i+=1
lmax=line2[i]

# all the parameters
# first is the contact, second is mirror
mparas=(('Pd','Al'),)

gparas=[(wl,a,b,t_etch,t_coat)
        for wl in np.array([3.8,])
        for a in np.array([8,])
        for b in a*np.array([0.125,])
        for t_etch in np.arange(0,3.0,0.3)
        for t_coat in np.arange(0,5.2,0.3) 
         ]

paras=list(itertools.product(mparas, gparas ))
paras=[ flatten.flattenArrayN(elem) for elem in paras ]

# paras=gparas
# [print(elem) for elem in paras]

njobs=len(paras)
print('{} jobs are generated'.format(njobs))
print('nslot={}, memory={}'.format(common.nslot,common.memory))
userInput = input('Any key please, Ctrl-c to quit')

count=0
for para in paras:
    os.chdir(path)
    
    mat_contact,mat_mirror,wl,a,b,t_etch,t_coat=para
    dname='G2D_C{}_M{}_wl{:0.2f}_A{:0.2f}B{:0.2f}_TE{:0.2f}C{:0.2f}'.format(
            mat_contact,mat_mirror,wl,a,b,t_etch,t_coat)

    print(dname)
    if not os.path.exists(dname):
        os.makedirs(dname)
        # write the new py file
        fname2=os.path.join(dname,fname)
        for i in range(lmin,lmax):
            line=mlines[i];
            # print(line)
            line=re.sub(r'^(\s*mat_contact=).*$', 
                        r"\g<1>'{}';".format(mat_contact), 
                        line)
            line=re.sub(r'^(\s*mat_mirror=).*$', 
                        r"\g<1>'{}';".format(mat_mirror), 
                        line)

            line=re.sub(r'^(\s*wl=).*$', r'\g<1>{:0.3f};'.format(wl), line)
            line=re.sub(r'^(\s*a=).*$', r'\g<1>{:0.3f};'.format(a), line)
            line=re.sub(r'^(\s*b=).*$', r'\g<1>{:0.3f};'.format(b), line)
            line=re.sub(r'^(\s*t_etch=).*$', 
                        r'\g<1>{:0.3f};'.format(t_etch), line)
            line=re.sub(r'^(\s*t_coat=).*$', 
                        r'\g<1>{:0.3f};'.format(t_coat), line)
            # print(line)
            mlines[i]=line

        with open(fname2, 'w') as f:
            for line in mlines:
                f.write("{}".format(line))
        print("Create {:s}".format(fname))

        os.chdir(path)
        copy_other_files(dname)
        os.chdir(path)
        create_job_file(dname)
        count+=1
        print('Generated {} jobs, {} to do'.format(count, njobs-count))
        print('-------------------')

print('{} jobs are generated'.format(count))

