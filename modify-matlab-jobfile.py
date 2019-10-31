#!/usr/bin/env python
import os,glob,sys
import re,argparse
from termcolor import colored, cprint
# I focus on one job file in this folder

def process_a_jobfile(jobfile,args):
    with open(jobfile, 'r') as f:
        mfile=f.readlines()

    ln1,ln2=None,None  # the begin and end of the parameter section
    newlines=[0]*5
    for idx, line in enumerate(mfile):
        if (ln1 is None) and re.search("^#+\$",line):
            ln1=idx
        if (ln1 is not None):
            # begin to find all job parameters
            if re.search("^#+\$",line):
                # get job name
                if re.search("^#\$\s+-N",line):
                    newlines[0]=line
                # current working directory
                newlines[1]='#$ -cwd\n'
                # Queue
                newlines[2]='#$ -q {}\n'.format(args.queue)
                # core number
                newlines[3]='#$ -pe smp {}\n'.format(args.core)
                # memory
                newlines[4]='#$ -l mt={}G\n'.format(args.memory)
            else:
                ln2=idx
                break

    for idx, line in enumerate(mfile):
        if re.search("^ncore=",line):
            mfile[idx]='ncore={}\n'.format(args.core)

    # print(ln1)
    # print(ln2)
    # print(newlines)

    mfile=mfile[0:ln1]+newlines+mfile[ln2:]
    with open(jobfile, 'w') as f:
        for line in mfile:
            f.write("{}".format(line))

        #     print('-'*30)


#########################
# main function
if __name__== "__main__":

    # parse args
    parser = argparse.ArgumentParser(description='Modify a meep job file.')

    parser.add_argument('jobfiles',
                        nargs='+',
                        help='py file to run')
    parser.add_argument('-q', '--queue', 
                        default='all.q',
                        help="Define a queue: (all.q|UI|INFORMATICS)")
    parser.add_argument('-c', '--core', type=int,
                        default=16,
                        help="Number of cores")
    parser.add_argument('-m', '--memory', type=int,
                        default=20,
                        help="Minimum of memory of the node")
    args = parser.parse_args()
    if re.search('^INF',args.queue):
        args.queue='INFORMATICS'
    print(args)
    for jobfile in args.jobfiles:
        if not re.search('\.job$',jobfile):
            text=colored('{}'.format(jobfile), 'red')
            userInput = input(text+' may not right! Continue? y/n   ');
            if userInput!='y':
                sys.exit('Get the right job file!')
        process_a_jobfile(jobfile,args)
