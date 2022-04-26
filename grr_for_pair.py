"""
Usage: python3 grr_for_pair.py genome1.faa genome2.faa output_folder
"""

import os
import subprocess
import sys

def command_print_run(c):
    print(c)
    subprocess.run(c,shell=True)

def extract_best(name,evalue=0.1):
    res={}
    try:
        for line in open(name):
            if line.startswith("#"):continue
            spl=line.split()
            if not spl[0] in res:
                if (float(spl[10])<evalue):
                    res[spl[0]]=spl[1]
        return res
    except FileNotFoundError:
        return res

def grr_for_pair(n1,n2,f):

    name1=n1[n1.rfind("/")+1:n1.rfind(".")]
    name2=n2[n2.rfind("/")+1:n2.rfind(".")]

    if name1==name2:
        sys.exit(0)

    elif name1>name2:
        n1,n2=n2,n1
        name1,name2=name2,name1

    name="{}/{}-{}".format(f,name1,name2)
    try:
        #check if the calculation done already
        open("{}.grr".format(name))
        sys.exit(0)
    except FileNotFoundError:
        pass

    command_print_run("blastp -subject {} -query {} -outfmt 7 -out {}/{}-{}.blast".format(n1,n2,f,name1,name2))
    command_print_run("blastp -subject {} -query {} -outfmt 7 -out {}/{}-{}.blast".format(n2,n1,f,name2,name1))

    res1=extract_best("{}/{}-{}.blast".format(f,name1,name2))
    res2=extract_best("{}/{}-{}.blast".format(f,name2,name1))

    outf=open('{}.bbh'.format(name),'w')
    for k in res1:
        if res2.get(res1[k],"")==k:
            outf.write("{}\t{}\n".format(res1[k],k))
    outf.close()

    command_print_run("cat {} {} >{}.faa".format(n1,n2,name))
    command_print_run("powerneedle -gapopen 10 -gapextend 0.5 -brief -pairs {}.bbh -identities {}.needle -alignment {}.needlealg {}.faa".format(name,name,name,name))
    os.remove("{}.faa".format(name))
    os.remove("{}.needlealg".format(name))

    outf=open("{}.grr".format(name),'w')
    g1=int(subprocess.run("grep '>' {}|wc -l".format(n1),shell=True, check=True, stdout=subprocess.PIPE).stdout.decode("ascii").split()[0])
    g2=int(subprocess.run("grep '>' {}|wc -l".format(n2),shell=True, check=True, stdout=subprocess.PIPE).stdout.decode("ascii").split()[0])
    bbh=int(subprocess.run("wc -l {}.bbh".format(name),shell=True, check=True, stdout=subprocess.PIPE).stdout.decode("ascii").split()[0])
    grr=0
    min80=0
    min35=0
    for line in open("{}.needle".format(name)):
        if not line.startswith("Sequence1"):
            ident=float(line.split()[2]) ### was line.split()[3] - similarity, changed when analyzing phaster
            grr+=ident/100
            if ident>=80:min80+=1
            if ident>=35:min35+=1

    try:
        grr2=grr/min(g1,g2)
    except ZeroDivisionError:
        grr2=0
    outf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(name1,name2,grr2,grr,g1,g2,bbh,min35,min80))

if __name__== "__main__":
    if len(sys.argv)!=4:
        print(__doc__)
    else:
        grr_for_pair(sys.argv[1],sys.argv[2],sys.argv[3])
