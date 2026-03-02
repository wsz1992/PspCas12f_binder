import os,sys
from scipy import stats
import numpy as np
import difflib


fasta={}
with open('Cas12f_all.fasta')as f:
    g=f.readline()
    while g:
        if g[0]=='>':
            g=g.replace('\n','')
            name=g.replace('>','')
            g=f.readline()
            g=g.replace('\n','')
            fasta[name]=g
        g=f.readline()
def findseq(l):
    site=int(l[0][1:])
    if l[2]in fasta:
        return(fasta[l[2]][(site-6):(site+5)])
    else:
        return('0')
def shortblast(mutseq,wtseq):
    seq1=mutseq[:5]
    seq2=mutseq[6:]
    maxdiff=0
    t=''
    for i in range(len(wtseq)-11):
        seq=wtseq[i:i+11]
        n=0
        for j in range(11):
            if mutseq[j]==seq[j]:
                n+=1
        diff1=n/11
        if(diff1>maxdiff):
            maxdiff=diff1
            t=seq
    
    return(t,maxdiff)

mutlist=[]
proteinlist={}
with open('mutlist.txt')as f:
    g=f.readline()
    g=f.readline()
    while g:
        g=g.replace('\n','')
        l=g.split('\t')
        l.append(findseq(l))
        mutlist.append(l)
        if l[2]not in proteinlist:
            proteinlist[l[2]]={}
        proteinlist[l[2]][l[0]]=l[1]
        g=f.readline()
sitelist={}
for i in mutlist:
    l=i
    site=int(i[0][1:])
    name=fasta[i[2]][site-1]+str(site)+i[0][0]+'_'+i[2]
    #print(name)
    if name not in sitelist:
        sitelist[name]=[]
        for j in fasta:
            if j!=i[2]:
                if len(i[3])<11:
                    for k in range(11-len(i[3])):
                        i[3]+='*'
                t,maxdiff=shortblast(i[3],fasta[j][(site-20):(site+20)])
                if (t)and(maxdiff>0.5):
                    seqsite=fasta[j].index(t)+5-1
                    q=fasta[j][seqsite+1]+str(seqsite+2)+i[0][0]
                    q1=i[0][0]+str(seqsite+2)
                    if q1 in proteinlist[j]:
                        confirm=str(proteinlist[j][q1])
                    else:
                        confirm=''
                    sitelist[name].append([q,j,str(maxdiff),i[3],str(i[1]),t,confirm])
with open('find_result.txt','w')as f:
    f.write('Site_old\tmut\tprotein\tmaxsame\tAsseq\tvalue\tpredictseq\tconfirm\n')
    for i in sitelist:
        for j in sitelist[i]:
            if j[0][0]!=j[0][-1]:
                f.write(str(i)+'\t'+'\t'.join(j)+'\n')

