from pymol import cmd
import os,sys
import re
aa_codes = {
 'A': 'ALA','C': 'CYS','D': 'ASP','E': 'GLU','F': 'PHE','G': 'GLY','H': 'HIS','K': 'LYS','I': 'ILE','L': 'LEU','M': 'MET','N': 'ASN','P': 'PRO','Q': 'GLN','R': 'ARG','S': 'SER','T': 'THR','V': 'VAL','Y': 'TYR','W': 'TRP'
}
aa_codes_reverse = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
    'PHE':'F','GLY':'G','HIS':'H','LYS':'K',
    'ILE':'I','LEU':'L','MET':'M','ASN':'N',
    'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TYR':'Y','TRP':'W'
}
Pos_aa=['R','K','H','M','L','Y','Q']
def getcode(pdbfile):
    seq = ''
    for line in open(pdbfile):
        if line[0:6] =="SEQRES":
            columns = line.split()
            for resname in columns[4:]:
                if resname in aa_codes_reverse:
                    seq = seq + aa_codes_reverse[resname]
    i = 0
    while i < len(seq) :
        print(seq[i:i+64])
        i =i+64
    return(seq)
def pymol_mutation(filepath,filename,asite,chain,mut,outpath):
    cmd.load(filepath+'/'+filename)
    cmd.remove('solvent')
    cmd.wizard('mutagenesis')
    for j in chain:
        chainID=j
        cmd.select('1','chain '+chainID+' and resi '+asite)
        cmd.indicate('1')
        cmd.get_wizard().set_mode(aa_codes[mut])
        cmd.get_wizard().do_select('chain '+chainID+' and resi '+asite)
        cmd.get_wizard().apply()
    cmd.set_wizard(None)
    cmd.save(outpath+'/'+filename.replace('.pdb','_'+asite+mut+'.pdb'))
    cmd.delete('all')
    result=[]
    with open(outpath+'/'+filename.replace('.pdb','_'+asite+mut+'.pdb'))as f:
        g=f.readline()
        while g:
            g=g.replace('\n','')
            result.append(g)
            g=f.readline()
    with open(outpath+'/'+filename.replace('.pdb','_'+asite+mut+'.pdb'),'w')as f:
        for k in filehead:
            f.write(k+'\n')
        for k in result:
            f.write(k+'\n')
        
mutlist=[]
with open(sys.argv[1])as f:
    g=f.readline()
    while g:
        g=g.replace('\n','')
        if g[-1].isalpha():
            mutlist.append(g)
        else:
            for i in Pos_aa:
                mutlist.append(g+i)
        g=f.readline()
filedir=os.listdir(sys.argv[2])
print(filedir)
if ('cas12f'in sys.argv[3])or('mle_trans'in sys.argv[3]):
    Castype='Cas12f'
else:
    Castype='Cas9'
if Castype=='Cas12f':
    chainselect=['A','B']
else:
    chainselect=['A']
change=sys.argv[3]
os.system('mkdir -p '+change)
for i in filedir:
    if '.pdb'in i:
        filehead=[]
        filename=sys.argv[2]+'/'+i
        with open(filename)as f:
            g=f.readline()
            while g:
                if('SEQRES'in g)or('REMARK'in g):
                    filehead.append(g.replace('\n',''))
                g=f.readline()
        protein=getcode(filename)
        for j in mutlist:
            n=int(j[:-1])
            if protein[n-1]!=j[-1]:
                pymol_mutation(sys.argv[2],i,j[:-1],chainselect,j[-1],sys.argv[3])
        #os.system('mv '+sys.argv[2]+'/'+i.replace('.pdb','_*.pdb')+' '+change+'/')
