from pymol import cmd

import os,sys
import re

sumlist = []
distance=[]

#####
aa_codes = {
     'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E',
     'PHE':'F', 'GLY':'G', 'HIS':'H', 'LYS':'K',
     'ILE':'I', 'LEU':'L', 'MET':'M', 'ASN':'N',
     'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
     'THR':'T', 'VAL':'V', 'TYR':'Y', 'TRP':'W'}#定义一个字典，将三字母氨基酸转换为单字母

DNA = ['DT','DG','DA','DC']##定义一个DNA列表

RNA = ['G', 'C', 'U', 'A']##定义一个RNA列表
#####
####定义函数来找配体
distance=6.0
disrange=[]
for j in range(int(distance*10)):
    disrange.append(round(1.5+j*0.1,1))
def search(lig,sym,chainID):
    
    cmd.load(file_path + '/' + i)
    sss = lig
    cmd.remove('solvent')
    cmd.select('1', ' resn %s in chain %s' % (sss,chainID))
    # cmd.set_name('')
    bindaa={}
    for j in disrange:
        cmd.select('2', 'byres 1 around '+str(round(distance+1.5-j,1)))#原子间的距离
        srlig =[]
        for a in cmd.get_model("2").atom:
            if a.resn in aa_codes:
                aazong = aa_codes[a.resn] + a.resi
                if aazong not in srlig:
                    srlig.append(aazong)
                #print(aazong)
        for k in srlig:
            bindaa[k]=(round(distance+1.5-j,1))
                

    cmd.delete('all')
    t=[]
    t.append(chainID)
    t.append(i.split('.')[0])
    for j in bindaa:
        t.append(j+':'+str(bindaa[j]))
    print(t)
    if len(t)>2 and t not in sumlist: #剔除没有配体的列表
        sumlist.append(t)
  
#######函数结束



file_path = sys.argv[1] ###这个地方是pdb文件的路径
print(file_path)
file = os.listdir(file_path)
for i in file:
    if '.pdb'not in i:
        file.remove(i)
    if 'clean.pdb'in i:
        file.remove(i)
for i in file:
    print(i)
    f = open(file_path+'/'+i)
    ls = []

    otherlig = []  # 定义从序列中找出的配体列表 otherlig 是小分子，金属离子
    DNAlig = [] ###与DNA相互作用
    RNAlig = []#与RNA相互作用
    proteinlig = []# 与蛋白质相互作用
    HET = []
    AAA=[]
    '''
    判断分辨率是否大于3.5A
    '''
    g = open(file_path + '/' + i)
    for reso in g:
         #if reso[10:15] =='X-RAY':
             #for reso in g:
        if 'REMARK   2 RESOLUTION.'in reso:
            AA = re.findall('\d+',reso[22:])
            print(AA)
            if not AA:
                AA=['2','60']
            Aa = ''.join(AA)
            AAA.append(Aa)

    if AAA:
        if int(AAA[0]):
            for j in f:
                if j[0:6]=='SEQRES':
                    colume = j.split()
                    if [colume[2],colume[3]] not in ls:
                        ls.append([colume[2],colume[3],colume[4]])
            f.close()
            #print(ls)
            minls = []
            for m in ls:
                minls.append(int(m[1]))
            minnum = min(minls)
            chainls = []
            chainssssssla = []
            for m in ls:
                if m[2]in DNA:
                    if m[0] not in chainls:
                        chainls.append(m[0])
            n = open(file_path+'/'+i)
            for j in n:
                if j[0:6]=='SEQRES':
                    colume = j.split()
                    if colume[2] in chainls and len(ls)>1:
                        #print(colume)
                        for col in colume[4:]:
                            if col in DNA and col not in DNAlig:
                                DNAlig.append(col)
                            elif col in RNA and col not in RNAlig:
                                RNAlig.append(col)
                            elif col in aa_codes and col not in proteinlig:
                                proteinlig.append(col)
                            elif col not in otherlig and col not in aa_codes and col not in RNA and col not in DNA:
                                otherlig.append(col)
                HETlist = ['HETNAM','HETSYN','HETATM']
                if j[0:3]=='HET' and j[0:6] not in HETlist:
                    colume = j.split()
                    if colume[1] not in HETlist:
                        HET.append(colume[1])
                    if colume[2] not in chainssssssla:
                        chainssssssla.append(colume[2])
            print(DNAlig)
            #print(RNAlig)
            if proteinlig:
                ligand = '+'.join(proteinlig)
                chainID = '+'.join(chainls)
                search(ligand,'III',chainID)
            if DNAlig:
                ligand = '+'.join(DNAlig)
                chainID = '+'.join(chainls)
                search(ligand,'DNA',chainID)
            if RNAlig:
                ligand = '+'.join(RNAlig)
                chainID = '+'.join(chainls)
                search(ligand,'RNA',chainID)
            if otherlig:
                chainID = '+'.join(chainls)
                for ligand in otherlig:
                    search(ligand,ligand,chainID)
            if HET:
                chainID = '+'.join(chainssssssla)
                for ligand in HET:
                    search(ligand,ligand,chainID)
        
with open('bind_result.xls','w')as f:
    for list in sumlist:
        f.write('\t'.join(list)+'\n')
