import os,sys
filename=os.listdir(sys.argv[1])
name=sys.argv[2]
hydrogen={}
with open(sys.argv[1]+'/'+name+'.hydrogen_bonds.txt')as f:
    g=f.readline()
    while g:
        g=g.replace('\n','')
        l=g.split(',')
        t=','.join(l[0:4])
        hydrogen[t]=g
        g=f.readline()
charge={}
with open(sys.argv[1]+'/'+name+'.charge-complementary.txt')as f:
    g=f.readline()
    while g:
        g=g.replace('\n','')
        l=g.split(',')
        t=','.join(l[0:4])
        charge[t]=g
        g=f.readline()
print(charge)
result_hydrogen={}
result_charge={}
hydrogenlist={}
chargelist={}
for i in filename:
    if(name in i)and('.hydrogen_bonds.txt'in i):
        mut=i.replace('.hydrogen_bonds.txt','').split('_')[-1]
        with open(sys.argv[1]+'/'+i)as f:
            g=f.readline()
            while g:
                g=g.replace('\n','')
                l=g.split(',')
                del(l[-3:-2])
                t=','.join(l[0:4])
                if t not in hydrogen:
                    if mut not in result_hydrogen:
                        result_hydrogen[mut]=[]
                    result_hydrogen[mut].append(mut+','+','.join(l))
                    if mut[:-1]==l[2].split(':')[1]:
                        if mut not in hydrogenlist:
                            hydrogenlist[mut]=[]
                        hydrogenlist[mut].append(mut+','+','.join(l))
                g=f.readline()
    if(name in i)and('.charge-complementary.txt'in i):
        mut=i.replace('.charge-complementary.txt','').split('_')[-1]
        with open(sys.argv[1]+'/'+i)as f:
            g=f.readline()
            while g:
                g=g.replace('\n','')
                l=g.split(',')
                q=l[-1]
                l[-1]=l[-2]
                l[-2]=q
                t=','.join(l[0:4])
                if t not in charge:
                    print(t)
                    if mut not in result_charge:
                        result_charge[mut]=[]
                    result_charge[mut].append(mut+','+','.join(l))
                    if mut[:-1]==l[2].split(':')[1]:
                        if mut not in chargelist:
                            chargelist[mut]=[]
                        chargelist[mut].append(mut+','+','.join(l))
                g=f.readline()

with open(name+'_select_hydrogen_bonds.csv','w')as f:
    for i in result_hydrogen:
        for j in result_hydrogen[i]:
            f.write(j+'\n')

with open(name+'_select_hydrogen_bonds_2.csv','w')as f:
    for i in hydrogenlist:
        for j in hydrogenlist[i]:
            f.write(j+'\n')

with open(name+'_select_charge-complementary.csv','w')as f:
    for i in result_charge:
        for j in result_charge[i]:
            f.write(j+'\n')
with open(name+'_select_charge-complementary_2.csv','w')as f:
    for i in chargelist:
        for j in chargelist[i]:
            f.write(j+'\n')
with open(name+'_select_both_all.csv','w')as f:
    for i in result_charge:
        if i in result_hydrogen:
            for j in result_charge[i]:
                f.write(j+'\n')
            for j in result_hydrogen[i]:
                f.write(j+'\n')
with open(name+'_select_both.csv','w')as f:
    for i in chargelist:
        if i in hydrogenlist:
            for j in chargelist[i]:
                f.write(j+'\n')
            for j in hydrogenlist[i]:
                f.write(j+'\n')
with open(name+'_mutlist.txt','w')as f:
    mutlist=[]
    for i in result_charge:
        if i in result_hydrogen:
            if result_charge[i][0]not in mutlist:
                mutlist.append(result_charge[i][0])
                f.write(result_charge[i][0]+'\n')
