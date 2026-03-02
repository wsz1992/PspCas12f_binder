from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from pyrosetta import ScoreFunction
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet,fill_hbond_set
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
import sys,os
mm = MoveMap()
mm.set_bb(False)
mm.set_chi(True)
init()
filename=os.listdir(sys.argv[1])
outdir=sys.argv[1]+'_bonds'
os.system('mkdir '+outdir)
for dname in filename:
    if ('.pdb'in dname)and('.clean'not in dname):
        
        cleanATOM(sys.argv[1]+'/'+dname)
        pose=pose_from_pdb(sys.argv[1]+'/'+dname.replace('.pdb','.clean.pdb'))
        def get_pdb_label(res_id):
            pdb_num = pose.pdb_info().number(res_id)
            icode = pose.pdb_info().icode(res_id)  # 插入码（如' '或'A'）
            chain = pose.pdb_info().chain(res_id)
            pdb_str = f"{pdb_num}{icode}" if icode != ' ' else f"{pdb_num}"
            return f"{pdb_str}{chain}"
        def calculate_net_charge(residue):
            """通过原子电荷累加计算净电荷"""
            total_charge = 0.0
            for atm in range(1, residue.natoms()+1):
                #print(residue.atomic_charge(atm))
                total_charge += residue.atomic_charge(atm)
            return round(total_charge)  # 四舍五入到整数
        scorefxn = create_score_function("ref2015.wts")  # 或者使用"dna_bsc1"
        scorefxn.set_weight(rosetta.core.scoring.fa_elec, 1.0)
        scorefxn(pose)
        min_mover = MinMover(mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, True)
        min_mover.apply(pose)
        protein_res = []
        dna_res = []
        for i in range(1, pose.total_residue()+1):
            if pose.residue(i).is_protein():
                protein_res.append(i)
            elif pose.residue(i).is_DNA():
                dna_res.append(i)
        print(protein_res)
        print(dna_res)
        hbond_set = pose.get_hbonds()
        hydrogen_bonds = []

        # 遍历氢键
        for hbond in hbond_set.hbonds():
            print(hbond)
            donor_res = hbond.don_res()
            acceptor_res = hbond.acc_res()
            donor_pdb = get_pdb_label(donor_res)
            acceptor_pdb = get_pdb_label(acceptor_res)
            # 判断是否为蛋白-DNA相互作用
            if (donor_res in protein_res and acceptor_res in dna_res) or \
            (donor_res in dna_res and acceptor_res in protein_res):
                #print(hbond)
                #print(donor_res,acceptor_res)
                donor_atom = pose.residue(donor_res).atom_name(hbond.don_hatm())
                acceptor_atom = pose.residue(acceptor_res).atom_name(hbond.acc_atm())
                # 获取残基类型
                donor_type = pose.residue(donor_res).name()  # 如 "LYS"
                acceptor_type = pose.residue(acceptor_res).name()  # 如 "DG"（DNA残基）
                donor_atom = pose.residue(donor_res).atom_name(hbond.don_hatm())
                acceptor_atom = pose.residue(acceptor_res).atom_name(hbond.acc_atm())
                hydrogen_bonds.append({
                    "donor_pdb": donor_pdb,
                    "acceptor_pdb": acceptor_pdb,
                    "donor_res": donor_res,
                    "acceptor_res": acceptor_res,
                    "donor_type": donor_type,  # 添加残基类型
                    "acceptor_type": acceptor_type,  # 添加残基类型
                    "donor_atom": donor_atom.strip(),
                    "acceptor_atom": acceptor_atom.strip(),
                    "energy": hbond.energy()
                })

        print(f"Found {len(hydrogen_bonds)} protein-DNA hydrogen bonds")
        with open(outdir+'/'+dname.replace('.pdb','.hydrogen_bonds.txt'),'w')as f:
            for i in hydrogen_bonds:
                t=[]
                for j in i:
                    t.append(str(j)+':'+str(i[j]))
                f.write(','.join(t)+'\n')
        # 配置打分函数（推荐使用包含DNA的力场）

        # 计算残基对能量
        res_pair_energy = {}

        for prot in protein_res:
            for dna in dna_res:
                try:
                    # 正确的参数顺序：
                    # eval_ci_2b(rsd1, rsd2, pose, emap)
                    emap = rosetta.core.scoring.EMapVector()
                    scorefxn.eval_ci_2b(
                        pose.residue(prot),
                        pose.residue(dna),
                        pose,  # 关键修正：传入完整的pose对象
                        emap
                    )
                    elec_energy = emap[rosetta.core.scoring.fa_elec]
                    #print(elec_energy)
                    if abs(elec_energy) > 0.5:
                        res_pair_energy[(prot, dna)] = elec_energy
                        
                except RuntimeError as e:
                    print(f"跳过残基对 {prot}-{dna}（无相互作用）: {str(e)}")
                    continue
        #print(res_pair_energy)
        # 分析电荷互补性
        charge_pairs = []
        for pair, energy in res_pair_energy.items():
            prot_charge = calculate_net_charge(pose.residue(pair[0]))
            dna_charge = calculate_net_charge(pose.residue(pair[1]))
            prot_pdb = get_pdb_label(pair[0])
            dna_pdb = get_pdb_label(pair[1])
            # 获取残基类型
            prot_type = pose.residue(pair[0]).name()  # 蛋白质残基类型
            dna_type = pose.residue(pair[1]).name()  # DNA残基类型
            if (prot_charge * dna_charge) < 0:  # 异号电荷
                charge_pairs.append({
                    "protein_pdb": prot_pdb,
                    "dna_pdb": dna_pdb,
                    "protein_res": pair[0],
                    "dna_res": pair[1],
                    "protein_type": prot_type,  
                    "dna_type": dna_type,  
                    "energy": energy,
                    "charge_complementary": True
                })

        print(f"Found {len(charge_pairs)} charge-complementary pairs")
        with open(outdir+'/'+dname.replace('.pdb','.charge-complementary.txt'),'w')as f:
            for i in charge_pairs:
                t=[]
                for j in i:
                    t.append(str(j)+':'+str(i[j]))
                f.write(','.join(t)+'\n')

