from pymol import cmd
import os
import sys

# 氨基酸三字母到单字母的映射
aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'
}


def get_residue_atoms_and_ca(chain):
    """
    获取指定链的所有原子，按残基分组，并记录每个残基的 Cα 坐标（如果存在）。
    返回字典：{残基号: {'atoms': 原子列表, 'ca_coord': (x,y,z) 或 None}}
    """
    atoms = cmd.get_model(f'chain {chain}').atom
    res_dict = {}
    for atom in atoms:
        resi = atom.resi
        if resi not in res_dict:
            res_dict[resi] = {'atoms': [], 'ca_coord': None}
        res_dict[resi]['atoms'].append(atom)
        if atom.name == 'CA':
            # 存储 Cα 坐标 (假设 atom.coord 是一个长度为3的列表/元组)
            res_dict[resi]['ca_coord'] = atom.coord
    return res_dict


def min_atom_distance(atomsA, atomsB):
    """计算两组原子之间的最小距离"""
    min_dist = float('inf')
    for a in atomsA:
        ax, ay, az = a.coord  # 解包坐标
        for b in atomsB:
            bx, by, bz = b.coord
            dx = ax - bx
            dy = ay - by
            dz = az - bz
            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
            if dist < min_dist:
                min_dist = dist
                if min_dist == 0:          # 原子重叠，直接返回0
                    return 0.0
    return min_dist


def find_interactions(pdb_file, chainA, chainB, threshold=5.0):
    """分析单个 PDB 文件中两条链之间的界面残基对"""
    print(f"Processing {pdb_file}")
    cmd.load(pdb_file)
    cmd.remove('solvent')          # 去除水分子，避免干扰

    # 检查链是否存在
    chains_present = cmd.get_chains()
    if chainA not in chains_present or chainB not in chains_present:
        print(f"Warning: chain {chainA} or {chainB} not found in {pdb_file}, skipping.")
        cmd.delete('all')
        return []

    # 获取原子数据
    resA = get_residue_atoms_and_ca(chainA)
    resB = get_residue_atoms_and_ca(chainB)

    if not resA or not resB:
        print(f"Warning: No residues found in chain {chainA} or {chainB} in {pdb_file}, skipping.")
        cmd.delete('all')
        return []

    # 只保留含有 Cα 的残基（确保可以计算 Cα 距离）
    resA_items = [(resi, data) for resi, data in resA.items() if data['ca_coord']]
    resB_items = [(resi, data) for resi, data in resB.items() if data['ca_coord']]

    ca_cutoff = threshold + 2.0      # Cα 筛选阈值，避免遗漏
    interactions = []

    for resiA, dataA in resA_items:
        caA = dataA['ca_coord']
        atomsA = dataA['atoms']
        for resiB, dataB in resB_items:
            caB = dataB['ca_coord']
            # 计算 Cα 距离，快速筛选
            dx = caA[0] - caB[0]
            dy = caA[1] - caB[1]
            dz = caA[2] - caB[2]
            ca_dist = (dx*dx + dy*dy + dz*dz) ** 0.5
            if ca_dist < ca_cutoff:
                atomsB = dataB['atoms']
                min_dist = min_atom_distance(atomsA, atomsB)
                if min_dist <= threshold:
                    resnA = atomsA[0].resn   # 同一残基所有原子残基名相同
                    resnB = atomsB[0].resn
                    interactions.append((chainA, resiA, resnA, chainB, resiB, resnB, round(min_dist, 1)))

    cmd.delete('all')
    return interactions


def main():
    if len(sys.argv) < 4:
        print("Usage: python find_dimer.py <pdb_directory> <chainA> <chainB> [threshold]")
        sys.exit(1)

    pdb_dir = sys.argv[1]
    chainA = sys.argv[2]
    chainB = sys.argv[3]
    threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0

    # 收集所有 PDB 文件（支持 .pdb 和 .ent）
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(('.pdb', '.ent'))]
    all_results = []

    for fname in pdb_files:
        pdb_path = os.path.join(pdb_dir, fname)
        interactions = find_interactions(pdb_path, chainA, chainB, threshold)
        for inter in interactions:
            chA, resiA, resnA, chB, resiB, resnB, dist = inter
            aaA = aa_codes.get(resnA, 'X')
            aaB = aa_codes.get(resnB, 'X')
            # 残基标识格式：链_残基号_单字母 (例如 A_123_A)
            resA_str = f"{chA}_{resiA}_{aaA}"
            resB_str = f"{chB}_{resiB}_{aaB}"
            all_results.append([fname, resA_str, resB_str, f"{dist:.1f}"])

    # 输出结果表格
    with open('dimer_bind_result.xls', 'w') as f:
        f.write("PDB\tResidue_A\tResidue_B\tDistance\n")
        for line in all_results:
            f.write("\t".join(line) + "\n")

    # 输出所有界面残基的标识（用于后续突变等）
    mutlist = set()
    for line in all_results:
        mutlist.add(line[1])   # Residue_A
        mutlist.add(line[2])   # Residue_B
    with open('mutlist_dimer.txt', 'w') as f:
        for item in sorted(mutlist):
            f.write(item + "\n")

    print(f"Done. Results written to dimer_bind_result.xls and mutlist.txt")


if __name__ == "__main__":
    main()