from pathlib import Path
import pandas as pd
from rdkit import Chem
import freesasa
import math
import re

FRAG_CSV = Path("fragment 정보 csv")     
LIG_DIR = Path("ligand 단독 PDB 폴더")                  
COMPLEX_DIR = Path("protein+ligand complex PDB 폴더")               
OUT_CSV = Path("output_sasa_results.csv")

# 단백질 ATOM 마지막 serial 번호 (complex PDB에서 ligand가 시작하는 atom index 기준)
N_REC = 2185  # receptor atom 수 (0-based index에서 ligand는 N_REC 이후)

def clean_smarts_to_pattern(smarts: str):
    if not isinstance(smarts, str):
        return None

    s = re.sub(r"\[\d+\*\]", "[*]", smarts)

    q = Chem.MolFromSmarts(s)
    if q is None:
        return None

    Chem.RemoveStereochemistry(q)

    return q


def fix_unknown_radii(structure):
    n = structure.nAtoms()
    radii = []

    for i in range(n):
        r = structure.radius(i)

        atom_name = structure.atomName(i).strip()
        el = atom_name.upper()

        if r <= 0 and el.startswith("CL"):
            r=1.75 # 문헌값
        elif r <= 0 and el.startswith("BR"):
            r=1.87 # 문헌값
        radii.append(r)

    structure.setRadii(radii)

    for i in range(structure.nAtoms()):
        r = structure.radius(i)
        atom_name = structure.atomName(i).strip()
        el = atom_name.upper()
        if r <= 0 and (el.startswith("CL") or el.startswith("BR")):
            rn = structure.atomResName(i).strip()
            raise RuntimeError(
                f"Cl/Br radius still <=0 after fix for atom {i} ({rn} {atom_name}), check mapping. (확인요망)"
            )
        
# 1. ligand 단독(unbound)에서 SASA 계산
def compute_unbound_sasa(ligand_pdb: Path):

    params = freesasa.Parameters()
    params.setAlgorithm(freesasa.LeeRichards)

    structure = freesasa.Structure(str(ligand_pdb))
    fix_unknown_radii(structure) 
    result = freesasa.calc(structure, params) 

    n_atoms = structure.nAtoms()
    atom_areas = [result.atomArea(i) for i in range(n_atoms)]
    lig_sasa_unbound = sum(atom_areas)

    return lig_sasa_unbound, atom_areas, n_atoms


# 2. complex(bound)에서 SASA 계산
def compute_bound_sasa(complex_pdb: Path):
    params = freesasa.Parameters()
    params.setAlgorithm(freesasa.LeeRichards)

    structure = freesasa.Structure(str(complex_pdb))
    fix_unknown_radii(structure) 
    result = freesasa.calc(structure, params) 

    total_sasa = result.totalArea()
    n_atoms = structure.nAtoms()
    atom_areas = [result.atomArea(i) for i in range(n_atoms)]

    if n_atoms <= N_REC:
        raise ValueError(f"{complex_pdb} : n_atoms={n_atoms}, N_REC={N_REC} 으로 ligand가 없는 것처럼 보임")

    ligand_indices = range(N_REC, n_atoms)
    lig_sasa_bound = sum(atom_areas[i] for i in ligand_indices)

    return total_sasa, lig_sasa_bound, atom_areas, n_atoms


# 3. 라벨링 함수
def classify_orientation(orientation_abs: float):

    if orientation_abs is None or (isinstance(orientation_abs, float) and math.isnan(orientation_abs)):
        return "NA"

    if orientation_abs >= 0.20:
        return "outward"
    elif orientation_abs <= 0.05:
        return "inward"
    else:
        return "middle"

def classify_rel_unbound(rel_unbound: float):

    if rel_unbound is None or (isinstance(rel_unbound, float) and math.isnan(rel_unbound)):
        return "NA"
    if rel_unbound >= 0.20:
        return "outer_like"
    elif rel_unbound <= 0.05:
        return "core_like"
    else:
        return "middle_region"


# 4. 메인 파이프라인
def main():
    df = pd.read_csv(FRAG_CSV)

    df["frag_sasa_unbound"] = float("nan")
    df["frag_sasa_bound"] = float("nan")
    df["lig_sasa_unbound"] = float("nan")
    df["lig_sasa_bound"] = float("nan")
    df["total_sasa_bound"] = float("nan")

    df["orientation_abs"] = float("nan")   # frag_bound / lig_bound
    df["rel_unbound"] = float("nan")       # frag_unbound / lig_unbound

    df["orientation_label"] = "NA"         # inward/middle/outward
    df["rel_unbound_label"] = "NA"         # core_like/middle_region/outer_like


    unbound_cache = {}  
    bound_cache = {}   

    for idx, row in df.iterrows():
        chembl_id = str(row["Molecule ChEMBL ID"])
        smiles = row["Smiles"]
        frag_smarts = row["fragment"]

        lig_pdb = LIG_DIR / f"{chembl_id}_docked.pdb"
        complex_pdb = COMPLEX_DIR / f"{chembl_id}_complex.pdb"

        if not lig_pdb.is_file():
            print(f"[SKIP] {chembl_id}: ligand PDB 없음 → {lig_pdb}")
            continue
        if not complex_pdb.is_file():
            print(f"[SKIP] {chembl_id}: complex PDB 없음 → {complex_pdb}")
            continue

        if chembl_id not in unbound_cache:
            try:
                unbound_cache[chembl_id] = compute_unbound_sasa(lig_pdb)
            except Exception as e:
                print(f"[SKIP] {chembl_id}: unbound SASA 계산 실패: {e}")
                continue

        lig_sasa_unbound, atom_areas_unbound, n_atoms_unbound = unbound_cache[chembl_id]

        if chembl_id not in bound_cache:
            try:
                bound_cache[chembl_id] = compute_bound_sasa(complex_pdb)
            except Exception as e:
                print(f"[SKIP] {chembl_id}: bound SASA 계산 실패: {e}")
                continue

        total_sasa_bound, lig_sasa_bound, atom_areas_bound, n_atoms_bound = bound_cache[chembl_id]

        df.at[idx, "lig_sasa_unbound"] = lig_sasa_unbound
        df.at[idx, "lig_sasa_bound"] = lig_sasa_bound
        df.at[idx, "total_sasa_bound"] = total_sasa_bound

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[SKIP] {chembl_id}: SMILES 파싱 실패")
            continue

        frag_query = clean_smarts_to_pattern(frag_smarts)
        if frag_query is None:
            print(f"[SKIP] {chembl_id}: fragment SMARTS 정리/파싱 실패: {frag_smarts}")
            continue

        matches = mol.GetSubstructMatches(frag_query)
        if not matches:
            print(f"[SKIP] {chembl_id}: fragment SMARTS 매칭 없음 (cleaned): {frag_smarts}")
            continue

        rdkit_indices = matches[0]  

        frag_indices_unbound = list(rdkit_indices)
        if any(i >= n_atoms_unbound for i in frag_indices_unbound):
            print(f"[SKIP] {chembl_id}: fragment index가 unbound n_atoms({n_atoms_unbound})를 넘김")
            continue

        frag_indices_bound = [N_REC + i for i in rdkit_indices]
        if any(i >= n_atoms_bound for i in frag_indices_bound):
            print(f"[SKIP] {chembl_id}: fragment index가 bound n_atoms({n_atoms_bound})를 넘김")
            continue

        frag_sasa_unbound = sum(atom_areas_unbound[i] for i in frag_indices_unbound)
        frag_sasa_bound = sum(atom_areas_bound[i] for i in frag_indices_bound)

        orientation_abs = frag_sasa_bound / lig_sasa_bound if lig_sasa_bound > 0 else float("nan")
        rel_unbound = (
            frag_sasa_unbound / lig_sasa_unbound
            if lig_sasa_unbound > 0
            else float("nan")
        )

        ori_label = classify_orientation(orientation_abs)
        rel_label = classify_rel_unbound(rel_unbound)

        df.at[idx, "frag_sasa_unbound"] = frag_sasa_unbound
        df.at[idx, "frag_sasa_bound"] = frag_sasa_bound
        df.at[idx, "orientation_abs"] = orientation_abs
        df.at[idx, "rel_unbound"] = rel_unbound
        df.at[idx, "orientation_label"] = ori_label
        df.at[idx, "rel_unbound_label"] = rel_label

        print(
            f"{idx}: {chembl_id}, frag={frag_smarts}, "
            f"ori={orientation_abs:.3f}, rel_unbound={rel_unbound:.3f}, "
            f"labels=({ori_label}, {rel_label})"
        )

    df.to_csv(OUT_CSV, index=False)
    print(f"\n저장 완료: {OUT_CSV}")


if __name__ == "__main__":
    main()
