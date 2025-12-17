from pathlib import Path

RECEPTOR_PDB = Path("리셉터 PDB 파일 경로")     
LIG_DIR = Path("리간드 PDB 파일 경로")        
OUT_DIR = Path("complex_pdb")            
OUT_DIR.mkdir(exist_ok=True)

# 1. 리셉터 ATOM/HETATM 줄 미리 읽어두기
with RECEPTOR_PDB.open() as f:
    receptor_lines = [
        line.rstrip("\n")
        for line in f
        if line.startswith(("ATOM", "HETATM"))
    ]

print(f"Receptor atoms: {len(receptor_lines)}개")

# 2. 각 리간드에 대해 복합체 PDB 생성
for lig_pdb in sorted(LIG_DIR.glob("*_docked.pdb")):
    with lig_pdb.open() as f:
        ligand_lines = [
            line.rstrip("\n")
            for line in f
            if line.startswith(("ATOM", "HETATM"))
        ]

    print(f"{lig_pdb.name}: ligand atoms {len(ligand_lines)}개")

    out_path = OUT_DIR / lig_pdb.name.replace("_docked", "_complex")

    with out_path.open("w") as out:
        all_lines = receptor_lines + ligand_lines
        serial = 1
        for line in all_lines:
            # PDB 포맷에서 atom serial number는 7~11열 (index 6:11)
            if line.startswith(("ATOM", "HETATM")):
                new_line = f"{line[:6]}{serial:5d}{line[11:]}"
                serial += 1
                out.write(new_line + "\n")
        out.write("END\n")

    print(f"→ {out_path} 저장 완료")
