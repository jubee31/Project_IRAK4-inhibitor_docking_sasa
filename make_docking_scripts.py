# make_docking_scripts.py
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
LIG_DIR = PROJECT_ROOT / "ligands_3d"
SCRIPT_DIR = PROJECT_ROOT / "dock_scripts"

SCRIPT_DIR.mkdir(exist_ok=True)

RECEPTOR = "리간드 PDBQT 파일 경로"

CENTER_X = -23.10
CENTER_Y = -16.47
CENTER_Z = -9.79

SIZE_X = 20
SIZE_Y = 20
SIZE_Z = 20

TEMPLATE = """#!/bin/bash
#SBATCH 설정 작성!!!!

SCRIPT_DIR="스크립트 폴더"
PROJECT_ROOT="기본 프로젝트 폴더"
cd "$PROJECT_ROOT"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate jol

LIG_NAME="{base}"

LIG_FILE="리간드 PDBQT 파일 경로"
OUT_FILE="도킹 결과 PDBQT 파일 경로"
LOG_FILE="도킹 로그 파일 경로"

if [ ! -f "$LIG_FILE" ]; then
  echo "[ERROR] Ligand file not found: $LIG_FILE"
  exit 1
fi

if [ -f "$OUT_FILE" ]; then
  echo "[SKIP] $LIG_NAME already docked"
  exit 0
fi

echo "[INFO] Start docking $LIG_NAME"
date

vina \\
  --receptor "{receptor}" \\
  --ligand "$LIG_FILE" \\
  --center_x {cx} --center_y {cy} --center_z {cz} \\
  --size_x {sx} --size_y {sy} --size_z {sz} \\
  --exhaustiveness 8 \\
  --num_modes 3 \\
  --cpu 6 \\
  --out "$OUT_FILE" \\
  > "$LOG_FILE" 2>&1

STATUS=$?

if [ $STATUS -eq 0 ]; then
  echo "[OK] Docking finished for $LIG_NAME"
else
  echo "[FAIL] Docking failed for $LIG_NAME (exit code $STATUS)"
fi

date
"""

def main():
    pdbqts = sorted(LIG_DIR.glob("*.pdbqt"))
    if not pdbqts:
        print(f"[ERROR] {LIG_DIR} 안에 .pdbqt가 없습니다. 먼저 SDF -> PDBQT 변환을 하세요.")
        return

    print(f"[INFO] Found {len(pdbqts)} ligands in {LIG_DIR}")

    for lig_file in pdbqts:
        base = lig_file.stem  # CHEMBLxxxx
        sh_path = SCRIPT_DIR / f"dock_{base}.sh"

        text = TEMPLATE.format(
            base=base,
            receptor=RECEPTOR,
            cx=CENTER_X,
            cy=CENTER_Y,
            cz=CENTER_Z,
            sx=SIZE_X,
            sy=SIZE_Y,
            sz=SIZE_Z,
        )

        sh_path.write_text(text)
        sh_path.chmod(0o755)
        print(f"[MAKE] {sh_path.name}")

    print("[DONE] 모든 dock_*.sh 생성 완료")

if __name__ == "__main__":
    main()
