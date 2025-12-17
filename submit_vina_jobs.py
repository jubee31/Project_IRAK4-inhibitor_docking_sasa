# submit_vina_jobs.py
import os
from pathlib import Path
from tqdm import tqdm

PROJECT_ROOT = Path(__file__).resolve().parent
SCRIPT_DIR = PROJECT_ROOT / "dock_scripts"
VINA_OUT_DIR = PROJECT_ROOT / "vina_out"


def is_done_from_script(sh: Path) -> bool:
    """
    dock_CHEMBLxxxx.sh 이름에서 CHEMBLxxxx를 뽑아서
    vina_out/CHEMBLxxxx_docked.pdbqt 존재 여부로 완료 여부 체크
    """
    name = sh.name
    if not name.startswith("dock_") or not name.endswith(".sh"):
        return False
    base = name[len("dock_"):-len(".sh")]  # CHEMBLxxxx
    out_file = VINA_OUT_DIR / f"{base}_docked.pdbqt"
    return out_file.is_file()

def main():
    scripts = sorted(SCRIPT_DIR.glob("dock_*.sh"))
    if not scripts:
        print(f"[ERROR] {SCRIPT_DIR} 안에 dock_*.sh가 없습니다. 먼저 make_docking_scripts.py를 실행하세요.")
        return

    print(f"[INFO] Found {len(scripts)} docking scripts in {SCRIPT_DIR}")

    for sh in tqdm(scripts, desc="Submitting docking jobs", unit="job"):
        if not sh.is_file():
            continue

        if is_done_from_script(sh):
            tqdm.write(f"[SKIP] {sh.name} (이미 도킹 완료)")
            continue

        os.system(f'while [ $(squeue | grep username | wc -l) -gt 50 ]; do sleep 5; done')
        os.system(f'cd "{sh.parent}" && sbatch "{sh.name}"')
        tqdm.write(f"[SUBMIT] {sh.name}")

    print("\n[SUMMARY] 제출 루프 종료")

if __name__ == "__main__":
    main()
