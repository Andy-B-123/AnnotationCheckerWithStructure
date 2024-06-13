import os
import subprocess
from pathlib import Path

def check_foldseek_availability(foldseek_path):
    if os.system("command -v foldseek > /dev/null 2>&1") == 0:
        return "foldseek"
    elif foldseek_path and os.path.isfile(foldseek_path) and os.access(foldseek_path, os.X_OK):
        return foldseek_path
    else:
        raise FileNotFoundError("foldseek is not available on the system path or the provided path is invalid.")

def run_foldseek(fasta_file, database_path, output_directory, foldseek_path=None, threads=1):
    foldseek_exec = check_foldseek_availability(foldseek_path)
    weights_path = r"/scratch3/bac050/AnnotationCheckerWithStructure/weights"
    output_file = str((output_directory / "filtered_proteins.foldseek.out"))
    tmp_dir = str((output_directory / "filtered_proteins.foldseek.tmp"))

    command = [
        foldseek_exec, "easy-search", str(fasta_file), database_path, output_file, tmp_dir,
        "--prostt5-model", weights_path, "--format-mode", "2", "--threads", str(threads)
    ]
    print(command)
    print(f"Executing command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
        print(f"Command executed successfully: {' '.join(command)}")
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
