# --*-- conding:utf-8 --*--
# @time:5/29/25 20:00
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File:create_msa.py

# msa_batcher.py

import subprocess
from pathlib import Path

class MsaBatcher:
    """
    Batch-generate MSAs (A3M) for each sequence in a multi-FASTA using Clustal Omega.

    Attributes
    ----------
    fasta_path : Path
        Path to the input multi-FASTA file.
    msa_dir : Path
        Directory where per-sequence .a3m files will be written.
    clustalo_exe : str
        The Clustal Omega executable name or full path.
    """

    def __init__(self,
                 fasta_path: str = "project/data/seqs/test_sequences.fasta",
                 msa_dir: str   = "project/data/msa",
                 clustalo_exe: str = "clustalo"):
        self.fasta_path = Path(fasta_path)
        self.msa_dir    = Path(msa_dir)
        self.clustalo_exe = clustalo_exe

        # ensure output directory exists
        self.msa_dir.mkdir(parents=True, exist_ok=True)

    def generate(self):
        """
        Parse the input FASTA and run Clustal Omega on each record to produce .a3m MSAs.
        """
        with self.fasta_path.open() as f:
            entries = f.read().strip().split('>')
            for entry in entries:
                if not entry:
                    continue

                lines = entry.splitlines()
                header = lines[0].split()[0]
                seq = ''.join(lines[1:])

                tmp_fasta = self.msa_dir / f"{header}.fasta"
                out_msa   = self.msa_dir / f"{header}.a3m"

                # write out the single-record FASTA
                with tmp_fasta.open('w') as tf:
                    tf.write(f">{header}\n")
                    for i in range(0, len(seq), 60):
                        tf.write(seq[i:i+60] + "\n")



                # run Clustal Omega
                cmd = [
                    self.clustalo_exe,
                    "-i", str(tmp_fasta),
                    "-o", str(out_msa),
                    "--outfmt", "a3m",
                    "--force"
                ]
                print(f"[INFO] Processing {header} -> {out_msa.name}")
                subprocess.run(cmd, check=True)

                # clean up temporary file
                tmp_fasta.unlink()

        print(f"[INFO] All MSAs generated in: {self.msa_dir}")
