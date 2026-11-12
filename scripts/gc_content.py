from Bio import SeqIO
from pathlib import Path
import pandas as pd
species = ["ecoli", "efergusonii", "shigella", "salmonella", "pseudomonas"]
rows = []
for sp in species:
    file_path = Path(f"data/{sp}/cds.fna")
    if not file_path.exists(): continue
    for record in SeqIO.parse(file_path, "fasta"):
        seq = str(record.seq).upper()
        gc = seq.count("G") + seq.count("C")
        valid = seq.count("A") + seq.count("T") + gc
        if valid == 0: continue
        rows.append({"species": sp, "gene_id": record.id.split()[0], "length": len(seq), "GC_percent": 100 * gc / valid})
df = pd.DataFrame(rows)
df["genome_mean_GC"] = df.groupby("species")["GC_percent"].transform("mean")
df["difference_from_mean"] = (df["GC_percent"] - df["genome_mean_GC"]).abs()
df["unusual_GC"] = df["difference_from_mean"] > 10
df.to_csv("results/gc_content.csv", index=False)
