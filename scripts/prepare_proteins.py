from Bio import SeqIO
from pathlib import Path
species = ["ecoli", "efergusonii", "shigella", "salmonella", "pseudomonas"]
Path("data/combined").mkdir(parents=True, exist_ok=True)
with open("data/combined/all_proteins.faa", "w") as out:
    for sp in species:
        file_path = Path(f"data/{sp}/proteins.faa")
        if not file_path.exists(): continue
        for record in SeqIO.parse(file_path, "fasta"):
            original_id = record.id.split()[0] 
            record.id = f"{sp}|{original_id}"
            record.description = "" 
            SeqIO.write(record, out, "fasta")
