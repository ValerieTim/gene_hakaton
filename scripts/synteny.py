from pathlib import Path
from collections import defaultdict
import pandas as pd
species = ["ecoli", "efergusonii", "shigella", "salmonella", "pseudomonas"]
def parse_attributes(text):
    result = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result
rows = []
for sp in species:
    gff_path = Path(f"data/{sp}/annotation.gff")
    if not gff_path.exists(): continue
    contigs = defaultdict(list)
    with open(gff_path) as file:
        for line in file:
            if line.startswith("#"): continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS": continue
            contig, _, _, start, end, _, strand, _, attr = fields
            attr = parse_attributes(attr)
            protein_id = attr.get("protein_id")
            if not protein_id: continue
            contigs[contig].append({"species": sp, "protein_id": protein_id, "start": int(start), "end": int(end), "strand": strand})
    for contig, genes in contigs.items():
        genes.sort(key=lambda x: x["start"])
        for i, gene in enumerate(genes):
            left = genes[max(0, i - 3):i]
            right = genes[i + 1:i + 4]
            rows.append({"species": sp, "protein_id": gene["protein_id"], "contig": contig, "left_neighbors": ";".join(x["protein_id"] for x in left), "right_neighbors": ";".join(x["protein_id"] for x in right)})
pd.DataFrame(rows).to_csv("results/gene_neighborhoods.csv", index=False)
