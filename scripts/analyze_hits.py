import pandas as pd
species = ["ecoli", "efergusonii", "shigella", "salmonella", "pseudomonas"]
columns = ["query", "subject", "identity", "length", "evalue", "bitscore", "qcov", "scov"]
hits = pd.read_csv("results/homology.tsv", sep="\t", names=columns)
hits["query_species"] = hits["query"].str.split("|").str[0]
hits["subject_species"] = hits["subject"].str.split("|").str[0]
hits = hits[hits["query"] != hits["subject"]]
hits = hits[(hits["evalue"] <= 1e-5) & (hits["qcov"] >= 50) & (hits["scov"] >= 50)].copy()
presence = pd.crosstab(hits["query"], hits["subject_species"])
for sp in species:
    if sp not in presence.columns: presence[sp] = 0
presence = (presence > 0).astype(int)
presence.to_csv("results/presence_matrix.csv")
same_species = hits[hits["query_species"] == hits["subject_species"]].copy()
same_species = same_species[(same_species["identity"] >= 30) & (same_species["qcov"] >= 50) & (same_species["scov"] >= 50)]
same_species.to_csv("results/possible_paralogs.csv", index=False)
