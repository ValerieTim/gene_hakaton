import pandas as pd

def main():
    # 1. Медиана по всему протеому E. coli
    gc = pd.read_csv('results/gc_content.csv')
    gc_ecoli = gc[gc.species=='ecoli'].copy()
    median_proteome = gc_ecoli.length.median()

    # 2. Поиск орфанов-одиночек
    m = pd.read_csv('results/presence_matrix_corrected.csv')
    orph = m[(m.ecoli==1) & (m.efergusonii==0) & (m.pseudomonas==0) & (m.salmonella==0) & (m.shigella==0)].copy()
    
    par = pd.read_csv('results/possible_paralogs.csv')
    pids = set(par['query']) | set(par['subject'])
    orph = orph[~orph['query'].isin(pids)]
    orph['acc'] = orph['query'].str.split('|').str[1]

    # 3. Слияние с GC-данными
    gc_ecoli['acc'] = gc_ecoli.gene_id.str.extract(r'(NP_\d+\.\d+|WP_\d+\.\d+|YP_\d+\.\d+)')
    o = orph.merge(gc_ecoli[['acc','length','GC_percent']], on='acc')

    # 4. Считаем отклонение от среднего (50.7929...)
    o['d'] = o.GC_percent - 50.7929621856939

    print(f"Всего орфанов-одиночек: {len(o)}")
    print(f"Медиана длины орфанов: {o.length.median()}")
    print(f"Медиана длины всего протеома: {median_proteome}")
    print(f"Сдвиг > +5 (плюс): {(o.d > 5).sum()}")
    print(f"Сдвиг < -5 (минус): {(o.d < -5).sum()}")

if __name__ == "__main__":
    main()
