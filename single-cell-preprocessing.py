import scanpy as sc
import numpy as np
import scrublet as scr
import sys,os
import pandas as pd


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80)

#'/Volumes/Xin-Y/matrices/BM/Gene/filtered/'
print("l'implementation de donnees control et traité ")
fichier_traite=sys.argv[1]
adata = sc.read_10x_mtx(fichier_traite, cache=True)
sc.pl.highest_expr_genes(adata, n_top=20)

print("matrice de comptage")
print(adata.to_df())

print("l'étape de filtrage des genes")

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(adata)

print("calcul de pourcentage de cellules mitochodriales")
mito_genes = adata.var_names.str.startswith('mt-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


print("visualiser les plots pour choisir les parametres à conserver")
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes',color="percent_mito")

print("filter des genes indésirés or ceux qui sont dans les cellules mitochondriales")
n_genes=input("choisir un nombre de genes que vous voulez conserver:")
percent_mito=input("choisir un percentage de mitochondrie que vous voulez conserver:")
adata = adata[adata.obs['n_genes'] < int(n_genes), :]
adata = adata[adata.obs['percent_mito'] < float(percent_mito), :]
print(adata)

print("prediction des doublets")
scrub = scr.Scrublet(adata.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
scrub.plot_histogram()

sum(adata.obs['predicted_doublets'])


print("filter des doublets")
adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
sc.pl.violin(adata, 'n_genes',jitter=0.4, groupby = 'doublet_info', rotation=45)
adata = adata[adata.obs['doublet_info'] =='False', :]
print(adata)

print("l'etape de normalisation")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("distinguer les genes differentiellement exprimés et conserver que des genes qui sont differentiel exprimé")
nomination_echantillon=input("inserer le nom d'échantillonage:")
adata.obs['sample'] = nomination_echantillon
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key = 'sample')
sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

print("fin de pre-processing")


if len(sys.argv)==2:
    print(adata)
elif len(sys.argv)==3:
    a=input("barcode? (yes or no)")
    if a =="yes":
        barcodes = sys.argv[2]
        df = pd.DataFrame(adata.obs)
        fichier = df.iloc()[:, 0]
        fichier.to_csv(os.path.join(barcodes), columns=[], header=False)
    else:
        save_file = sys.argv[2]
        adata.write_h5ad(save_file)
    print(adata)
elif len(sys.argv)==4:
    print("Le fichier de result et barcodes filtés sont générés")
    print(adata)
    save_file = sys.argv[2]
    adata.write_h5ad(save_file)
    barcodes = sys.argv[3]
    df = pd.DataFrame(adata.obs)
    fichier = df.iloc()[:, 0]
    fichier.to_csv(os.path.join(barcodes), columns=[], header=False)


else:
    raise ImportError ("Inserer au moins d'un fichier ou votre fromat de fichier n'est pas correct.")
