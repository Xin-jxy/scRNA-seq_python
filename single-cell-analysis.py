import scanpy as sc
import anndata as ad
import sys
import pandas as pd
import numpy as np


def concatenation_meme_echantillon(fichier_control,fichier_trait1,fichier_trait2):
    print("l'étape d'implementation de donnees traites et controles et la concatenation")
    adata_ctl = ad.read_h5ad(fichier_control)
    adata_trait1 = ad.read_h5ad(fichier_trait1)
    adata_trait2 = ad.read_h5ad(fichier_trait2)
    trait1=input("saisir le nom d'echantillon traite 1:")
    trait2=input("saisir le nom d'echantillon traite 2:")
    ctl=input("saisir le nom d'echantillon controle:")
    trait=input("saisir le nom d'echantillon traite:")
    adadas_trait=adata_trait1.concatenate(adata_trait2,batch_categories=[trait1,trait2])
    adatas=adata_ctl.concatenate(adadas_trait,batch_categories=[ctl,trait])

    return adatas

def concatenation_entre_echantillon(echantillon1, echantillon2):
    print('la concatenation entre deux differents exemples')
    adata_trait1 = ad.read_h5ad(echantillon1)
    adata_trait2 = ad.read_h5ad(echantillon2)
    e1=input("saisir le nom d'echantillon 1:")
    e2=input("saisir le nom d'echantillon 2:")
    adadas_trait=adata_trait1.concatenate(adata_trait2,batch_categories=[e1,e2])

    return adadas_trait


def un_seul_fichier_calcul(fichier_trait):
    adata_trait= ad.read_h5ad(fichier_trait)
    try:
        sc.pp.combat(adata_trait)
        print("filtrage d'effet batch")
    except:
        print("Attention,vous etes en train d'analyser un seul echantillon ou verifier s'il y a des erreurs!!")
    print(adata_trait)
    sc.tl.pca(adata_trait, svd_solver='arpack')
    sc.pp.neighbors(adata_trait, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_trait)
    sc.tl.leiden(adata_trait)
    sc.tl.louvain(adata_trait)
    sc.tl.rank_genes_groups(adata_trait, "louvain", use_raw=True,method='t-test')

    test=input("est-ce que vous voulez visualiser les listes de gene avant l'annotation? (yes or no)")
    if test=="yes":
        try:
            nb_g=int(input("inserer le nombre de groupe exist pour avoir les donnees de tous les groups clusterises:"))
            g=range(0,nb_g)
            n_g=0
            for i in g:
                cluster=sc.get.rank_genes_groups_df(adata_trait,group=str(i))
                df=pd.DataFrame(cluster)
                sort_by_logfc=df.sort_values(by="logfoldchanges",ascending=False)
                print("G"+str(n_g))
                print(sort_by_logfc.iloc[0:20])
                n_g+=1
        except:
            print("vous saisez un faux nombre de clusters!!")
    else:
        print("cet étape est sauté")

    return adata_trait

def anno_cell_D2(adata):
    print("annotation cellulaire")
    cluster = [adata.obs.louvain.isin(['0', '2', '6', '9', '24', '10', '15']),
               adata.obs.louvain.isin(['7' ,'11', '14', '22']),
               adata.obs.louvain.isin(['4']),
               adata.obs.louvain.isin(['12', '16', '17']),
                adata.obs.louvain.isin(['1', '3', '5', '8', '18', '20', '23', '25']),
               adata.obs.louvain.isin(['13','19']),
               adata.obs.louvain.isin(['26']),
               adata.obs.louvain.isin(['21'])]
    cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Myonuclei','Schwann cells']
    adata.obs.louvain=np.select( cluster, cell_group )

    return adata

def anno_cell_B6(adata):
    print("annotation cellulaire")
    cluster = [adata.obs.louvain.isin(['1', '6', '10', '16', '19']),
               adata.obs.louvain.isin(['0' ,'11', '20', '21']),
               adata.obs.louvain.isin(['3','9']),
               adata.obs.louvain.isin(['4', '14','15']),
               adata.obs.louvain.isin(['2', '5', '7', '8', '13']),
               adata.obs.louvain.isin(['12','18']),
               adata.obs.louvain.isin(['17'])]
    cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Schwann cells']
    adata.obs.louvain=np.select( cluster, cell_group )

    return adata

def sort_marker_gene(adata):
    sc.tl.rank_genes_groups(adata, "louvain", use_raw=True,method='t-test')
    distingtion=input("Quel échantilon que tu veux analyser comme cluster?")
    if distingtion=="D2":
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Myonuclei','Schwann cells']
    elif distingtion=="B6":
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Schwann cells']
    else:
        print("error")
    for i in cell_group:
        cluster=sc.get.rank_genes_groups_df(adata,group=i)
        df=pd.DataFrame(cluster)
        sort_by_logfc=df.sort_values(by="logfoldchanges",ascending=False)
        print(i)
        print(sort_by_logfc.iloc[0:10])
    print("sort table de barcode")
    barcode_tous=sc.get.obs_df(adata,keys=['louvain'])
    df_barcode=pd.DataFrame(barcode_tous)
    Muscle_stem_myo=df_barcode.loc[(df_barcode['louvain']=="Myonuclei")|(df_barcode['louvain']=="Muscle stem cells")]
    fibro=df_barcode.loc[df_barcode['louvain']=="Fibroblasts"]
    a=input("est-ce que tu veux generer un nouveau fichier contenant tous les barcodes? (yes or no)")
    if a =="yes":
        df_barcode.to_csv("barcode.csv")
        Muscle_stem_myo.to_csv("barcode_Muscle_stem_myo.csv")
        fibro.to_csv("barcode_fibro.csv")
    else:
        print("no fichier généré")


def trace_plot(adatas):
    #tracer le plot de pca
    sc.pl.pca(adatas, color=["sample"])
    sc.pl.pca_variance_ratio(adatas, log=True)

    #tracer le plot umap

    sc.pl.umap(adatas, color=['sample'])

    #l'analyse differentiel
    #le clustering de gènes avec l'algo louvain
    sc.pl.umap(adatas, color=['louvain'],legend_loc='on data', frameon=False)

    #l'analyse avec t-test
    sc.pl.rank_genes_groups(adatas, n_genes=15, sharey=False)


    #visualisation avec l'analyse précise
    #sc.tl.paga(adatas, groups='louvain')
    #sc.pl.paga(adata, color=['leiden', 'MEGF6', 'MT-CO2', 'AKR7A3'])

    #visualiser un group contre un autre
    #sc.tl.rank_genes_groups(adatas, 'louvain', groups=['0'], reference='1', method='wilcoxon')
    #sc.pl.rank_genes_groups(adatas, group=["0"],n_genes=20)
    #print(sc.get.rank_genes_groups_df(adatas, group="0"))


if __name__ == '__main__':
    if len(sys.argv)==2:
        fichier_trait=sys.argv[1]
        result_sort=un_seul_fichier_calcul(fichier_trait)
        print(result_sort)
        distingtion=input("Quel échantilon que tu veux analyser?")
        if distingtion=="D2":
            print("C'est D2")
            annot=anno_cell_D2(result_sort)
            mg=sort_marker_gene(annot)
            plot=input("Est-ce que tu veux affichir plot? (yes or no)")
            if plot=="yes":
                trace_plot(annot)
            else:
                print("l'affichage de plot est sauté!!")
        elif distingtion=="B6":
            print("C'est B6")
            annot=anno_cell_B6(result_sort)
            mg=sort_marker_gene(annot)
            plot=input("Est-ce que tu veux affichir plot? (yes or no)")
            if plot=="yes":
                trace_plot(annot)
            else:
                print("l'affichage de plot est sauté!!")
        else:
            trace_plot(result_sort)
    elif len(sys.argv)==4:
        echantillon1=sys.argv[1]
        echantillon2=sys.argv[2]
        res=concatenation_entre_echantillon(echantillon1, echantillon2)
        save_file = sys.argv[3]
        res.write_h5ad(save_file)
    elif len(sys.argv)==5:
        fichier_control=sys.argv[1]
        fichier_trait1=sys.argv[2]
        fichier_trait2=sys.argv[3]
        save_file = sys.argv[4]
        result_sort=concatenation_meme_echantillon(fichier_control,fichier_trait1,fichier_trait2)
        result_sort.write_h5ad(save_file)

    else :
        raise ImportError("Verifier le nombre de fichiers insérés est correcte")
