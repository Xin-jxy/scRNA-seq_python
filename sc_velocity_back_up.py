import cellrank as cr
import scvelo as scv
import scanpy as sc
import sys,os
import numpy as np
import scrublet as scr
from scipy import sparse
import pandas as pd


#sys.path.insert('/Users/xinyue/Draft-python/venv/lib/python3.7/site-packages/cellrank','/usr/local/Caskroom/miniconda/base/bin/python3.8/site-packages')

def creation_newadata(load_adata,spliced_input,unspliced_input,ambiguous_input):
    print("la création de new anndata")

    adata = sc.read_10x_mtx(load_adata)

    spliced=np.loadtxt(spliced_input, skiprows=3, delimiter=' ')
    shape = np.loadtxt(spliced_input, skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

    adata.layers['spliced']=sparse.csr_matrix((spliced[:,2], (spliced[:,0]-1, spliced[:,1]-1)), shape = (shape)).tocsr().T

    unspliced=np.loadtxt(unspliced_input, skiprows=3, delimiter=' ')
    adata.layers['unspliced']=sparse.csr_matrix((unspliced[:,2], (unspliced[:,0]-1, unspliced[:,1]-1)), shape = (shape)).tocsr().T

    ambiguous= np.loadtxt(ambiguous_input, skiprows=3, delimiter=' ')
    adata.layers['ambiguous']=sparse.csr_matrix((ambiguous[:,2], (ambiguous[:,0]-1, ambiguous[:,1]-1)), shape = (shape)).tocsr().T

    return adata


def data_preprocess(adata,result):
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80)

    #'/Volumes/Xin-Y/matrices/BM/Gene/filtered/'
    print("l'implementation de donnees control et traité ")
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
    adata.write_h5ad(filename=result)
    return adata


def calculate_scanpy(adata):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.tl.louvain(adata)
    sc.tl.rank_genes_groups(adata, "louvain", use_raw=True,method='t-test')
    #sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False)
    #print(sc.get.rank_genes_groups_df(adata,group="louvain"))
    return adata

def annotation_cell_type(adata):
    distingtion=input("Quel échantilon que tu veux analyser?")
    if distingtion=="D2":
        print("C'est D2")
        cluster = [adata.obs.louvain.isin(['0', '8', '10', '6', '22']),
                   adata.obs.louvain.isin(['4' ,'19', '21']),
                   adata.obs.louvain.isin(['5']),
                   adata.obs.louvain.isin(['13', '17']),
                   adata.obs.louvain.isin(['1', '2', '3','7', '9', '12','14','20']),
                   adata.obs.louvain.isin(['11','15']),
                   adata.obs.louvain.isin(['16']),
                   adata.obs.louvain.isin(['18'])]
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Schwann cells','Myonuclei']

    elif distingtion=="DM":
        print("C'est DM")
        cluster = [adata.obs.louvain.isin(['2', '8']),
                   adata.obs.louvain.isin(['4' ,'6', '14']),
                   adata.obs.louvain.isin(['11']),
                   adata.obs.louvain.isin(['5', '10']),
                   adata.obs.louvain.isin(['0', '1', '3','7', '9', '13','15']),
                   adata.obs.louvain.isin(['12']),
                   adata.obs.louvain.isin(['16'])]
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Myonuclei']

    elif distingtion=="DW":
        print("C'est DW")
        cluster = [adata.obs.louvain.isin(['0', '1', '3', '5', '7','12','18','21']),
                   adata.obs.louvain.isin(['2' ,'10', '13','19','22','24']),
                   adata.obs.louvain.isin(['6','9']),
                   adata.obs.louvain.isin(['4', '15']),
                   adata.obs.louvain.isin(['11']),
                   adata.obs.louvain.isin(['8','14','17','23']),
                   adata.obs.louvain.isin(['16','20']),]
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Schwann cells']

    elif distingtion=="B6":
        print("C'est B6")
        cluster = [adata.obs.louvain.isin(['0', '5', '10', '11', '13']),
                   adata.obs.louvain.isin(['1']),
                   adata.obs.louvain.isin(['4']),
                   adata.obs.louvain.isin(['6']),
                   adata.obs.louvain.isin(['2', '3', '7', '8']),
                   adata.obs.louvain.isin(['9']),
                   adata.obs.louvain.isin(['12'])]
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Schwann cells']

    elif distingtion=="BM":
        print("C'est BM")
        cluster = [adata.obs.louvain.isin(['2', '10', '17']),
                   adata.obs.louvain.isin(['0' ,'15', '18']),
                   adata.obs.louvain.isin(['6']),
                   adata.obs.louvain.isin(['3', '8','9']),
                   adata.obs.louvain.isin(['1', '4', '5','7', '12','13']),
                   adata.obs.louvain.isin(['11','14']),
                   adata.obs.louvain.isin(['16'])]
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells','Schwann cells']

    elif distingtion=="BW":
        print("C'est BW")
        cluster = [adata.obs.louvain.isin(['0', '3', '5', '8', '11']),
                   adata.obs.louvain.isin(['1']),
                   adata.obs.louvain.isin(['2','6','7','10','13','14']),
                   adata.obs.louvain.isin(['12', '15']),
                   adata.obs.louvain.isin(['4']),
                   adata.obs.louvain.isin(['9'])]
        cell_group=['Endothelial cells','Fibroblasts','Tenocytes','Muscle stem cells','Immune cells','Mural cells']

    adata.obs.louvain=np.select( cluster, cell_group )

    return adata

def cluster_analyse(adata):
    #adata.obs['clusters']=scv.load(barcodes)
    distingtion=input("Quel échantilon que tu veux analyser?")
    cluster_dis=input("Quel cluster que tu veux analyser?(MuSC,Endo,C.T,Immune,MC)")
    if distingtion=="DW":
        print(f"C'est {distingtion}")
        if cluster_dis=="MuSC":
            adata_a = adata[adata.obs['louvain'].isin(['4','15'])]
        elif cluster_dis=="Endo":
            adata_a = adata[adata.obs['louvain'].isin(['0', '1', '3', '5', '7','12','18','21'])]
        elif cluster_dis=="C.T":
            adata_a = adata[adata.obs['louvain'].isin(['2' ,'10', '13','19','22','24','6','9'])]
        elif cluster_dis=="Immune":
            adata_a = adata[adata.obs.louvain.isin(['11'])]
        elif cluster_dis=="MC":
            adata_a = adata[adata.obs.louvain.isin(['8','14','17','23'])]
        else:
            raise ValueError(f"Saisir a correct name for our analyse{distingtion}!!")
    elif distingtion=="DM":
        print(f"C'est {distingtion}")
        if cluster_dis=="MuSC":
            adata_a = adata[adata.obs['louvain'].isin(['5', '10'])]
        elif cluster_dis=="Endo":
            adata_a = adata[adata.obs['louvain'].isin(['2', '8'])]
        elif cluster_dis=="C.T":
            adata_a = adata[adata.obs['louvain'].isin(['4' ,'6', '14','11'])]
        elif cluster_dis=="Immune":
            adata_a = adata[adata.obs.louvain.isin(['0', '1', '3','7', '9', '13','15'])]
        elif cluster_dis=="MC":
            adata_a = adata[adata.obs.louvain.isin(['12'])]
        else:
            raise ValueError(f"Saisir a correct name for our analyse{distingtion}!!")

    elif distingtion=="BM":
        print(f"C'est {distingtion}")
        if cluster_dis=="MuSC":
            adata_a = adata[adata.obs['louvain'].isin(['3', '8','9'])]
        elif cluster_dis=="Endo":
            adata_a = adata[adata.obs['louvain'].isin(['2', '10', '17'])]
        elif cluster_dis=="C.T":
            adata_a = adata[adata.obs['louvain'].isin(['0' ,'15', '18','6'])]
        elif cluster_dis=="Immune":
            adata_a = adata[adata.obs.louvain.isin(['1', '4', '5','7', '12','13'])]
        elif cluster_dis=="MC":
            adata_a = adata[adata.obs.louvain.isin(['11','14'])]
        else:
            raise ValueError(f"Saisir a correct name for our analyse{distingtion}!!")

    elif distingtion=="BW":
        print(f"C'est {distingtion}")
        if cluster_dis=="MuSC":
            adata_a = adata[adata.obs['louvain'].isin(['12', '15'])]
        elif cluster_dis=="Endo":
            adata_a = adata[adata.obs['louvain'].isin(['0', '3', '5', '8', '11'])]
        elif cluster_dis=="C.T":
            adata_a = adata[adata.obs['louvain'].isin(['1','2','6','7','10','13','14'])]
        elif cluster_dis=="Immune":
            adata_a = adata[adata.obs.louvain.isin(['4'])]
        elif cluster_dis=="MC":
            adata_a = adata[adata.obs.louvain.isin(['9'])]
        else:
            raise ValueError(f"Saisir a correct name for our analyse{distingtion}!!")

    return adata_a


def scvelo_analyse(adata):

    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print("l'estimation d'ARN velocité")
    scv.tl.velocity(adata, mode='stochastic')
    #scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)

    #dynamical modeling
    scv.tl.recover_dynamics(adata)
    #scv.tl.differential_kinetic_test(adata) #dynamical

    scv.tl.velocity_clusters(adata)
    scv.tl.velocity_confidence(adata)

    scv.tl.rank_velocity_genes(adata, groupby='louvain') #ici le groupby est cencé fait en fonction de groups dans le var,
    # je sais plus si ça peut colorer en obsm, ce dont on a besoin (X_fitsne_by_harmony), pareil pour autres fonctions.

    #scv.tl.rank_dynamical_genes(adata) #dynamical

    scv.tl.terminal_states(adata)

    scv.tl.velocity_pseudotime(adata)
    scv.tl.latent_time(adata)
    scv.tl.paga(adata,groups='louvain')


    return adata


def plot(adata):
    print("visaulaisation d'ARN vélocité")

    scv.pl.velocity_embedding_grid(adata, basis='X_umap',color='louvain')

    scv.pl.velocity_embedding_stream(adata, basis='X_umap',color='louvain')

    #scv.pl.velocity(adata, basis='X_umap',color='louvain')
    scv.pl.velocity_graph(adata,basis='X_umap',color='louvain')

    #scv.pl.paga(adata,basis='umap',color='velocity_clusters')

    scv.pl.scatter(adata, color=['root_cells', 'end_points'])

    scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot')

    scv.pl.scatter(adata, color='louvain')

    scv.pl.scatter(adata, color='velocity_confidence')

    #driver genes
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='louvain', n_convolve=100)


def cellrank_analyse(adata):
    scv.tl.recover_dynamics(adata, n_jobs=8)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)
    scv.pl.velocity_embedding_stream(adata, basis="X_umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4)

    cr.tl.terminal_states(adata, cluster_key="louvain", weight_connectivities=0.2)
    cr.pl.terminal_states(adata)

    cr.tl.initial_states(adata, cluster_key="louvain")
    cr.pl.initial_states(adata, discrete=True)

    cr.tl.lineages(adata)
    cr.pl.lineages(adata, same_plot=False)
    cr.pl.lineages(adata, same_plot=True)

    scv.tl.recover_latent_time(adata, root_key="initial_states_probs", end_key="terminal_states_probs")



if __name__ == '__main__':
    if len(sys.argv)==6:
        load_adata=sys.argv[1]
        spliced_input=sys.argv[2]
        unspliced_input=sys.argv[3]
        ambiguous_input=sys.argv[4]
        store=sys.argv[5]
        result=creation_newadata(load_adata,spliced_input,unspliced_input,ambiguous_input)
        adata=data_preprocess(result,store)

    elif len(sys.argv)==2:
        r=sys.argv[1]
        result=scv.read(filename=r,cache=True)
        print("c'est le dataset de départ")
        print(result)
        annot_dis=input("est-ce que vous voulez annoter les clusters pour les plots?(yes or no)")

        if annot_dis=="yes":
            velo=calculate_scanpy(result)
            annot=annotation_cell_type(velo)
            print(annot)
            analysis=scvelo_analyse(annot)
            plot(analysis)
            #cellrank_analyse(annot)

        elif annot_dis=="no":
            methode=input("Quel objet de resultat que t'as initialment? Seurat? (yes or no)")
            if methode=="yes":
                analyse=scvelo_analyse(result)
                plot(analyse)
            elif methode=="no":
                velo=calculate_scanpy(result)
                cluster=cluster_analyse(velo)
                #cellrank_analyse(cluster)
                analysis=scvelo_analyse(cluster)
                plot(analysis)

        else:
            raise ValueError("inserer les valeurs reconnus par la commande")

    elif len(sys.argv)==3:
        r=sys.argv[1]
        barcodes=sys.argv[2]
        result=scv.read(r, cache=True)
        velo=calculate_scanpy(result)
        print(velo)
        print("print les barcodes qui sont associés avec cluster")
        df = pd.DataFrame(velo.obs["louvain"])
        print(df)
        df.to_csv(os.path.join(barcodes))
        print("Attention, il n'y a pas d'annotation cellulaire pour un seul echantillon")

    else:
        raise ImportError ("Inserer au moins d'un fichier ou votre fromat de fichier n'est pas correct.")

