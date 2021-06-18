import scvelo as scv
import scanpy as sc
import sys
import numpy as np
import pandas as pd
from scipy import sparse

def creation_loom(load_adata,spliced_input,unspliced_input,ambiguous_input):
    print("la création de new anndata")

    adata = sc.read_10x_mtx(load_adata)

    spliced=np.loadtxt(spliced_input, skiprows=3, delimiter=' ')
    f = sc.read(spliced_input)
    print(f)
    rowattr=input("inserer l'attribut de row:")
    colattr=input("inserer l'attribut de col:")
    adata.layers['spliced']=sparse.csr_matrix((spliced[:,2], (spliced[:,0]-1, spliced[:,1]-1)), shape = (int(rowattr),int(colattr))).tocsr().T

    unspliced=np.loadtxt(unspliced_input, skiprows=3, delimiter=' ')
    adata.layers['unspliced']=sparse.csr_matrix((unspliced[:,2], (unspliced[:,0]-1, unspliced[:,1]-1)), shape = (int(rowattr),int(colattr))).tocsr().T

    ambiguous= np.loadtxt(ambiguous_input, skiprows=3, delimiter=' ')
    adata.layers['ambiguous']=sparse.csr_matrix((ambiguous[:,2], (ambiguous[:,0]-1, ambiguous[:,1]-1)), shape = (int(rowattr),int(colattr))).tocsr().T

    # Subset Cells based on STAR filtering
    #selected_barcodes = pd.read_csv(barcode_filtered, header = None)
    #adata = adata[selected_barcodes[0]]

    return adata
    #adata.write_loom(loom,write_obsm_varm=True)

def combine_result(r1,r2,r3):
    ctl=scv.read(r1, cache=True)
    trait_ctl=data_preprocess(ctl)
    trait1=scv.read(r2, cache=True)
    trait_trait1=data_preprocess(trait1)
    trait2=scv.read(r3, cache=True)
    trait_trait2=data_preprocess(trait2)
    trait_coll = scv.utils.merge(trait_trait1, trait_trait2)
    all= scv.utils.merge(trait_coll, trait_ctl)

    return all

def data_preprocess(adata):
    scv.set_figure_params()

    print("pre-processing de data")
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)

    print("l'étape de normalisation")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    print(adata)

    return adata

def analyse_velocite(adata):
    print("le début d'analyse vélocité")
    scv.pl.proportions(adata)

    print("API")
    scv.pp.pca(adata)
    scv.pp.neighbors(adata)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.umap(adata)
    scv.tl.louvain(adata)


    print("l'estimation d'ARN velocité")
    scv.tl.velocity(adata, mode='stochastic')
    #scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)

    print("visaulaisation d'ARN vélocité")
    #dynamical modeling
    scv.tl.recover_dynamics(adata)
    #scv.tl.differential_kinetic_test(adata)

    scv.tl.velocity_clusters(adata)
    scv.tl.velocity_confidence(adata)

    scv.tl.rank_velocity_genes(adata, groupby='velocity_clusters')
    #scv.tl.rank_dynamical_genes(adata)

    scv.tl.terminal_states(adata)
    scv.tl.velocity_pseudotime(adata)
    scv.tl.latent_time(adata)
    #scv.tl.paga(adata,groups='velocity_clusters')


    #tirer une liste de marker gene
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    print(df.head())
    #print(scv.get_df(adata,'rank_velocity_genes'))
    #print(scv.get_df(adata,'rank_dynamical_genes'))


    return adata

def plot(adata):
    #print(adata)
    scv.pl.velocity_embedding(adata, basis='X_umap',color='velocity_clusters')

    scv.pl.velocity_embedding_grid(adata, basis='X_umap',color='velocity_clusters')

    scv.pl.velocity_embedding_stream(adata, basis='X_umap',color='velocity_clusters')

    #scv.pl.velocity(adata, basis='X_umap',color='louvain')
    scv.pl.velocity_graph(adata,basis='X_umap',color='velocity_clusters')

    #scv.pl.paga(adata,basis='umap',color='velocity_clusters')

    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='velocity_clusters', n_convolve=100)


    scv.pl.scatter(adata, color=['root_cells', 'end_points'])

    scv.pl.scatter(adata, color='velocity_pseudotime', color_map='gnuplot')

    scv.pl.scatter(adata, color='velocity_clusters')

    scv.pl.scatter(adata, color='velocity_confidence')

def analyse_partiel(adata):
    df = adata.var
    df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

    kwargs = dict(xscale='log', fontsize=16)
    with scv.GridSpec(ncols=3) as pl:
        pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
        pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
        pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

    scv.get_df(adata, 'fit*', dropna=True).head()


if __name__ == '__main__':
    if len(sys.argv)==5:
        load_adata=sys.argv[1]
        spliced_input=sys.argv[2]
        unspliced_input=sys.argv[3]
        ambiguous_input=sys.argv[4]
        #barcode_input=sys.argv[5]
        result=creation_loom(load_adata,spliced_input,unspliced_input,ambiguous_input)
        adata=data_preprocess(result)
        velo=analyse_velocite(adata)
        plot(velo)

    elif len(sys.argv)==4:
        r1=sys.argv[1]
        r2=sys.argv[2]
        r3=sys.argv[3]
        r=combine_result(r1,r2,r3)
        print(r)
        velo=analyse_velocite(r)
        a=input("est-ce que tu veux visualiser le plot globale (yes or no)")
        if a =="yes":
            plot(velo)
        else:
            analyse_partiel(velo)

    elif len(sys.argv)==2:
        r=sys.argv[1]
        result=scv.read(r, cache=True)
        adata=data_preprocess(result)
        velo=analyse_velocite(adata)
        print(velo)
        print("Attention, il n'y a pas d'annotation cellulaire pour un seul echantillon")
        plot(velo)
    else:
        raise ImportError ("Inserer au moins d'un fichier ou votre fromat de fichier n'est pas correct.")

'''if len(sys.argv)>5:
        load_adata=sys.argv[1]
        spliced_input=sys.argv[2]
        unspliced_input=sys.argv[3]
        ambiguous_input=sys.argv[4]
        loom=sys.argv[5]
        creation_loom(load_adata,spliced_input,unspliced_input,ambiguous_input,loom)
        if len(sys.argv)==5:
            result=creation_loom(load_adata,spliced_input,unspliced_input,ambiguous_input,loom)
            print(result)'''