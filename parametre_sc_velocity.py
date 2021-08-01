import sc_velocity
import scvelo as scv
import os,sys
import pandas as pd

def add_and_create_newadata(load_adata,spliced_input,unspliced_input,ambiguous_input,stocage):
    if os.path.isdir(load_adata):
        if os.path.splitext(stocage)[-1]=="h5ad":
            result=sc_velocity.creation_newadata(load_adata,spliced_input,unspliced_input,ambiguous_input)
            sc_velocity.data_preprocess(result,stocage)
        else:
            raise TypeError("Vérifiez la nommination d'extension de resultat sortie, il faut nomminer comme h5ad")
    else:
        raise TypeError("Vérifiez le premier input et rassurez c'est un repertoire!!")

def analyse_annotation_velocity(result):
    print("c'est le dataset de départ")
    print(result)
    annot_dis=input("est-ce que vous voulez annoter les clusters pour les plots?(yes or no)")

    if annot_dis=="yes":
        velo=sc_velocity.calculate_scanpy(result)
        annot=sc_velocity.annotation_cell_type(velo)
        print(annot)
        analysis=sc_velocity.scvelo_analyse(annot)
        sc_velocity.plot(analysis)
        #cellrank_analyse(annot)

    elif annot_dis=="no":
        methode=input("Quel objet de resultat que t'as initialment? Seurat? (yes or no)")
        if methode=="yes":
            analyse=sc_velocity.scvelo_analyse(result)
            sc_velocity.plot(analyse)
        elif methode=="no":
            velo=sc_velocity.calculate_scanpy(result)
            print(velo.layers["spliced"])
            cluster=sc_velocity.cluster_analyse(velo)
            #cellrank_analyse(cluster)
            analysis=sc_velocity.scvelo_analyse(cluster)
            sc_velocity.plot(analysis)

    else:
        raise ValueError("inserer les valeurs reconnus par la commande")

def barcodes_file_associated_cluster(result,barcode):

    velo=sc_velocity.calculate_scanpy(result)
    print(velo)
    print("print les barcodes qui sont associés avec cluster")
    df = pd.DataFrame(velo.obs["louvain"])
    print(df)
    df.to_csv(os.path.join(barcode))
    print("Attention, il n'y a pas d'annotation cellulaire pour un seul echantillon")

if len(sys.argv)==6:
    load_adata=sys.argv[1]
    spliced_input=sys.argv[2]
    unspliced_input=sys.argv[3]
    ambiguous_input=sys.argv[4]
    stocage=sys.argv[5]
    add_and_create_newadata(load_adata,spliced_input,unspliced_input,ambiguous_input,stocage)

elif len(sys.argv)==2:
    r=sys.argv[1]
    result=scv.read(filename=r,cache=True)
    analyse_annotation_velocity(result)

elif len(sys.argv)==3:
    r=sys.argv[1]
    barcodes=sys.argv[2]
    result=scv.read(r, cache=True)
    barcodes_file_associated_cluster(result,barcodes)

else:
    raise ImportError ("Inserer au moins d'un fichier ou votre fromat de fichier n'est pas correct.")
