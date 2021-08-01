import sys

def cluster_barcodes(file):
    l=[]
    with open(file,"r") as f:
        rd=f.readlines()
        for line in rd[1:]:
            d={}
            sp_l=line.strip().split(",")
            d["barcodes"]=sp_l[0]
            d["cluster"]=sp_l[1]
            l.append(d)
        return l

def barcodes_celltype(file):
    l=[]
    with open(file,"r") as f:
        rd=f.readlines()
        distingtion=input("Quel échantilon que tu veux analyser?")
        for line in rd[1:]:
            d={}
            sp_l=line.strip().split(",")
            sp_e=sp_l[0].split("-")
            try:
                if sp_e[2]==distingtion or sp_e[1]==distingtion:
                    d["sample"]=sp_e[2]
                    d["barcodes"]=sp_e[0]
                    d["cell_type"]=sp_l[1]
                    l.append(d)
            except:
                d["barcodes"]=sp_e[0]
                d["cell_type"]=sp_l[1]
                l.append(d)
        return l

def res(file_scanpy,file_scvelo):
    l=[]
    l_cluster=[]
    d_cell_group={}
    a=cluster_barcodes(file_scanpy)
    b=barcodes_celltype(file_scvelo)
    for i in a:
        for ii in b:
            if i["barcodes"]==ii["barcodes"]:
                annot={}
                annot["cell_type"]=ii["cell_type"]
                annot["cluster"]=i["cluster"]
                l.append(annot)

    for ii in l:
        if ii["cluster"] not in l_cluster:
            l_cluster.append(ii["cluster"])
    for iii in l_cluster:
        d_cell_group[iii]=[]
    for i3 in l:
        if i3["cell_type"] not in d_cell_group[i3["cluster"]]:
            d_cell_group[i3["cluster"]].append(i3["cell_type"])

    new_dico_count={}
    new_annot=[]
    for i1 in l:
        new_annot.append((i1["cluster"],i1["cell_type"]))

    for aaa in new_annot:
        if aaa[0] not in new_dico_count.keys():
            new_dico_count[aaa[0]]={aaa[1]:1}
        else:
            if aaa[1] not in new_dico_count[aaa[0]].keys():
                new_dico_count[aaa[0]][aaa[1]]=1
            else:
                new_dico_count[aaa[0]][aaa[1]]+=1
    dico_count={int(bb):new_dico_count[bb] for bb in new_dico_count}
    sorted_value=sorted(dico_count)
    sorted_dico={cc:dico_count[cc] for cc in sorted_value}
    return sorted_dico


print(res(sys.argv[1],sys.argv[2]))

