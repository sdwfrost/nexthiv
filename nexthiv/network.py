from collections import Counter
import networkx as nx
from BioExt.tn93 import tn93

def tn93_network(query,threshold,matchMode,minOverlap):
    nodes=Counter()
    nnodes=len(query)
    for q in query:
        nodes[q.id]=0
    edges=[]
    L = len(str(query[0].seq))
    for i in range(nnodes-1):
        q = query[i]
        for j in range(i+1,nnodes):
            r = query[j]
            newd = tn93(str(q.seq),str(r.seq),L,matchMode,minOverlap)
            if newd < threshold:
                nodes[q.id]+=1
                edges.append([q.id,r.id])
    return(nodes,edges)

# Not working for repeated sequences yet
def tn93_nx(query,threshold,matchMode,minOverlap):
    G=nx.Graph()
    nnodes=len(query)
    for q in query:
        G.add_node(q.id)
    L = len(str(query[0].seq))
    for i in range(nnodes-1):
        q = query[i]
        for j in range(i+1,nnodes):
            r = query[j]
            newd = tn93(str(q.seq),str(r.seq),L,matchMode,minOverlap)
            if newd < threshold:
                G.add_edge(q.id,r.id)
    return(G)

def cluster_ids(nodes):
    result = [x for x in nodes if nodes[x]>0]
    return(result)

def el_to_nx(el):
    lines=[" ".join(x) for x in el]
    G = nx.parse_edgelist(lines, nodetype = str)
    return(G)

def csd(el):
    c=Counter()
    G=el_to_nx(el)
    cc=list(nx.connected_components(G))
    for x in cc:
        c[x]+=1
    return(c)
