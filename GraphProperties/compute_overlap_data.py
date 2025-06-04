import numpy as np
from tqdm import tqdm, trange
from itertools import combinations
from collections import defaultdict
from math import comb
import pickle
import os

from load_graph import *


def make_icdfs(endpoints, nodes_per_com, coms_per_node, n_most_to_cut):
    print("Two Overlaps")
    two_overlaps = defaultdict(set)
    for i, node_coms in tqdm(enumerate(coms_per_node), total=len(coms_per_node), mininterval=0.5):
        for x in combinations(node_coms,2):
            two_overlaps[tuple(sorted(x))].add(i)
    
    icdf2 = np.zeros(len(endpoints), dtype="float64")
    for overlap in two_overlaps.values():
        icdf2[endpoints <= len(overlap)] += 1
    icdf2 = icdf2 / icdf2[0]

    print("Three Overlaps")
    icdf3 = np.zeros(len(endpoints), dtype="float64")
    with tqdm(total=sum(comb(len(i),3) for i in coms_per_node), mininterval=0.5) as pbar:
        for node_coms in coms_per_node:
            for coms in combinations(node_coms, 3):
                size = len(set.intersection(*[two_overlaps[tuple(sorted(x))] for x in combinations(coms,2)]))
                icdf3[endpoints <= size] += 1/size
                pbar.update()
    icdf3 = icdf3 / icdf3[0]

    print("Four Overlaps")
    icdf4 = np.zeros(len(endpoints), dtype="float64")
    sizes = np.array([len(i) for i in coms_per_node])
    most_coms = np.argsort(sizes)[::-1][:n_most_to_cut]
    most_coms_size_one_overlaps = {int(i):float(comb(sizes[i], 4)) for i in most_coms}
    most_coms = frozenset(most_coms)
    print("Nodes without most coms")
    n_intersections_to_compute = 0
    for i in range(len(coms_per_node)):
        if i not in most_coms:
            n_intersections_to_compute += comb(len(coms_per_node[i]),4)
    with tqdm(total=n_intersections_to_compute, mininterval=0.5) as pbar:
        for i in range(len(coms_per_node)):
            if i in most_coms:
                continue
            for coms in combinations(coms_per_node[i], 4):
                intersection = set.intersection(*[two_overlaps[tuple(sorted(x))] for x in combinations(coms,2)])
                size = len(intersection)
                icdf4[endpoints <= size] += 1/size
                # Add a bit if the nodes with the most coms are in this overlap
                # For each node in the overlap and not in most coms, add 1/(size * number of nodes in the overlap we will see in the iteration)
                these_four_and_most_coms = intersection.intersection(most_coms)
                for j in these_four_and_most_coms:
                    icdf4[endpoints <= size] += 1/(size*(len(intersection) - len(these_four_and_most_coms)))
                    most_coms_size_one_overlaps[j] -= 1/(len(intersection) - len(these_four_and_most_coms))
                pbar.update()

    print("Nodes with most coms")
    # for pair of nodes in most_coms
    shared_coms = dict()
    for i,j in combinations(most_coms, 2):
        # get communities in common
        shared_coms[(i,j)] = frozenset.intersection(coms_per_node[i], coms_per_node[j])
    n_intersections_to_compute = 0
    for i in shared_coms.values():
        n_intersections_to_compute += comb(len(i), 4)
    with tqdm(total=n_intersections_to_compute, mininterval=0.5) as pbar:
        for (i,j), shared in shared_coms.items():
            # for every 4 community overlaps
            for x in combinations(shared, 4):
                # get the intersection of these 4 communities
                intersection = set.intersection(*[two_overlaps[tuple(sorted(y))] for y in combinations(x,2)])
                # if there are any nodes not in most_coms we've already counted this overlap
                pbar.update()
                if not intersection.issubset(most_coms): 
                    continue
                size = len(intersection)
                icdf4[endpoints <= size] += 1/comb(size, 2) # we see this overlap (size choose 2) times

                most_coms_size_one_overlaps[i] -= 1/(size - 1) # we do this for i once for each other node in the overlap
                most_coms_size_one_overlaps[j] -= 1/(size - 1)

    s = sum(round(x) for x in most_coms_size_one_overlaps.values())
    icdf4[0] += s
    icdf4 = icdf4 / icdf4[0]

    return icdf2, icdf3, icdf4

# Right endpoints for all graphs
endpoints = np.logspace(0, 5, 200) # 10^5 is larger than the largest overlap
endpoints = np.unique(np.round(endpoints)) # Only need integer endpoints


# Youtube
file = "data/youtube_overlap_data.pkl"
if os.path.isfile(file):
    print(f"Overlap data already exists. If you want to regenereate please delete/rename the file {file}")
else:
    data = dict()
    data["endpoints"] = endpoints

    print("Youtube")
    g, coms = load_snap("data/com-youtube.ungraph.txt", "data/com-youtube.all.cmty.txt", drop_outliers=True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, coms, g.vs["comms"], 300)
    data["Youtube"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("d=2")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_youtube_d2_com.dat", False)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 225)
    data["d=2"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("d=8")
    nodes_per_com, coms_per_node  = load_coms("data/abcdoo_youtube_d8_com.dat", False)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 500)
    data["d=8"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("CBK")
    nodes_per_com, coms_per_node  = load_coms("data/cbk_youtube.dat", False)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 600)
    data["CBK"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    with open(file, "wb") as f:
        pickle.dump(data, f)
    print("Saved everything except d=64")

    print("d=64")
    nodes_per_com, coms_per_node  = load_coms("data/abcdoo_youtube_d64_com.dat", False)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 600)
    data["d=64"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    with open(file, "wb") as f:
        pickle.dump(data, f)


# DBLP
file = "data/dblp_overlap_data.pkl"
if os.path.isfile(file):
    print(f"Overlap data already exists. If you want to regenereate please delete/rename the file {file}")
else:
    data = dict()
    data["endpoints"] = endpoints

    print("DBLP")
    g, coms = load_snap("data/com-dblp.ungraph.txt", "data/com-dblp.all.cmty.txt")
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, coms, g.vs["comms"], 400)
    data["DBLP"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("d=2")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_dblp_d2_com.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 225)
    data["d=2"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("d=8")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_dblp_d8_com.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 50)
    data["d=8"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("CBK")
    nodes_per_com, coms_per_node = load_coms("data/cbk_dblp.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 600)
    data["CBK"] = {"two":icdf2, "three":icdf3, "four":icdf4}
    
    with open(file, "wb") as f:
        pickle.dump(data, f)
    print("Saved everything except d=64")

    print("d=64")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_dblp_d64_com.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 325)
    data["d=64"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    with open(file, "wb") as f:
        pickle.dump(data, f)


# Amazon
file = "data/amazon_overlap_data.pkl"
if os.path.isfile(file):
    print(f"Overlap data already exists. If you want to regenereate please delete/rename the file {file}")
else:
    data = dict()
    data["endpoints"] = endpoints

    print("Amazon")
    g, coms = load_snap("data/com-amazon.ungraph.txt", "data/com-amazon.all.dedup.cmty.txt")
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, coms, g.vs["comms"], 25)
    data["Amazon"] = {"two":icdf2, "three":icdf3, "four":icdf4}
    print("d=2")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_amazon_d2_com.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 25)
    data["d=2"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("d=8")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_amazon_d8_com.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 50)
    data["d=8"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    print("d=CBK")
    nodes_per_com, coms_per_node = load_coms("data/cbk_amazon.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 600)
    data["CBK"] = {"two":icdf2, "three":icdf3, "four":icdf4}
    
    with open(file, "wb") as f:
        pickle.dump(data, f)
    print("Saved everything except d=64")

    print("d=64")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_amazon_d64_com.dat", True)
    icdf2, icdf3, icdf4 = make_icdfs(endpoints, nodes_per_com, coms_per_node, 500)
    data["d=64"] = {"two":icdf2, "three":icdf3, "four":icdf4}

    with open(file, "wb") as f:
        pickle.dump(data, f)