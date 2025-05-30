import numpy as np
from tqdm import tqdm, trange
from itertools import combinations
from collections import defaultdict
from math import comb
import pickle
import os

from load_graph import *


def advanced_method(endpoints, nodes_per_com, coms_per_node):
    print("Advanced 4 overlap method")
    icdf = np.zeros(len(endpoints), dtype="float64")
    sizes = np.array([len(i) for i in coms_per_node])
    most_coms = np.argsort(sizes)[::-1][:1000]
    most_coms_size_one_overlaps = {int(i):float(comb(sizes[i], 4)) for i in most_coms}
    most_coms = frozenset(most_coms)

    print("Not most coms nodes")
    for i in trange(len(coms_per_node)):
        if i in most_coms:
            continue
        for coms in combinations(coms_per_node[i], 4):
            intersection = frozenset.intersection(*[nodes_per_com[j] for j in coms])
            size = len(intersection)
            icdf[endpoints <= size] += 1/size
            # Add a bit if the nodes with the most coms are in this overlap
            # For each node in the overlap and not in most coms, add 1/(size * number of nodes in the overlap we will see in the iteration)
            these_four_and_most_coms = intersection.intersection(most_coms)
            for j in these_four_and_most_coms:
                icdf[endpoints <= size] += 1/(size*(len(intersection) - len(these_four_and_most_coms)))
                most_coms_size_one_overlaps[j] -= 1/(len(intersection) - len(these_four_and_most_coms))

    print("Nodes with most coms")
    # for pair of nodes in most_coms
    for i,j in tqdm(combinations(most_coms, 2), total=comb(len(most_coms), 2)):
        # get communities in common
        shared_coms = frozenset.intersection(coms_per_node[i], coms_per_node[j])
        # for every 4 community overlaps
        for x in combinations(shared_coms, 4):
            # get the overlap of these 4 communities
            overlap = frozenset.intersection(*[nodes_per_com[y] for y in x])
            # if there are any nodes not in most_coms we've already counted this overlap
            if not overlap.issubset(most_coms): 
                continue
            overlap_size = len(overlap)
            icdf[endpoints <= overlap_size] += 1/comb(overlap_size, 2) # we see this overlap (size choose 2) times

            most_coms_size_one_overlaps[i] -= 1/(overlap_size - 1) # we do this for i once for each other node in the overlap
            most_coms_size_one_overlaps[j] -= 1/(overlap_size - 1)

    s = sum(round(x) for x in most_coms_size_one_overlaps.values())
    icdf[0] += s
    return icdf / icdf[0]


def intersection_icdf(endpoints, nodes_per_com, coms_per_node, n_to_overlap):
    icdf = np.zeros(len(endpoints), dtype="float64")
    for node_coms in tqdm(coms_per_node):
        for coms in combinations(node_coms, n_to_overlap):
            size = len(frozenset.intersection(*[nodes_per_com[i] for i in coms]))
            icdf[endpoints <= size] += 1/size
    return icdf / icdf[0]


def overlaps_dict(endpoints, nodes_per_com, coms_per_node, use_advanced=True):
    data = dict()
    data["two_overlaps"] = intersection_icdf(endpoints, nodes_per_com, coms_per_node, 2)
    data["three_overlaps"] = intersection_icdf(endpoints, nodes_per_com, coms_per_node, 3)
    if use_advanced:
        data["four_overlaps"] = advanced_method(endpoints, nodes_per_com, coms_per_node)
    else:
        data["four_overlaps"] = intersection_icdf(endpoints, nodes_per_com, coms_per_node, 4)
    return data


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
    data["Youtube"] = overlaps_dict(endpoints, coms, g.vs["comms"])
    print("d=2")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_youtube_d2_com.dat", False)
    data["d=2"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=8")
    nodes_per_com, coms_per_node  = load_coms("data/abcdoo_youtube_d8_com.dat", False)
    data["d=8"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=64")
    nodes_per_com, coms_per_node  = load_coms("data/abcdoo_youtube_d64_com.dat", False)
    data["d=64"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("CBK")
    nodes_per_com, coms_per_node  = load_coms("data/cbk_youtube.dat", False)
    data["CBK"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
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
    data["DBLP"] = overlaps_dict(endpoints, coms, g.vs["comms"])
    print("d=2")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_dblp_d2_com.dat", True)
    data["d=2"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=8")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_dblp_d8_com.dat", True)
    data["d=8"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=64")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_dblp_d64_com.dat", True)
    data["d=64"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=CBK")
    nodes_per_com, coms_per_node = load_coms("data/cbk_dblp.dat", True)
    data["CBK"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
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
    data["Amazon"] = overlaps_dict(endpoints, coms, g.vs["comms"])
    print("d=2")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_amazon_d2_com.dat", True)
    data["d=2"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=8")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_amazon_d8_com.dat", True)
    data["d=8"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=64")
    nodes_per_com, coms_per_node = load_coms("data/abcdoo_amazon_d64_com.dat", True)
    data["d=64"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    print("d=CBK")
    nodes_per_com, coms_per_node = load_coms("data/cbk_amazon.dat", True)
    data["CBK"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    with open(file, "wb") as f:
        pickle.dump(data, f)