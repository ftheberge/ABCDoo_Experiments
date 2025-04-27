import numpy as np
from tqdm import tqdm
from itertools import combinations
import pickle
import os

from load_graph import *


def intersection_icdf(endpoints, nodes_per_com, coms_per_node, n_to_overlap):
    icdf = np.zeros(len(endpoints), dtype="float64")
    for node_coms in tqdm(coms_per_node):
        for coms in combinations(node_coms, n_to_overlap):
            size = len(frozenset.intersection(*[nodes_per_com[i] for i in coms]))
            icdf[endpoints <= size] += 1/size
    return icdf / np.sum(icdf)


def overlaps_dict(endpoints, nodes_per_com, coms_per_node, do_four=True):
    data = dict()
    data["two_overlaps"] = intersection_icdf(endpoints, nodes_per_com, coms_per_node, 2)
    data["three_overlaps"] = intersection_icdf(endpoints, nodes_per_com, coms_per_node, 3)
    if do_four:
        data["four_overlaps"] = intersection_icdf(endpoints, nodes_per_com, coms_per_node, 4)
    return data


# Right endpoints for all graphs
endpoints = np.logspace(0, 5, 200) # 10^5 is larger than the largest overlap

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
    data["d=64"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node, do_four=False)
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
    data["d=64"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node, do_four=False)
    print("d=CBK")
    nodes_per_com, coms_per_node = load_coms("data/cbk_amazon.dat", True)
    data["CBK"] = overlaps_dict(endpoints, nodes_per_com, coms_per_node)
    with open(file, "wb") as f:
        pickle.dump(data, f)