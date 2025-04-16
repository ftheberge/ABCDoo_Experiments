# ABCDoo_Experiments
Experiments to generate plots for ABCDoo paper(s)

The experiments are found in various directories:

## GraphProperties
To get started, ensure there is a folder titled ```GraphProperties/data```.
In this folder put the DBLP, AMAZON and YOUTUBE graph and community files from the [SNAP](https://snap.stanford.edu/data/#communities) website.
The folder should look like:

data/
- com-amazon.all.dedup.cmty.txt
- com-amazon.ungraph.txt
- com-dblp.all.cmty.txt
- com-dblp.ungraph.txt
- com-youtube.all.cmty.txt
- com-youtube.ungraph.txt

To ensure local file paths work properly, be sure to run each script from within the GraphProperties folder.
Run the julia file ```abcdoo_snap_graph_sampler.jl``` and the python file ```CKB.py```.
This will generate ABCD+O^2 graphs and CKB communities using parameters measured from the real graphs and save them in the data folder.
The parameters are hard-coded into the files.
We measured the parameters using ```measure_params.ipynb```.

At this point the majority of the experiments can be run.
Below are the notebooks listed in the order that the figures they produce appear.

Experiments in Section 3.1
- ```community_size.ipynb```
- ```communities_per_node.ipynb```
- ```community_overlap.ipynb```

Experiments in Section 3.2
- ```rho.ipynb```
- ```degree.ipynb``` (this may get / have been cut)
- ```overlap_density.ipynb```
- ```IEF.ipynb```

The experiment compaing rho to overlap density requires new data.
Run ```abcdoo_rho_sampler.jl``` (again the parameters have been hard-coded).
Then use the notebook ```rho_vs_overlap_density.ipynb```.

## CompareAlgorithms

Code to compare clustering algorithms over several ABCDoo graphs with varying parameters.