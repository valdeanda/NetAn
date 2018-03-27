# Network analysis
We provide a single script  that use the output derived from [Metamis sofware](https://www.ncbi.nlm.nih.gov/pubmed/27887570) sofware, which is a directed and weighted adyacence list , and compute several statistics 
```
Firmicutes      Ascomycota      -0.000272753971115013
Firmicutes      Chlorobi        -0.000242343947634859
Firmicutes      Ignavibacteriae 0.000178797125901064
Firmicutes      Dikarya unclassified    -0.00017738367107812
Firmicutes      Candidatus Beckwithbacteria     -0.000173505314274337
Firmicutes      Candidatus Woesebacteria        -0.00015914042503068
Firmicutes      Bacteroidetes/Chlorobi group unclassified       -0.000153790507711771
Firmicutes      Armatimonadetes 0.000145294003838928
Firmicutes      Cloacimonetes   -0.000130612871481698
Firmicutes      Candidatus Shapirobacteria      -0.00013061259484347
```

The main script :  [NetworkAnalysis.py](./scripts/NetworkAnalysis.py) recieves the above mentioned list and compute several topological features such as density, mean degree, hubs, connected components, clustering coefficient modularity etc.  

# Dependencies 

In order to use the script, the following dependencies must be installed first.
1. [ Python3  ](https://www.python.org/)
2. [Networkx](https://networkx.github.io/) 
3. [Community library](https://github.com/taynaud/python-louvain/) 

```
sudo -H pip install python-louvain
```

Some extra python libreries 

Numpy and matplotlib
```
sudo pip2 install numpy matplotlib
```
Librarpy tk 

```
sudo apt-get install python-tk
```

### Command line for directed networks   
```
python Networks/scripts/NetworkAnalysis.py -d Networks/data/a_phylum_consensus.txt
```
### Command line for undirected networks 
In the case of having and adyacency list with no directionality, you can use the option -u for undirected networks
```
python Networks/scripts/NetworkAnalysis.py -u Networks/data/a_phylum_consensus.txt
```

### Network output 

Once the program has finished, the following message will appears in the terminal
```
Analysing Network:    Networks/data/a_phylum_consensus.txt    type:   Directed
        - Parsing file to network ...
        - Obtaining order ...
        - Obtaining size ...
        - Obtaining diameter (for undir network version) ...
        - Obtaining radius (for undir network version) ...
        - Obtaining density ...
        - Obtaining mean degree ...
        - Obtaining max degree ...
        - Obtaining mean clusttering coefficient ...
        - Obtaining hubs ...
        - Obtaining mincut vertex set (for undir network version) ...
        - Obtaining mincut edge set (for undir network version) ...
        - Obtaining connected components ...
        - Obtaining strongly connected components (just for directed case) ...
        - Obtaining maximal cliques (for undir network version) ...
        - Obtaining cycle basis (for undir network version) ,,,
        - Obtaining maximal independent set ...
        - Obtaining degree distribution ...
        - Drawing in degree distribution ...
        - Drawing out degree distribution ...
        - Obtaining communities (for undir network version)...
        - Drawing networks ...
        - Running Random Network Analysis ...
        - Finished analysis for this network
```

## Analyzing the results  

If you did not have any problems with the dependencies, the following output files should have been generated

```
a_phylum_consensus_directed.txt
a_phylum_consensus_random_results.txt
a_phylum_consensus_distribution_indegree.png
a_phylum_consensus_distribution_outdegree.png
a_phylum_consensus_community_network.png
a_phylum_consensus_directed_network.png
```
Two  files contains the statistical analysis, one real and one generating 1000 random networks with the same number of nodes and edges.  
### Analyzing the random test  
```
less a_phylum_consensus_random_results.txt
```

```
Random Analysis Results, Mean Measures:                 runs(100)
 - [/] order:   61
 - [/] size:    3560
 - [u] diameter:        1.79
 - [u] radius:  1.0
 - [d] density: 0.972677595628
 - [d] mean degree:     58.0
 - [d] clustering coefficient:  0.499653314753
 - [d] maximum in degree:       60.0
 - [d] maximum out degree:      60.0
 - [d] hubs with max in degree: 11.4
 - [d] hubs with max out degree:        11.41
 - [u] modularity:      0.0
 - [u] number of communities:   1.0
 - [/] number of connected components:  1.0
 - [d] number of strongly connected components: 1.0
 - [u] number of cycles in cycle basis: 1768.73

[u] undirected associated graph
[d] directed graph
[/] both undirected and directed

```







