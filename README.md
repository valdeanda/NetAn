# Networks
We provide a workflow that use the output derived from [Metamis sofware](https://www.ncbi.nlm.nih.gov/pubmed/27887570) sofware, which is a directed and weighted adyacence list: 
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


#Getting the statistics by using the NetworkAnalysis.py script 

```
python Networks/scripts/NetworkAnalysis.py -d Networks/data/a_phylum_consensus.txt
```

