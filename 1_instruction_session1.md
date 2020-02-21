I, use ARACNE 
1.	Install geWorkbench, create a shortcut on your desktop.
2.	Click on the shortcut to evoke the geWorkbench. 
3.	Now you should see the CommendsAnalysisARACNE analysis.
4.	Open data file of Bcell-100.exp, it will ask an annotation file, choose the HG…csv file. This way the gene symbols will be added to the probe list, with a “()” separating it from other information. 
5.	Now select a few important genes that you want to find out their targets. Myc, Ccnd1, E2F1 etc. 
6.	Use ARACNE analysis hub markers (From Sets selection, the probes will appear in the box). Mode (Discovery), Algorithm (Adaptive Partitioning), Threshold (P-Value, 1.e-7, no Correction). DPI Tolerance (Do not apply), Bootstrap number (1). Click Analyze. This will result in a network, that will be open in the Cytoscape inside the window of geWorkbench. 
7.	From the Cytoscape export the graph as “MycCcnd1E2F1.xgmml”.
8.	Now you can start Cytoscape and import network from file, using the file named “MycCcnd1E2F1.xgmml”, and explore it.


II, Run Julia and NetworkInference
1.	### Click Julia icon to start up
2.	### Start up Jupyter notebook for Julia
using IJulia 
notebook()
### a terminal will pop out and start a URL with Jupyter notebook. 
3.	### find the 2_infer_network_ipynb from notebook folder and click it to open a worksheet. 
4.	### Run the steps in the worksheet. It will export the network and genes file to a csv file called 2_graphJulia.csv and genes list to 2_genesJulia.csv.
5.	### Run 2_JuliaGraphToCSV.R to give gene name to the network file and save it as 2_JuliaNet.csv
6.	### Import 2_JuliaNet.csv to the Cytoscape, by 
File Import  network from file 2_JuliaNet.csv
Use ““” and “,” as separate, unselect all and reselect “from” and “to” from the file. 
7.	### Change layout and use Mcode to explore the different levels of threshold on the connectivity.

