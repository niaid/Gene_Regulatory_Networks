# Programs to install,
geWorkbench. 
Download and install. (it needs email, name, affiliation for a free download)
http://wiki.c2b2.columbia.edu/workbench/index.php/Download_and_Installation#geWorkbench_2.6.0.3_.28released_December_21.2C_2015.29

make an alias of the “launch_geworkbench_macosx_16G.command” file from the unzipped folder to the desktop for easy deployment.

# Julia language installation,
https://julialang.org/downloads/
once installed, evoke Julia, 
install a few programs for the tutorial.

using Pkg

Pkg.add(“GraphPlot”)

Pkg.add(“LightGraphs”)

Pkg.add(“NetworkInference”)

Pkg.add(“Printf”)  
Pkg.add(“CSV”)

#### After installation, test to make sure there is no error messages.
using GraphPlot 

using LightGraphs

using NetworkInference

using Printf

using CSV

#### If everything is fine. Then exit using
quit()

In Cytoscape, install the iRegulon app.
Apps-->Apps Manager --> search for iRegulon -->click it to install it.


# Install the SCENIC package following this website bellow. 
This will take a while because there are 4 major pieces in it (Genie3, cisTarget, AUCell, and SCENIC) and lot of independencies.
Mainly just follow the “Installation” section. You can ignore the “Species-specific database”. 

I have downloaded them into the folder on the server.
https://rawcdn.githack.com/aertslab/SCENIC/701cc7cc4ac762b91479b3bd2eaf5ad5661dd8c2/inst/doc/SCENIC_Setup.html

I have transferred some data files to the server.\n  /home/bcbb_teaching_files/


