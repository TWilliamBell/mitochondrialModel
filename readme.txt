This repository includes a public-facing version of the Layton Lab's mitochondrial model, as developed for hepatocytes, proximal tubule cells, and medullary thick ascending limb cells.

There are two subdirectories, one for hepatocytes, and another for the nephron cell types.  Inside each directory is the necessary directories for any script to run smoothly, including a directory to output plots to (dataVis), a directory to output data to (results), a directory containing all the scripts used in studying the given tissue (modelScripts), and a directory containing the bare minimum necessary files for running the model (minimumWorkingModel).

I (Will) suggest starting with minimumWorkingModel, which includes better commenting and allows one to see the crucial components of the model efficiently.  The directory modelScripts includes more scripts, many performing very specific functions, such as outputting the results necessary for a single plot.  The results are designed to be a reproducible as possible, requiring no more than perhaps uncommenting a certain section of code in order to reproduce any work done for this project.

This model is adopted in large part from a model developed by Wu et al. in which Beard Lab for cardiac tissue.  
