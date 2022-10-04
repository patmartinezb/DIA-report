## A tool for generating reports from DIA proteomics data

This tool launches a local Shiny-based web application that enables an overall analysis of a DIA proteomics data matrix (from FragPipe DIA-NN). It provides an interactive interface, easy to use for those users with little or no experience in the R programming language, and outputs an HTML report, a .pptx file with the figures, and an Excel file with different types of data.


## **Requirements for the data matrix file**

An Excel (.xlsx) or .tsv file.

The data matrix should be the output of FragPipe's DIA-NN, which means that it has to include the following annotation columns: "Protein.Group", "Protein.Ids", "Protein.Names", "Genes" and "First.Protein.Description"; and the columns with the **actual data** (with their names matching the `BioReplicate` column from the Metadata file).

## **Requirements for the metadata file**

An Excel (.xlsx) file.

For the app to work, the following variables, with those *exact* same names (capital letters included), have to be present:

-   **BioReplicate**: each sample within each condition - it should be copy-pasted from the column names of the raw data matrix file

-   **Group**: which group the sample belongs to - the names should be informative and always the same (i.e., GAPDH in all instances, not variations like Gapdh)


## **Requirements for the comparison file**

An Excel (.xlsx) file.

For the app to work, the following variables, with those *exact* same names (capital letters included), have to be present:

-   **Condition**: name of the group that is acting as treatment/mutated, or simply one of the groups

-   **Control**: name of the group that is acting as control/wildtype, or simply one of the groups

**The names given to the Condition and Control must be written exactly the same as in the Group column of the Metadata file**


> For the app to work, it is **mandatory** that all headers of all files are written correctly, that the BioReplicates are written (copy-pasted) as in the data matrix file, and that de Control-Condition names are the same as in the Group variable of the Metadata file.