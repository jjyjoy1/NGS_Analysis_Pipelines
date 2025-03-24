# Install the package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("curatedMetagenomicData")

# Load the package
library(curatedMetagenomicData)

# See available datasets (includes HMP)
availableDatasets()

# Get HMP data specifically
hmp_datasets <- availableDatasets()[grep("HMP", availableDatasets())]
print(hmp_datasets)

# Load a specific dataset
hmp_data <- curatedMetagenomicData("HMP_2012.metaphlan_bugs_list.stool", dryrun = FALSE)



