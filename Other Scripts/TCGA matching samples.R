#Vamos a ver cuál es la muestra control de todas las muestras que tenemos


######Debido  a que desconocemos cual es la muestra control de nuestros datos accederemos a los datos del TCGA para matxh su metadata con nuestros read counts.
metadata_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/RNAseq/PAAD count and metadata.xlsx"

counts_data <- read_excel(metadata_path, sheet = "COUNTS ")

counts_data <- read.csv(metadata_file, stringsAsFactors = FALSE)

#Observa la estrcutura del metadata
print(head(counts_data))

# Check the column names of your count matrix and the sample identifiers in metadata
print(colnames(count_matrix))

print(counts_data$Column1)
