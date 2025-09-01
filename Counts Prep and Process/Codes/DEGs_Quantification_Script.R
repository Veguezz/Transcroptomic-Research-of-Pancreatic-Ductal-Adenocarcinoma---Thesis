




df <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/res_df_PanCuRx_V1.rds")

UpGPanC <- df[df$log2FoldChange > 0 & df$log2FoldChange < 1 & df$pvalue < 0.05, ]

DownGPanC <- df[df$log2FoldChange > -1 & df$log2FoldChange < 0 & df$pvalue < 0.05, ]


dfT <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.rds") 

UpGTCGA <- dfT[dfT$log2FoldChange > 0 & dfT$log2FoldChange < 1 & dfT$pvalue < 0.05, ]

DownGTCGA <- dfT[dfT$log2FoldChange > -1 & dfT$log2FoldChange < 0 & dfT$pvalue < 0.05, ]






