#librerías
library(tidyverse) #incluye ggplot2 para visualización de datos, y dplyr para manipulación de los datos.
library (RColorBrewer) #Para ampliar nuestra gama de colores
library(ggrepel) #Para lindas anotaciones

res_dfT <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.2.rds") 

genes_interest <- c(
  "CD3D", "GRAP", "ITGAL", "LTB", "LYVE1", "NCF1", 
  "PLA2G2A", "PLA2G2D", "TNFRSF9", "GBP6", 
  "IL20RB", "MME", "WNT11")

# Subset de expresión sólo para esos genes
Genes <- res_dfT[genes_interest, , drop = FALSE]






#Primero ponemos de menor a mayor en función al p-value
res_dfT <- res_dfT[order(res_dfT$pvalue), ]

#Tema
theme_set(theme_classic(base_size = 20)+
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), colour = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), colour = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

#Graficar el esqueleto principal
ggplot(data = res_dfT,aes(x=log2FoldChange, y= -log10(pvalue))) +
  geom_vline(xintercept = c(-0.6,0.6), col='gray', linetype='dashed') + #para crear treshold lines que me separe cuando considero up and down regulated
  geom_hline(yintercept = 1.3, col='gray', linetype='dashed') +
  geom_point()

#Para poder diferenciar por color los Up y Down regulated primero debemos clasificarlos creando una nueva columna
res_dfT$expresiondiferencial <- 'No' 
res_dfT$expresiondiferencial[res_dfT$log2FoldChange > 1 & res_dfT$pvalue < 0.05] <- 'UP'
res_dfT$expresiondiferencial[res_dfT$log2FoldChange < -1 & res_dfT$pvalue < 0.05] <- 'DOWN'

#Graficar el esqueleto principal costumizado
ggplot(data = res_dfT,aes(x=log2FoldChange, y= -log10(pvalue), colour = expresiondiferencial)) +
  geom_vline(xintercept = c(-1,1), col='gray', linetype='dashed') + #para crear treshold lines que me separe cuando considero up and down regulated
  geom_hline(yintercept = 1.3, col='gray', linetype='dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("royalblue","grey","darkred"),
                     labels = c("Genes Subexpresados", "No Significativo", "Genes Sobrexpresados")) +
  coord_cartesian(ylim = c(0,7), xlim = c(-5,5)) + #Ahora los ejes y los limites
  scale_x_continuous(breaks = seq(-10,10,2)) +
  labs(color='Expresión Diferencial', #Titulo leyenda
       x= expression("log"[2]*"FC"), y= expression("-log"[10]*"p-value")) + #Ponemos el nombre de los ejes
  ggtitle('Expresión Diferencial entre Pacientes Respondedores y No Respondedores')  #Ponemos el titulo



#Añadir etiquetas para los genes significativos usando 'delabel'

#Seleccionamos los primeros 20 genes según el p-value 

#top50seg <- head(res_dfT[order(res_dfT$pvalue, na.last = NA), "gene"], 50)

#top100seg <- head(res_dfT[order(res_dfT$pvalue, na.last = NA), "gene"], 100)

#res_dfT$delabel <- ifelse(res_dfT$gene %in% top100seg, res_dfT$gene, NA)

genes_interes <- c("LTB", "CD3D", "GRAP", "TNFRSF9", "PLA2G2A", "PLA2G2D","NCF1", "ITGAL", "LYVE1", 
                   "GBP6", "MME", "IL20RB", "WNT11")

res_dfT$delabel <- ifelse(res_dfT$hgnc_symbol %in% genes_interes, res_dfT$hgnc_symbol, NA)




#Graficar el esqueleto principal costumizado nuevamente pero ahora incluyendo "delabel"
volcanoplot <- ggplot(data = res_dfT,aes(x=log2FoldChange, y= -log10(pvalue), colour = expresiondiferencial, label= delabel)) +
  geom_vline(xintercept = c(-1,1), col='gray', linetype='dashed') + #para crear treshold lines que me separe cuando considero up and down regulated
  geom_hline(yintercept = 1.3, col = 'gray', linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("royalblue","grey","darkred"),
                     labels = c("Genes Subexpresados", "Genes No Significativos", "Genes Sobrexpresados")) +
  coord_cartesian(ylim = c(0,10), xlim = c(-10,10)) + #Ahora los ejes y los limites
  scale_x_continuous(breaks = seq(-8,8,2)) +
  labs(color='Expresión Diferencial', #Titulo leyenda
       x= expression("log"[2]*"FC"), y= expression("-log"[10]*"p-value")) +  #Ponemos el nombre de los ejes
  ggtitle('Expresión Diferencial entre pacientes Respondedores y No Respondedores TCGA')  #Ponemos el titulo

#Agregamos etiquetas
library(extrafont)

volcano_plot <- volcanoplot +
  geom_text_repel(
    data = subset(res_dfT, !is.na(delabel)),
    aes(label = delabel),
    color = "black",
    size = 3.5,
    max.overlaps = Inf,
    force = 1,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0.1
  ) +
  theme(text = element_text(family = "Times New Roman"))

dev.new()
print(volcano_plot)




#Si quisiera un listado de los genes upregulated y downregulated
library(openxlsx)

UpregulatedT <- res_dfT$hgnc_symbol[res_dfT$expresiondiferencial == "UP"]

Upregulated_dfT <- data.frame(Gene = UpregulatedT)

write.xlsx(Upregulated_dfT, "Upregulated_genesT.xlsx")

DownregulatedT <- res_dfT$hgnc_symbol[res_dfT$expresiondiferencial == "DOWN"]

Downregulated_dfT <- data.frame(Gene = DownregulatedT)

write.xlsx(Downregulated_dfT, "Downregulated_genesT.xlsx")

