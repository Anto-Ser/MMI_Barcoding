## --- Stacked Histograms and bubble plots----

rm(list=setdiff(ls(), c("df_MI", "MI_info", "Number_barcode", "Bubble_plot", "Stacked_histo", "Select_and_order_Organs")))

## --- Figure 1 c) Bubble plot for MDA tumours ----

Col.MDA.bbp = read.csv("./data/Y.Position_Colours_MDA_bubbles.csv")$Colours
Y.MDA.bbp = read.csv("./data/Y.Position_Colours_MDA_bubbles.csv")$Y.Position

MDA.tum = Select_and_order_Organs(df_MI, "MDA-231", "Tum", info_dataframe = Number_barcode)
names(MDA.tum) =str_split(names(MDA.tum), "_", simplify = TRUE)[,2] # select only mice ID
MDA.tum = MDA.tum %>% select(sort(names(MDA.tum)),everything()) # Sort by ID

MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x20

## --- Supp Figure 2 d) Bubble plot MDA tumours by experiments ----
MDA.tum = Select_and_order_Organs(df_MI, "MDA-231", "Tum", info_dataframe = Number_barcode)
MDA.tum.exp1= MDA.tum[,grep("Exp38_", names(MDA.tum))]
MDA.tum.exp2= MDA.tum[,grep("Exp1009_", names(MDA.tum))]
MDA.tum.exp3= MDA.tum[,grep("Exp1011_", names(MDA.tum))]

MDA.tum.bbp_1 = Bubble_plot(MDA.tum.exp1, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp_2 = Bubble_plot(MDA.tum.exp2, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp_3 = Bubble_plot(MDA.tum.exp3, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)

MDA.tum.bbp_1$BBP #5x8
MDA.tum.bbp_2$BBP #5x11
MDA.tum.bbp_3$BBP#5x8

## --- Supp Figure 1 i) Stack histogram PDX-412 tumours ----

PDX.tum = Select_and_order_Organs(df_MI, "PDX-T412", "Tum", info_dataframe = Number_barcode)
PDX.color = read.csv("./data/Colour_PDX_Histo.csv")

names(PDX.tum) =str_split(names(PDX.tum), "_", simplify = TRUE)[,2] # select only mice ID
PDX.tum = PDX.tum %>% select(sort(names(PDX.tum)),everything()) # Sort by ID
PDX.tum.histo = Stacked_histo(PDX.tum, angle.x = 90, my_col = PDX.color$x)

PDX.tum.histo$Stack.plot #5x12

## --- Supp Figure 7 Bubble plot MDA-231 organs and PDX-412 organs ----

MDA.tum = Select_and_order_Organs(df_MI, "MDA-231", "CTC", info_dataframe = Number_barcode)
names(MDA.tum) =str_split(names(MDA.tum), "_", simplify = TRUE)[,2] # select only mice ID
MDA.tum = MDA.tum %>% select(sort(names(MDA.tum)),everything()) # Sort by ID

MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x20

# Code rerun with "Lung", "Liver", "CTCs", could be transform in loop or function to generate everything in one page.


## --- Bubble plot for all organs grouped for PDX ----

Col.MDA.bbp = read.csv("./data/Y.Position_Colours_PDX_bubbles_v6.csv")$Colours
Y.MDA.bbp = read.csv("./data/Y.Position_Colours_PDX_bubbles_v6.csv")$Y.Position

MDA.tum = Select_and_order_Organs(df_MI, Model = "PDX-T412", Organ = "Lung", info_dataframe = Number_barcode)
names(MDA.tum) =str_split(names(MDA.tum), "_", simplify = TRUE)[,2] # select only mice ID
MDA.tum = MDA.tum %>% select(sort(names(MDA.tum)),everything()) # reorder by mice ID

MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x15

MDA.tum = Select_and_order_Organs(df_MI, "PDX-T412", "CTC", info_dataframe = Number_barcode)
names(MDA.tum) =str_split(names(MDA.tum), "_", simplify = TRUE)[,2] # select only mice ID
MDA.tum = MDA.tum %>% select(sort(names(MDA.tum)),everything()) # reorder by mice ID
MDA.tum.bbp = Bubble_plot(MDA.tum, angle.x = 90, Y = "input", data2=Y.MDA.bbp, my_col = Col.MDA.bbp)
MDA.tum.bbp$BBP #6x8


## --- Specific examples stacked histograms ----

Examples_MDA =c("022","005","030", "036") # Figure 3
Examples_PDX =c("056","047","066", "077") # Figure 3

pdf(file= "./Examples_histo.pdf", height = 5, width = 5)
for (i in Examples_MDA) {

  Ms_Tum_exp = df_MI[,grep(paste("_", i, "_", sep = ""), names(df_MI)), drop=F]

  colnames(Ms_Tum_exp) <- str_split(names(Ms_Tum_exp), "_", simplify = T)[,3]
  Ms_Tum_exp = Ms_Tum_exp[,c("Tum", "Lung", "Liver", "CTC")]

  p2 = Stacked_histo(Ms_Tum_exp, my_col = MDA.color$x)

  print(p2$Stack.plot + labs(title = paste("Mouse",i)))#4x4(MDA) #4x3(PDX)
}
dev.off()

Examples_MDA_Organs =c("019","008","030", "036","083") # Figure 4

pdf(file= "./Examples_MDA_Organs_histo.pdf", height = 5, width = 5)
for (i in Examples_MDA_Organs) {

  Ms_Tum_exp = df_MI[,grep(paste("_", i, "_", sep = ""), names(df_MI)), drop=F]

  colnames(Ms_Tum_exp) <- str_split(names(Ms_Tum_exp), "_", simplify = T)[,3]
  Ms_Tum_exp = Ms_Tum_exp[,c("Lung", "CTC", "Liver")]

  p2 = Stacked_histo(Ms_Tum_exp, my_col = MDA.color$x)

  print(p2$Stack.plot + labs(title = paste("Mouse",i)))#4x4(MDA) #4x3(PDX)
}
dev.off()

