## --- Data processing for scatter plots -----

rm(list=setdiff(ls(), c("df_MI", "MI_info", "Number_barcode", "Bubble_plot", "Stacked_histo", "Select_and_order_Organs")))

df_scatter = data.frame() #Empty df to paste all the results

Mouse.ID = unique(str_split(names(df_MI), "_", simplify = T)[,2]) # Mouse ID for the loop
# Filtering out all invitro data for scatter plots on organs:
Mouse.ID = Mouse.ID[!Mouse.ID %in% c("Ctrl", "Orga",   "vitroA", "vitroB", "vitroC", "vitroD",
                                     "Amp1","Amp2","Amp3","Amp4","cltrinj","cltrvitro","25k")]

for (i in Mouse.ID) { # Loop to extract barcodes in Tum, lung, liver, CTC per mouse when barcode present

  TEST_name = paste("_", i, "_", sep="") # Create the name of the mouse: _XXX_
  Ms_x_Tum = df_MI[,grep(TEST_name, names(df_MI)),drop=F] # Select all samples for this mouse in data frame
  colnames(Ms_x_Tum) <- str_split(names(Ms_x_Tum), "_", simplify = T)[,3] # Simplify name to only organs

  Orga_MI = c("Tum", "Lung", "Liver", "CTC")
  Ms_df_temp = data.frame(ID=1:2608, Tum=NA, Lung=NA, Liver=NA, CTC=NA) #Create a temporary data frame for the loop with first column as ID barcode

  for (j in Orga_MI) {
    if(sum(names(Ms_x_Tum) == j)==0){ # if no barcodes in organs column (no colunm = FALSE => sum() == 0) => values = NA
      Ms_df_temp[,j] = NA
    } else { # if barcodes in organ column all values pasted
      Ms_df_temp[,j] = Ms_x_Tum[,j]
    }
  }

  columns_to_check = 1:5

  # Filtering for rows containing barcodes (rowSums >0) before merging, rowSums on all rows minus ID row and NA
  Ms_df_temp = Ms_df_temp[rowSums(Ms_df_temp[,columns_to_check[-c(1,which(is.na(colSums(Ms_df_temp))==TRUE))],drop = FALSE])>0,]
  Ms_df_temp$Mouse = i # add the mouse ID

  df_scatter = rbind(df_scatter, Ms_df_temp) #combine each iteration of the loop

}

df_scatter= left_join(df_scatter,MI_info, by="Mouse")

df_scatter$facet = factor(df_scatter$MI, levels = c("IMFP", "IMFPc", "SC", "ID", "IV"))


## --- Scatter Tum vs Lung Log -----
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df_scatter,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=Mouse), alpha=0.9)+
  #geom_quasirandom(data= subset(df_scatter,Tum==0&Lung>0),aes(x=0.0001, y=Lung/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df_scatter,Tum>0&Lung==0),aes(x=Tum/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df_scatter,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_grid(Cells~facet)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  labs(title = "Scatter plot Tumour vs Lung")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Lung (%)")#5x11

df_scatter %>% filter(Tum>0 & Lung>0) %>% select(Mouse, Exp, Cells, facet) %>% group_by(Cells, facet) %>%
  dplyr::summarise(n_Ms = length(unique(Mouse)), n_Exp =length(unique(Exp)))

## --- Scatter MDA Tumour VS Lung, Liver, CTC Log scale -----

#Tum vs Lung:

ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= subset(df_scatter,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=Mouse), alpha=0.9)+
  #geom_quasirandom(data= subset(df_scatter,Tum==0&CTC>0),aes(x=0.0001, y=CTC/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= subset(df_scatter,Tum>0&CTC==0),aes(x=Tum/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= subset(df_scatter,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_grid(Cells~facet)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  labs(title = "Scatter plot Tumour vs CTC")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in CTC (%)")#5x11


df_scatter %>% filter(Cells == "MDA-231") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= ~subset(.,Tum==0&Lung>0),aes(x=0.0001, y=Lung/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= ~subset(.,Tum>0&Lung==0),aes(x=Tum/10000 , y=0.00005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=10),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Lung")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Lung (%)")# 4x14

#Tum vs Liver
df_scatter %>% filter(Cells == "MDA-231") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Tum>0&Liver>0),aes(x=Tum/10000, y=Liver/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= ~subset(.,Tum==0&Liver>0),aes(x=0.0001, y=Liver/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= ~subset(.,Tum>0&Liver==0),aes(x=Tum/10000 , y=0.00005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Tum>0&Liver>0),aes(x=Tum/10000, y=Liver/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=10),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Liver")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Liver (%)")#4x14

#Tum vs CTC
df_scatter %>% filter(Cells == "MDA-231") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= ~subset(.,Tum==0&CTC>0),aes(x=0.0001, y=CTC/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= ~subset(.,Tum>0&CTC==0),aes(x=Tum/10000 , y=0.00005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x, n=10),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs CTC")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in CTC (%)")#4x14


## --- MDA Scatter plots IV and IMFP individual mouse (No Log) ----

scatter.IV = df_scatter %>% filter(Cells == "MDA-231" & MI == "IV") %>%
  ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse), alpha=0.9)+
  #geom_quasirandom(data= ~subset(.,Lung==0&Liver>0),aes(x=0.0005, y=Liver/10000, color=Mouse ),width=0.1,groupOnX=TRUE,alpha=0.7)+
  #geom_quasirandom(data= ~subset(.,Lung>0&Liver==0),aes(x=Lung/10000 , y=0.0005, color=Mouse ),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse),method = "lm", se=F)+
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #facet_grid(Cells~MI)+
  #scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  #scale_color_grey()+
  scale_x_continuous(limits = c(0,100))+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  labs(title = "Scatter plot Lung vs Liver")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in Liver (%)")#5x6

scatter.IMFP = df_scatter %>% filter(Cells == "MDA-231" & MI == "IMFP") %>%
  ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse), alpha=0.9)+
  #geom_quasirandom(data= ~subset(.,Lung==0&Liver>0),aes(x=0.0005, y=Liver/10000, color=Mouse ),width=0.1,groupOnX=TRUE,alpha=0.7)+
  #geom_quasirandom(data= ~subset(.,Lung>0&Liver==0),aes(x=Lung/10000 , y=0.0005, color=Mouse ),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse),method = "lm", se=F)+
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  #facet_grid(Cells~MI)+
  theme_classic()+
  #scale_color_grey()+
  #scale_color_grey()+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  scale_x_continuous(limits = c(0,100))+
  labs(title = "Scatter plot Lung vs Liver")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in Liver (%)")#5x6

ggarrange(scatter.IMFP, scatter.IV, nrow=1)

## --- Scatter PDX Tumour VS Lung, Liver, CTC Log scale -----

#Tum vs Lung
df_scatter %>% filter(Cells == "PDX-T412") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= ~subset(.,Tum==0&Lung>0),aes(x=0.0001, y=Lung/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= ~subset(.,Tum>0&Lung==0),aes(x=Tum/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Tum>0&Lung>0),aes(x=Tum/10000, y=Lung/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs Lung")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in Lung (%)")#4x14

# Tum vs CTC
df_scatter %>% filter(Cells == "PDX-T412") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI), alpha=0.9)+
  #geom_quasirandom(data= ~subset(.,Tum==0&CTC>0),aes(x=0.0001, y=CTC/10000, color=MI ),width=0.5,groupOnX=TRUE,alpha=0.7)
  geom_quasirandom(data= ~subset(.,Tum>0&CTC==0),aes(x=Tum/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Tum>0&CTC>0),aes(x=Tum/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Tumour vs CTC")+
  xlab("Barcode frequencies in Tumour (%)")+ylab("Barcode frequencies in CTC (%)")#4x10

#Lung vs CTC
df_scatter %>% filter(Cells == "PDX-T412") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI), alpha=0.9)+
  geom_quasirandom(data= ~subset(.,Lung==0&CTC>0),aes(x=0.001, y=CTC/10000),width=0.5,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= ~subset(.,Lung>0&CTC==0),aes(x=Lung/10000 , y=0.05),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow=1)+
  theme_classic()+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Lung vs CTC")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in CTC (%)")#4x10


## --- Scatter MDA Organs VS Lung, Liver, CTC Log scale -----


# Lung vs CTC
df_scatter %>% filter(Cells == "MDA-231") %>%
  ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=Mouse), alpha=0.9)+
  geom_quasirandom(data= ~subset(.,Lung==0&CTC>0),aes(x=0.0005, y=CTC/10000),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= ~subset(.,Lung>0&CTC==0),aes(x=Lung/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Lung vs CTC")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in CTC (%)")#4x16

# Lung vs CTC
df_scatter %>% filter(Cells == "PDX-T412") %>%
  ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  geom_point(data= ~subset(.,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=Mouse), alpha=0.9)+
  geom_quasirandom(data= ~subset(.,Lung==0&CTC>0),aes(x=0.001, y=CTC/10000),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= ~subset(.,Lung>0&CTC==0),aes(x=Lung/10000 , y=0.05),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Lung>0&CTC>0),aes(x=Lung/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Lung vs CTC")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in CTC (%)")#4x16


# Lung vs Liver

# R square Lung vs Liver:
ddply(filter(df_scatter, Cells == "MDA-231"),.(MI),function(x) round(summary(lm(x$Liver ~ x$Lung))$r.squared, digit = 3))

df_scatter %>% filter(Cells == "MDA-231") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  #stat_cor(data= ~subset(.,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=MI), alpha=0.9)+
  #geom_text(data=r2,aes(x=0.001,y=100, label = paste("R^2: ", round(r2,digit=4),sep="")),parse=T)+
  geom_point(data= ~subset(.,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=Mouse), alpha=0.9)+
  geom_quasirandom(data= ~subset(.,Lung==0&Liver>0),aes(x=0.0005, y=Liver/10000),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= ~subset(.,Lung>0&Liver==0),aes(x=Lung/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Lung>0&Liver>0),aes(x=Lung/10000, y=Liver/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Lung vs Liver")+
  xlab("Barcode frequencies in Lung (%)")+ylab("Barcode frequencies in Liver (%)")#4x16


# Liver vs CTC
df_scatter %>% filter(Cells == "MDA-231") %>%
ggplot()+
  geom_hline(yintercept=1, linetype = 'dotted')+ geom_vline(xintercept = 1, linetype = 'dotted')+
  #stat_cor(data= ~subset(.,Liver>0&CTC>0),aes(x=Liver/10000, y=CTC/10000, color=MI), alpha=0.9)+
  #geom_text(data=r2,aes(x=0.001,y=100, label = paste("R^2: ", round(r2,digit=4),sep="")),parse=T)+
  geom_point(data= ~subset(.,Liver>0&CTC>0),aes(x=Liver/10000, y=CTC/10000, color=Mouse), alpha=0.9)+
  geom_quasirandom(data= ~subset(.,Liver==0&CTC>0),aes(x=0.0005, y=CTC/10000),width=0.1,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= ~subset(.,Liver>0&CTC==0),aes(x=Liver/10000 , y=0.0005),width=0.3,groupOnX=F,alpha=0.7)+
  stat_smooth(data= ~subset(.,Liver>0&CTC>0),aes(x=Liver/10000, y=CTC/10000, color=MI),method = "lm", se=F)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(~facet, nrow = 1)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = sample(pals::tableau20(), 77, replace = TRUE))+
  #scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Liver vs CTC")+
  xlab("Barcode frequencies in Liver (%)")+ylab("Barcode frequencies in CTC (%)")#4x16



