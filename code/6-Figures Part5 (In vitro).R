## -- In Vitro analysis --

rm(list=setdiff(ls(), c("df_MI", "MI_info", "Number_barcode", "df_scatter", "Bubble_plot", "Stacked_histo", "Select_and_order_Organs")))

## ---  Passage 1 controls ----
P1 =  df_MI %>% select(c("BC09_25k_P0.1M_",
                          "Exp1009_cltrvitro_500k1_", "Exp1009_cltrvitro_500k2_",
                          "Exp1011_Ctrl_Ctrl3_", "Exp1011_Ctrl_Ctrl4_",
                          "Exp1009_Amp1_P1_", "Exp1009_Amp2_P1_", "Exp1009_Amp3_P1_", "Exp1009_Amp4_P1_",
                          "Exp1011_vitroA_P1_","Exp1011_vitroB_P1_","Exp1011_vitroC_P1_","Exp1011_vitroD_P1_"))
names(P1) = c("P0", "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4",
              "Flask1 P1", "Flask2 P1", "Flask3 P1", "Flask4 P1", "Flask5 P1", "Flask6 P1", "Flask7 P1", "Flask8 P1")

In.Vitro.Colours = read.csv("./data/Colour_InVitro_Histo.csv")
vitro.col = In.Vitro.Colours$x

Vitro.p1 = Stacked_histo(P1, angle.x = 90, my_col = vitro.col)
Vitro.p1$Stack.plot #6x10
#vitro.col = Vitro.p1$Colour
#write.csv(vitro.col, file="Colour_InVitro_Histo.csv")



corrplot(as.matrix(cor(log(P1+1))),is.corr = F,
         tl.col = 'black',col = viridis(50),
         #type="upper",
         method = "color") #7x7

## ---  Passages overtime ----

vitro_raw = df_MI[,grep("BC09|Amp|vitro|ctrl|Ctrl|cltr", names(df_MI)),drop=F]

df_vitro = data.frame(nam=names(vitro_raw), bc =colSums(vitro_raw>0), shannon= diversity(t(vitro_raw)))
df_vitro = cbind(df_vitro, str_split(df_vitro$nam, "_", simplify = T))
df_vitro = df_vitro[,-7]
names(df_vitro) = c("Full_name","bc","Shannon", "Exp","Flask","Passage")

#write.csv(df.sum.vitro, file="Vitro_barcode_Shannon.csv")

old_names = c("Amp1", "Amp2", "Amp3", "Amp4", "vitroA", "vitroB", "vitroC", "vitroD", "25k", "cltrvitro")
new_names = c(paste0("Flask", 1:8), "Population", "Control")

df_vitro = df_vitro %>%
  #select(-last(names(.))) %>%
  mutate(New = stringi::stri_replace_all_regex(Flask, pattern = old_names, replacement = new_names, vectorize = FALSE)) %>%
  filter(New %in% new_names & New !="Control") %>%
  mutate(Passage_n = as.numeric(str_replace(str_remove(str_remove(Passage, "P"), "bis"), "0.1M", "0"))) %>%
  arrange(New, Passage_n) %>%
  mutate(new_name = paste0(New, "-P", Passage_n))


## --- Comparison Tumours vs Passage 10 in vitro -----

vitro_clean = vitro_raw[,df_vitro$Full_name]
names(vitro_clean) = df_vitro$new_name
vitro_clean = vitro_clean %>% select("Population-P0", everything()) # Move last column to first position

P10 = vitro_clean[,grep("P10", names(vitro_clean))]

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Tum", info_dataframe = Number_barcode)

MDA.tum.vitro = MDA.tum[,grep(paste("Exp1009", "Exp1011", sep="|"), names(MDA.tum))]

MDA.tum.vitro.P10 = cbind(MDA.tum.vitro, P10)

corrplot(as.matrix(cor(log(MDA.tum.vitro.P10+1))),is.corr = F,cor.lim =c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') #7x7


## --- Correlation by experiment and by MI -----

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Tum", info_dataframe = Number_barcode)

corrplot(as.matrix(cor(log(MDA.tum+1))),is.corr = F,cor.lim =c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') #10x10

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Tum", info_dataframe = Number_barcode)
MDA.tum = MDA.tum[,grep("Exp38", names(MDA.tum))]

corrplot(as.matrix(cor(log(MDA.tum+1))),is.corr = F, col.lim = c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') # 6x6

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Tum", info_dataframe = Number_barcode)
MDA.tum = MDA.tum[,grep("Exp1009", names(MDA.tum))]

corrplot(as.matrix(cor(log(MDA.tum+1))),is.corr = F, col.lim = c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') # 6x6

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Tum", info_dataframe = Number_barcode)
MDA.tum = MDA.tum[,grep("Exp1011", names(MDA.tum))]

corrplot(as.matrix(cor(log(MDA.tum+1))),is.corr = F, col.lim = c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') # 6x6


MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Lung", info_dataframe = Number_barcode)
MDA.tum = MDA.tum[,grep("Exp38", names(MDA.tum))]



corrplot(as.matrix(cor(log(MDA.tum+1))),is.corr = F, col.lim = c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') # 6x6

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Lung", info_dataframe = Number_barcode)
MDA.tum = MDA.tum[,grep("Exp1009", names(MDA.tum))]

corrplot(as.matrix(cor(log(MDA.tum+1))),is.corr = F, col.lim = c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') # 6x6

MDA.tum = Select_and_order_Organs(df_MI, Model = "MDA-231", Organ = "Lung", info_dataframe = Number_barcode)
MDA.tum = MDA.tum[,grep("Exp1011", names(MDA.tum))]

corrplot(as.matrix(cor(log(MDA.tum+1))),is.corr = F, col.lim = c(-1,1),
         tl.col = 'black', order = "hclust",
         #type="upper",col = viridis(100),
         method = "color",
         hclust.method = 'ward.D2') # 6x6


## --- Bubble plots tumours and passage 10 -----

MDA.tum.vitro.P10$Pop_P0 = vitro_clean$`Population-P0`
MDA.tum.vitro.P10 = MDA.tum.vitro.P10 %>% select("Pop_P0", everything())

bubble.plot.MDA.Tumour.and.P10 = Bubble_plot(MDA.tum.vitro.P10, angle.x = 90, my_col = vitro.col, Y = "input", data2 = (MDA.tum.vitro.P10$Pop_P0/10000)) #6x15
bubble.plot.MDA.Tumour.and.P10$BBP +
  scale_y_continuous(trans = 'log10',limits=c(3e-3,2),
                     breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
                     labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))

## --- Simulation for frequencies in tumours based on P0 frq -----

df_for_sampling = data.frame(Var1=as.character(1:2608))

n = colSums(MDA.tum.vitro>0) # number of barcode per tumour for sampling same number of barcodes

set.seed(101)

for(i in 1:length(n)){
  t = table(sample(1:2608, n[i], prob = MDA.tum.vitro.P10[,"Pop_P0"], replace = TRUE)) %>% as.data.frame()
  names(t) = c("Var1", paste0("Sim_Tum_",i))
  df_for_sampling = df_for_sampling %>% full_join(t, by="Var1")
}

df_for_sampling[is.na(df_for_sampling)]=0

df_Tum_Sim = data.frame(Tum_n_bc = rowSums(MDA.tum.vitro>0),
                         Sim_n_bc = rowSums(df_for_sampling[,-1]>0),
                         P0 = MDA.tum.vitro.P10[,"Pop_P0"])

ggplot()+
  geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)+
  geom_point(data = df_Tum_Sim,aes(Tum_n_bc, P0/10000, colour="Tumours", size = Tum_n_bc/Sim_n_bc), position = "jitter", alpha=0.4)+
  geom_point(data= df_Tum_Sim,aes(Sim_n_bc, P0/10000, colour="Simulation"), position = "jitter", alpha =0.4)+
  labs(subtitle = "Barcode appearance tumour vs simulation based on frequency at P0",
       x= "Number of appearence in tumour", y="Barcode frequency in initial population")+
  theme_bw()+
  theme(panel.grid.minor.y = element_blank())+
  scale_y_continuous(trans = 'log10',limits=c(3e-3,2),
                     breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
                     labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))

## --- Extract top 10 barcode per flask at passage 10 -----


P10 = vitro_clean[,grep("P10", names(vitro_clean))]
head(df_vitro)
unique(df_vitro$Passage)[unique(df_vitro$Passage) != "P0.1M"]

top10_all_flask_P10 = P10 %>%
  dplyr::mutate(ID=1:2608) %>%
  pivot_longer(!ID,names_to = "name", values_to = "freq_bc") %>%
  group_by(name)%>%
  arrange(desc(freq_bc)) %>%
  slice(1:10) %>%
  pull(ID)

# Biomass of total barcode top10 barcode per flask at P10 in MDA tumours:

df_Biomass_vitro_Tum = data.frame(matrix(nrow=0, ncol=3))
colnames(df_Biomass_vitro_Tum) = c("names", "biomass", "passage")
length_bc_in_top = vector()
vitro_clean2 = vitro_clean
names(vitro_clean2) = paste0(names(vitro_clean), "_")

for (i in 1:31) {
  #print(i)
  #i=1
  c = paste0("P",i, "_")
  top10_all_flask = vitro_clean2[,grep(c, names(vitro_clean2))] %>%
    dplyr::mutate(ID=1:2608) %>%
    pivot_longer(!ID,names_to = "name", values_to = "freq_bc") %>%
    group_by(name) %>%
    arrange(desc(freq_bc)) %>%
    #slice(1:10) %>% pull(ID) # For top 10 barcode per flask remove code below
    dplyr::mutate(cumul_frq = cumsum(freq_bc)) %>% # For top 95% barcode
    filter(cumul_frq < 950000) %>% # For top 95% barcode
    pull(ID) # For top 95% barcode
  length_bc_in_top = c(length_bc_in_top, length(unique(top10_all_flask)))
  vitro_temp = vitro_clean2[,grep(c, names(vitro_clean2))]
  df_Biomass_vitro_Tum = rbind(df_Biomass_vitro_Tum, data.frame(names = names(MDA.tum.vitro),
                                                                biomass = colSums(MDA.tum.vitro[unique(top10_all_flask),])/10000,
                                                                passage_dis = c,
                                                                passage_con = i))
}


df_Biomass_vitro_Tum2 = df_Biomass_vitro_Tum %>% filter(names %in%names(MDA.tum.vitro)[1:16])

df_Biomass_vitro_Tum2 %>% filter(passage_con %in%c(1,seq(5,15,5)))%>%
  group_by(passage_con) %>%
  dplyr::mutate(median_bio= median(biomass), sd_bio = sd(biomass))%>%
  ungroup()%>%
  ggplot(aes(x=factor(passage_con), y=biomass))+
  geom_bar(stat = "summary", fun = "median")+
  geom_errorbar(aes(ymin=median_bio-sd_bio, ymax=median_bio+sd_bio), width=0.1)+
  geom_quasirandom(size=2, width = 0.03)+
  theme_classic()

df_Biomass_vitro_Tum2 %>% filter(passage_con %in%c(1,seq(5,15,5)))%>%
  group_by(passage_con) %>% dplyr::summarise(Median_bio = median(biomass))

stat_biom = df_Biomass_vitro_Tum2 %>% filter(passage_con %in%c(1,seq(5,15,5)))
model <- aov(biomass~factor(passage_con), data=stat_biom)
summary(model)
TukeyHSD(model, conf.level=.95)

