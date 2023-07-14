# Section0: Libraries ----------------------------------------

Packages = c("dplyr", "stringr", "readxl", "tidyr", "ggplot2",
             "reshape2", "corrplot", "viridis", "grid", "gridExtra", "gridGraphics",
             "janitor", "RColorBrewer","ggbeeswarm", "ggpubr", "vegan", "ggcorrplot", "randomcoloR",
             "scales", "networkD3", "gridGraphics", "grid", "plyr", "tibble","ggthemes", "circlize",
             "writexl", "mdthemes", "ggpmisc", "pals")

new_packages <- Packages[!(Packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(Packages, library, character.only = TRUE)

# Section1: Functions -----------------------------------------------------

# --- Merge PCR to merge replicate together ---

merge_PCR = function(data,min.rep=2,CPM=TRUE,show.corr=FALSE, cor.filter = 0.6){
  #Code to generate filtered data by replicate
  #min.rep: is the minimal number of replicates
  #CPM: normalize to counts per million (CPM) prior and after merging
  #note:remove samples with zero counts before running this
  samples <- substr(names(data),1,nchar(names(data))-1)
  merged = data.frame(row.names = row.names(data))

  for(s in unique(samples)){

    # Loop status:
    if(which(s==unique(samples)) %in% round(quantile(1:length(unique(samples)), probs = seq(0.1,1,0.1)))){
      print(paste0("Merging replicates: ",
                   round(which(s==unique(samples))/length(unique(samples))*100),
                   "% samples completed"))}

    # Merging replicates:
    d = data[,grep(s,names(data))]
    d = as.data.frame(d)
    if(CPM==TRUE & ncol(d)>=2 & nrow(d)>0) d=sweep(d,2,colSums(d)/1e6,`/`)
    if(ncol(d)>=2 & show.corr==TRUE & nrow(d)>1) {
      corr=cor(d,use="pair",method="p")
      print(corr[1,2])
      diag(corr)=NA #replace diago
    }

    if(ncol(d)<=1){merged=cbind(merged,0)}
    else {merged=cbind(merged,(rowSums(d>0,na.rm = TRUE)>=min.rep)*rowSums(d,na.rm = TRUE))}
  }
  if(CPM==TRUE) merged=sweep(merged,2,colSums(merged)/1e6,`/`)
  names(merged)=unique(samples)
  merged = merged[,!is.na(colSums(merged))]
  return(merged)
}

# --- Merge Tumour samples ---

merge_tumor = function(data,min.rep=1,CPM=TRUE, Remove_NA_col=FALSE){
  tum.pieces = names(select(data,contains("Tum")))
  tum.samples <- paste0(apply(str_split_fixed(tum.pieces,"_",3)[,1:2],1, paste,collapse="_"),"_Tum")
  tum.merged = data.frame(row.names = row.names(data))

  for(s in unique(tum.samples)){

    # Loop status:
    Progress_loop_quantiles = round(quantile(1:length(unique(tum.samples)), probs = seq(0.25,1,0.25)))
    if(which(s==unique(tum.samples)) %in% Progress_loop_quantiles){
      print(paste0("Merging tumours: ",
                   names(Progress_loop_quantiles[Progress_loop_quantiles == which(s==unique(tum.samples))]),
                   " samples completed"))}

    d = data[,grep(s,names(data))]
    if(Remove_NA_col==TRUE){d=d[!is.na(colSums(d))]}
    if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`)
    if(ncol(d)<=1){tum.merged=cbind(tum.merged,0)}
    else {tum.merged=cbind(tum.merged,(rowSums(d>0)>=min.rep)*rowSums(d))}
  }
  names(tum.merged)=unique(tum.samples)
  #tum.merged=tum.merged[,colSums(tum.merged)>0]
  if(CPM==TRUE) tum.merged=sweep(tum.merged,2,colSums(tum.merged)/1e6,`/`)

  ##now merge this with the non.tumor samples
  return(merge(tum.merged,select(data,!contains("Tum")),by="row.names",all=TRUE)%>%tibble::column_to_rownames(var="Row.names") )
}


#--- Function to filter data per mouse by presence of barcode in primary tumour:

Filter_by_Tum = function(data, CPM=TRUE, filter=TRUE){
  exp = str_split_fixed(names(data),"_",3)[,1]
  mice = str_split_fixed(names(data),"_",3)[,2]
  mice = mice[nchar(mice)>0]
  mice = paste(exp, mice, sep = "_")
  df.merged = data.frame(row.names = row.names(data))
  for(s in unique(mice)){
    #print(s) #s = "339"
    d = as.data.frame(data[,grepl(paste0(s,"_"),names(data))])
    DA = which(d[,grep("Tum", names(d)), drop=F]>0)
    if(filter==TRUE) d[-DA,]=0
    if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`)
    df.merged = cbind(df.merged,d)
  }
  return(df.merged)
}

#--- Merge samples labeled nameI (ex:Lung LungI) as a single sample

merge_sampleI = function(data, CPM=TRUE){
  #data=df.PCR.Tum.Filt
  #CPM=TRUE
  Mouse.nb = str_split(names(data), "_", simplify = T)[,2]
  Mega.merged = data.frame(row.names = row.names(data))
  for (i in unique(Mouse.nb)) {
    #i="333"
    #print(i)
    TEST.name = paste("_", i, "_", sep="")
    Ms.x = data[,grep(TEST.name, names(data)), drop=F]
    id.name <- paste(str_split(names(Ms.x[1]), "_", simplify = T)[,c(1,2)], collapse = "_")
    colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
    names(Ms.x)[which(endsWith(names(Ms.x), 'I'))] = str_sub(names(Ms.x)[which(endsWith(names(Ms.x), 'I'))],1,nchar(names(Ms.x)[which(endsWith(names(Ms.x), 'I'))])-1)

    samples <- names(Ms.x)
    merged = data.frame(row.names = row.names(Ms.x))
    for(s in unique(samples)){
      #print(s)
      #s="Lung"
      d = Ms.x[,grep(s,names(Ms.x)), drop=F]
      if(ncol(d)==1){merged=cbind(merged,d);next}
      d = d[,!is.na(colSums(d)), drop=F]
      if(ncol(d)==1){merged=cbind(merged,d);next}
      if(ncol(d)==2){d=sweep(d,2,colSums(d)/1e6,`/`)}
      merged=cbind(merged,(rowSums(d>0)*rowSums(d)))
    }
    if(CPM==TRUE) merged=sweep(merged,2,colSums(merged)/1e6,`/`)
    names(merged) = paste(id.name,unique(samples), sep = "_")
    Mega.merged = cbind(Mega.merged,merged)
  }
  return(Mega.merged)
}

#--- Stacked histogram function: needs a dataframe with one column ID

Stacked_histo = function(data, PCT=TRUE, angle.x=0, my_col = "Default"){
  if(PCT==TRUE) {data= data/10000}
  if(sum(names(data)=="ID")==0){data$ID = 1:nrow(data)}

  combo_color = c(brewer.pal(11, "Spectral"), brewer.pal(11, "RdYlGn"), brewer.pal(11, "RdBu"), brewer.pal(11, "PuOr"), brewer.pal(11, "PiYG"))
  getPalette = colorRampPalette(combo_color)
  if(any(my_col== "Default")){COLOR = sample(getPalette(2608), 2608)} else {
    COLOR = my_col
  }

  data$COLOR = COLOR
  COLOR.to.save = COLOR
  data = data[rowSums(data[1:(ncol(data)-2)])>0,] #remove row w/o bc (ignoring last 2 col (ID, colours))
  colnames(data) =gsub("^[^_]*_","",names(data))

  COLOR= data$COLOR
  data = data[,-ncol(data)]

  mA = melt(data, id=c("ID"))
  mA$ID = as.factor(mA$ID) #transform barcode number in factor

  mA = mA %>% group_by(variable) %>% mutate(pos = sum(value)-(cumsum(value)-0.5*value))
  #Add new colunm in mA table for position of each barcode (to later on be labeled on graph)

  insert_minor <- function(major_labs, n_minor) {labs <-
    c( sapply( major_labs, function(x) c(x, rep("", n_minor) ) ) )
  labs[1:(length(labs)-n_minor)]}

  #PLot:
  p123 = ggplot(data =mA, aes(y = value, x = variable, fill = ID, label= ID)) +
    geom_bar(stat="identity",alpha=0.9) + #pos="fill"
    #geom_text(data = subset(mA, value>5), aes(y = pos, label = ID), size = 3) +
    theme_bw()+
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "grey10", face="bold",size=10,angle = 0),
          axis.text.x = element_text(color = "black", face="bold",size=10,angle = angle.x,
                                     vjust = 0.3),
          #hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "grey10", size = 15))+
    xlab("")+ylab("Number of reads(%)")+ ggtitle("")+
    scale_fill_manual(values = COLOR )+
    scale_y_continuous(breaks = seq(0,100, 10),
                       labels = insert_minor (seq(0, 100, by=20),1),
                       limits = c(0,100.001),
                       expand = c(0,0))
  my_list = list(Stack.plot = p123, Colour= COLOR.to.save)
  return(my_list)
}

#--- Bubble function:
Bubble_plot = function(data, PCT=TRUE, Y="avg", data2, angle.x=0, my_col="Default"){

  if(PCT==TRUE){data = data/10000} #CPM to percent
  if(Y=="avg"){data$Y = rowSums(data)/ncol(data)} else if (Y=="Random"){#Average all sample for Y order or Random order
    data$Y= sample(length(data[,1]), length(data[,1]))
  } else if(Y=="input") {data$Y = data2 #manual input of previously saved Y
  } else{data$Y = data[,1]}#Take first column as Y ref = Tum collapsed
  colnames(data) =gsub("^[^_]*_","",names(data))
  data$ID = as.numeric(rownames(data))

  data1 = data[rowSums(data[1:(ncol(data)-2)])>0,1:(ncol(data)-2)] #remove empty rows
  data1$Y = data[rownames(data1),"Y"]
  data1$ID = data[rownames(data1),"ID"]

  mG =melt(data1, id=c("ID","Y"))
  mG$ID = as.factor(mG$ID)
  mG[mG==0] = NA

  Pallett = unname(distinctColorPalette(500))
  Max.color = sample(Pallett, length(unique(data$ID)),replace = TRUE)


  if(any(my_col== "Default")){col.bbp = Max.color[data1$ID]} else {
    col.bbp = my_col[data1$ID]
  }
  Plot.bbp = ggplot(data=mG)+
    #geom_point(aes(x=variable,y=Y,color=ID,size=value), data=subset(mG,variable=!"X25k_P0_1M", stat="identity"),alpha=0.5)+
    geom_quasirandom(aes(x=variable, y=Y, color=ID, size=value, group=1),alpha=0.8, width = 0.2)+
    scale_size_continuous(range = c(1.5,15))+
    #geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)+
    scale_color_manual(values = col.bbp,guide="none")+ #breaks=length(unique(mG$ID)) removed
    xlab("")+ylab("Frequency in Tumours")+ ggtitle("")+labs(size="Frq %")+
    #scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
    #                  breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
    #                 labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    theme_minimal()+
    theme( panel.grid.major.x = element_blank() ,
           panel.grid.minor.y = element_blank(),
           axis.text.x = element_text(angle = angle.x))

  if(Y %in% c("avg")){Plot.bbp = Plot.bbp +
    scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
                       breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
                       labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)} else if(Y %in% c("input", "Random")){

      Plot.bbp = Plot.bbp +
        scale_y_continuous(limits = c(0,length(data[,1])))}
  Plot.bbp


  my_list = list(BBP = Plot.bbp, Colour = Max.color, Y.order = data$Y)
  return(my_list)
}

#Select the specific organ and order the sample by mode of injection
Select_and_order_Organs = function(data, Model, Organ, info_dataframe){
  MI.tum = data[,grep(paste("_", Organ, sep = ""), names(data)), drop=F]

  order.tum.MDA= vector()
  for (i in unique(info_dataframe$MI)) {
    print(i)
    #i="IMFPc"
    order.temp = info_dataframe %>%
      filter(Org == Organ) %>%
      filter(MI == i & Cells == Model) %>% arrange(Mouse) %>% .$Full_name
    order.tum.MDA = c(order.tum.MDA, order.temp)
  }
  Tum.order = MI.tum[,order.tum.MDA]
  newdf <- t(na.omit(t(Tum.order)))
  newdf = as.data.frame(newdf)
  Tum.order = newdf
}


