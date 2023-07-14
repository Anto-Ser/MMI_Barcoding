# Section2: Pre - Figures analysis -----------------------------------------------------

raw_data = read.csv("./data/RAW_DATA.csv")

# Set to 0 all barcode with 10 reads or less:
raw_data[raw_data<=10]=0

##   Tresholding: Sample with more than 10k reads
raw_data= raw_data[,colSums(raw_data)>10000]

## Merging PCR replicates:
PCR_replicates = merge_PCR(raw_data)

## Merging Tumour samples: (combine all pieces of each tumour as one)
PCR_tum_merged = merge_tumor(PCR_replicates)

PCR_tum_merged$id <- as.integer(row.names(PCR_tum_merged))#Reorder by id name
PCR_tum_merged = PCR_tum_merged[order(PCR_tum_merged$id), ]#Reorder by id name
PCR_tum_merged = PCR_tum_merged[,-ncol(PCR_tum_merged)]#Reorder by id name


## Filter data for barcodes only found in tumours:
PCR_tum_filter = Filter_by_Tum(PCR_tum_merged)
names(PCR_tum_filter)[grep("data", names(PCR_tum_filter))] = setdiff(names(PCR_tum_merged), names(PCR_tum_filter))

## Separate vitro and vivo samples before next function:
vitro = PCR_tum_filter[,grep("BC09|Amp|vitro|ctrl|Ctrl|cltr", names(PCR_tum_filter)),drop=F]
vivo = select(PCR_tum_filter, -all_of(names(vitro)))

## Merge samples duplicate presence of "I" after names for in vivo samples:
PCR_tum_filter_final = merge_sampleI(vivo)

## Final data frame for analysis
df_MI = PCR_tum_filter_final[,grep(paste(c("Exp38_","Exp1009", "Exp1011", "Exp72", "Exp1012", "BC09"), collapse="|"), names(PCR_tum_filter_final)), drop=F]

df_MI = cbind(df_MI, vitro)

## -- Data cleaning in working data frame:

names(df_MI) = str_replace(names(df_MI), "_Lung1", "_Lung")
names(df_MI) = str_replace(names(df_MI), "_Liver1", "_Liver")
names(df_MI) = str_replace(names(df_MI), "_Blood", "_CTC")

## -- Info data frame for experiment and mice:

MI_info = read.csv("./data/Info_mice.csv")

MI_info$Mouse = str_remove(MI_info$Mouse, "'")
