# Remove bioactivity repetitions: some of the compounds are annotated on the same target more than one time
# We take the mean of the observed values as the bioactivity value for a given compound-target pair

bio_repetitions <- readRDS("bioactivity_repetitions.rds")

# Read the information of the dataset
df1 <- readRDS("COX_dataset_info.rds") #read.table(file="COX_dataset_info.csv",sep=',',header=T)

##################################################################
#### Preprocessing of the fingerprints 
##################################################################
df2 <- readRDS("fps_COX_512.rds") #read.table(file="fps_COX_counts_512.csv",sep=",",header=FALSE)
nzv.columns <- nearZeroVar(df2, freqCut = 30/1)
nb_descs = ncol(df2) - length(nzv.columns)
nzv.names <- names(df2)[nzv.columns]
df2 <- df2[, -nzv.columns]
is.highly.correlated <- findCorrelation(cor(df2), 0.95)
df2 <- df2[,-is.highly.correlated]


##################################################################
#### Amino Acid Descriptors
#################################################################
df3 <- readRDS("Z3_COX.rds")
nzv.columns <- nearZeroVar(df3, freqCut = 30/1)
nb_descs = ncol(df3) - length(nzv.columns)
nzv.names <- names(df3)[nzv.columns]
df3 <- df3[, -nzv.columns]
#is.highly.correlated <- findCorrelation(cor(df3), 0.95)
#df3 <- df3[,-is.highly.correlated]

#################################################################
### Physicochemical (PaDEL) descriptors
################################################################
df4 <- readRDS("Padel_COX.rds")

rows = seq(1,nrow(df1))
df1=cbind(rows,df1)
df3_not_repeated = c()
df2_not_repeated = c()
df4_not_repeated = c()
new_data = c();new_data2 = c();new_data3 = c();new_data4 = c()

stdDevs = c()
Bioactivities=c()
ranges=c()
rows_repeated= c()

for (i in 1:nrow(bio_repetitions)){
bioactivity=c()
stdDev=c()
datos=c();datos2=c();datos3=c();datos4=c();datos5=c()
subset1=c();subset1_df2=c();subset1_df3=c();subset1_df4=c()
subset2=c();subset2_df2=c();subset2_df3=c();subset2_df4=c()

#first filter
subset1 = df1[which(df1$accession == bio_repetitions$V2[i]),]
subset1_df2 = df2[which(df1$accession == bio_repetitions$V2[i]),]
subset1_df3 = df3[which(df1$accession == bio_repetitions$V2[i]),]
subset1_df4 = df4[which(df1$accession == bio_repetitions$V2[i]),]

#second filter
subset2 = subset1[which(subset1$chembl_id.1 == bio_repetitions$V3[i]),]
subset2_df2 = subset1_df2[which(subset1$chembl_id.1 == bio_repetitions$V3[i]),]
subset2_df3 = subset1_df3[which(subset1$chembl_id.1 == bio_repetitions$V3[i]),]
subset2_df4 = subset1_df4[which(subset1$chembl_id.1 == bio_repetitions$V3[i]),]

values = (subset2$standard_value)*10^-9
values = -log(values,base=10)
bioactivity = mean(values)
rangenow=range(values)[2] - range(values)[1]
ranges=c(ranges,rangenow)
stdDev = sd(values)
datos = subset2[1,]
datos2 = subset2_df2[1,]
datos3 = subset2_df3[1,]
datos4 = subset2_df4[1,]

datos$standard_value = bioactivity
stdDevs = c(stdDevs,stdDev)
Bioactivities=c(Bioactivities,bioactivity)
# We add to the data.frame contatining all the mean values
new_data = rbind(new_data, datos)
new_data2 = rbind(new_data2, datos2)
new_data3 = rbind(new_data3, datos3)
new_data4 = rbind(new_data4, datos4)
# We get the rows that which compounds repeated
rows_repeated = c(rows_repeated, subset2$rows)
}

# Data not repeated
df_not_repeated= df1[-rows_repeated,]
df3_not_repeated = df3[-rows_repeated,]
df2_not_repeated = df2[-rows_repeated,]
df4_not_repeated = df4[-rows_repeated,]

hh= as.vector(df_not_repeated$standard_value)*10^-9
hh= -log(hh,base=10)
df_not_repeated$standard_value = hh

data_all1 = rbind(df_not_repeated,new_data)
data_all2 = rbind(df2_not_repeated,new_data2)
data_all3 = rbind(df3_not_repeated,new_data3)
data_all4 = rbind(df4_not_repeated,new_data4)

# We save the dataset
final_data = data.frame(data_all1,data_all2,data_all3,data_all4)
saveRDS(final_data, file="Whole_dataset.rds")

# We save the datasets for QSAR and QSAM comparisons
final_AA_descs = data.frame(data_all1,data_all3)
final_Compound_descs = data.frame(data_all1,data_all2,data_all4)
saveRDS(final_AA_descs, file="Whole_dataset_aa_descriptors.rds")
saveRDS(final_Compound_descs, file="Whole_dataset_compound_descriptors.rds")

save(list=ls(), file="data_processing_repetitions_COX.RData")

