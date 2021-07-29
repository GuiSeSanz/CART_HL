
# library(reticulate)
# reticulate::use_python("/usr/bin/python3")
# library(Seurat)
# library(cowplot)
# library(dplyr)
# library(future)
# library(viridis)
# library(ggplot2)
# library(plyr)
# library(reshape2)
# library(gridExtra)
# library(ggridges)

# weigths <- py_load_object('/home/sevastopol/data/gserranos/SimiC/Simic_Validation/Data/simic_val_Granja_filtered_L10.01_L20.001_Ws.pickle') #Weigths

# # files <- c('Bees_filtered_L10.01_L20.0001_Ws.pickle', 'ClonalKinetics_filtered_L10.01_L20.01_Ws.pickle', 'simic_val_Granja_filtered_L10.01_L20.01_Ws.pickle', 'simic_val_Kalsotra_QTH_L10.01_L20.01_AUCs.pickle', )

# for (name_w in names(weigths[['weight_dic']])){
#     print(name_w)
#     bottom <- weigths$weight_dic[[name_w]][101,]
#     tmp <- as.data.frame(weigths$weight_dic[[name_w]][-nrow(weigths$weight_dic[[name_w]]),])
#     colnames(tmp) <- weigths[['query_targets']]
#     rownames(tmp) <- weigths[['TF_ids']]
#     tmp <- scale(tmp, center=FALSE)
#     for( target in weigths[['query_targets']]){
#         tf_2_keep <- rownames(tmp)[order(abs(tmp[,target]), decreasing=T)[1:15]]
#         tmp[!rownames(tmp) %in% tf_2_keep, target] <- 0
#         print(sum(tmp[, target]!=0)) 
#     }
#     tmp <- rbind(tmp, bottom)
#     tmp <- as.matrix(tmp)
#     weigths$weight_dic[[name_w]] <- tmp
# }
# reticulate::py_save_object(weigths, filename = '/home/sevastopol/data/gserranos/SimiC/Simic_Validation/Data/simic_val_Granja_filtered_L10.01_L20.001_Ws_Filtered_15.pickle')





library(reticulate)
reticulate::use_python("/usr/bin/python3")
library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(ggridges)
library(stringr)





data_folder <- '/home/sevastopol/data/gserranos/CART_HL/SimiC/Data/'

# files <- c( 'Bees_filtered_L10.01_L20.0001_Ws.pickle', 
#             'simic_val_Granja_filtered_L10.01_L20.01_Ws.pickle', 
#             'ClonalKinetics_filtered_L10.01_L20.01_Ws.pickle', 
#             'simic_val_Kalsotra_QTH_L10.01_L20.01_Ws.pickle', 
#             'simic_val_Kalsotra_QTP_L10.001_L20.01_Ws.pickle' ,
#             'simic_val_Kalsotra_Time_L10.001_L20.001_Ws.pickle'
# )

files <- c('CART_HighLow_L10.01_L20.1_Ws.pickle')
for (fl in files){
    ifelse(file.exists(paste0(data_folder, fl)), print('ok'), print(fl))
}
 

for(fl in files){
    name <- str_extract( fl, '^[\\w]+(?=_L10)')
    pdf(paste0(data_folder, name, '_Weigths_BIC.pdf'))
    weigths <- py_load_object(paste0(data_folder, fl)) #Weigths
    print(fl)
    for (name_w in names(weigths[['weight_dic']])){
        print(name_w)
        max_l <- c()
        bottom <- weigths$weight_dic[[name_w]][101,]
        tmp <- as.data.frame(weigths$weight_dic[[name_w]][-nrow(weigths$weight_dic[[name_w]]),])
        colnames(tmp) <- weigths[['query_targets']]
        rownames(tmp) <- weigths[['TF_ids']]
        tmp <- scale(tmp, center=FALSE)
        for( target in weigths[['query_targets']]){
            all_tf <- rownames(tmp)[order(abs(tmp[,target]), decreasing=T)]
            l <- 1
            tf_2_keep <- all_tf[1:l]
            while ( sum(tmp[tf_2_keep, target]^2) /sum(tmp[,target]^2)  < 0.9 ){
                l <- l +1
                tf_2_keep <- all_tf[1:l]
            }
            max_l <- c(max_l, l)
            tmp[!rownames(tmp) %in% tf_2_keep, target] <- 0
        }
        tmp <- rbind(tmp, bottom)
        tmp <- as.matrix(tmp)
        weigths$weight_dic[[name_w]] <- tmp
        hist(max_l, breaks = 20, main = paste0(name, '_Weights_',name_w ), col = 'darkgrey')
    }
    new_fl_name = paste0(data_folder,sub('.pickle', '_filtered_BIC.pickle', fl))
    reticulate::py_save_object(weigths, filename = new_fl_name)
    dev.off()
    print(new_fl_name)
}
