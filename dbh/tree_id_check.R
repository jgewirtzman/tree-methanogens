master_tree_list<-read.csv("/Users/jongewirtzman/Downloads/Data Inventory & Library - Tree ID Dictionary.csv")

dataset<-read.csv("/Users/jongewirtzman/Downloads/Gas_Chromatograph_Data.csv")
dataset<-dataset[which(dataset$Tree.ID!="BLANK"),]

dataset_tree_list<-dataset$Tree.ID

colSums(sapply(dataset_tree_list, `==`, master_tree_list$Tree.ID..No.Separator.))==1


