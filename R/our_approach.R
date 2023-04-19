rm(list=ls())
#setwd("/home/behnaz/Desktop/Our_approach/")
library(rlist)
source("./R/Functions.R")
#Data has 3 columns: Gene	Position	Quantification

Our_approach=function(Data,length_of_path=1){
  prepared_data=data.frame(Phosphosite=paste(Data$Gene,Data$Position,sep="_"),
                           Quantification=Data$Quantification)
  if(sum(duplicated(prepared_data$Phosphosite))>0){
    prepared_data=aggregate(.~Phosphosite,data = prepared_data,FUN = mean)
    print("each phosphosite should give once, but you have duplicated row, we consider mean of Quantification for these rows")
  }
  list_of_target=prepared_data$Quantification
  names(list_of_target)=prepared_data$Phosphosite
  all_path_nodes=find_all_path_nodes(list_of_target=names(list_of_target),length_of_path=length_of_path,Edges=Edges,end_point_gene=NaN,only_end_point=F)
  all_adjacency_matrices=create_all_adjacency_matrices(all_path_nodes)
  scores=compute_scores(all_adjacency_matrices=all_adjacency_matrices,
                        all_path_nodes=all_path_nodes,IF_of_target=list_of_target)
  data_scores=data.frame(Score=unlist(scores),Kinase_Gene = names(scores))
  scores_old_forward=compute_scores_with_old_forward(all_adjacency_matrices=all_adjacency_matrices,
                        all_path_nodes=all_path_nodes,IF_of_target=list_of_target)
  data_scores_old_forward=data.frame(Score=unlist(scores_old_forward),Kinase_Gene = names(scores_old_forward))
  scores_without_Forward=compute_scores_without_Forward(all_adjacency_matrices=all_adjacency_matrices,
                        all_path_nodes=all_path_nodes,IF_of_target=list_of_target)
  data_scores_without_Forward=data.frame(Score=unlist(scores_without_Forward),Kinase_Gene = names(scores_without_Forward))
  results=list()
  results[["scores"]]=data_scores
  results[["scores_without_forward"]]=data_scores_without_Forward
  results[["scores_with_old_forward"]]=data_scores_old_forward
  return(results)
}


####sample

#Data=read.csv("Sample_data.csv")
#result=Our_approach(Data)
