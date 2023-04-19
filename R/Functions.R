#' @importFrom rlist list.sort
#' @import rlist
#' @export
Edges=readRDS(file = "./data/Edges.rds")

find_all_path_nodes=function(list_of_target,length_of_path,Edges,end_point_gene=NaN,only_end_point=F){
  Edges=unique(Edges)
  #forbidden_nodes=gsub("_.*","",(list_of_target[list_of_target %in% Edges$GeneSymbol_b]))
  all_path=list()
  co_path=list()
  all_nodes=list()
  all_par_nodes=list()
  all_co_nodes=list()
  all_nodes[[1]]=list_of_target[list_of_target %in% Edges$GeneSymbol_b]
  #all_previous_nodes=c(all_nodes[[1]],forbidden_nodes)
  all_previous_nodes=c(all_nodes[[1]])
  i=1
  while(i<(length_of_path+1)){
    print(i)
    new_edge=unique(Edges[,c("GeneSymbol_a","GeneSymbol_b","Interaction")])
    colnames(new_edge)=c(paste('node_',as.character(i+1),sep=''),paste('node_',as.character(i),sep=''),paste('w_e',as.character(i),sep=''))
    new_weighted_edge=new_edge[(new_edge[,paste('node_',as.character(i),sep='')] %in% all_nodes[[i]]),]
    all_new_nodes=new_weighted_edge[,paste('node_',as.character(i+1),sep='')]
    new_nodes=unique(all_new_nodes[!all_new_nodes %in% all_previous_nodes])
    #new_nodes=unique(all_new_nodes)
    if(length(new_nodes)==0){
      break
    }else{
      all_previous_nodes=c(all_previous_nodes,new_nodes)
      all_path[[i]]=new_weighted_edge[new_weighted_edge[,paste('node_',as.character(i+1),sep='')] %in% new_nodes,]
      all_nodes[[i+1]]=new_nodes
      co_path[[i+1]]=Edges[Edges$GeneSymbol_a %in% new_nodes & Edges$GeneSymbol_b %in% new_nodes,]
      i=i+1
    }
  }
  return(list(all_path=all_path,all_nodes=all_nodes,co_path=co_path))
}

##this function create adjacency matrices for all layers
create_all_adjacency_matrices=function(all_path_nodes){
  child_adjacency_matrices=c()
  co_adjacency_matrices=c()
  length_of_path=length(all_path_nodes$all_path)
  child_adjacency_matrices[[1]]=create_child_adjacency_matrix(all_path_nodes$all_path[[1]],all_path_nodes$all_nodes[[2]],all_path_nodes$all_nodes[[1]])
  if(length_of_path>1){
    for(i in (2:length_of_path)){
      child_adjacency_matrices[[i]]=create_child_adjacency_matrix(all_path_nodes$all_path[[i]],all_path_nodes$all_nodes[[i+1]],all_path_nodes$all_nodes[[i]])
      co_adjacency_matrices[[i]]=create_co_adjacency_matrix(all_path_nodes$co_path[[i]],all_path_nodes$all_nodes[[i]])
    }
  }
  if(dim(all_path_nodes$co_path[[length_of_path+1]])[1]>0){
    co_adjacency_matrices[[length_of_path+1]]=create_co_adjacency_matrix(all_path_nodes$co_path[[length_of_path+1]],all_path_nodes$all_nodes[[length_of_path+1]])
  }
  return(list(child_adjacency_matrices=child_adjacency_matrices,co_adjacency_matrices=co_adjacency_matrices))
}

##this function create adjacency matrix for one layer
create_child_adjacency_matrix=function(path,par_nodes,child_nodes){
  path_matrix=as.matrix(unique(path))
  unique_row=unique(par_nodes)
  unique_col=unique(child_nodes)
  Adjacency_matrix=matrix(0,nrow = length(unique_row),ncol = length(unique_col))
  colnames(Adjacency_matrix)=unique_col
  row.names(Adjacency_matrix)=unique_row
  for(row in (1:dim(path_matrix)[1])){
    Adjacency_matrix[path_matrix[row,1],path_matrix[row,2]]=as.numeric(path_matrix[row,3])
  }
  return(Adjacency_matrix)
}

###if nodes were in one layer we create adjacency matrix with this
create_co_adjacency_matrix=function(path,nodes){
  path_matrix=as.matrix(unique(path))
  unique_row=unique(nodes)
  unique_col=unique(nodes)
  co_Adjacency_matrix=matrix(0,nrow = length(unique_row),ncol = length(unique_col))
  colnames(co_Adjacency_matrix)=unique_col
  row.names(co_Adjacency_matrix)=unique_row
  for(row in (1:dim(path_matrix)[1])){
    co_Adjacency_matrix[path_matrix[row,1],path_matrix[row,2]]=as.numeric(path_matrix[row,3])
  }
  return(co_Adjacency_matrix)
}
NewForwardAlg=function(all_adjacency_matrices,all_path_nodes){
  P=list()
  final_P=list()
  length_of_path=length(all_adjacency_matrices$child_adjacency_matrices)
  #P[[length_of_path+1]]=rep(1/length(all_path_nodes$all_nodes[[length_of_path+1]]),length(all_path_nodes$all_nodes[[length_of_path+1]]))
  P[[length_of_path+1]]=rep(1,length(all_path_nodes$all_nodes[[length_of_path+1]]))
  
  names(P[[length_of_path+1]])=all_path_nodes$all_nodes[[length_of_path+1]]
  final_P[[length_of_path+1]]=P[[length_of_path+1]]
  #compute p of other layers iteratively
  if(length_of_path>1){
    for(i in (length_of_path:2)){
      num_child=rowSums(all_adjacency_matrices$child_adjacency_matrices[[i]])
      #P[[i]]=rep(1/length(all_path_nodes$all_nodes[[i]]),length(all_path_nodes$all_nodes[[i]]))
      P[[i]]=rep(1,length(all_path_nodes$all_nodes[[i]]))
      names(P[[i]])=all_path_nodes$all_nodes[[i]]
      final_P[[i]]=(P[[i]][names(P[[i]])]+as.matrix(t(all_adjacency_matrices$child_adjacency_matrices[[i]])[names(P[[i]]),names(final_P[[i+1]])]) %*%
                (final_P[[i+1]]/num_child[names(final_P[[i+1]])]))[,1]
    }
  }
  return(unlist(final_P))
  }


ForwardAlg=function(all_adjacency_matrices,all_path_nodes){
  P=list()
  final_P=list()
  length_of_path=length(all_adjacency_matrices$child_adjacency_matrices)
  #start from last layer
  P[[length_of_path+1]]=rep(1/length(all_path_nodes$all_nodes[[length_of_path+1]]),length(all_path_nodes$all_nodes[[length_of_path+1]]))
  names(P[[length_of_path+1]])=all_path_nodes$all_nodes[[length_of_path+1]]
  final_P[[length_of_path+1]]=P[[length_of_path+1]]
  #compute p of other layers iteratively
  if(length_of_path>1){
  for(i in (length_of_path:2)){
    num_child=rowSums(all_adjacency_matrices$child_adjacency_matrices[[i]])
    num_co=colSums(all_adjacency_matrices$co_adjacency_matrices[[i]])
    num_child_of_co=rowSums(all_adjacency_matrices$child_adjacency_matrices[[i-1]])
    P[[i]]=(as.matrix(t(all_adjacency_matrices$child_adjacency_matrices[[i]])[,names(final_P[[i+1]])]) %*%
              (final_P[[i+1]]/num_child[names(final_P[[i+1]])]))[,1]
    final_P[[i]]=(P[[i]]+t(all_adjacency_matrices$co_adjacency_matrices[[i]][names(P[[i]]),names(P[[i]])])%*%
      (P[[i]]/pmax(1,num_co[names(P[[i]])]+num_child_of_co[names(P[[i]])])))[,1]
    #final_P[[i]]=final_P[[i]]/sum(final_P[[i]])
    }
  }
  return(unlist(final_P))
}

OldForwardAlg=function(all_adjacency_matrices,all_path_nodes){
  P=list()
  final_P=list()
  length_of_path=length(all_adjacency_matrices$child_adjacency_matrices)
  P[[length_of_path+1]]=rep(1/length(all_path_nodes$all_nodes[[length_of_path+1]]),length(all_path_nodes$all_nodes[[length_of_path+1]]))
  #P[[length_of_path+1]]=rep(1,length(all_path_nodes$all_nodes[[length_of_path+1]]))
  
  names(P[[length_of_path+1]])=all_path_nodes$all_nodes[[length_of_path+1]]
  final_P[[length_of_path+1]]=P[[length_of_path+1]]
  #compute p of other layers iteratively
  if(length_of_path>1){
    for(i in (length_of_path:2)){
      num_child=rowSums(all_adjacency_matrices$child_adjacency_matrices[[i]])
      P[[i]]=rep(1/length(all_path_nodes$all_nodes[[i]]),length(all_path_nodes$all_nodes[[i]]))
      #P[[i]]=rep(1,length(all_path_nodes$all_nodes[[i]]))
      names(P[[i]])=all_path_nodes$all_nodes[[i]]
      final_P[[i]]=(P[[i]][names(P[[i]])]+as.matrix(t(all_adjacency_matrices$child_adjacency_matrices[[i]])[names(P[[i]]),names(final_P[[i+1]])]) %*%
                      (final_P[[i+1]]/num_child[names(final_P[[i+1]])]))[,1]
    }
  }
  return(unlist(final_P))
}

BackwardAlg=function(all_adjacency_matrices,all_path_nodes,IF_of_target){
  IF=list()
  final_IF=list()
  length_of_path=length(all_adjacency_matrices$child_adjacency_matrices)
  #compute IF of downstream layer(first layer)
  IF[[1]]=IF_of_target
  final_IF[[1]]=IF[[1]]
  #compute IF for other layer iteratively
  for(i in (1:(length_of_path))){
    num_par=colSums(all_adjacency_matrices$child_adjacency_matrices[[i]])
    if( (i+1)< length(all_adjacency_matrices$co_adjacency_matrices)){
      num_co=colSums(all_adjacency_matrices$co_adjacency_matrices[[i+1]])
    }
    if(i!=length_of_path){
      num_par_of_co=colSums(all_adjacency_matrices$child_adjacency_matrices[[i+1]])
    }  
    IF[[i+1]]=((all_adjacency_matrices$child_adjacency_matrices[[i]][,names(num_par)])%*%
                         ((final_IF[[i]][names(num_par)])/pmax(1,num_par)))[,1] ### barresi konam ke chera
    if( (i+1)< length(all_adjacency_matrices$co_adjacency_matrices) ){
        final_IF[[i+1]]=IF[[i+1]]+(((all_adjacency_matrices$co_adjacency_matrices[[i+1]][names(IF[[i+1]]),names(IF[[i+1]])])%*%
          (IF[[i+1]]/pmax(1,num_co[names(IF[[i+1]])]+num_par_of_co[names(IF[[i+1]])]))))[,1]
     }
    else{
        final_IF[[i+1]]=IF[[i+1]]
    }
    if(length(final_IF[[i+1]])==1){
      names(final_IF[[i+1]])=row.names(all_adjacency_matrices$child_adjacency_matrices[[i]])
    }
  }
  

  #kinase_target=names(final_IF[[1]])[gsub("_.*","",names(unlist(final_IF[[1]]))) %in% c(Edges$GeneSymbol_a,Edges$GeneSymbol_b)]
  #kinase_target_scores=final_IF[[1]][kinase_target]
  #names(kinase_target_scores)=gsub("_.*","",names(kinase_target_scores))
  #return(c(unlist(final_IF)[!names(unlist(final_IF)) %in% names(unlist(final_IF[[1]]))],kinase_target_scores))
  return(unlist(final_IF[2:length(final_IF)]))
  
}

compute_scores=function(all_adjacency_matrices,all_path_nodes,IF_of_target){
  #all_p=ForwardAlg(all_adjacency_matrices,all_path_nodes)
  all_p=NewForwardAlg(all_adjacency_matrices,all_path_nodes)
  all_IF=BackwardAlg(all_adjacency_matrices,all_path_nodes,IF_of_target=IF_of_target)
  final_scores=abs(all_p) * all_IF[names(all_p)]
  sorted_scores=final_scores[names(list.sort(abs(final_scores)))]
  return(sorted_scores)
}

compute_scores_with_old_forward=function(all_adjacency_matrices,all_path_nodes,IF_of_target){
  all_p=OldForwardAlg(all_adjacency_matrices,all_path_nodes)
  all_IF=BackwardAlg(all_adjacency_matrices,all_path_nodes,IF_of_target=IF_of_target)
  final_scores=abs(all_p) * all_IF[names(all_p)]
  sorted_scores=final_scores[names(list.sort(abs(final_scores)))]
  return(sorted_scores)
}

compute_scores_without_Forward=function(all_adjacency_matrices,all_path_nodes,IF_of_target){
  all_IF=BackwardAlg(all_adjacency_matrices,all_path_nodes,IF_of_target=IF_of_target)
  final_scores=all_IF
  sorted_scores=final_scores[names(list.sort(abs(final_scores)))]
  return(sorted_scores)
}

NumDownstream=function(all_adjacency_matrices){
  
}