#' Process DET output
#'
#' Process a fitted density estimation tree; return the modified data frame
#' and necessary components to plot the tree with igraph.
#'
#' @param x a fitted tree (output of \code{\link[det11]{ml_det}})
#' @param dataset rows must be variables, columns are observations
#'
#' @export
DET_maketree <- function(x, dataset){

  ids <- rank(unique(x$tree$ID))
  orig_ids <- unique(x$tree$ID)

  sort_id <- rbind(cbind(ids,orig_ids),c(0,0))

  tree <- x$tree
  tree$ID <- ids
  tree$parent_ID <- sort_id[match(tree$parent_ID,sort_id[,2]),1]

  # flip right and left node
  tree <- tree[order(tree$ID),]
  
  nodes <- data.frame(name = tree$ID,stringsAsFactors = F)
  relations <- data.frame(from = tree$parent_ID[-1], to = tree$ID[-1],stringsAsFactors = F)
  
  

  # # need to find edge split
  # edge_labels <- mapply(FUN = function(from,to){
  #   row <- which(tree$parent_ID == from & tree$ID == to)
  #   dir <- as.character(tree[row,"parent_split_edge"])
  #   var <- tree[which(tree$ID == from),"variable"]
  #   ineq <- ifelse(dir == "left"," <= "," > ")
  #   value <- tree[which(tree$ID == from),which(names(tree) == dir)]
  #   paste0(from,"->",to,": x",var,ineq,value)
  # },from=relations$from,to=relations$to)
  
  # need to find edge split
  split.equ <- mapply(FUN = function(from,to){
    row <- which(tree$parent_ID == from & tree$ID == to)
    dir <- as.character(tree[row,"parent_split_edge"])
    var <- tree[which(tree$ID == from),"variable"]
    ineq <- ifelse(dir == "left"," <= "," > ")
    value <- tree[which(tree$ID == from),which(names(tree) == dir)]
    paste0(from,"->",to,": x",var,ineq,value)
  },from=relations$from,to=relations$to)
  

  node_shape <- ifelse(tree$leaf == 0,"circle","square")
  node_col <- ifelse(tree$leaf == 0,"#D1E8E2","#EFE2BA")
  g <- igraph::graph.data.frame(relations, directed = FALSE, vertices = nodes)

  plot(g,
       layout = igraph::layout_as_tree(g, root = 1),
       vertex.label = round(tree$below,4),
       vertex.shape = node_shape,
       vertex.size = 12,
       vertex.color = node_col,
       vertex.label.color = "black",
       # edge.label = edge_labels,
       edge.label = round(merge(relations, tree, by.x = "to", by.y = "ID")$below*ncol(dataset),0),
       edge.label.color = "black",
       edge.label.cex = 0.7)
  
  # add the threshold and area
  # get the coordinates for each box
  coordinates <- igraph::layout_as_tree(g, root = 1)
  coordinates[,1] <- scales::rescale(coordinates[,1], to = c(-1,1))
  coordinates[,2] <- scales::rescale(coordinates[,2], to = c(-1,1))
  
  variable.name <- if(!is.null(rownames(dataset))) rownames(dataset) else paste("V",1:nrow(dataset), sep = "")
  
  # annotate the split criteria for intermediate node
  for(i in 1:nrow(coordinates))
  {
    if(tree$leaf[i] == 0)
    {
      label = paste(variable.name[tree$variable[i]], "<=", round(tree$right[i],3), sep = " ")
      text(coordinates[i,1], coordinates[i,2]-0.1, label, cex = 0.8)
    } 
  }
  

  return(
    list(
      tree=tree,
      nodes=nodes,
      relations=relations,
      edge_labels=round(merge(relations, tree, by.x = "to", by.y = "ID")$below*ncol(dataset),0),
      node_shape=node_shape,
      node_color=node_col,
      g=g
    )
  )
}
