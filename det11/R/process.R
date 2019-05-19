#' Process DET output
#'
#' Process a fitted density estimation tree; return the modified data frame
#' and necessary components to plot the tree with igraph.
#'
#' @param x a fitted tree (output of \code{\link[det11]{ml_det}})
#'
#' @export
DET_maketree <- function(x){

  ids <- rank(unique(x$tree$ID))
  orig_ids <- unique(x$tree$ID)

  sort_id <- rbind(cbind(ids,orig_ids),c(0,0))

  tree <- x$tree
  tree$ID <- ids
  tree$parent_ID <- sort_id[match(tree$parent_ID,sort_id[,2]),1]

  nodes <- data.frame(name = tree$ID,stringsAsFactors = F)
  relations <- data.frame(from = tree$parent_ID[-1], to = tree$ID[-1],stringsAsFactors = F)

  # need to find edge split
  edge_labels <- mapply(FUN = function(from,to){
    row <- which(tree$parent_ID == from & tree$ID == to)
    dir <- as.character(tree[row,"parent_split_edge"])
    var <- tree[which(tree$ID == from),"variable"]
    ineq <- ifelse(dir == "left"," <= "," > ")
    value <- tree[which(tree$ID == from),which(names(tree) == dir)]
    paste0(from,"->",to,": x",var,ineq,value)
  },from=relations$from,to=relations$to)

  node_shape <- ifelse(tree$leaf == 0,"circle","rectangle")
  node_col <- ifelse(tree$leaf == 0,"#D1E8E2","#EFE2BA")
  g <- igraph::graph.data.frame(relations, directed = FALSE, vertices = nodes)

  plot(g,
       layout = igraph::layout_as_tree(g, root = 1),
       vertex.label = paste0(tree$ID,"\n",round(tree$below,4)),
       vertex.shape = node_shape,
       vertex.size = 12,
       vertex.color = node_col,
       vertex.label.color = "black",
       edge.label = edge_labels,
       edge.label.color = "black",
       edge.label.dist = 20)

  return(
    list(
      tree=tree,
      nodes=nodes,
      relations=relations,
      edge_labels=edge_labels,
      node_shape=node_shape,
      node_color=node_col,
      g=g
    )
  )
}
