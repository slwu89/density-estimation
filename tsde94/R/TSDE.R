################################################################################
#   TSDE algorithm
#   R package by: Qu Cheng, Sean Wu, Robert Spear
#   Original TSDE algorithm written by Nong Shang, see:
#   Spear, R. C., Grieb, T. M., & Shang, N. (1994). Parameter uncertainty and interaction in complex environmental models. Water Resources Research. https://doi.org/10.1029/94WR01732
################################################################################

#' Generate Tree
#'
#' Using the output of \code{\link{TSDE}}, generate a tree
#' in igraph format
#'
#' @param df output of \code{\link{TSDE}} (a \code{\link{data.frame}} object)
#' @param pass a \code{\link{data.frame}} of the training data
#' @param varnames a character vector of variable names (columns of pass)
#' @param plot produce a plot of the tree?
#'
#' @export
TSDE_maketree <- function(df,pass,varnames,plot=FALSE){

  pass <- as.data.frame(pass)
  colnames(pass) <- varnames

  df$sleft <- as.character(df$sleft)
  df$srigt <- as.character(df$srigt)

  # get parent IDs
  df$parentID <- ifelse(df$NodeID == 1L, NA, df$NodeID %/% 2L)

  # 0s are leaf nodes; they are rectangles
  df$shape <- ifelse(df$var == 0L, "rectangle", "circle")
  df$color <- ifelse(df$var == 0L, "#EFE2BA", "#D1E8E2")

  # just a vector of IDs for the rows
  df$rowID <- 1:nrow(df)

  # replace the x's with variable names
  replace <- grep(pattern = "x",df$sleft)
  df$sleft[replace] <- sapply(df$sleft[replace],function(x,varnames){
    var <- regexpr(pattern = "[[:alpha:]]{1}[[:digit:]]{1,}",text = x)
    start <- 1L
    end <- attributes(var)$match.length
    num <- as.integer(substr(x = x,start = start+1L,stop = end))
    paste0(varnames[num],substr(x = x,start = end+1L,stop = nchar(x)))
  },USE.NAMES = FALSE,varnames=varnames)

  replace <- grep(pattern = "x",df$srigt)
  df$srigt[replace] <- sapply(df$srigt[replace],function(x,varnames){
    var <- regexpr(pattern = "[[:alpha:]]{1}[[:digit:]]{1,}",text = x)
    start <- 1L
    end <- attributes(var)$match.length
    num <- as.integer(substr(x = x,start = start+1L,stop = end))
    paste0(varnames[num],substr(x = x,start = end+1L,stop = nchar(x)))
  },USE.NAMES = FALSE,varnames=varnames)

  # list of passes
  pass_list <- list()
  pass_list[[1]] <- pass
  for(i in 2:nrow(df))
  {
    current_row <- df[i,]
    parent_row <- df[df$NodeID == current_row$parentID,]
    which_child <- ifelse(current_row$NodeID %% 2L == 0L, "left", "right")
    pass_parent_list <- pass_list[[parent_row$rowID]]
    if(which_child == "left")
    {
      pass_list[[i]] <- pass_parent_list[eval(parse(text = paste("pass_parent_list$", parent_row$sleft, sep = ""))), ]
    } else
    {
      pass_list[[i]] <- pass_parent_list[!eval(parse(text = paste("pass_parent_list$", parent_row$sleft, sep = ""))), ]
    }
  }

  df$PassNo <- unlist(lapply(pass_list, nrow))

  nodes <- data.frame(name = df$NodeID)
  relations <- data.frame(from = df$parentID[-1], to = df$NodeID[-1])
  g <- igraph::graph.data.frame(relations, directed = FALSE, vertices = nodes)
  
  if(plot){
    plot(g, 
         layout = igraph::layout_as_tree(g, root = 1), 
         vertex.label = round(df$inside,3), 
         vertex.shape = df$shape, 
         vertex.size = 16, 
         vertex.color = df$color, 
         vertex.label.color = "black", 
         edge.label = df$PassNo[-1], 
         edge.label.color = "black", 
         edge.label.dist = 20)
  }

  return(g)
}
