#' ontoNW
#'
#' \strong{Ontology Network}. Build gene-ontology relationship static and dynamic network.
#' @param onto data.frame. Includes an 'id' column with term id (e.g. G0:0008150), a 'type' colum (e.g. BP, KEGG), a 'name' column (e.g. biological_process), a 'genes' column with string of genes of the query belonging to the ontology term (e.g. IL1RL1, IL5RA, TP53), and a 'fdr' column with FDR values.
#' @param cutoff numeric. A fdr value as cutoff threshold to plot edges/nodes below cutoff
#' @param dynamic boolean. Produce a dynamic network.
#' @param pdf boolean. Saves the static network as pdf.
#' @param html boolean. Saves the dynamic network as html. Applies only if dynamic=T.
#' 
#' @return nodes. nodesCut. Data frames with nodes id, pval, name, and type. Cutoff applies to nodesCut.
#' @return edges. edgesCut. Data frames with edges from node, to node, and weight (p-value). Cutoff applies to edgesCut.
#' @return net. netCut. igraph object. Cutoff applies to netCut.
#' @return DynNetwork. visNetwork object and htmlwidget. Cutoff applies.
#' 
#' 
#' @details Layout are lgl for static network and mds for dynamic network. Read igraph and visNetwork documentation.
#' 
#' @import igraph
#' @import scales
#' @import visNetwork
#' 
#' @author Nathan Lemonnier \email{nathanael.lemonnier@@gmail.com}
#' 
#' @examples
#' onto=data.frame(id  = c("GO:0002376", "GO:0032634")
#'                 type= c("BP","BP")
#'                 name= c("immune system process","interleukin-5 production")
#'                 gene= c(c(IL5RA,SMPD3,CLC,IL1RL1,CYSLTR2,ALOX15,PIK3R6), c(IL5RA,IL1RL1))
#'                 fdr = c(0.005, 0.00005)
#'                 )
#' ontoNW(onto, cutoff=0.05, dynamic=T, pdf=F, html=T)       

#function
ontoNW <- function(onto, cutoff=0.05, pdf=T, html=T, dynamic=T){
  cat("loading required packages...\n")
  cat("   igraph.\n")
  cat("   scales...\n")

  #required packages  
  library(igraph) #for function sumlog
  library(scales) #for function barplot2

  cat("...creating directories...\n")
  #create directories
  
  cat("...building nodes table...\n")
  
  #list genes
  genelist <- unique(unlist(strsplit(as.character(onto$genes), ",")))
  
  colnames(onto) <- tolower(colnames(onto))
  
  #nodes table
  nodes <- data.frame(
    id =   c(
      genelist
    , as.character(onto$id)
    )
  , pval = c(
      rep(NA,length(genelist))
    , as.character(onto$BH.FDR.p.value)
  )
  , name = c(
      genelist
    , as.character(onto$name)
  )
  , type = c(
      rep("gene",length(genelist))
    , as.character(onto$type)
  )
    )
  
  cat("...building edges table...\n")
  #edges table
  edges <- t(read.table(text=as.character(onto$genes[1]), sep=",", header=F, colClasses = "character"))
  for (i in 2:nrow(onto)){
    edges <- rbind(edges,
                       t(read.table(text=as.character(onto$genes[i]), sep=",", header=F, colClasses = "character")))
  }
  edges <- cbind(edges,
                     matrix(,nrow=nrow(edges),ncol=2))
  colnames(edges) <- c("from","to","weight")
  edges <- as.data.frame(edges)
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$weight <- as.numeric(edges$weight)
  
  ngenes <- ncol(read.table(text=as.character(onto$genes[1]), sep=","))
  for (i in 2: nrow(onto)){
    ngenes <- rbind(ngenes,
                    ncol(read.table(text=as.character(onto$genes[i]), sep=",")))
  }
  onto <- cbind(onto, ngenes)
  
  edges$to[1:sum(onto$ngenes[1])] <- as.character(onto$id[1])
  for (i in 2:nrow(onto)){
    edges$to[(sum(onto$ngenes[1:i-1])+1):sum(onto$ngenes[1:i])] <- as.character(onto$id[i])
  }
  edges$weight[1:sum(onto$ngenes[1])] <- as.character(onto$fdr[1])
  for (i in 2:nrow(onto)){
    edges$weight[(sum(onto$ngenes[1:i-1])+1):sum(onto$ngenes[1:i])] <- as.character(onto$fdr[i])
  }
  rownames(edges) <- NULL
  class(edges$weight) <- "numeric"
  
  # build network
  net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
  V(net)$label <- names(V(net))
  l <- layout_with_lgl(net)
  l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  nodesgenes <- which(V(net)$label %in% nodes$name[which(nodes$type=='gene')]==T)
  nodesonto  <- which(V(net)$label %in% nodes$name[which(nodes$type!='gene')]==T)
  V(net)$size[nodesonto ] <- -log2(as.numeric(V(net)$pval[nodesonto]))/2
  V(net)$size[nodesgenes] <- 10
  

  #network and plot with FDR cutoff
  net.cut <- delete_edges(net, E(net)[weight>cutoff])
  V(net.cut)$degree[1:17] <- as.data.frame(table(as_edgelist(net.cut)[,1]))[match(V(net.cut)$label[1:17],as.character(as.data.frame(table(as_edgelist(net.cut)[,1]))[,1])),2]
  V(net.cut)$degree[18:nrow(nodes)] <- as.data.frame(table(as_edgelist(net.cut)[,2]))[match(V(net.cut)$label[18:nrow(nodes)],as.character(as.data.frame(table(as_edgelist(net.cut)[,2]))[,1])),2]
  colrs <- cbind(as.character(as.data.frame(table(V(net.cut)$type))[,1]), rbind("skyblue","skyblue2","skyblue4","yellow3","lightgreen","red"))
  V(net.cut)$col <- colrs[match(V(net.cut)$type,colrs[,1]),2]
  net.cut <- delete_vertices(net.cut, V(net.cut)[is.na(V(net.cut)$degree)])
  
  
  l.cut <- layout_with_lgl(net.cut)
  l.cut <- norm_coords(l.cut, ymin=-1, ymax=1, xmin=-1, xmax=1)
  if(pdf==T){
    pdf("StaticNetwork_lgl.pdf", width = 10, height = 10, compress = F)
    plot(net.cut,
         edge.arrow.size=.4,
         edge.curved=0,
         edge.color=alpha("black",0.1),
         vertex.color=V(net.cut)$col,
         vertex.frame.color=alpha("white",0.3),
         vertex.label=V(net.cut)$description,
         vertex.label.color="black",
         vertex.label.cex=.7,
         vertex.size=V(net.cut)$degree/1.5,
         rescale=T,
         layout=l.cut*1.2)
    title(main=paste("Functional enrichment with FDR < ",cutoff, sep=""))
    legend("topleft",
           legend=colrs[,1], pch=21,
           col=alpha("black",0.3),
           pt.bg=colrs[,2],
           pt.cex=2.5,
           bty="n",
           ncol=1)
    legend("bottomleft",
           legend=c(paste("max degree",max(V(net.cut)$degree),sep="=")
                    ,paste("min degree",min(V(net.cut)$degree),sep="=")),
           pch="",
           col="black",
           pt.bg=alpha("white",0.1),
           #pt.cex=c(max(V(net.cut)$degree/2),min(V(net.cut)$degree/2)),
           bty="n",
           ncol=1)
        dev.off()
    
  }else{
    plot(net.cut,
         edge.arrow.size=.4,
         edge.curved=0,
         edge.color=alpha("black",0.1),
         vertex.color=V(net.cut)$col,
         vertex.frame.color=alpha("white",0.3),
         vertex.label=V(net.cut)$description,
         vertex.label.color="black",
         vertex.label.cex=.7,
         vertex.size=V(net.cut)$degree/1.5,
         rescale=T,
         layout=l.cut*1.2)
    title(main=paste("Functional enrichment with FDR < ",cutoff, sep=""))
    legend("topleft",
           legend=colrs[,1], pch=21,
           col=alpha("black",0.3),
           pt.bg=colrs[,2],
           pt.cex=2.5,
           bty="n",
           ncol=1)
    legend("bottomleft",
           legend=c(paste("max degree",max(V(net.cut)$degree),sep="=")
                    ,paste("min degree",min(V(net.cut)$degree),sep="=")),
           pch="",
           col="black",
           pt.bg=alpha("white",0.1),
           #pt.cex=c(max(V(net.cut)$degree/2),min(V(net.cut)$degree/2)),
           bty="n",
           ncol=1)
  }

  
  #dynamic network and plot with FDR cutoff  
  if(dynamic==T){
  library("visNetwork")
  colnames(nodes)[1] <- "id"
  
  edgesCut <- edges[which(edges$weight<cutoff),]
  nodesCut <- nodes[match(unique(edgesCut$from),nodes$id),]
  nodesCut <- rbind(nodesCut, nodes[match(unique(edgesCut$to),nodes$id),])
  nodesCut$shape <- "dot"  
  nodesCut$shadow <- TRUE # Nodes will drop shadow
  nodesCut$title <- nodesCut$description # Text on click
  nodesCut$label <- nodesCut$id # Node label
  
  nodesgenesCut <- which(nodesCut$name %in% nodes$name[which(nodes$type=='gene')]==T)
  nodesontoCut  <- which(nodesCut$name %in% nodes$name[which(nodes$type!='gene')]==T)
  
  nodesCut$size <- c(
    as.data.frame(table(edgesCut$from)[match(nodesCut$name[nodesgenesCut]
                                            ,as.character(as.data.frame(table(edgesCut$from))$Var1)
                                            )
                                     ]
                         )$Freq
  , as.data.frame(table(edgesCut$to)[match(nodesCut$id[nodesontoCut]
                                             ,as.character(as.data.frame(table(edgesCut$to))$Var1)
                                            )
                                       ]
    )$Freq
  )
  nodesCut$borderWidth <- 0 # Node border width
  
  nodesCut$color.background <- colrs[match(nodesCut$type,colrs[,1]),2]
  nodesCut$color.border <- alpha("white",0.1)
  nodesCut$color.highlight.background <- alpha("white",0.1)
  nodesCut$color.highlight.border <- alpha("white",0.1)
  
  edgesCut$color <- alpha("black",0.3)    # line color  
  
  DynNetwork <- visNetwork(nodesCut, edgesCut, main=paste("Functional enrichment with BH FDR < ",cutoff, " (mds layout)", sep="")) %>%
    visOptions(highlightNearest = TRUE, 
               selectedBy = "type") %>%
    visIgraphLayout(layout = "layout_with_lgl", physics = F,
                    smooth = F, type = "full", randomSeed = NULL,
                    layoutMatrix = NULL)
  if(html==T){
  visSave(DynNetwork
          , file="DynamicNetwork_mds.html")
  }else{DynNetwork}
  }else{}
  
}