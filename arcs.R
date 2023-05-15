#install.packages("BiocManager")
#BiocManager::install("ggtreeExtra")
library("ggtree")
library("ggtreeExtra")
library("treeio")
library("ggplot2")
library(tibble)
library(caper)
library(dplyr)
#library(RColorBrewer)


args = commandArgs(trailingOnly=TRUE)

#arc.pos=0.0011
#label.pos=0.0012
arc.pos=0.0010
label.pos=0.0011

col_vector=c(
 "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
 "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
 "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
 "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
 "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3",
 "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
 "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9",
 "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8")
# "#4DAF4A",
# "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
# "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
# "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
# "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
# "#CCEBC5", "#FFED6F")

genera.ordering=c(
 "Abactochromis",     "Alticorpus",        "Aristochromis",    
 "Astatotilapia",     "Aulonocara",        "Buccochromis",     
 "Champsochromis",    "Cheilochromis",     "Chilotilapia",     
 "Chindongo",         "Copadichromis",     "Corematodus",      
 "Ctenopharynx",      "Cyathochromis",     "Cynotilapia",      
 "Cyrtocara",         "Dimidiochromis",    "Diplotaxodon",     
 "Fossorochromis",    "Genyochromis",      "Gephyrochromis",   
 "Hemitaeniochromis", "Hemitilapia",       "Labeotropheus",    
 "Labidochromis",     "Lethrinops",        "Maylandia",        
 "Mchenga",           "Melanochromis",     "Mylochromis",      
 "Naevochromis",      "Nimbochromis",      "Otopharynx",       
 "Pallidochromis",    "Petrotilapia",      "Placidochromis",   
 "Protomelas",        "Pseudotropheus",    "Rhamphochromis",   
 "Sciaenochromis",    "Stigmatochromis",   "Taeniochromis",    
 "Taeniolethrinops",  "Tramitichromis",    "Trematocranus",    
 "Tropheops",         "Tyrannochromis")

names(col_vector)=genera.ordering

chr.name=args[1]
# filename for meta data
meta="./inversion_man_mt_2022-11-16.tsv"
# filename for tree
tree.file=paste0("./inversion_samples_noninverted_pwd_nj_tree_ancestral_rooted.newick.nwk")
#pdf.file=paste0("C:\\work\\Fishes\\ggtree\\arcs_", chr.name, ".pdf")

#tree.file="C:\\work\\Fishes\\ggtree\\OneDrive_1_19.12.2022\\inversion_haplotype_jn_tree_chr2.newick"
#pdf.file="C:\\work\\Fishes\\ggtree\\arcs.pdf"

# read tree data
tree=read.newick(tree.file)
#    
p=ggtree(tree,lwd=0.1)

tips.df=p$data
tip.labels=p$data$label
tip.labels=tip.labels[!is.na(tip.labels)]

individuals=unique(unlist(lapply(tip.labels,function(x) substr(x,0,nchar(x)-2))))

#arc.pos=0.00098
tip.coord=lapply(
  individuals,
  function(x){
    first=tips.df[tips.df$label==paste0(x,"00"),][1,]
    second=tips.df[tips.df$label==paste0(x,"_1"),][1,]
    if (sum(is.na(first))>0 || sum(is.na(second))>0){
      return(NA)
    }
    if (first$y<second$y){
      buffer=first
      first=second
      second=buffer
      }
    return(list(  
                f.x= arc.pos,
                f.y=first$y,
                s.x= arc.pos,
                s.y=second$y,
                f.node=first$node,
                s.node=second$node,
                label=substr(second$label,0,nchar(second$label)-2)
                ))
    
  }
)

arc.df=data.frame()
for(tip.pair in tip.coord){
  if(is.na(tip.pair)){
    next
  }
  arc.df=rbind(arc.df,tip.pair)
} 

#meta.df=read.delim(meta,sep="\t",header=T,as.is=T)
meta.df.full=read.delim(meta,sep="\t",header=T,as.is=T)
meta.df=meta.df.full[c("X","clade","simple_id")]
meta.df$clade=factor(meta.df$clade)
merged.df=merge(arc.df, meta.df, by.x="label", by.y="X")

p$data$label=unlist(lapply(p$data$label,function(x) {
  parts=strsplit(x,"_")[[1]]
  ind=parts[1]
  ht=parts[2]
  return(paste0(merged.df[merged.df$label==ind,'simple_id'][1],"_",ht))
}))



#
# mind the layers order! Tips should be layers[[3]] for patch to work - see further
#

#pic = p + geom_tiplab(size=1, geom="text", align=T, lwd=1) + 
#  geom_curve(
#    aes(x = f.x, y = f.y, xend = s.x, yend = s.y, colour= clade),
#    curvature = -1, data = merged.df) +
#  geom_hilight(data=mark.nodes, mapping=aes(node=node,fill=genus)) +
#  xlim(c(0,0.005)) +
#  scale_colour_manual(values = c("Shallow"= "#ff6247",  "Deep" = "#4876ff",  "Utaka" = "#006400")) +
#  theme(legend.position = "bottom")

############################ HERE WE ARE PATCHING ggtree ON THE FLY PROBABLY NOT NECESSARY FOR YOUR VERSION!###
# fixes this:                                                                                                 #
# 'gpar' element 'lwd' must not be length 0                                                                   #
###############################################################################################################

library(grid)

empty <- getFromNamespace("empty", "ggplot2")
`%||%` <- getFromNamespace("%||%", "ggplot2")
draw_panel_deb = function(data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
                          lineend = "butt", linejoin = "round", na.rm = FALSE, nudge_x = 0) {
  
  data$x <- data$x + nudge_x
  
  data <- ggplot2::remove_missing(data, na.rm = na.rm, c("x", "y", "xend",
                                                         "yend", "linetype", "linewidth", "shape"), name = "geom_segment")
  if (empty(data))
    return(zeroGrob())
  if (!coord$is_linear()) {
    tmpgroup <- data$group
    starts <- subset(data, select = c(-xend, -yend))
    starts$group <- 1
    ends <- rename(subset(data, select = c(-x, -y)), c("x" = "xend", "y" = "yend"))
    ends$group <- 2
    pieces <- rbind(starts, ends)
    
    trans <- coord$transform(pieces, panel_params)
    starts <- trans[trans$group==1, ,drop=FALSE]
    ends <- trans[trans$group==2, ,drop=FALSE]
    ends <- rename(subset(ends, select=c(x, y)), c("xend"="x", "yend"="y"))
    data <- cbind(starts, ends)
    data$group <- tmpgroup
  }else{
    data <- coord$transform(data, panel_params)
  }
  
  arrow.fill <- arrow.fill %||% data$colour
  return(grid::segmentsGrob(data$x, data$y, data$xend, data$yend,
                            default.units = "native", gp = gpar(col = alpha(data$colour,
                                                                            data$alpha), fill = alpha(arrow.fill, data$alpha),
                                                                lwd = data$linewidth * ggplot2::.pt, lty = data$linetype,
                                                             lineend = lineend, linejoin = linejoin), arrow = arrow)
  )   
  
  
  ## data$x <- data$x - sapply(data$label, function(x) convertWidth(grobWidth(textGrob(x, gp=gpar(fontsize=.04* .pt))), "native", TRUE))
  ##GeomSegment$draw_panel(data = data, panel_params = panel_params, coord = coord,
  ##                       arrow = arrow, arrow.fill = arrow.fill,
  ##                       lineend = lineend, linejoin = linejoin, na.rm = na.rm)
}

#pic$layers[[3]]$geom$draw_panel=draw_panel_deb

###################################################### END OF PATCH! ##############################################

pic2 = p + 
  geom_curve(
    aes(x = f.x, y = f.y, xend = s.x, yend = s.y, colour= clade),
    curvature = -1, data = merged.df) +
  xlim(c(0,0.005)) +
  scale_colour_manual(values = c("Shallow"= "#ff6247",  "Deep" = "#4876ff",  "Utaka" = "#006400")) +
  theme(legend.position = "bottom")


depth.table=node.depth(tree, method = 1)
int.nodes=seq(length(tree$tip.label):length(tree$tip.label)+tree$Nnode)

monophyletic.group<-function(n){
  tip.list=clade.members(n, tree, tip.labels = T, include.nodes=FALSE)
  if(length(tip.list)<2){
    return(NA)
  }
  ids=unique(unlist(lapply(tip.list,function(x) substr(x,0,nchar(x)-2))))
  genera=meta.df.full[meta.df.full$X %in% ids,"genus"]
  if(length(unique(genera))==1){
    return(c(n,genera[1],depth.table[n]))
  }else{
    return(NA)
  }
}

mono.nodes.list=lapply(int.nodes,monophyletic.group)
mono.nodes.list=mono.nodes.list[!is.na(mono.nodes.list)]
mono.nodes.df=as.data.frame(do.call(rbind, mono.nodes.list))
colnames(mono.nodes.df)=c('node','genus','depth')
mono.nodes.df$node=as.numeric(mono.nodes.df$node)


desc.nodes=function(x){
  node.descendants=clade.members(x, tree, tip.labels = T, include.nodes=T)
  if(length(node.descendants$nodes)>1){
    return(tail(node.descendants$nodes,-1))
  }else{
    return(NA)
  }
}

non.nested.nodes=function(genus){
  rows=mono.nodes.df[mono.nodes.df$genus==genus,]
  internal.nodes=unique(unlist(lapply(rows$node,desc.nodes)))
  internal.nodes=internal.nodes[!is.na(internal.nodes)]
  return(rows[!rows$node %in% internal.nodes,])
}

mark.nodes=data.frame()
for (genus in unique(mono.nodes.df$genus)){
  mark.nodes=rbind(mark.nodes, non.nested.nodes(genus))
}

mark.nodes$node=as.integer(mark.nodes$node)
#mark.nodes$genus=as.factor(mark.nodes$genus)

mark.coords=p$data[mark.nodes$node,c('x','y','node')]
mark.merged=merge(mark.coords,mark.nodes,by='node')

#n <- length(unique(mark.nodes$genus))
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

for(row.n in 1:nrow(mark.nodes)){

  node.row=mark.nodes[row.n,]
  pic2 <- collapse(
    pic2,
    node = node.row$node ,
    mode = "max",
    clade_name = node.row$genus,
    alpha = 0.8,
    color = col_vector[node.row$genus],
    fill = col_vector[node.row$genus]) 
#    geom_cladelabel(node=node.row$node,
#                    label=node.row$genus,
#                    align=T, 
#                    color='black',
#                    fontsize=3)
}




#pdf(pdf.file,40,60)
pic2 + geom_text(aes(label=genus,color=genus, x=label.pos),data = mark.merged,size=3)
dev.off()

#?geom_text
#pic2 
#mark.nodes$node=as.integer(mark.nodes$node)
#mark.nodes.temp=mark.nodes[1:10,]
#mark.nodes.temp$node=sample(seq(1,n.fake)+n.fake,10)
#mark.nodes.temp
#p=ggtree(tree,lwd=0.1)
#p + geom_hilight(data=mark.nodes, mapping=aes(node=node,fill=genus)) +
#  theme(legend.position = "bottom")
