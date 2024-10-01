library(readr)
library(tidyr)
library(dplyr)
library(igraph)
library(plyr)
library(ComplexHeatmap)
library(ggplotify)
library(grid)
library(ggpubr)
library(geosphere)
library(lubridate)
#library(circlize)

## concordance heatmap function from Katie Evans
concordance_heatmap <- function(strains, version) {
  
  if (length(strains) > 2) {
    
    if (version == "new") {
      gtcheck_table <- pwcc
      color_discordance = "royalblue1"
    }
    
    if (version == "old") {
      gtcheck_table <- gtcheck_old
      color_discordance = "yellow"
    }
    
    df1 <- gtcheck_table %>%
      dplyr::select(i, j, concordance) %>% 
      dplyr::filter(i %in% strains & j %in% strains) %>% 
      unique()
    
    df2 <- df1 %>% 
      dplyr::rename(i=j, j=i)
    
    df <- dplyr::bind_rows(df1, df2) %>% 
      unique() %>% 
      tidyr::spread(j, concordance) 
    
    if (nrow(df) > 2) {   
      
      if (nrow(df) <= 10) { font_size = 12}
      if (nrow(df) > 10 & nrow(df) < 15) { font_size = 10}
      if (nrow(df) >= 15 & nrow(df) < 30) { font_size = 8}
      if (nrow(df) >= 30 ) { font_size = 6}
      
      # make heatmap of concordance
      m = select(df, -i) %>% as.matrix()
      rownames(m) = df$i
      
      col_fun = circlize::colorRamp2(c(min(m, na.rm=T), max(m, na.rm=T)), c(color_discordance, "indianred1"))
      
      p  = ComplexHeatmap::Heatmap(m, col=col_fun, rect_gp = gpar(col = "white", lwd = 2), cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.5f", m[i, j]), x, y, gp = gpar(fontsize = font_size)) 
      })
      
      print(p)
    } else { print("Not more than 2 strains in group in this dataset")}
    
  } else { print("Not more than 2 strains in group in this dataset")}
}

## adaptation of Katie's function to plot discordant sites
discordance_heatmap <- function(strains, version) {
  
  if (length(strains) > 2) {
    
    if (version == "new") {
      gtcheck_table <- pwcc
      color_discordance = "royalblue1"
    }
    
    if (version == "old") {
      gtcheck_table <- gtcheck_old
      color_discordance = "yellow"
    }
    
    df1 <- gtcheck_table %>%
      dplyr::select(i, j, discordance) %>% 
      dplyr::filter(i %in% strains & j %in% strains) %>% 
      unique()
    
    df2 <- df1 %>% 
      dplyr::rename(i=j, j=i)
    
    df <- dplyr::bind_rows(df1, df2) %>% 
      unique() %>% 
      tidyr::spread(j, discordance) 
    
    if (nrow(df) > 2) {   
      
      if (nrow(df) <= 10) { font_size = 12}
      if (nrow(df) > 10 & nrow(df) < 15) { font_size = 10}
      if (nrow(df) >= 15 & nrow(df) < 30) { font_size = 8}
      if (nrow(df) >= 30 ) { font_size = 6}
      
      # make heatmap of concordance
      m = select(df, -i) %>% as.matrix()
      rownames(m) = df$i
      
      col_fun = circlize::colorRamp2(c(min(m, na.rm=T), max(m, na.rm=T)), c(color_discordance, "indianred1"))
      
      p  = ComplexHeatmap::Heatmap(m, col=col_fun, rect_gp = gpar(col = "white", lwd = 2), cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.5f", m[i, j]), x, y, gp = gpar(fontsize = font_size)) 
      })
      
      print(p)
    } else { print("Not more than 2 strains in group in this dataset")}
    
  } else { print("Not more than 2 strains in group in this dataset")}
}

#this function takes a graph with detected subcommunities and splits them into individual graphs
induceNS <- function(graphL) {
  mems <- length(unique(graphL[[2]]$membership))
  subgraphs <- list()
  for (mem in 1:mems) {
    sub <- induced_subgraph(graphL[[1]], which(membership(graphL[[2]]) == mem))
    subgraphs[[mem]] <- sub
  }
  return(subgraphs)
}

setwd("/vast/eande106/projects/Nicolas/github/concordance-nf/graph_dev/")
#read gtcheck (pairwise concordance) file

#briggsae
#pwcc <- readr::read_tsv("/vast/eande106/projects/Nicolas/collabs/forNikita/graph_analysis/gtcheck.tsv") %>%
#elegans
pwcc <- readr::read_tsv("raw_data/elegans/gtcheck.tsv") %>%
  dplyr::mutate(concordance=(sites-discordance)/sites) %>%
  dplyr::mutate(discordant=discordance) %>%
  dplyr::mutate(discordance=discordant/sites)

#get sampling metadata from species sheet
#briggsae
#winfo <- readr::read_tsv("/vast/eande106/projects/Nikita/c_briggsae/cbriggsae_popgen/raw_data/c_briggsae_species_sheet_2024-01-29.tsv") %>%
#elegans
winfo <- readr::read_tsv("raw_data/elegans/c_elegans_species_sheet_2024-02-17.tsv") %>%
  dplyr::select(strain,latitude,longitude,source_lab,sampling_date)

#join gtcheck and metadata & estimate pairwise haversine distance
reldat <- pwcc %>% 
  dplyr::left_join(winfo,by = c("i"="strain")) %>%
  dplyr::rename("lat.i"=latitude,"lon.i"=longitude,"source.i"=source_lab,"date.i"=sampling_date) %>%
  dplyr::left_join(winfo,by = c("j"="strain")) %>%
  dplyr::rename("lat.j"=latitude,"lon.j"=longitude,"source.j"=source_lab,"date.j"=sampling_date) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(dist=distHaversine(c(lon.i,lat.i),c(lon.j,lat.j))) %>%
  dplyr::ungroup() 

#estimate pairwise time difference in days and sampling lab matches; keep only relevant columns
reldat_adj <- reldat %>%
  dplyr::mutate(timedif=abs(as.numeric(difftime(ymd(date.i),ymd(date.j),units="days")))) %>%
  dplyr::mutate(labdiff=ifelse(source.i==source.j,T,F)) %>%
  dplyr::select(-lat.i,-lon.i,-lat.j,-lon.j,-date.i,-date.j,-source.i,-source.j,-avg_min_depth,-sites) 

#set base edge threshold
thrsh <- 0.9995

#filter edges by thresh
pass <- reldat_adj %>% dplyr::filter(concordance > thrsh)

#extract vertices
verts <- unique(c(pwcc$i,pwcc$j))

#get 2nd minimum distance between two samples
min2dist <- nth(sort(unique(pass$dist)),2)

#select for relevant columns, set concordance as the weight variable, replace distances of 0 with the 2nd minimum distance
rels <- pass %>% 
  dplyr::select(i,j,concordance,discordant,dist,timedif,labdiff) %>%
  dplyr::rename(from=i,to=j,weight=concordance) %>%
  #dplyr::mutate(weight=((weight-min(weight))/(max(weight)-min(weight)))) %>%
  #dplyr::mutate(weight=ifelse(weight==0,nth(weight,2,order_by = weight),weight)) %>%
  dplyr::mutate(dist=ifelse(dist==0,min2dist,dist))

#build and plot graph from edges and vertices
g <- graph_from_data_frame(rels, directed=F, vertices=verts)
plot(g,layout=layout_with_mds,vertex.size=1,vertex.label=NA)

#use cluster_fg to draw subcommunities from presence/absence of edge
c <- cluster_fast_greedy(g)
plot(c,g,layout=layout_with_mds,mark.border="black",mark.col="pink",vertex.size=1,vertex.label=NA)

#decompose subcomunities into list of graphs
gd <- decompose(g)
#length(gd)

#empty list and counter
dfl <- list()
rell <- list()
#loop through list of subcommunities, normalise weight, replace weights of 0 with 2nd minimum weight
#subcommunity data frames are stored in 'dfl'
#subcommunity graphs are stored in 'rell' (1-vertex subcommunities are NULL)
for (i in 1:length(gd)) {
  print(i)
  names <- vertex.attributes(gd[[i]])$name
  count <- length(names)
  tmpdf <- data.frame(members=paste(names,collapse = ","),nmem=count,gnum=i)
  dfl[[i]] <- tmpdf
  if (count>1) {
    print(count)
    tmptbl <- rels %>% 
      dplyr::filter((from %in% names) & (to %in% names)) %>%
      dplyr::mutate(concordance=weight) %>%
      dplyr::mutate(weight=((weight-min(weight))/(max(weight)-min(weight)))) %>%
      dplyr::mutate(weight=ifelse(weight==0,nth(weight,2,order_by = weight),weight))
    tmpverts<-unique(c(tmptbl$from,tmptbl$to)) 
    rell[[i]] <- graph_from_data_frame(tmptbl, directed=F, vertices=tmpverts)
  } else {
    rell[[i]] <- NA
  }
}

dfl_all <- ldply(dfl,data.frame)
singles <- dfl_all %>% dplyr::filter(nmem==1)
# groups <- list()
# for (i in 1:length(gd)) {
#   groups[[i]] <-data.frame(group=i,members= V(gd[[i]])$name)
# }
# groups_tbl <- ldply(groups,data.frame)
# 
# write.table(groups_tbl,'/vast/eande106/projects/Nicolas/github/concordance-nf/graph_dev/processed_data/c_briggsae_groups_pass1.tsv', quote = FALSE)
#save.image("/vast/eande106/projects/Nicolas/github/concordance-nf/graph_dev/processed_data/c_briggsae_env.Rda")


#search for substructure across isotypes using cluster_edge_betweeness()
glist_fg<- list()
glist_eb<- list()
cc <- list()
##iterate through graphs 
for (i in 1:length(rell)) {
  #skip single-strain groups
  if (any(is.na(rell[[i]]))) { 
  #if (length(rell[[i]])<2) { #length appears to work for both graph objects and NA
    next
  }
  
  #temporary
  sg <- rell[[i]]
  
  #temporarily skip large graphs for development testing
  # if(length(V(sg))>80) {
  #   next
  # }
  
  #plot(sg,layout=layout_with_fr,axes = T,vertex.size=2)
  sc_fg <- cluster_fast_greedy(sg)
  sc_eb <- cluster_edge_betweenness(sg,directed = F)
  V(sg)$label.cex=1.5
  #cc stores information about graphs that I thought could be useful
  # ngraph = initial graph index before NULL removal
  # mems_x = number of members according to community search algorithm x (eb or fg)
  # nmem = all members in the graph
  # ccrange = the min and max weights
  # drange = the min and max number of discordant sites among all pairs of trains
  # min/maxddiff = the values that make drange^
  cc[[i]] <- data.frame(ngraph=i,
                        mems_eb=length(unique(sc_eb$membership)),
                        mems_fg=length(unique(sc_fg$membership)),
                        nmem=length(sc_eb$membership),
                        ccrange=paste0(round(min(edge.attributes(sg)$concordance),digits = 5),"-",round(max(edge.attributes(sg)$concordance),digits = 5)),
                        drange=paste0(min(edge.attributes(sg)$discordant),"-",max(edge.attributes(sg)$discordant)),
                        minddiff=min(edge.attributes(sg)$discordant),
                        maxddiff=max(edge.attributes(sg)$discordant))
  glist_fg[[i]] <- as.grob(expression(plot(sc_fg,sg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  glist_eb[[i]] <- as.grob(expression(plot(sc_eb,sg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
}


memberships <- ldply(cc,data.frame) %>%
  dplyr::filter(!is.na(ngraph))

#subset to isotypes with detected substructure using edge_betweeness and have more than 2 vertices
subcomm <- memberships %>% 
  dplyr::filter(mems_eb > 1 & nmem> 2) %>%
  dplyr::mutate(size=ifelse(nmem >10,"Large","Small")) #arbitrary description of size


#store isotypes with no substructure
simples <- memberships %>% 
  dplyr::filter(!(ngraph %in% subcomm$ngraph)) 


#nrow(singles) + nrow(simples) + nrow(subcomm)
#yes, good

#write info about isotypes with substructure
#readr::write_delim(subcomm,"processed_data/elegans/isotypes_wSubcomm_localnorm.tsv",delim = "\t",quote = 'none')

#generate list of plots that contain representations of graphs and subgraphs of isotypes with detected subcommunities across genetic and environmental parameters as edge weights
graphCC <- list() #concordance graphs
graphDI <- list() #distance graphs
graphTI <- list() #time graphs
graphLB <- list() #lab graphs
subgrList2<- list() #plots
k=1 #plot index
for (j in subcomm$ngraph) {

  #store concordance graph
  sg <- rell[[j]]
  
  #temporary table to generate distance graph
  names <- vertex.attributes(sg)$name
  count <- length(names)
  tmptbl <- rels %>% 
    dplyr::filter(!is.na(dist)) %>%
    dplyr::filter((from %in% names) & (to %in% names)) %>%
    dplyr::mutate(weight=dist) %>%
    dplyr::mutate(weight=(abs(weight-max(weight)))+1)
  tmpverts<-unique(c(tmptbl$from,tmptbl$to)) 
  #store distance graph
  dg <- graph_from_data_frame(tmptbl, directed=F, vertices=tmpverts)
  
  #temporary table to generate time graph
  tmptbl <- rels %>% 
    dplyr::filter(!is.na(timedif)) %>%
    dplyr::filter((from %in% names) & (to %in% names)) %>%
    dplyr::mutate(weight=timedif) %>%
    dplyr::mutate(weight=(abs(weight-max(weight)))+1)
  tmpverts<-unique(c(tmptbl$from,tmptbl$to)) 
  #store time graph
  tg <- graph_from_data_frame(tmptbl, directed=F, vertices=tmpverts)
  
  #labdiff is boolean, replace weights of zero
  lg <- sg
  E(lg)$weight <- E(lg)$labdiff + 0.00000001 
  #E(tg)$weight <- E(tg)$dist

  #cluster subcommunities - fg seems to polarize intermediate strains to most genetically similar group
  dc_fg <- cluster_fast_greedy(dg)
  tc_fg <- cluster_fast_greedy(tg)
  sc_fg <- cluster_fast_greedy(sg)
  lc_eb <- cluster_edge_betweenness(lg,directed = F) #graph is probably not needed for lab, but easier comparison than a matrix
  
  #eb graphs not used
  #dc_eb <- cluster_edge_betweenness(dg,directed = F)
  #tc_eb <- cluster_edge_betweenness(tg,directed = F)
  #sc_eb <- cluster_edge_betweenness(sg,directed = F)

  #store graphs for each parameter
  graphCC[[k]] <- list(sg,sc_fg)
  graphDI[[k]] <- list(dg,dc_fg)
  graphTI[[k]] <- list(tg,tc_fg)
  graphLB[[k]] <- list(lg,lc_eb)
  
  #store graph plots
  p1 <- as.grob(expression(plot(dc_fg,dg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  p3 <- as.grob(expression(plot(sc_fg,sg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  p5 <- as.grob(expression(plot(tc_fg,tg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  p7 <- as.grob(expression(plot(lc_eb,lg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  
  #eb graphs not used
  #p2 <- as.grob(expression(plot(dc_eb,dg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  #p4 <- as.grob(expression(plot(sc_eb,sg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  #p6 <- as.grob(expression(plot(tc_eb,tg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  #p8 <- as.grob(expression(plot(lc_eb,lg,layout=layout_with_fr,mark.border="black",mark.col="lightpink",vertex.size=3.5,edge.width=0.2)))
  
  #store composite plots
  plot <- ggarrange(ggarrange(p3,hjust=-1,vjust=3),
                    ggarrange(p1,hjust=-1,vjust=3),
                    ggarrange(p5,hjust=-1,vjust=3),
                    ggarrange(p7,hjust=-1,vjust=3),
                    ncol=2,nrow=2,labels=c("By concordance:","By distance:","By time:","By lab:"),hjust=-0.2)
   
  
  subgrList2[[k]] <-annotate_figure(plot, top = text_grob(paste0("Group ID:",j), color = "red", face = "bold", size = 14))
  k=k+1
}

#perform all pairwise comparisons between members of subgraphs across environmental parameters
all_comps <- list()
for (k in 1:length(graphCC)) {
  
  #partition subcommunities of concordance graph
  target <- induceNS(graphCC[[k]])
  #partition subcommunities of distance graph
  dist <- induceNS(graphDI[[k]])

  set_comp <- list()
  for (i in 1:length(target)) {
    comm_comp <- list()
    for (j in 1:length(dist)) {
      ntar <- V(target[[i]])$name
      nparam <- V(dist[[j]])$name
      comp<- as.data.frame.list(c(k,i,j,length(ntar),length(nparam),length(ntar[(ntar %in% nparam)]),length(nparam[!(nparam %in% ntar)])))
      colnames(comp) <- c("cgraph","cgraph_subID","pgraph_subID","cgraph_NE","pgraph_NE","cgraph_MA","pgraph_MI")
      comm_comp[[j]] <- comp
      #print(comps[[k]])
    }
    set_comp[[i]] <- ldply(comm_comp,data.frame)
  }
  all_comps[[k]] <- ldply(set_comp,data.frame) 
  time <- induceNS(graphTI[[k]])
  lab <- induceNS(graphLB[[k]])
}

#estimate shared membership between subgraphs across environmental parameters and classify as SEP, KEEP, or UNK
#NE = number of elements
#MA = matches
#MI = missing
#pp = proportion
comp_subgr <- ldply(all_comps,data.frame) %>%
  dplyr::group_by(cgraph) %>%
  dplyr::mutate(max_subC=max(cgraph_subID)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!(cgraph_MA==0)) %>%
  dplyr::mutate(pp_MA=cgraph_MA/cgraph_NE,pp_MI=pgraph_MI/pgraph_NE) %>%
  dplyr::group_by(cgraph,cgraph_subID) %>%
  dplyr::mutate(groupID=cur_group_id()) %>%
  dplyr::filter(pp_MA==max(pp_MA) & pp_MI==(min(pp_MI))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(perf_MA=ifelse(pp_MA == 1 & pp_MI ==0,T,F)) %>%
  dplyr::group_by(cgraph) %>%
  dplyr::mutate(sep_iso=ifelse(any(perf_MA==T),"CORR","MAYBE")) %>%
  dplyr::mutate(partial=ifelse(pp_MA < 1 & pp_MI >0,F,T)) %>%
  dplyr::mutate(keep_iso=ifelse(any(partial)==T,"MAYBE","UNCORR")) %>%
  dplyr::mutate(outcome=ifelse(sep_iso=="CORR","SEP",ifelse(keep_iso=="UNCORR","KEEP","UNK"))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(targetUID=paste0(cgraph,"_",cgraph_subID),queryUID=paste0(cgraph,"_",pgraph_subID))

#isotypes that must be separated (complete shared membership)
matches <- comp_subgr %>%dplyr::filter(outcome=="SEP")
#isotypes that must be kept together (no shared membership)
keeps <- comp_subgr %>% dplyr::filter(outcome== "KEEP")
#unknowns
unk <- comp_subgr %>% dplyr::filter(outcome == "UNK") 
#classify unknowns by abitrary likelihood of being separated
class_unk <- unk %>%
  dplyr::mutate(maybe=ifelse(((pp_MI==0 & pp_MA >= 1/cgraph_NE & pp_MA > 0.5) | (pp_MI <= 1/pgraph_NE & pp_MA==1)),T,F)) %>%
  dplyr::group_by(cgraph) %>%
  dplyr::mutate(maybe_sep=ifelse(any(maybe)==T,"MSEP","UNK")) %>%
  dplyr::ungroup()
#isotypes that COULD be separated (partial matches with high shared membership)
maybs <- class_unk %>% dplyr::filter(maybe_sep=="MSEP")
#isotypes that are unlikely to be separated (partial matches with low shared membership)
unk2 <- class_unk %>% dplyr::filter(maybe_sep=="UNK")

#split `matches` into individual isotypes
splitList <- list()
for (i in unique(matches$cgraph)) {
  split <- induceNS(graphCC[[i]])
  
  subgList <- list()
  for (j in 1:length(graphCC[[i]])) {
    mem <- paste0(V(split[[j]])$name,collapse = ",")
    subgList[[j]] <- data.frame(members=mem,
                                nmem=length(V(split[[j]])$name),
                                gnum=paste0(subcomm$ngraph[i],"_",letters[j]))
  }
  splitList[[i]] <- ldply(subgList,data.frame)
}
split_df <- ldply(splitList,data.frame)

#get isotypes to keep
keep_subcomm <- subcomm[keeps$cgraph,] %>% dplyr::select(-size)
#merge simples and keepers
simp_keep <- rbind(simples,keep_subcomm)

#make df of simple+keep cases
simpleList<-list()
for (i in 1:length(simp_keep$ngraph)) {
  simpleList[[i]] <- data.frame(members=paste(V(rell[[simp_keep$ngraph[i]]])$name,collapse=","),
                                nmem=simp_keep$nmem[i],
                                gnum=simp_keep$ngraph[i])
}
simple_df <- ldply(simpleList,data.frame)



#make final isotype assignment (excluding unresolved graphs for manual curation)
resolved_groups <- rbind(singles,simple_df,split_df) %>%
  tidyr::separate(gnum,into = c("og_group","partition"),sep="_",remove = F) %>%
  dplyr::mutate(og_group=as.numeric(og_group)) %>%
  dplyr::rename(ID=gnum) %>%
  dplyr::arrange(og_group) %>%
  dplyr::select(-partition,-nmem) %>%
  dplyr::mutate(members = strsplit(as.character(members), ",")) %>% 
  tidyr::unnest(members)

#original CE isotype count 604
#594 are singles, simples, or keeps
#7 are split into 14
#17 are unknown and require manual curation 

#get unknowns
unk_subcomm <- subcomm[unique(unk$cgraph),] %>% dplyr::select(-size)
unk_graphs <- subgrList2[unique(unk$cgraph)]
