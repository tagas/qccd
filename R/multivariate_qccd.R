rm(list = ls())
library(pcalg)
library(readr)
library(CAM)
library(ggraph)
library(igraph)
library(naglr)
source("./cd_utils.R")

# 1. load dataset
octet_ds <-
  read_delim("./../data/octet.tab",
             "\t",
             escape_double = FALSE,
             trim_ws = TRUE)

octet_ds = data.frame(octet_ds)

# 2. get  get CPDAG from CAM
cam_adj <- CAM(octet_ds)
cam_adj = cam_adj$Adj
g_cam <- graph_from_adjacency_matrix(cam_adj)
cam_cpdag = dag2cpdag(igraph.to.graphNEL(g_cam))
cam_adj = as(cam_cpdag, "matrix")


d = ncol(octet_ds)
dir_mat = matrix(0, nrow = d, ncol = d)
conf_mat = matrix(0, nrow = d, ncol = d)
dir_mat_dag = matrix(0, nrow = d, ncol = d)
dir_mat_fin = matrix(0, nrow = d, ncol = d)


# 3. orent edges in CPDAG with QCCD
for (i in 1:d) {
  for (j in 1:d) {
    if (i < j && (cam_adj[i, j] == 1 || cam_adj[j, i] == 1)) {
      res = QCCD(cbind(octet_ds[, i], octet_ds[, j]))
      cd = res$cd
      conf = res$eps
      if (cd == T) {
        dir_mat[i, j] = 1
        dir_mat[j, i] = 0
        conf_mat[i, j] = conf
      } else{
        dir_mat[i, j] = 0
        dir_mat[j, i] = 1
        conf_mat[j, i] = conf
      }
    }
  }
}


# 4. Rank edges by confidence and insert into final Adj matrix
conf_tmp = conf_mat
conf_mat[conf_mat < 0.5] = 1 - conf_mat[conf_mat < 0.5]
conf_mat[conf_mat == 1] = 0

for (i in 1:sum(cam_adj)) {
  c_max = which(conf_mat == max(conf_mat), arr.ind = TRUE)
  conf_mat[c_max[1], c_max[2]] = 0
  dir_mat_dag[c_max[1], c_max[2]] = 1
  
  m2 <- matrix(dir_mat_dag, d, d)
  if (isValidGraph(m2, type = "dag")) {
    print(paste("Added: "))
    dir_mat_fin[c_max[1], c_max[2]] = 1
  } else{
    dir_mat_dag[c_max[1], c_max[2]] = 0
    print(paste("Skipped: "))
  }
}


g_cam <- graph_from_adjacency_matrix(cam_adj)
g_cam_qqcd <- graph_from_adjacency_matrix(dir_mat_fin)

# 5. is a Valid DAG
m1 <- matrix(cam_adj, d, d)
isValidGraph(m1, type = "dag")
# FALSE because m1 is CPDAG

m2 <- matrix(dir_mat_fin, d, d)
isValidGraph(m2, type = "dag")
# TRUE

pdf("./../results/octet_cam_qccd.pdf", width = 10, height = 6)


par(mfrow = c(1, 2))
set.seed(0)
names <- c(
  "rpA",
  "rsA",
  "r_sigma",
  "|rpA-rpB|",
  "(ipA-eaB) \n /rpA",
  "|rsA-rsB|",
  "(eaB-ipB) \n /rpA_sq",
  "(ipA-eaA)\n /rpA",
  "(ipA-eaB) \n /rsA",
  "hA-lA",
  "energy"
)

# 6. plotting results 
par(mfrow = c(1,2))
l1 <- layout_with_sugiyama(g_cam)
l2 <- layout_with_sugiyama(g_cam_qqcd)

plot(
  g_cam,
  vertex.size =  25,
  layout = l1$layout,
  edge.arrow.size = 1 / 3,
  vertex.color = "white",
  vertex.label = names,
  label.dist = 5,
  label.cex  = 5,
  main = "CAM"
)


plot(
  g_cam_qqcd,
  vertex.size =  25,
  layout = l1$layout,
  edge.arrow.size = 1 / 3,
  vertex.color = "white",
  vertex.label = names,
  vertex.color = "blue",
  label.dist = 5,
  label.cex  = 5,
  main = "CAM QQCD"
)
dev.off()


library(ggraph)

set.seed(0)
gg_cam <- g_cam %>%
  set_vertex_attr("label", value = names)
vertex_attr(gg_cam, "label")
vertex_attr(gg_cam)

names <- 1:11
V(gg_cam)$name <- names
# plot using ggraph
ggraph(gg_cam, layout = "igraph", algorithm = 'kk') +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                 end_cap = circle(6, 'mm'), colour = "black") +
   geom_point(aes(x=x, y=y), colour = "darkgray", size = 15) +
  geom_node_text(aes(label = name), size = 4) + theme_naglr(grid = FALSE, ticks = FALSE, axis = FALSE,plot_margin = margin(30, 30, 30, 30) ) + labs(x = " ", y = " ")



set.seed(0)
g <- g_cam_qqcd %>%
  set_vertex_attr("label", value = names)
vertex_attr(g, "label")
vertex_attr(g)
names <- 1:11

V(g)$name <- names
# plot using ggraph
ggraph(g, layout = 'igraph', algorithm = 'graphopt') + 
  geom_edge_link(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                 end_cap = circle(6, 'mm'), colour = "black") +
  geom_point(aes(x=x, y=y), colour = "darkgray", size = 15) +
  geom_node_text(aes(label = name), size = 4) + theme_naglr(grid = FALSE, ticks = FALSE, axis = FALSE,plot_margin = margin(30, 30, 30, 30) ) + labs(x = " ", y = " ")



