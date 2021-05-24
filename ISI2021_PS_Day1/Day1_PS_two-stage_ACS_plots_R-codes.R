# ----------------------------------------------------------------------------------------------------------
# Population graph and BIG representation of Thompsons data (1991,ACS: designs with primary and secondary units,
# Biometrics, p.1115)
# ----------------------------------------------------------------------------------------------------------
# Load igraph
library(igraph)

# --------------------------
# Plot 1: grids with cases
# --------------------------
skthPSUACS <- function(){
x.lim <- 100
y.lim <- 105
grid.dim <- 5

plot(0,xaxt="n",yaxt="n",type="l",ylab="",xlab="",xlim=c(0,x.lim),ylim=c(0,y.lim),bty="n")
lines(c(0,0),c(0,x.lim)); lines(c(0,x.lim),c(0,0)); lines(c(0,x.lim),c(x.lim,x.lim)); lines(c(x.lim,x.lim),c(0,x.lim))

for(k in 1:20){
  tmp <- k*grid.dim
  lines(c(0,x.lim),c(tmp,tmp))
  lines(c(tmp,tmp),c(0,x.lim))
}

# Strip labels
for(k in 0:19){
  tmp <- grid.dim/2 + k*grid.dim
  tmp.label <- LETTERS[(k+1)]
  text(tmp,103,label=tmp.label,cex=1.2)
}


# First case-network
points(x=runif(2,0.5,4.5),y=runif(2,60.5,64.5),pch=20,cex=0.8)
points(x=runif(2,0.5,4.5),y=runif(2,55.5,59.5),pch=20,cex=0.8)
points(x=runif(1,5.5,9.5),y=runif(1,65.5,69.5),pch=20,cex=0.8)
points(x=runif(11,5.5,9.5),y=runif(11,60.5,64.5),pch=20,cex=0.8)
points(x=runif(22,5.5,9.5),y=runif(22,55.5,59.5),pch=20,cex=0.8)
points(x=runif(3,5.5,9.5),y=runif(3,50.5,54.5),pch=20,cex=0.8)
points(x=runif(1,10.5,14.5),y=runif(1,65.5,69.5),pch=20,cex=0.8)
points(x=runif(26,10.5,14.5),y=runif(26,60.5,64.5),pch=20,cex=0.8)
points(x=runif(19,10.5,14.5),y=runif(19,55.5,59.5),pch=20,cex=0.8)
points(x=runif(5,10.5,14.5),y=runif(5,50.5,54.5),pch=20,cex=0.8)
points(x=runif(1,15.5,19.5),y=runif(1,65.5,69.5),pch=20,cex=0.8)
points(x=runif(5,15.5,19.5),y=runif(5,60.5,64.5),pch=20,cex=0.8)
points(x=runif(8,15.5,19.5),y=runif(8,55.5,59.5),pch=20,cex=0.8)

# Second case-network
points(x=runif(1,10.5,14.5),y=runif(1,15.5,19.5),pch=20,cex=0.8)
points(x=runif(5,15.5,19.5),y=runif(5,10.5,14.5),pch=20,cex=0.8)
points(x=runif(10,15.5,19.5),y=runif(10,15.5,19.5),pch=20,cex=0.8)
points(x=runif(17,15.5,19.5),y=runif(17,20.5,24.5),pch=20,cex=0.8)
points(x=runif(2,20.5,24.5),y=runif(2,10.5,14.5),pch=20,cex=0.8)
points(x=runif(26,20.5,24.5),y=runif(26,15.5,19.5),pch=20,cex=0.8)
points(x=runif(26,20.5,24.5),y=runif(26,20.5,24.5),pch=20,cex=0.8)
points(x=runif(1,20.5,24.5),y=runif(1,25.5,29.5),pch=20,cex=0.8)
points(x=runif(1,25.5,29.5),y=runif(1,10.5,14.5),pch=20,cex=0.8)
points(x=runif(6,25.5,29.5),y=runif(6,15.5,19.5),pch=20,cex=0.8)
points(x=runif(9,25.5,29.5),y=runif(9,20.5,24.5),pch=20,cex=0.8)
points(x=runif(1,25.5,29.5),y=runif(1,25.5,29.5),pch=20,cex=0.8)

# Third case-network
points(x=runif(2,65.5,69.5),y=runif(2,75.5,79.5),pch=20,cex=0.8)
points(x=runif(2,70.5,74.5),y=runif(2,75.5,79.5),pch=20,cex=0.8)
points(x=runif(22,70.5,74.5),y=runif(22,80.5,84.5),pch=20,cex=0.8)
points(x=runif(14,70.5,74.5),y=runif(14,85.5,89.5),pch=20,cex=0.8)
points(x=runif(6,74.5,79.5),y=runif(6,75.5,79.5),pch=20,cex=0.8)
points(x=runif(38,74.5,79.5),y=runif(38,80.5,84.5),pch=20,cex=0.8)
points(x=runif(25,74.5,79.5),y=runif(25,85.5,89.5),pch=20,cex=0.8)
points(x=runif(1,80.5,84.5),y=runif(1,75.5,79.5),pch=20,cex=0.8)
points(x=runif(3,80.5,84.5),y=runif(3,80.5,84.5),pch=20,cex=0.8)
points(x=runif(2,80.5,84.5),y=runif(2,85.5,89.5),pch=20,cex=0.8)
}

# ----------------------------
# Plot 2: BIG representation
# -----------------------------
skthPSUACS_BIG <- function(){
a <- make_empty_graph(n=400)
vertex_attr(a) <- list(name = c(1:400))
vertex_attr(a, "label") <- V(a)$name


b <- make_empty_graph()
b <- add.vertices(b,nv=20,attr=list(name=LETTERS[1:20],type=rep(TRUE,20)))
b <- add.vertices(b,nv=400,attr=list(name=V(a)$name,type=rep(FALSE,length(V(a)))))


b <- add_edges(b,edges=as.character(c(8,9,8,28,9,29,28,29,28,27,29,30,28,48,29,49,27,47,47,48,48,49,49,50,47,67,48,68,49,69,67,68,68,69)),
               attr=list(type=rep(FALSE,17)))
b <- add_edges(b,edges=as.character(c(57,77,76,77,77,78,78,98,77,97,76,96,95,96,96,97,97,98,98,118,97,117,96,116,115,116,116,117,117,118)),
               attr=list(type=rep(FALSE,15)))
b <- add_edges(b,edges=as.character(c(265,285,285,284,284,283,283,303,284,304,285,305,303,304,304,305,303,323,304,324,305,325,323,
                                      324,324,325)),attr=list(type=rep(FALSE,13)))



k <- 0
temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- 0
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),
                         temp.letter,27,temp.letter,28,temp.letter,29,temp.letter,30,temp.letter,47,temp.letter,48,temp.letter,49,temp.letter,50,temp.letter,67,temp.letter,68,temp.letter,69),
               attr=list(name=paste(temp.letter,c(1:31),sep=''),type=rep(TRUE,31)))

temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),
                         temp.letter,8,temp.letter,9,temp.letter,47,temp.letter,48,temp.letter,49,temp.letter,50,temp.letter,67,temp.letter,68,temp.letter,69),
               attr=list(name=paste(temp.letter,c(1:29),sep=''),type=rep(TRUE,29)))

temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),
                         temp.letter,8,temp.letter,9,temp.letter,27,temp.letter,28,temp.letter,29,temp.letter,30,temp.letter,67,temp.letter,68,temp.letter,69,temp.letter,76,temp.letter,77,temp.letter,78,
                         temp.letter,95,temp.letter,96,temp.letter,97,temp.letter,98,temp.letter,115,temp.letter,116,temp.letter,117,temp.letter,118),
               attr=list(name=paste(temp.letter,c(1:40),sep=''),type=rep(TRUE,40)))

temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),
                         temp.letter,8,temp.letter,9,temp.letter,27,temp.letter,28,temp.letter,29,temp.letter,30,
                         temp.letter,47,temp.letter,48,temp.letter,49,temp.letter,50,temp.letter,57,
                         temp.letter,95,temp.letter,96,temp.letter,97,temp.letter,98,
                         temp.letter,115,temp.letter,116,temp.letter,117,temp.letter,118),
               attr=list(name=paste(temp.letter,c(1:39),sep=''),type=rep(TRUE,39)))

temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),
                         temp.letter,57,temp.letter,76,temp.letter,77,temp.letter,78,
                         temp.letter,115,temp.letter,116,temp.letter,117,temp.letter,118),
               attr=list(name=paste(temp.letter,c(1:28),sep=''),type=rep(TRUE,28)))

temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),
                         temp.letter,57,temp.letter,76,temp.letter,77,temp.letter,78,
                         temp.letter,95,temp.letter,96,temp.letter,97,temp.letter,98),
               attr=list(name=paste(temp.letter,c(1:28),sep=''),type=rep(TRUE,28)))


for(j in 1:7){
temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                          temp.letter,(start.node+6),temp.letter,(start.node+7),
                          temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                          temp.letter,(start.node+13),temp.letter,(start.node+14),
                          temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                          temp.letter,(start.node+20)),
                attr=list(name=paste(temp.letter,c(1:20),sep=''),type=rep(TRUE,20)))
}


temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),temp.letter,283,temp.letter,284,temp.letter,285,temp.letter,303,
                         temp.letter,304,temp.letter,305,temp.letter,323,temp.letter,324,temp.letter,325),
               attr=list(name=paste(temp.letter,c(1:29),sep=''),type=rep(TRUE,29)))

temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),temp.letter,265,temp.letter,303,
                         temp.letter,304,temp.letter,305,temp.letter,323,temp.letter,324,temp.letter,325),
               attr=list(name=paste(temp.letter,c(1:27),sep=''),type=rep(TRUE,27)))


temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),temp.letter,265,temp.letter,283,temp.letter,284,
                         temp.letter,285,temp.letter,323,temp.letter,324,temp.letter,325),
               attr=list(name=paste(temp.letter,c(1:27),sep=''),type=rep(TRUE,27)))



temp.letter <- LETTERS[k+1]
k <- k + 1
start.node <- start.node + 20
b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                         temp.letter,(start.node+6),temp.letter,(start.node+7),
                         temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                         temp.letter,(start.node+13),temp.letter,(start.node+14),
                         temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                         temp.letter,(start.node+20),temp.letter,265,temp.letter,283,temp.letter,284,
                         temp.letter,285,temp.letter,303,temp.letter,304,temp.letter,305),
               attr=list(name=paste(temp.letter,c(1:27),sep=''),type=rep(TRUE,27)))


for(j in 1:3){
  temp.letter <- LETTERS[k+1]
  k <- k + 1
  start.node <- start.node + 20
  b <- add.edges(b,edges=c(temp.letter,(start.node+1),temp.letter,(start.node+2),temp.letter,(start.node+3),temp.letter,(start.node+4),temp.letter,(start.node+5),
                           temp.letter,(start.node+6),temp.letter,(start.node+7),
                           temp.letter,(start.node+8),temp.letter,(start.node+9),temp.letter,(start.node+10),temp.letter,(start.node+11),temp.letter,(start.node+12),
                           temp.letter,(start.node+13),temp.letter,(start.node+14),
                           temp.letter,(start.node+15),temp.letter,(start.node+16),temp.letter,(start.node+17),temp.letter,(start.node+18),temp.letter,(start.node+19),
                           temp.letter,(start.node+20)),
                 attr=list(name=paste(temp.letter,c(1:20),sep=''),type=rep(TRUE,20)))
}


b <- as.undirected(b,mode='each')


g <- b
set.seed(25052021)
plot(g,vertex.label=c(NA,V(g)$name[V(g)$type==TRUE])[c(c(1:20)+1,(V(b)$type[c(21:420)]+1))],
     vertex.size=c(2,10)[V(g)$type+1],
     vertex.color=c('orange','yellow')[V(g)$type+1],
     edge.color=c('green','grey')[E(g)$type+1],vertex.label.color='black',
     vertex.label.cex=0.8)
}

skthPSUACS()
skthPSUACS_BIG()
