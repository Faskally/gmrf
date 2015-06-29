
#' @export
plotgraph <- function(graph) {
  plot(graph,       #the graph to be plotted
  #layout=layout.fruchterman.reingold, # the layout method. see the igraph documentation for details
  #layout=layout.kamada.kawai,
  #vertex.frame.color='blue',    #the color of the border of the dots 
  #vertex.label.dist=0.5,
  #vertex.label.color='black',   #the color of the name labels
  #vertex.label.cex=.7,      #specifies the size of the font of the labels. can also be made to vary
  #vertex.size = 2,
  #edge.label = E(graph)$name,
  #edge.arrow.size = 0.5

  vertex.label=NA,
  vertex.frame.color = NA,
  edge.color = NA,
  vertex.size = 5
)

}

#' Simulate from an inhomogenous GMRF
#'
#' Details
#'
#' This function does stuff.
#'
#' @param Q a symmetric positive semi-definate matrix corresponding
#'     to an inhomogenous GMRF
#' @return a single draw from with the appropriate covariance structure
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
simplify.graph <- function(g) {
  # choose a vertex with degree 2
  i <- which(degree(g) == 2)
  while (length(i)) {
    i <- i[1]
    Vs <- as.numeric(V(g)[nei(i)])
    g[Vs[1], Vs[2]] <- sum(g[i,Vs])
    g <- delete.vertices(g, i)
    i <- which(degree(g) == 2)
  }
  g
}



#' Simulate from an inhomogenous GMRF
#'
#' Details
#'
#' This function does stuff.
#'
#' @param Q a symmetric positive semi-definate matrix corresponding
#'     to an inhomogenous GMRF
#' @return a single draw from with the appropriate covariance structure
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' wk <- ctm[ctm $ CATCHMENT %in% c(3:13, 15:50),]
#' # work out neighbours
#' wk_nb <- poly2nb(wk, queen = TRUE)
#' # create a rw2 GMRF precision matrix and simulate some spatial
#' # structure
#' Q <- getQMat(wk_nb)
#' x <- simQ(exp(-4)*Q)
#' plot(wk, col = grey( pmax(0, pmin(1, (x + 5) / 10))))
simQ <- function(Q, tol = 1e-9, rank = NULL) {
  eQ <- eigen(Q)
  if (is.null(rank)) rank <- sum(eQ $ values > tol)
  colSums(t(eQ $ vectors[,1:rank]) * rnorm(rank, 0, 1/sqrt(eQ $ values[1:rank])))
}



#' Simulate from an inhomogenous GMRF
#'
#' Details
#'
#' This function does stuff.
#'
#' @param eQ an eigen decomposition of a symmetric positive semi-definate matrix corresponding
#'     to an inhomogenous GMRF
#' @return a single draw from with the appropriate covariance structure
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' wk <- ctm[ctm $ CATCHMENT %in% c(3:13, 15:50),]
#' # work out neighbours
#' wk_nb <- poly2nb(wk, queen = TRUE)
#' # create a rw2 GMRF precision matrix and simulate some spatial
#' # structure
#' Q <- getQMat(wk_nb)
#' x <- simQ(exp(-4)*Q)
#' plot(wk, col = grey( pmax(0, pmin(1, (x + 5) / 10))))
simeQ <- function(eQ, k = 1, tol = 1e-9, rank = NULL) {
    if (is.null(rank)) rank <- sum(eQ $ values > tol)
  colSums(t(eQ $ vectors[,1:rank]) * rnorm(rank, 0, 1/sqrt(k * eQ $ values[1:rank])))
}

#' Title
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' dim(ctm)
getQMat <- function(g) {

  # is this a connected graph?
  if (!is.connected(g)) {
    stop("there are some unconnected vertices...")
  }

  Q <- g[]
  if (any(Q @ x == 0)) stop("something went wrong")
  Q @ x[] <- -1/Q @ x
  diag(Q) <- degree(g)
  Q
}


#' Compute RW1 matrix
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
# compute RW1 matrix
getRW1Mat <- function(g) {
  # is this a connected graph?
  if (!is.connected(g)) {
    stop("there are some unconnected vertices...")
  }

  Q <- g[]
  if (any(Q @ x == 0)) stop("something went wrong")
  Q @ x[] <- -1
  diag(Q) <- apply(g[], 1, function(x) sum(x>0))
  Q
}




#' Compute Weighted RW1 matrix
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
getWRW1Mat <- function(g) {
  # is this a connected graph?
  if (!is.connected(g)) {
    stop("there are some unconnected vertices...")
  }

  Q <- g[]
  if (any(Q @ x == 0)) stop("something went wrong")
  Q @ x[] <- -1/Q @ x
  diag(Q) <- apply(g[], 1, function(x) sum(1/x[x>0]))
  Q
}


#' Title
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' plotcatch(4)
plotcatch <- function(area, add = FALSE) {
    y <- rivs[rivs $ CATCH_ == area,]

    HC <- unique(floor(y $ HYDRO_CODE))
    HC <- HC[HC > 0]

    plot(y, col = grey(0.8), add = add)
    plot(y[y $ HYDRO_CODE == 0,], col = "lightgreen", add = TRUE)
    if (sum(y $ HYDRO_CODE > 0)) {
        plot(y[floor(y $ HYDRO_CODE) %in% HC,], add = TRUE, col = "red")
        if (sum(y $ HYDRO_CODE %in% HC))  
          plot(y[y $ HYDRO_CODE %in% HC,], add = TRUE, col = "blue")
    }
    title(main = y $ CATCH_NAME[1])
    plot(ctm[ctm $ CATCHMENT != area, ], col = grey(0.9), border = grey(0.8), add = TRUE)
    if (sum(ctm $ CATCHMENT == area)) plot(ctm[ctm $ CATCHMENT == area,], add = TRUE)
    points(gis[gis $ CATCHMENT == area,], pch = 16, cex = 0.8)
}



#' Title
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' plotcatchtile(90)
plotcatchtile <- function(area) {
    # get the osm tile
    z <- tiles[[paste(area)]]

    # get the river network and convert projection
    y <- rivs[rivs $ CATCH_ == area,]
    y <- spTransform(y, z $ tiles[[1]] $ projection)

    # get the sites and convert projection
    x <- gis[gis $ CATCHMENT == area,]
    x <- spTransform(x, z $ tiles[[1]] $ projection)

    # convert local catchmentboundaries
    #w <-     
    plot(z)
    plot(y, add = TRUE, col = "darkblue", lwd = 2)
    plot(x, add = TRUE, pch = 16, col = "red")
    #plot(ctm, border = grey(0.3), add = TRUE)
    title(main = y $ CATCH_NAME[1])
} 



#' Title
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' dim(ctm)
plotrivers <- function(area, osm = FALSE) {
  y <- rivs[rivs $ CATCH_ == area,]

  bbox <- as.data.frame(t(y @ bbox))

  wgs84 = '+proj=longlat +datum=WGS84'
  bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
  merc = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"

  y @ proj4string <- CRS(bng)

  ## assign original coordinate system
  ## Then to transform it to WGS84
  bboxll <- spTransform(SpatialPoints(cbind(bbox$x, bbox$y), CRS(bng)), 
                        CRS(wgs84)) @ coords

  # get map
  topright <- c(lat = bboxll[2,2], lon = bboxll[1,1])
  botleft <- c(lat = bboxll[1,2], lon = bboxll[2,1])

  if (osm) {
    gm <- openmap(topright,botleft, type = "bing")
  } else {
    gm <- get_map(c(t(bboxll)), maptype = "terrain")
  }

  if (osm) {
    # convert y to tiles projection
    y <- spTransform(y, CRS(merc))
  } else {
    # convert y to tiles projection
    y <- spTransform(y, CRS(wgs84))
  }

  # Then extract line data
  xys <- do.call(rbind, lapply(1:length(y), function(i) y @ lines[[i]] @ Lines[[1]] @ coords))
  colnames(xys) <- c("x", "y")
  xys <- as.data.frame(xys)

  np <- sapply(1:length(y), function(i) nrow(y @ lines[[i]] @ Lines[[1]] @ coords))
  id <- rep(1:length(y), np)
  xys <- cbind(xys, y @ data[id,])
  xys $ seg <- id
  xys $ grp <- with(xys, ifelse(HYDRO_CODE == 0, "No Code", 
                           ifelse(HYDRO_CODE == 10, "Main", 
                              ifelse(floor(HYDRO_CODE) == HYDRO_CODE, "Other Baseline",
                                  "Other river"))))

  mytheme <- 
    theme(axis.line        = element_blank(),
          axis.text.x      = element_blank(),
          axis.text.y      = element_blank(),
          axis.ticks       = element_blank(),
          axis.title.x     = element_blank(),
          axis.title.y     = element_blank(),
          legend.position  = "bottom",
          panel.background = element_blank(),
          panel.border     = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background  = element_blank())


  gl <- geom_line(aes(x = x, y = y, 
                      group = OBJECTID,
                      colour = factor(grp)
                      ), 
                  size = 0.5,
                  data = xys)
  if (osm) {
    gg <- autoplot(gm) + mytheme + gl
  } else {
    gg <- ggmap(gm, darken = 0.2) + mytheme + gl
  }

  list(gg = gg, om = gm)
}



#' Title
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' y <- rivs[rivs $ CATCH_NAME == "River Irvine",]
#' g <- buildTopo(y)
#' plot(g, vertex.label=NA, vertex.size=2,vertex.size2=2)
buildTopo = function(lines) {

    g = gIntersection(lines, lines)
    edges = do.call(rbind, lapply(g@lines[[1]]@Lines, function(ls) {
        as.vector(t(ls@coords))
    }))
    lengths = sqrt((edges[, 1] - edges[, 3])^2 + (edges[, 
        2] - edges[, 4])^2)

    froms = paste(edges[, 1], edges[, 2])
    tos = paste(edges[, 3], edges[, 4])

    graph = graph.edgelist(cbind(froms, tos), directed = FALSE)
    E(graph)$weight = lengths

    xy = do.call(rbind, strsplit(V(graph)$name, " "))

    V(graph)$x = as.numeric(xy[, 1])
    V(graph)$y = as.numeric(xy[, 2])

    # better names?
    V(graph)$name <- paste(1:length(V(graph)$name))

    return(graph)
}



#' Title
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' y <- rivs[rivs $ CATCH_NAME == "River Irvine",]
#' g <- buildTopo(y)
#' 
#' from = cbind(291867.1, 933646.4)
#' to = cbind(312009.1, 922430.1)
#' pp = routePoints(g, from, to)
#' plot(y)
#' plotRoute(g, pp, col = "blue", lwd = 2, asp = 1,  xlab = "", ylab = "", add = TRUE)
#' points(rbind(from, to))
plotRoute <- function(graph, nodes, add = FALSE, ...) {
    if (add) {
        lines(V(graph)[nodes]$x, V(graph)[nodes]$y, ...)
    } else {
        plot(V(graph)[nodes]$x, V(graph)[nodes]$y, type = "l", 
            ...)
    }
}



#' Title
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' y <- rivs[rivs $ CATCH_NAME == "River Irvine",]
#' g <- buildTopo(y)
#' 
#' from = cbind(291867.1, 933646.4)
#' to = cbind(312009.1, 922430.1)
#' pp = routePoints(g, from, to)
routePoints <- function(graph, from, to) {
    xyg = cbind(V(graph)$x, V(graph)$y)

    ifrom = get.knnx(xyg, from, 1)$nn.index[1, 1]
    ito = get.knnx(xyg, to, 1)$nn.index[1, 1]

    p = get.shortest.paths(graph, ifrom, ito, output = "vpath")
    p $ vpath [[1]]
}



#' Finds and plots the shortest route between two points
#' on a river network
#'
#' Details
#'
#' Description - This function does stuff.
#'
#' @param nb an input parameter
#' @return what does it return
#' @seealso \code{\link{getQMat}} for creating a GMRF from a graph
#' @export
#' @examples
#' y <- rivs[rivs $ CATCH_NAME == "River Irvine",]
#' g <- buildTopo(y)

#' 
#' from = cbind(291867.1, 933646.4)
#' to = cbind(312009.1, 922430.1)
#' pp = routePoints(g, from, to)
#' plot(y)
#' plotRoute(g, pp, col = "blue", lwd = 2, asp = 1,  xlab = "", ylab = "", add = TRUE)
#' points(rbind(from, to))
#' 
#' # go auto
#' plotMyRoute(g, y)
plotMyRoute <- function(graph, river) {
  # one of mine!
  plotcatch(river @ data $ CATCH_[1])
  from = do.call(cbind, locator(1))
  to = do.call(cbind, locator(1))
  points(rbind(from, to))
  pp = routePoints(graph, from, to)
  plotRoute(graph, pp, col = "blue", lwd = 2, asp = 1,  xlab = "", ylab = "", add = TRUE)
  #plot(river, add = TRUE)    
}





# some useful functions
riverNext <- function(ID, by.distance = FALSE) {
  L1 <- rivs[rivs $ OBJECTID == ID,]
  if (L1 $ TNODE_ == 0) {
    cat("unconnected node\n")
    return(NA)
  }
  if (sum(rivs $ FNODE_ == L1 $ TNODE_)==0) {
    cat("river mouth\n")
    return(NA)
  }
  L2s <- rivs[rivs $ FNODE_ == L1 $ TNODE_,]

  flocs <- sapply(L2s @ lines, function(ll) ll @ Lines[[1]] @ coords[1,])
  tloc <- tail(L1 @ lines[[1]] @ Lines[[1]] @ coords, 1)
  dist <- sqrt(colSums((flocs - c(tloc))^2))
  if (all(dist > 10))  {
    cat("river mouth\n")
    return(NA)
  }
  L2ind <- which.min(dist) 
  L2 <- L2s[L2ind,]

  rbind(L1, L2)
}

riverEnd <- function(ID, add = FALSE, old = NULL, plot = FALSE) {
  route <- rep(NA, 100)
  route[1] <- ID
  rnext <- riverNext(ID)
  i <- 1
  while(!identical(rnext, NA) & !(route[i] %in% old)) {
    route[i <- i + 1] <- rnext $ OBJECTID[2]
    plot(rivs[rivs $ OBJECTID %in% route,], col = "blue", lwd = 2)
    rnext <- riverNext(route[i])
  }
  route[!is.na(route)]
}





