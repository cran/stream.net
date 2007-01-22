############################################################
# stream.net.R
#   version 1.0.5
#
# tree network functions, d white, 22 January 2007
#
#   net.arcinput         read Arc/Info AAT and gen format lines
#   net.qmodel           create a random network
#   net.orders           compute Strahler and Shreve orders
#   net.lengths          compute segment lengths
#   net.shortlengths     compute segment straight line lengths
#   net.addsegs          add segments to links
#   net.addatt           add attribute to segments
#   net.dist             compute segment up/down distance matrix
#   net.total.dist       compute segment total distance matrix
#   net.dir              compute segment direction matrix
#   net.prox             compute neighborhood of a segment
#   net.maxupslope       compute segment max upstream slope
#   net.correlogram      compute network correlogram
#   net.autocorr.onelag  compute autocorrelation on network
#   net.autocorr.att     add autocorrelated attribute to segments
#   net.interp           interpolate sparse attribute to network
#   net.group            classify attribute for net.map
#   net.map              draw network geography with labels
#   net.map.key          draw key for net.map
#   net.read             read table format dump of network
#   net.write            write out table format of network
#   read.arcgenlin       read Arc/Info generate format lines
#
# The network data structure, as returned from net.qmodel 
# and net.arcinput, is a list with three components:
#
# $links    a data frame of the topological structure with
#           the following components
#   $lid      identifier for the link
#   $parent   index of parent link (or 0 if root)
#   $left     index of left child (or 0 if terminal)
#   $right    index of right child (or 0 if terminal)
#   $depth    topological depth of link (# links from root)
#   $first    index of first segment of this link
#   $last     index of last segment of this link
#   $strahler Strahler order
#   $shreve   Shreve order
#   ...       additional attributes
#
# $segs     a data frame of the segment structure (pieces of
#           links) with the following components
#   $sid      identifier for the segment (for external ref)
#   $link     index to parent link
#   $nxt      index of the next segment for this link
#   $prev     index of the previous segment for this link
#   $up       if segment (and coordinates) oriented upstream
#             = 1, else = 0
#   $length   segment coordinate length
#   ...       additional attributes
#
# $cords    a list of the x and y coordinates for each 
#           segment
#   $x        x coordinates with a sublist for each segment
#   $y        y coordinates with a sublist for each segment
# (The names for the sublists are the respective sids)
#
#
############################################################

net.arcinput <- function (aatname, linname, lineformat="R")
#
# Read in a network from Arc/Info data.  Arguments:
#
# aatname     name of AAT object; if "character", then
#             character string of the AAT file name,
#             else R object with AAT data
# linname     name of lines object; if "character", then
#             character string of the lines file name,
#             else R object with lines data
# lineformat  "arcgen" for Arc/Info generate format lines;
#             "R" for S/R format lines; only applicable
#             to file input
#
# Arc/Info input consists of two objects, either from 
# files in text format, or existing R objects.
#
# First, information from the coverage AAT, unloaded from 
# the tables command, and containing at least the user-id 
# for each arc and the from# and to# node numbers.  Any 
# number of arc attributes can follow.  After unloading, 
# the file should be edited so that the first record, the 
# "header", has the names of the fields in the file, in
# order, separated by commas.  The user-id is expected to 
# have the name "arcid", the from# node "from", and the to#
# node "to".  Otherwise the names can be any legal R names
# (see the R reference manual).
#
# Second, the arc coordinates from an ungenerate command.
# The ungenerate command should use the "fixed" option (and
# the "line" option, of course).  If the lines have already  
# been converted to S/R format, then argument "lineformat" 
# should be set to "R" (much faster read).
#
# It is essential that each arc have a unique user-id ("arcid"),
# that there be no nodes of degree greater than 2, that there 
# be no cycles in the network (loops), and that the topology be 
# clean for this process to work correctly.
#
{
    follow <- function (startnode, parent, depth)
    # sub-function that recursively follows links 
    # in the network
    {
        paths <- ((aat$from == startnode) | 
                (aat$to == startnode)) &
                ! aat$done
        while (sum (paths) == 1)
        {
            start <- which (paths)[1]
            aat$done[start] <<- TRUE
            if (aat$from[start] == startnode) up <- 1
            else up <- 0
            seg <<- seg + 1
            links$last[link] <<- seg
            segs[seg, 1:5] <<- 
                c(aat$arcid[start], link, seg, seg-1, up)
            if (maxaat > 3) 
                segs[seg, 6:(maxaat+2)] <<- aat[start, 4:maxaat]
            segs[seg-1, "nxt"] <<- seg
            if (up == 1) startnode <- aat$to[start]
            else startnode <- aat$from[start]
            paths <- ((aat$from == startnode) | 
                    (aat$to == startnode)) &
                    ! aat$done
        }
        if (sum (paths) > 2)
        {
            string <- paste ("node", startnode, "degree > 2")
            stop (string)
        }
        else if (sum (paths) == 2)
        {
            start <- which (paths)[1]
            aat$done[start] <<- TRUE
            if (aat$from[start] == startnode) up <- 1
            else up <- 0
            plink <- link <<- link + 1
            seg <<- seg + 1
            links[link,] <<- 
                c(link, parent, 0, 0, depth+1, seg, seg)
            links$left[parent] <<- link
            segs[seg, 1:5] <<- 
                c(aat$arcid[start], link, seg, seg, up)
            if (maxaat > 3) 
                segs[seg, 6:(maxaat+2)] <<- aat[start, 4:maxaat]
            if (up == 1) follow (aat$to[start], plink, depth+1)
            else follow (aat$from[start], plink, depth+1)

            start <- which (paths)[2]
            aat$done[start] <<- TRUE
            if (aat$from[start] == startnode) up <- 1
            else up <- 0
            plink <- link <<- link + 1
            seg <<- seg + 1
            links[link,] <<- 
                c(link, parent, 0, 0, depth+1, seg, seg)
            links$right[parent] <<- link
            segs[seg, 1:5] <<- 
                c(aat$arcid[start], link, seg, seg, up)
            if (maxaat > 3) 
                segs[seg, 6:(maxaat+2)] <<- aat[start, 4:maxaat]
            if (up == 1) follow (aat$to[start], plink, depth+1)
            else follow (aat$from[start], plink, depth+1)
        }
        invisible ()
    } # end of function follow

    if (mode (aatname) == "character") aat <- 
        read.table (aatname, header=TRUE, sep=",", as.is=TRUE)
    else aat <- aatname
    links <- data.frame (matrix (c(1,0,0,0,1,1,1),
        nrow=1, ncol=7, dimnames=list(1, 
        c("lid","parent","left","right","depth","first","last"))))
    maxaat <- ncol (aat)
    froms <- unique (sort (aat$from))
    tos <- unique (sort (aat$to))
    startnode <- tos[which (! tos %in% froms)]
    if (length (startnode) == 1) 
    {
        start <- aat$to == startnode
        up <- 0
    }
    else
    {
        startnode <- froms[which (! froms %in% tos)]
        if (length (startnode) != 1) stop (
            "arcinput: unique mouth not found")
        start <- aat$from == startnode
        up <- 1
    }
    segs <- data.frame (matrix (c(0,1,1,1,0), nrow=1, dimnames=
        list (1, c("sid","link", "nxt","prev","up"))))
    if (maxaat > 3) segs <- cbind (segs, aat[start, 4:maxaat])
    segs$sid[1] <- aat$arcid[start]
    segs$up[1] <- up
    aat$done <- rep (FALSE, nrow(aat))
    aat$done[start] <- TRUE
    link <- 1
    seg <- 1
    if (up == 1) follow (aat$to[start], 1, 1)
    else follow (aat$from[start], 1, 1)

    if (mode (linname) == "character") 
    {
        if (lineformat == "arcgen") arcs <- 
            read.arcgenlin (linname)
        else arcs <- read.table (linname, header=TRUE)
    }
    else arcs <- linname
    firsts <- is.na (arcs$y)
    pos <- which (c(firsts, TRUE))
    arcids <- arcs$x[firsts]
    arcs$x[firsts] <- NA
    cords <- list (x=list(), y=list())
    for (i in 1:length(arcids))
    {
        cords$x[[i]] <- arcs$x[(pos[i]+1):pos[i+1]]
        cords$y[[i]] <- arcs$y[(pos[i]+1):pos[i+1]]
    }
    names (cords$x) <- names (cords$y) <- arcids
    net <- list (links=net.orders(links), segs=segs, cords=cords)
    net$segs$length <- net.lengths (net)
    return (net)
}

############################################################

net.qmodel <- function (size=100, Q=0.5)
#
# Generate a stochastic network according to the Q model
# (as described in Costa-Cabral & Burges, 1997, Water
# Resources Research 33:2179-2197).  Arguments:
#
# size  magnitude of the network in number of terminal branches
# Q     probability of an internal branch rather than a 
#       terminal branch (0 <= Q <= 1)
#
# Any value of Q will generate all possible topologies of the 
# given size, however, the likelihood of any topology is 
# controlled by the value provided.  If Q = 0.01, then networks 
# with short lengths to terminal branches are more likely 
# (because branching from existing terminal branches is much 
# more probable than from existing internal branches).  If
# Q = 0.99, then networks with a few long sequences of branches 
# are more likely.
#
{
    incdepth <- function (link)
    # increment depths recursively
    {
        links$depth[link] <<- links$depth[link] + 1
        if (links$left[link] > 0) 
        {
            incdepth (links$left[link])
            incdepth (links$right[link])
        }
    }

    getx <- function (start)
    # x coordinates number end nodes left to right
    {
        left <- links$left[start]
        right <- links$right[start]
        if (left == 0)
        {
            x <<- x + 1
            cords$x[[start]][2] <<- x
        }
        else 
        {
            getx (left)
            getx (right)
            x1 <- cords$x[[left]][2]
            x2 <- cords$x[[right]][2]
            hx <- 0.5 * (x1 + x2)
            cords$x[[start]][2] <<- hx
            cords$x[[left]][1] <<- hx
            cords$x[[right]][1] <<- hx
        }
    }

    add <- size - 2
    links <- data.frame (matrix (c(1,2,3,0,1,1,2,0,0,3,0,0,1,2,2),
        nrow=3, ncol=5, dimnames=list(seq(3), 
        c("lid","parent","left","right","depth"))))
    if (add > 0) 
    {
        for (i in 1:add)
        {
            lids <- links$lid
            ends <- ! as.logical (links$left)
            s <- sum (ends)
            n <- nrow (links)
            r <- runif (1)
            if (r > Q) 
            {
                k <- ceiling ((r - Q) / (1 - Q) * s)
                j <- seq(n)[ends][k]
                links$left[j] <- n + 1
                links$right[j] <- n + 2
                links[n+1,] <- c(n+1,j,0,0,links$depth[j]+1)
                links[n+2,] <- c(n+2,j,0,0,links$depth[j]+1)
            }
            else
            {
                k <- ceiling (r / Q * (n - s))
                j <- seq(n)[!ends][k]
                lr <- round (runif (1))
                if (lr == 0) 
                {
                    incdepth (links$left[j])
                    incdepth (links$right[j])
                    links[n+1,] <- c(n+1,j,0,0,links$depth[j]+1)
                    links[n+2,] <- c(n+2,j,links$left[j],
                        links$right[j],links$depth[j]+1)
                    links$parent[links$left[j]] <- n + 2
                    links$parent[links$right[j]] <- n + 2
                    links$left[j] <- n + 1
                    links$right[j] <- n + 2
                }
                else
                {
                    incdepth (links$left[j])
                    incdepth (links$right[j])
                    links[n+1,] <- c(n+1,j,links$left[j],
                        links$right[j],links$depth[j]+1)
                    links[n+2,] <- c(n+2,j,0,0,links$depth[j]+1)
                    links$parent[links$left[j]] <- n + 1
                    links$parent[links$right[j]] <- n + 1
                    links$left[j] <- n + 1
                    links$right[j] <- n + 2
                }
            }
        }
    }
    links$last <- links$first <- links$lid
    segs <- data.frame (matrix (0, 
        nrow=nrow(links), ncol=5, dimnames=
        list (seq(nrow(links)), c("sid","link","nxt","prev","up"))))
    for (i in 1:nrow(segs)) segs[i,] <- c(i, i, i, i, 1)

    cords <- list (x=rep(list(c(0, 0, NA)), nrow(links)), 
                   y=rep(list(c(0, 0, NA)), nrow(links)))
    x <- 0
    getx (1)
    cords$x[[1]][1] <- cords$x[[1]][2]
    for (i in 1:nrow(links))
        cords$y[[links$lid[i]]] <- c(links$depth[i]-1, links$depth[i], NA)
    names (cords$x) <- names (cords$y) <- seq (nrow (links))

    net <- list (links=net.orders (links), segs=segs, cords=cords)
    net$segs$length <- net.lengths (net)
    return (net)
}

############################################################

net.orders <- function (links)
#
# Calculate the Strahler and Shreve orders for the links
# of a network.  Argument is the link table for a network.
# Returns new link table.
#
{
    ends <- ! as.logical (links$left)
    left <- links$left
    right <- links$right
    st <- sh <- rep (0, nrow (links))
    st[ends] <- 1
    sh[ends] <- 1
    while (any (st == 0))
    {
        for (i in 1:nrow(links))
        {
            a <- st[left[i]]
            b <- st[right[i]]
            if (st[i] == 0 && a > 0 && b > 0)
            {
                st[i] <- max (((a + b) %/% 2) + 1, max (a, b))
                sh[i] <- sh[left[i]] + sh[right[i]]
            }
        }
    }
    links$strahler <- as.integer (st)
    links$shreve <- as.integer (sh)
    return (links)
}

############################################################

net.lengths <- function (net)
#
# Compute segment lengths based on coordinate distance.  
# Argument is a network object.  Returns lengths.
#
{
    getlength <- function (sid)
    {
        x <- net$cords$x[[as.character(sid)]]
        y <- net$cords$y[[as.character(sid)]]
        z <- 0
        for (i in 2:(length(x)-1))
        {
            z <- z + sqrt ((x[i]-x[i-1])^2 + (y[i]-y[i-1])^2)
        }
        z
    }
    depthfirst <- function (seg)
    { # depth first traversal of net segments
        lengths[seg] <<- getlength (net$segs$sid[seg])
        if (net$segs$up[seg] == 1) path <- "nxt" 
        else path <- "prev"
        if (net$segs[seg, path] == seg) 
            { # go to links
            link <- net$segs$link[seg]
            left <- net$links$left[link]
            if (left == 0) return ()
            if (net$segs$up[seg] == 1) seg <- net$links$first[left]
            else seg <- net$links$last[left]
            depthfirst (seg)
            right <- net$links$right[link]
            if (net$segs$up[seg] == 1) seg <- net$links$first[right]
            else seg <- net$links$last[right]
            depthfirst (seg)
        }
        else depthfirst (net$segs[seg, path])
    } # end of recursive function depthfirst
    lengths <- rep (0, nrow (net$segs))
    # depthfirst (1)
    for (seg in 1:nrow(net$segs)) lengths[seg] <- 
        getlength (net$segs$sid[seg])
    return (lengths)
}

############################################################

net.shortlengths <- function (net)
#
# Calculate straight line distance for segments from first to 
# last nodes ignoring intermediate nodes.  Return vector of 
# these lengths.
#
{
    getlength <- function (sid)
    {
        x <- net$cords$x[[as.character(sid)]]
        y <- net$cords$y[[as.character(sid)]]
        n <- length(x) - 1
        sqrt ((x[n]-x[1])^2 + (y[n]-y[1])^2)
    }
    lengths <- rep (0, nrow (net$segs))
    for (seg in 1:nrow(net$segs)) lengths[seg] <- 
        getlength (net$segs$sid[seg])
    return (lengths)
}

############################################################

net.addsegs <- function (net, dist="Uniform", min=1, max=2, 
    shape=1.35, scale=2.3)
#
# Add segments to (random) network.  Number of segments to 
# add per link is random draw from uniform distribution  
# between min and max, inclusive.  Does not handle segment 
# attributes, so should be used before attributes are added.
#
# 3/2/06 scott leibowitz
# Revised version allows for Gamma distribution as well as
# Uniform (default).  The "min" and "max" arguments only
# apply to "Uniform", and "shape and "scale" only apply to
# "Gamma".  This version also allows for adding zero segments
# by breaking out of the loop after the assignment of n if n=0.
#
{
    if (dist == "Gamma") 
    {
        ndist <- parse (text = 
        "floor (rgamma (1, shape=shape, scale=scale))")
    }
    else # Uniform
    {
        ndist <- parse (text = 
        "floor (runif (1, min=min, max=(max + 1 - 1e-10)))")
    }
    for (i in 1:nrow (net$links))
    {
        send <- nrow (net$segs)
        first <- net$links$first[i]
        last <- prev <- net$links$last[i]
        n <- eval(ndist)
        if (n == 0) next
        x0 <- net$cords$x[[first]][1]
        dx <- (net$cords$x[[first]][2] - x0) / (n+1)
        net$cords$x[[first]][2] <- x0 + dx
        y0 <- net$cords$y[[first]][1]
        dy <- (net$cords$y[[first]][2] - y0) / (n+1)
        net$cords$y[[first]][2] <- y0 + dy
        for (j in 1:n)
        {
            net$segs[send+j, 1:5] <- c(send+j,i,send+j,prev,1)
            net$segs[last, "nxt"] <- send + j
            last <- prev <- net$links$last[i] <- send + j
            net$cords$x[[send+j]] <- c(x0+j*dx, x0+(j+1)*dx, NA)
            net$cords$y[[send+j]] <- c(y0+j*dy, y0+(j+1)*dy, NA)
        }
    }
    names (net$cords$x) <- names (net$cords$y) <- net$segs$sid
    net$segs$length <- net.lengths (net)
    return (net)
}

############################################################

net.addatt <- function (net, name=NULL, func=NULL, dist=NULL,
    ind="shreve", sd=0.02, arg1=NULL, arg2=NULL, 
    boundscaling=FALSE, outscaling=FALSE, prob=0.01, 
    min=0, max=1, vector=TRUE)
#
# Add a segment attribute generated by either a functional
# relationship with optional noise added, or a distributional 
# property.  See code below for valid names of functions and 
# distributions, and arguments to the distribution functions.
# Arguments:
#
# net           a network object
# name          name of the attribute to be added (if ! vector)
# func          name of the functional form
# dist          name of the distribution function
# ind           the independent variable for functional 
#                 values
# sd            standard deviation for added noise to 
#                 func output.  Used as the proportion of 
#                 the range of the produced values.
# arg1          first argument to distribution functions
#                 e.g., shape1, mean, etc.  See below.
# arg2          second argument to distribution functions
# boundscaling  true to scale distribution values to either 
#                 absolute bounds of distribution or to the 
#                 prob and/or 1-prob quantile; implies
#                 outscaling is true
# outscaling    true to scale return values to [min, max]
# prob          probability value for scaling distributions
# min           minimum value for scaling output values
# max           maximum value for scaling output values
# vector        true to return a vector else a field in $segs
#
# Returns network with attribute added to $segs else vector
{
    inscale <- function (x, ymin, ymax)
    {
        yxscale <- (ymax - ymin) / diff (range (x))
        y <- (x - min (x)) * yxscale + ymin
        y
    }
    outscale <- function (x, xbounds, ybounds)
    {
        yxscale <- diff (ybounds) / diff (xbounds)
        y <- (x - xbounds[1]) * yxscale + ybounds[1]
        y
    }
    perturb <- function (x, sd=0.02) 
    {
        r <- sd * diff (range (x))
        e <- rnorm (length (x), sd=r)
        y <- x + e
        y
    }
    if (is.null (func) && is.null (dist)) dist <- "Uniform"
    n <- nrow (net$segs)
    if (! is.null (func))
    {
        if (ind == "strahler" || ind == "shreve")
        {
            x <- rep (0, n)
            for (i in 1:nrow (net$links))
            {
                nxt <- net$links$first[i]
                last <- net$links$last[i]
                xi <- net$links[i, ind]
                x[nxt] <- xi
                while (nxt != last)
                {
                    nxt <- net$segs[nxt, "nxt"]
                    x[nxt] <- xi
                }
            }
        }
        else x <- net$segs[ , ind]
        y <- switch (func,
            Linear = perturb (x, sd=sd),
            NegLinear = perturb (-x, sd=sd),
            Exponential = perturb (exp (inscale (x, 1, 10)), sd=sd),
            NegExponential = perturb (exp (-inscale (x, 1, 10)), sd=sd),
            Sigmoid = perturb (1/(1+exp(-inscale (x, -5, 5))), sd=sd),
            RevSigmoid = perturb (rev (1/(1+exp(-inscale (x, -5, 5)))), sd=sd),
            Logarithmic = perturb (log (inscale (x, 1, 100)), sd=sd),
            RevLogarithmic = perturb (log (rev (inscale (x, 1, 100))), sd=sd),
            Sinusoidal = perturb (cos (inscale (x, -pi, pi)), sd=sd),
            Quadratic = perturb (-(inscale (x, -10, 10))^2, sd=sd))
        bounds <- range (y)
    }
    else
    {
        args <- switch (dist,
            Uniform = c(0, 1),
            Normal = c(0, 1),
            Ramp = c(2, 1),
            NegRamp = c(1, 2),
            Exponential = c(100, 1),
            NegExponential = c(1, 100),
            Lognormal = c(0, 0.4),
            Weibull = c(100, 100),
            Unimodal = c(2, 2),
            Logistic = c(0, 1))
        if (is.null (arg1)) arg1 <- args[1]
        if (is.null (arg2)) arg2 <- args[2]
        y <- switch (dist,
            Uniform = runif (n, min=arg1, max=arg2),
            Normal = rnorm (n, mean=arg1, sd=arg2),
            Ramp = rbeta (n, shape1=arg1, shape2=arg2),
            NegRamp = rbeta (n, shape1=arg1, shape2=arg2),
            Exponential = rbeta (n, shape1=arg1, shape2=arg2),
            NegExponential = rbeta (n, shape1=arg1, shape2=arg2),
            Lognormal = rlnorm (n, meanlog=arg1, sdlog=arg2),
            Weibull = rweibull (n, shape=arg1, scale=arg2),
            Unimodal = rbeta (n, shape1=arg1, shape2=arg2),
            Logistic = rlogis (n, location=arg1, scale=arg2))
        if (boundscaling) bounds <- switch (dist,
            Uniform = c(0, 1),
            Normal = c(qnorm (prob, mean=arg1, sd=arg2), 
                       qnorm (1-prob, mean=arg1, sd=arg2)),
            Ramp = c(0, 1),
            NegRamp = c(0, 1),
            Exponential = c(0, 1),
            NegExponential = c(0, 1),
            Lognormal = c(0, qlnorm (1-prob, meanlog=arg1, sdlog=arg2)),
            Weibull = c(0, qweibull (1-prob, shape=arg1, scale=arg2)),
            Unimodal = c(0, 1),
            Logistic = c(qlogis (prob, location=arg1, scale=arg2),
                         qlogis (1-prob, location=arg1, scale=arg2)))
        else bounds <- range (y)
    }
    if (outscaling || boundscaling) 
        y <- outscale (y, bounds, c(min, max))
    if (vector) 
    {
        names (y) <- net$segs$sid
        return (y)
    }
    else
    {
        present <- which (names (net$segs) == name)
        if (length (present) > 0) net$segs <- net$segs[,- present]
        net$segs$y <- y
        names (net$segs) <- c (names (net$segs)[-ncol(net$segs)], name)
        return (net)
    }
}

############################################################

net.dist <- function (net, ends=0.5, method="coordinate",
    digits=10)
#
# Creates network segment pairwise distance matrix
#
# net     a network object
# ends    how to handle from and to segments
#         if ends=0, do not use from and to lengths
#         if ends=1, use total of from and to lengths
#         if ends=0.5, use half of from and to lengths
# method  "coordinate" means network coordinate distance
#         "segment" means distance in # of segments
# digits  if NULL, do not round output matrix, 
#         else precision of rounding
#
# Returns distance matrix.
#
# Upstream distances are found by indexing the source segment 
# by its row and the destination segment by its column. 
# Downstream distances are the opposite; the source is the 
# column and the destination is the row.
#
# Algorithm adapted from that of SG Leibowitz.
# This algorithm does produce roundoff residue and thus the 
# matrix can be purged of very small non-zero values with 
# rounding.
{
    getsegs <- function ()
    {
        segs <- vector ("list", length=nrow(net$links))
        for (i in 1:nrow(net$segs))
        {
            link <- net$segs$link[i]
            segs[[link]] <- c(i, segs[[link]])
        }
        segs
    }

    allsegs <- getsegs ()
    if (method == "segment") lengths <- rep (1, nrow (net$segs))
    else lengths <- net$segs$length
    dist <- matrix (0, nrow(net$segs), nrow(net$segs))
    droot <- rep (0, nrow (net$segs))
    link <- seg <- 1
    if (method == "segment") droot[seg] <- 1
    else droot[seg] <- lengths[seg]
    leftsibs <- list ()
    parents <- NULL
    nodes <- NULL
    links <- NULL
    lastsegs <- NULL
    fini <- FALSE
    repeat
    {
        parent <- seg
        parents <- c(parent, parents)
        seg <- net$segs$nxt[parent]
        if (seg == parent)
        {
            left <- net$links$left[link]
            if (left > 0)
            {
                # branch to new left link
                links <- c(left, net$links$right[link], NA, links)
                nodes <- c(seg, nodes)
                link <- left
                seg <- net$links$first[link]
            }
            else
            {
                # pop queues
                once <- TRUE 
                repeat
                {
                    segs <- allsegs[[links[1]]]
                    links <- links[-1]
                    link <- links[1]
                    if (is.na (link)) 
                    {
                        if (length (links) > 1)
                        {
                            links <- links[-1] # pop NA
                            nodes <- nodes[-1]
                            parent <- parents[1]
                            if (once) 
                            {
                                once <- FALSE
                                leftsibs[[1]] <- c(segs, leftsibs[[1]])
                                lastsegs <- lastsegs[-1]
                            }
                            else
                            {
                                leftsibs[[2]] <- c(segs, leftsibs[[2]], 
                                    leftsibs[[1]])
                                leftsibs <- leftsibs[-1]
                                lastsegs <- lastsegs[-1]
                            }
                            i <- match (nodes[1], parents)
                            parents <- parents[i:length(parents)]
                            link <- links[1]
                        }
                        else
                        {
                            fini <- TRUE
                            break
                        }
                    }
                    else
                    {
                        parent <- parents[1]
                        if (once) 
                        {
                            l <- list (segs)
                            leftsibs <- c(l, leftsibs)
                        }
                        else leftsibs[[1]] <- c(segs, leftsibs[[1]])
                        lastsegs <- c(nodes[1], lastsegs)
                        i <- match (nodes[1], parents)
                        parents <- parents[i:length(parents)]
                        break
                    }
                }
                if (fini) break
                parent <- net$links$last[net$links$parent[link]]
                seg <- net$links$first[link]
            }
        }

        droot[seg] <- droot[parent] + lengths[seg]
        dist[parents, seg] <- droot[seg] - droot[parents] + 
            ends * lengths[parents] - (1 - ends) * lengths[seg]

        if (length (leftsibs) > 0)
        {
            for (i in 1:length(lastsegs))
            {
                l <- leftsibs[[i]]
                n <- lastsegs[i]
                dist[l, seg] <- droot[seg] - droot[n] -
                    (1 - ends) * lengths[seg]
                dist[seg, l] <- droot[l] - droot[n] -
                    (1 - ends) * lengths[l]
            }
        }
    }
    if (! is.null (digits)) dist <- round (dist, digits)
    rownames (dist) <- colnames (dist) <- net$segs$sid
    return (dist)
}

############################################################

net.total.dist <- function (dist)
#
# Create symmetric total distance matrix between pairs of segments,
# ignoring upstream and downstream
#
# dist  an upstream/downstream distance matrix from net.dist
#
# Returns the total distance matrix.
{
    total <- dist + t(dist)
    return (total)
}

############################################################

net.dir <- function (dist)
#
# Create direction matrix between pairs of segments
#
# dist  an upstream/downstream distance matrix
#
# direction is direction between all pairs of segments,
#   where 0 = sib net, -1 = down, and 1 = up
# Returns the direction matrix.
{
    dir <- matrix (0, nrow=nrow(dist), ncol=nrow(dist))
    dir[dist == 0] <- -1
    dir[t(dist == 0)] <- 1
    diag (dir) <- 0
    return (dir)
}

############################################################

net.prox <- function (dist, seg, lag=1, direction="both")
#
# Find neighborhood of a segment.  
#
# dist       a distance matrix from net.dist
# seg        index to a segment
# lag        distance within which neighborhood is computed
# direction  "up", "down", or "both"
#
# This returns a vector of segments in proximity of seg in 
# the direction "up", "down", or "both".  Siblings that
# are in a combination of upstream and downstream directions 
# are not included in "up" or "down".  The distance matrix 
# determines whether the neighborhood is defined by number 
# of segments or by coordinates. 
#
{
    x <- dist[seg, ]
    y <- dist[, seg]
    if (direction == "up") 
    {
        xs <- x != 0
        ys <- y == 0
        d <- x[xs & ys]
        e <- which (xs & ys)
        proxvec <- e[d <= lag] 
    }
    else if (direction == "down")
    {
        xs <- x == 0
        ys <- y != 0
        d <- y[xs & ys]
        e <- which (xs & ys)
        proxvec <- e[d <= lag] 
    }
    else
    {
        d <- net.total.dist (dist)
        d <- d[seg, ]
        proxvec <- which ((d != 0) & (d <= lag))
    }
    return (proxvec)
}

############################################################

net.maxupslope <- function (net, slopes)
#
# Creates network segment pairwise matrix to record maximum 
# upstream slope
#
# net     a network object
# slopes  vector of "slopes" in segment order 
#
# Uses the distance matrix algorithm to calculate the maximum 
# upstream slope between two segments.  The maximum is the 
# maximum of the slopes of all segments between the two 
# segments.  Returns matrix with these max upstream slopes.
{
    getsegs <- function ()
    {
        segs <- vector ("list", length=nrow(net$links))
        for (i in 1:nrow(net$segs))
        {
            link <- net$segs$link[i]
            segs[[link]] <- c(i, segs[[link]])
        }
        segs
    }

    allsegs <- getsegs ()
    maxupslope <- matrix (0, nrow(net$segs), nrow(net$segs))
    link <- seg <- 1
    leftsibs <- list ()
    parents <- NULL
    nodes <- NULL
    links <- NULL
    lastsegs <- NULL
    fini <- FALSE
    repeat
    {
        parent <- seg
        parents <- c(parent, parents)
        seg <- net$segs$nxt[parent]
        if (seg == parent)
        {
            left <- net$links$left[link]
            if (left > 0)
            {
                # branch to new left link
                links <- c(left, net$links$right[link], NA, links)
                nodes <- c(seg, nodes)
                link <- left
                seg <- net$links$first[link]
            }
            else
            {
                # pop queues
                once <- TRUE 
                repeat
                {
                    segs <- allsegs[[links[1]]]
                    links <- links[-1]
                    link <- links[1]
                    if (is.na (link)) 
                    {
                        if (length (links) > 1)
                        {
                            links <- links[-1] # pop NA
                            nodes <- nodes[-1]
                            parent <- parents[1]
                            if (once) 
                            {
                                once <- FALSE
                                leftsibs[[1]] <- c(segs, leftsibs[[1]])
                                lastsegs <- lastsegs[-1]
                            }
                            else
                            {
                                leftsibs[[2]] <- c(segs, leftsibs[[2]], 
                                    leftsibs[[1]])
                                leftsibs <- leftsibs[-1]
                                lastsegs <- lastsegs[-1]
                            }
                            i <- match (nodes[1], parents)
                            parents <- parents[i:length(parents)]
                            link <- links[1]
                        }
                        else
                        {
                            fini <- TRUE
                            break
                        }
                    }
                    else
                    {
                        parent <- parents[1]
                        if (once) 
                        {
                            l <- list (segs)
                            leftsibs <- c(l, leftsibs)
                        }
                        else leftsibs[[1]] <- c(segs, leftsibs[[1]])
                        lastsegs <- c(nodes[1], lastsegs)
                        i <- match (nodes[1], parents)
                        parents <- parents[i:length(parents)]
                        break
                    }
                }
                if (fini) break
                parent <- net$links$last[net$links$parent[link]]
                seg <- net$links$first[link]
            }
        }

        path <- parents
        while (length(path) > 1)
        {
            s <- path[length(path)]
            path <- path[- length(path)]
            maxupslope[s, seg] <- max (slopes[path])
        }
    }
    return (maxupslope)
}

############################################################

net.correlogram <- function (net, dist, segatt, nlags=10)
#
# Calculate autocorrelation function on segment attribute 
# using network distances
#
# net     a network object
# dist    a distance matrix from net.dist
#         (probably using ends=0.5)
# segatt  attribute associated with each segment.
#         If segatt is a character vector of length one,
#         then get attribute from net$segs$segatt, 
#         else assume segatt is a numeric vector of length  
#         equal to number of segments and in correct order.
# nlags   number of lags in the output function
#
# Uses total distances ignoring upstream/downstream.
# Returns list with autocorrelation coefficient, lag 
# distance, and number of pairs for each lag.
#
# Reference: Henley S. 1975. Autocorrelation coefficients 
# from irregularly spaced areal data. Computers and Geosciences 
# 2(4):437-438.
{
    nsegs <- nrow (net$segs)
    dist <- net.total.dist (dist)
    maxdist <- max (dist)
    if (is.character (segatt) && length (segatt) == 1) 
        z <- net$segs[[segatt]]
    else z <- segatt
    z <- z - mean (z)
    index <- t (combn (seq (nsegs), 2))
    scp <- ssq <- num <- lag <- rep (0, nlags)
    r <- rep (NA, nlags)
    for (k in 1:nrow(index))
    {
        i <- index[k, 1]
        j <- index[k, 2]
        d <- round (dist[i, j] / maxdist * nlags)
        scp[d] <- scp[d] + z[i] * z[j] * 2
        ssq[d] <- ssq[d] + (z[i] * z[i]) + (z[j] * z[j])
        num[d] <- num[d] + 1
    }
    for (k in 1:nlags) if (num[k] != 0)
    {
        r[k] <- scp[k] / ssq[k]
        lag[k] <- k * (maxdist / nlags)
    }
    return (list (r=r, lag=lag, num=num))
}

############################################################

net.autocorr.onelag <- function (net, dist, segatt, lag=1, 
    eps=1e-6)
#
# Calculate one lag autocorrelation for given lag on segment 
# attribute using network distances
#
# net     a network object
# dist    a distance matrix from net.dist,
#         probably using ends=0.5, and possibly using 
#         method="segment"
# segatt  attribute associated with each segment.
#         If segatt is a character vector of length one,
#         then get attribute from net$segs$segatt, 
#         else assume segatt is a numeric vector of length  
#         equal to number of segments and in correct order.
# lag     lag at which to produce autocorrelation
# eps     precision of calculating lag distance
#
# Uses total distances ignoring upstream/downstream.
# Returns list with autocorrelation coefficient and 
# number of pairs used.
{
    dist <- net.total.dist (dist)
    if (is.character (segatt) && length (segatt) == 1) 
        z <- net$segs[[segatt]]
    else z <- segatt
    nsegs <- nrow (net$segs)
    index <- t (combn (seq (nsegs), 2))
    keep <- NULL
    for (k in 1:nrow(index))
    {
        i <- index[k, 1]
        j <- index[k, 2]
        if (abs (dist[i, j] - lag) < eps) keep <- c(keep, k)
    }
    index <- index[keep, ]
    if (nrow (index) < 1) stop (
        "autocorr.onelag: no pairs at specified lag")
    z <- z - mean (z)
    scp <- ssq <- 0
    for (k in 1:nrow(index))
    {
        i <- index[k, 1]
        j <- index[k, 2]
        scp <- scp + z[i] * z[j] * 2
        ssq <- ssq + (z[i] * z[i]) + (z[j] * z[j])
    }
    return (list (r=scp / ssq, num=nrow(index)))
}

############################################################

net.autocorr.att <- function (net, dist, target=0.5, lag=1,
    outscaling=TRUE, min=0, max=1, eps=1e-6,
    vector=TRUE, name=NULL)
#
# Create a network autocorrelated attribute for segments
#
# net         a network object
# dist        a distance matrix from net.dist
#             (probably using ends=0.5)
# target      goal for autocorrelation
# lag         lag at which to produce autocorrelation
# outscaling  true to scale return values to [min, max]
# min         minimum value for scaling output values
# max         maximum value for scaling output values
# eps         precision of calculating lag distance
# vector      true to return a vector else a field in $segs
# name        name of the attribute to be added (if ! vector)
#
# Returns network with attribute added to $segs else vector
{
    outscale <- function (x, xbounds, ybounds)
    {
        yxscale <- diff (ybounds) / diff (xbounds)
        y <- (x - xbounds[1]) * yxscale + ybounds[1]
        y
    }
    onelag <- function (z)
    {
        z <- z - mean (z)
        scp <- ssq <- 0
        for (k in 1:nrow(index))
        {
            i <- index[k, 1]
            j <- index[k, 2]
            scp <- scp + z[i] * z[j] * 2
            ssq <- ssq + (z[i] * z[i]) + (z[j] * z[j])
        }
        return (scp / ssq)
    }
    tdist <- net.total.dist (dist)
    nsegs <- nrow (net$segs)
    index <- t (combn (seq (nsegs), 2))
    keep <- NULL
    for (k in 1:nrow(index))
    {
        i <- index[k, 1]
        j <- index[k, 2]
        if (abs (tdist[i, j] - lag) < eps) keep <- c(keep, k)
    }
    index <- index[keep, ]
    if (nrow (index) < 1) stop (
        "autocorr.att: no pairs at specified lag")
    att <- runif (nsegs)
    decay <- max (1, min (10 * lag, round (0.2 * nsegs)))
    while (onelag (att) < target)
    {
        start <- sample (nsegs, 1)
        done <- NULL
        for (i in 1:decay)
        {
            near <- net.prox (dist, seg=start, lag=i)
            near <- near[which (! near %in% done)]
            vals <- att[start] - 
                rnorm (length(near), i*0.05, sd=0.02)
            att[near] <- vals
            done <- c(done, near)
        }
    }
    if (outscaling) att <- outscale (att, range (att), c(min, max))
    if (vector)
    {
        names (att) <- net$segs$sid
        return (att)
    }
    else
    {
        present <- which (names (net$segs) == name)
        if (length (present) > 0) net$segs <- net$segs[,- present]
        net$segs$att <- att
        names (net$segs) <- c (names (net$segs)[-ncol(net$segs)], name)
        return (net)
    }
}

############################################################

net.interp <- function (net, dist, samples, predict=NULL, 
    maxdist=1e32, method="inverseDistance", power=2, 
    vector=TRUE, name=NULL)
#
# Interpolate sparse attributes to whole network
#
# net       a network object
# dist      a distance matrix from net.dist
#           (probably using ends=0.5)
# samples   two column matrix or data frame with first column 
#           sids to segments of sampled data and second
#           column sampled data values
# predict   vector of segments to which to predict, unless null 
#           in which case predict to all other segments
# maxdist   maximum distance for neighborhood
# method    "inverseDistance"
# power     exponent for inverse distance function
# vector    true to return a vector, else a field in $segs
# name      name of the attribute to be added (if ! vector)
#
# Returns network with interpolated attribute added to $segs 
# else vector
#
{
    nsegs <- nrow (net$segs)
    j <- match (samples[,1], net$segs$sid)
    if (any (is.na (j))) stop (
        "interp: unmatched sample sid")
    x <- rep (NA, nsegs)
    x[j] <- samples[,2]
    # samples[,1] <- j
    if (is.null (predict))
    {
        predict <- seq (nsegs) 
        predict <- predict[ ! (predict %in% j)]
    }
    tdist <- net.total.dist (dist)
    for (i in 1:length (predict))
    {
        sampdist <- tdist[predict[i], j]
        k <- which (sampdist <= maxdist)
        if (length (k) > 0)
        {
            v <- samples[k, 2]
            w <- sampdist[k] ^ power
            num <- sum (v / w)
            den <- sum (w)
            x[predict[i]] <- num / den
        }
    }
    if (vector) 
    {
        names (x) <- net$segs$sid
        return (x)
    }
    else
    {
        present <- which (names (net$segs) == name)
        if (length (present) > 0) net$segs <- net$segs[,- present]
        net$segs$x <- x
        names (net$segs) <- c (names (net$segs)[-ncol(net$segs)], name)
        return (net)
    }
}

############################################################

net.group <- function (net, segatt=NULL, linkatt=NULL, 
    ngroups=5, method="quantile", spread=NULL)
#
# Make attribute groups for mapping.  Currently, attribute 
# must be numerical.  Arguments:
#
# net      a network object
# segatt   if not NULL, attribute associated with each segment.
#          If segatt is a character vector of length one,
#          then get attribute from net$segs$segatt, 
#          else assume segatt is a numeric vector of length  
#          equal to number of segments and in correct order.
# segatt   if not NULL, attribute associated with each link.
#          If linkatt is a character vector of length one,
#          then get attribute from net$links$linkatt, 
#          else assume linkatt is a numeric vector of length  
#          equal to number of links and in correct order.
# ngroups  the number of groups to create
# method   either "quantile" or "equalInterval".  "quantile" 
#          means quartiles, quintiles, etc., depending on 
#          ngroups.  "equalInterval" means dividing the 
#          range of the attribute into ngroups equal 
#          intervals.
# spread   if not NULL, then for "equalInterval" classification, 
#          set the range of intervals to be the two element 
#          spread vector
#
# One of segatt or linkatt must be non-NULL.
#
# Returns a list with components $group, a valid grouping to use 
# with net.map, and $cuts, the break points that were generated 
# for the classification
{
    if (! is.null (segatt))
    {
        if (is.character (segatt) && length (segatt) == 1) 
            x <- net$segs[[segatt]]
        else x <- segatt
    }
    else if (! is.null (linkatt))
    {
        if (is.character (linkatt) && length (linkatt) == 1) 
            vec <- FALSE
        else vec <- TRUE
        x <- rep (0, nrow (net$segs))
        for (i in 1:nrow (net$links))
        {
            nxt <- net$links$first[i]
            last <- net$links$last[i]
            if (vec) xi <- linkatt[i]
            else xi <- net$links[i, linkatt]
            x[nxt] <- xi
            while (nxt != last)
            {
                nxt <- net$segs[nxt, "nxt"]
                x[nxt] <- xi
            }
        }
    }
    else stop ("group: no attribute for grouping")
    if (method == "quantile")
    {
        n <- ngroups
        while (n > 0)
        {
            cuts <- quantile (x, probs=seq(0,1,1/n), na.rm=TRUE)
            if (any (diff (cuts) == 0)) n <- n - 1
            else break
        }
        if (n < ngroups) warning ("quantile cardinality reduced")
    }
    else if (method == "equalInterval")
    {
        if (is.null (spread)) rang <- range (x, na.rm=TRUE)
        else rang <- spread
        cuts <- seq (from=rang[1], to=rang[2], length=(ngroups+1))
    }
    group <- cut (x, breaks=cuts, labels=FALSE, include.lowest=TRUE)
    names (group) <- net$segs$sid
    return (list (group=group, cuts=cuts))
}

############################################################

net.map <- function (net, group=NULL, linkatt=NULL, segatt=NULL, 
    col=NULL, lwd=NULL, cex=par("cex"), new=TRUE, outline=NULL,
    uniquegroup=FALSE)
#
# draw map of network.  Arguments:
#
# net      a network object
# group    if not NULL, then a vector of color codes for 
#          each segment.  names(group) must match sids.
#          (see below)
# linkatt  if not NULL, then write the attribute value for
#          the link attribute in text format
#          at the average of first and last points
# segatt   if not NULL, then write the attribute value for
#          the segment attribute in text format
#          at the midpoint of the segment
#          If segatt is a character vector of length one,
#          then get attribute from net$segs$segatt, 
#          else assume segatt is a numeric vector of length  
#          equal to number of segments and in correct order.
# col      par parameter for color codes for segments
#          default ramp is yellow-red-brown
# lwd      par parameter for line widths of segments
# cex      par parameter for size of text
# new      if TRUE, then create a new plot
# outline  if not NULL, then shade outline in light gray
#          first before drawing network; outline format is
#          S/R polygon format
# uniquegroup  if TRUE, group categories are not members of
#              the integers 1:n
#
# This function is patterned after map.groups in contributed
# library maptree.  See that help file for more details.
#
{
    if (new) 
    {
        if (! is.null (outline))
        {
            xor <- range (outline$x[!is.na(outline$y)], 
                na.rm=TRUE)
            yor <- range (outline$y[!is.na(outline$x)], 
                na.rm=TRUE)
        }
        xcr <- range (as.vector (
            sapply (net$cords$x, range, na.rm=TRUE)))
        ycr <- range (as.vector (
            sapply (net$cords$y, range, na.rm=TRUE)))
        if (! is.null (outline))
        {
            xrange <- ifelse (xor[1] < xcr[1], xor[1], xcr[1])
            xrange <- c(xrange, 
                ifelse (xor[2] > xcr[2], xor[2], xcr[2]))
            yrange <- ifelse (yor[1] < ycr[1], yor[1], ycr[1])
            yrange <- c(yrange, 
                ifelse (yor[2] > ycr[2], yor[2], ycr[2]))
        }
        else
        {
            xrange <- xcr
            yrange <- ycr
        }
        plot.new ()
        plot.window (xrange, yrange, asp=1)
    }
    if (! is.null (outline))
        polygon (outline, col="grey90", border="black")
    if (! is.null (group))
    {
        dense <- sort (unique (na.omit (group)))
        if (uniquegroup) nc <- length (dense) 
        else nc <- max(dense)
        if (is.null (col)) 
        {
            hue <- seq (0.167, 0.0, length=nc)
            sat <- rep (1, nc)
            bri <- c(1,1,1,1,1,0.9,0.8,0.7,0.6)[round 
                   (seq(1,9,length=nc))]
            kol <- hsv(hue, sat, bri)
        }
        else kol <- rep (col, nc)
        if (! is.null (lwd))
        {
            if (lwd[1] == "seq") lwd <- seq (1, nc)
            else if (length (lwd) < nc) lwd <- rep (lwd, nc)
        }
        else lwd <- par ("lwd")
    }
    else 
    {
        if (is.null (col)) kol <- "black" else kol <- col
        if (is.null (lwd)) lwd <- par ("lwd") else lwd <- lwd
    }
    # draw in ascending shreve order
    x <- rep (0, nrow (net$segs))
    for (i in 1:nrow (net$links))
    {
        nxt <- net$links$first[i]
        last <- net$links$last[i]
        xi <- net$links[i, "shreve"]
        x[nxt] <- xi
        while (nxt != last)
        {
            nxt <- net$segs[nxt, "nxt"]
            x[nxt] <- xi
        }
    }
    o <- order (x)
    sids <- net$segs$sid[o]
    for (i in 1:nrow (net$segs))
    {
        j <- match (sids[i], as.integer (names (net$cords$x)))
        xt <- net$cords$x[[j]]
        yt <- net$cords$y[[j]]
        if (! is.null (group))
        {
            j <- match (sids[i], as.integer (names (group)))
            if (uniquegroup) k <- match (group[j], dense)
            else k <- group[j]
            if (is.na (k)) lines (xt, yt, col="black", lwd=1)
            else lines (xt, yt, col=kol[k], lwd=lwd[k])
        }
        else lines (xt, yt, col=kol, lwd=lwd)
        n <- length (xt)
        if (n == 3)
        {
            xt <- 0.5 * (xt[1] + xt[2])
            yt <- 0.5 * (yt[1] + yt[2])
        }
        else
        {
            j <- ceiling (n %/% 2)
            xt <- xt[j]
            yt <- yt[j]
        }
        if (! is.null (segatt))
        {
            if (is.character (segatt) && length (segatt) == 1) 
                z <- net$segs[[segatt]][o[i]]
            else z <- segatt[o[i]]
            text (xt, yt, z, cex=cex, font=2)
        }
    }
    if (! is.null (linkatt)) for (i in 1:nrow (net$links))
    { # this could be improved by adding a midpoint attribute
      # to each link
        first <- net$segs$sid[net$links$first[i]]
        last <- net$segs$sid[net$links$last[i]]
        j <- match (first, as.integer (names (net$cords$x)))
        x1 <- net$cords$x[[j]][1]
        y1 <- net$cords$y[[j]][1]
        j <- match (last, as.integer (names (net$cords$x)))
        l <- length (net$cords$x[[j]])
        x2 <- net$cords$x[[j]][l-1]
        y2 <- net$cords$y[[j]][l-1]
        xt <- 0.5 * (x1 + x2)
        yt <- 0.5 * (y1 + y2)
        if (is.character (linkatt) && length (linkatt) == 1) 
            z <- net$links[[linkatt]][i]
        else z <- linkatt[i]
        text (xt, yt, z, cex=cex, font=2)
    }
    if (is.null (group)) invisible () else invisible (kol)
}

############################################################

net.map.key <- function (x, y, labels=NULL, cex=par("cex"), 
    pch=par("pch"), size=2.5*cex, col=NULL, head="", 
    sep=0.25*cex, horizontal=FALSE, new=FALSE)
#
# Draw key to (presumably) match map
#
# x,y         lower left coordinates of key 
#             in proportional units [0, 1]
# labels      vector of labels for classes
#             if NULL, integers 1:length(col), or "1"
# cex         par parameter, size of text
# pch         par parameter for type of symbols
# size        size of key symbols in cex units
# col         par parameter for color codes for segments
#             default ramp is yellow-red-brown
# head        text heading for key
# sep         separation in cex units between adjacent symbols
#             if sep=0 assume continuous scale and use gon=4
#             and put lables at breaks between squares
# new         if TRUE, call plot
# horizontal  if TRUE, key runs horizontal
#
# returned value is col or generated colors
#
{
    if (is.null (labels))
        if (is.null (col)) labels <- as.vector ("1")
        else labels <- as.character (seq (length (col) - 1))
    nsym <- length (labels)
    if (sep == 0) nsym <- nsym - 1
    if (is.null (col)) 
    {
        hue <- seq (0.167, 0.0, length=nsym)
        sat <- rep (1, nsym)
        bri <- c(1,1,1,1,1,0.9,0.8,0.7,0.6)[round 
            (seq(1,9,length=nsym))]
        kol <- hsv(hue, sat, bri)
    }
    else kol <- rep (col, nsym)
    if (new)
        plot(c(0,1), c(0,1), type="n", axes=FALSE, 
            xlab="", ylab="")
    oldadj <- par ("adj")
    par (adj=0)
    u <- par ("usr")
    ex <- par ("pin")[1]
    ey <- par ("pin")[2]
    uxr <- u[2] - u[1]
    uyr <- u[4] - u[3]
    px <- x * uxr + u[1]
    py <- y * uyr + u[3]
    if (horizontal)
    {
        halfy <- (size * par("cin")[2]) / 3
        xstep <- (size + sep) * par("cin")[1] / 2.0
        ystep <- halfy + 0.05
        hy <- halfy * uyr / ey
        dx <- xstep * uxr / ex
        dy <- ystep * uyr / ey
        qx <- px - dx
        qy <- py
    }
    else # vertical
    {
        halfx <- (size * par("cin")[1]) / 3
        xstep <- halfx + 0.05
        ystep <- (size + sep) * par("cin")[2] / 3.0
        hx <- halfx * uxr / ex
        dx <- xstep * uxr / ex
        dy <- ystep * uyr / ey
        qx <- px
        qy <- py - dy
    }
    if (horizontal)
    {
        if (sep == 0) 
        {
            for (i in 1:nsym) 
            {
                qx <- qx + dx
                points (qx, qy, pch=15, cex=size*1.1, col=kol[i])
                text (qx - dx, qy - dy, labels[i], cex=cex) 
            }
            text (qx, qy - dy, labels[nsym+1], cex=cex)
        }
        else # separated
        {
            for (i in 1:nsym) 
            {
                qx <- qx + dx
                points (qx, qy, col=kol[i], pch=pch, cex=size)
                text (qx - dx/2, qy - dy, labels[i], cex=cex) 
            }
        }
        if (length (head) > 0)
        {
            qx <- qx - (nsym * dx) / 2
            qy <- qy + (dy * length (grep ("$", head)))
            # if (sep == 0) qy <- qy + 0.5 * dy
            text (qx, qy, head, cex=cex, adj=0.5)
        }
    }
    else # vertical
    {
        if (sep == 0) 
        {
            for (i in 1:nsym) 
            {
                qy <- qy + dy
                points (qx, qy, pch=15, cex=size*1.1, col=kol[i])
                text (qx+dx, qy - dy/2, labels[i], cex=cex) 
            }
            text (qx+dx, qy + dy/2, labels[nsym+1], cex=cex)
        }
        else # separated
        {
            for (i in 1:nsym) 
            {
                qy <- qy + dy
                points (qx, qy, col=kol[i], pch=pch, cex=size)
                text (qx+dx, qy, labels[i], cex=cex) 
            }
        }
        if (length (head) > 0)
        {
            qy <- qy + (dy * length (grep ("$", head)))
            if (sep == 0) qy <- qy + 0.5 * dy
            text (qx-hx, qy, head, cex=cex)
        }
    }
    par (adj=oldadj)
    invisible (kol)
}

############################################################

net.write <- function (net, prefix="somenet")
#
# Write out a network in tables format.  Arguments:
#
# net     a network object
# prefix  a filename prefix, see below
#
# Writes out the links, segs, and cords of net in R table 
# format to three separate files.  The prefix is pasted 
# onto suffixes to make full filenames.
#
{
    filename <- paste (prefix, ".links.dat", sep="")
    write.table (net$links, file=filename, row.names=FALSE, 
        quote=FALSE)
    filename <- paste (prefix, ".segs.dat", sep="")
    write.table (net$segs, file=filename, row.names=FALSE,
        quote=FALSE)
    arcs <- matrix (0, nrow=0, ncol=2)
    colnames (arcs) <- c("x", "y")
    filename <- paste (prefix, ".cords.dat", sep="")
    for (i in 1:nrow (net$segs))
    {
        id <- as.integer (names (net$cords$x)[i])
        arcs <- rbind (arcs, c(id, NA))
        n <- length (net$cords$x[[i]]) - 1
        for (j in 1:n)
            arcs <- rbind (arcs, c(net$cords$x[[i]][j],
                net$cords$y[[i]][j]))
    }
    write.table (arcs, file=filename, row.names=FALSE, 
        quote=FALSE)
    invisible ()
}

############################################################

net.read <- function (prefix="somenet")
#
# Read in a network in tables format that has been written 
# by net.write.  Arguments:
#
# prefix  a filename prefix, see below
#
# Reads in the links, segs, and cords of a network in R table 
# format from three separate files.  The prefix is pasted 
# onto suffixes to make full filenames.
#
# Returns the net object created.
#
{
    filename <- paste (prefix, ".links.dat", sep="")
    links <- read.table (file=filename, header=TRUE)
    filename <- paste (prefix, ".segs.dat", sep="")
    segs <- read.table (file=filename, header=TRUE)    
    filename <- paste (prefix, ".cords.dat", sep="")
    arcs <- read.table (filename, header=TRUE)
    firsts <- is.na (arcs$y)
    pos <- which (c(firsts, TRUE))
    arcids <- arcs$x[firsts]
    arcs$x[firsts] <- NA
    cords <- list (x=list(), y=list())
    for (i in 1:length(arcids))
    {
        cords$x[[i]] <- arcs$x[(pos[i]+1):pos[i+1]]
        cords$y[[i]] <- arcs$y[(pos[i]+1):pos[i+1]]
    }
    names (cords$x) <- names (cords$y) <- arcids
    net <- list (links=net.orders(links), segs=segs, cords=cords)
    net$segs$length <- net.lengths (net)
    return (net)
}

############################################################

read.arcgenlin <- function (filename, coord=c("x","y")) 
#
# Read in an Arc/Info generate format lines file.
# Arguments:
#
# filename  name of file with generate data
# coord     the names to be given to the columns for the 
#           coordinates.  For stream.net, use the default.
#
# Returns the S/R format lines file with the Arc/Info user-id 
# placed in the first record of each line in the "x" coordinate 
# field accompanied by NA in the "y" coordinate field.
{
    string <- readLines (filename)
    sizes <- sapply (sapply (string, strwrap, width=1), length)
    x <- NULL
    for (i in 1:length(string)) 
    {
        if (sizes[i] == 1) 
        {
            if (string[i] == "END") x[i] <- "NA NA" 
            else x[i] <- paste(string[i], NA)
        }
        else x[i] <- string[i]
    }
    x <- strwrap (x, width=1)
    x[x == "NA"] <- NA
    x <- matrix (as.numeric(x), nrow=length(x)/2, ncol=2, byrow=TRUE)
    colnames (x) <- coord
    x <- x[! (is.na (x[,"x"]) & is.na (x[,"y"])),]
    data.frame (x)
}
