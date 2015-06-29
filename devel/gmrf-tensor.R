#' @export
tegmrf <- function (..., k = NA, bs = "cr", m = NA, d = NA, by = NA, fx = FALSE, 
    mp = TRUE, np = TRUE, xt = NULL, id = NULL, sp = NULL) {
    vars <- as.list(substitute(list(...)))[-1]
    dim <- length(vars)
    by.var <- deparse(substitute(by), backtick = TRUE)
    term <- deparse(vars[[1]], backtick = TRUE)
    if (dim > 1) 
        for (i in 2:dim) term[i] <- deparse(vars[[i]], backtick = TRUE)
    for (i in 1:dim) term[i] <- attr(terms(reformulate(term[i])), 
        "term.labels")
    if (sum(is.na(d)) || is.null(d)) {
        n.bases <- dim
        d <- rep(1, dim)
    }
    else {
        d <- round(d)
        ok <- TRUE
        if (sum(d <= 0)) 
            ok <- FALSE
        if (sum(d) != dim) 
            ok <- FALSE
        if (ok) 
            n.bases <- length(d)
        else {
            warning("something wrong with argument d.")
            n.bases <- dim
            d <- rep(1, dim)
        }
    }
    if (sum(is.na(k)) || is.null(k)) 
        k <- 5^d
    else {
        k <- round(k)
        ok <- TRUE
        if (sum(k < 3)) {
            ok <- FALSE
            warning("one or more supplied k too small - reset to default")
        }
        if (length(k) == 1 && ok) 
            k <- rep(k, n.bases)
        else if (length(k) != n.bases) 
            ok <- FALSE
        if (!ok) 
            k <- 5^d
    }
    if (sum(is.na(fx)) || is.null(fx)) 
        fx <- rep(FALSE, n.bases)
    else if (length(fx) == 1) 
        fx <- rep(fx, n.bases)
    else if (length(fx) != n.bases) {
        warning("dimension of fx is wrong")
        fx <- rep(FALSE, n.bases)
    }
    xtra <- list()
    if (is.null(xt) || length(xt) == 1) 
        for (i in 1:n.bases) xtra[[i]] <- xt
    else if (length(xt) == n.bases) 
        xtra <- xt
    else stop("xt argument is faulty.")
    if (length(bs) == 1) 
        bs <- rep(bs, n.bases)
    if (length(bs) != n.bases) {
        warning("bs wrong length and ignored.")
        bs <- rep("cr", n.bases)
    }
    bs[d > 1 & (bs == "cr" | bs == "cs" | bs == "ps" | bs == 
        "cp")] <- "tp"
    if (!is.list(m) && length(m) == 1) 
        m <- rep(m, n.bases)
    if (length(m) != n.bases) {
        warning("m wrong length and ignored.")
        m <- rep(0, n.bases)
    }
    if (!is.list(m)) 
        m[m < 0] <- 0
    if (length(unique(term)) != dim) 
        stop("Repeated variables as arguments of a smooth are not permitted")
    j <- 1
    margin <- list()
    for (i in 1:n.bases) {
        j1 <- j + d[i] - 1
        if (is.null(xt)) 
            xt1 <- NULL
        else xt1 <- xtra[[i]]
        stxt <- "s("
        for (l in j:j1) stxt <- paste(stxt, term[l], ",", sep = "")
        stxt <- paste(stxt, "k=", deparse(k[i], backtick = TRUE), 
            ",bs=", deparse(bs[i], backtick = TRUE), ",m=", deparse(m[[i]], 
                backtick = TRUE), ",xt=xt1", ")")
        margin[[i]] <- eval(parse(text = stxt))
        j <- j1 + 1
    }
    if (mp) 
        mp <- TRUE
    else mp <- FALSE
    if (np) 
        np <- TRUE
    else np <- FALSE
    full.call <- paste("tegmrf(", term[1], sep = "")
    if (dim > 1) 
        for (i in 2:dim) full.call <- paste(full.call, ",", term[i], 
            sep = "")
    label <- paste(full.call, ")", sep = "")
    if (!is.null(id)) {
        if (length(id) > 1) {
            id <- id[1]
            warning("only first element of `id' used")
        }
        id <- as.character(id)
    }
    ret <- list(margin = margin, term = term, by = by.var, fx = fx, 
        label = label, dim = dim, mp = mp, np = np, id = id, 
        sp = sp, inter = FALSE)
    class(ret) <- "gmrftensor.smooth.spec"
    ret
}

#' @export
smooth.construct.gmrftensor.smooth.spec <- function (object, data, knots) {
    inter <- object$inter
    m <- length(object$margin)
    if (inter) {
        object$mc <- if (is.null(object$mc)) 
            rep(TRUE, m)
        else as.logical(object$mc)
    }
    else {
        object$mc <- rep(FALSE, m)
    }
    Xm <- list()
    Sm <- list()
    nr <- r <- d <- array(0, m)
    C <- NULL
    object$plot.me <- TRUE
    for (i in 1:m) {
        knt <- dat <- list()
        term <- object$margin[[i]]$term
        for (j in 1:length(term)) {
            dat[[term[j]]] <- data[[term[j]]]
            knt[[term[j]]] <- knots[[term[j]]]
        }
        object$margin[[i]] <- if (object$mc[i]) 
            smoothCon(object$margin[[i]], dat, knt, absorb.cons = TRUE, 
                n = length(dat[[1]]))[[1]]
        else smooth.construct(object$margin[[i]], dat, knt)
        Xm[[i]] <- object$margin[[i]]$X
        if (!is.null(object$margin[[i]]$te.ok)) {
            if (object$margin[[i]]$te.ok == 0) 
                stop("attempt to use unsuitable marginal smooth class")
            if (object$margin[[i]]$te.ok == 2) 
                object$plot.me <- FALSE
        }
        if (length(object$margin[[i]]$S) > 1) 
            stop("Sorry, tensor products of smooths with multiple penalties are not supported.")
        Sm[[i]] <- object$margin[[i]]$S[[1]]
        d[i] <- nrow(Sm[[i]])
        r[i] <- object$margin[[i]]$rank
        nr[i] <- object$margin[[i]]$null.space.dim
        if (!inter && !is.null(object$margin[[i]]$C) && nrow(object$margin[[i]]$C) == 
            0) 
            C <- matrix(0, 0, 0)
    }
    XP <- list()
    if (object$np) 
        for (i in 1:m) {
            if (object$margin[[i]]$dim == 1) {
                if (!inherits(object$margin[[i]], c("cs.smooth", 
                  "cr.smooth", "cyclic.smooth", "random.effect"))) {
                  x <- get.var(object$margin[[i]]$term, data)
                  np <- ncol(object$margin[[i]]$X)
                  knt <- if (is.factor(x)) {
                    unique(x)
                  }
                  else {
                    seq(min(x), max(x), length = np)
                  }
                  pd <- data.frame(knt)
                  names(pd) <- object$margin[[i]]$term
                  sv <- if (object$mc[i]) 
                    svd(PredictMat(object$margin[[i]], pd))
                  else svd(Predict.matrix(object$margin[[i]], 
                    pd))
                  if (sv$d[np]/sv$d[1] < .Machine$double.eps^0.66) {
                    XP[[i]] <- NULL
                    warning("reparameterization unstable for margin: not done")
                  }
                  else {
                    XP[[i]] <- sv$v %*% (t(sv$u)/sv$d)
                    Xm[[i]] <- Xm[[i]] %*% XP[[i]]
                    Sm[[i]] <- t(XP[[i]]) %*% Sm[[i]] %*% XP[[i]]
                  }
                }
                else XP[[i]] <- NULL
            }
            else XP[[i]] <- NULL
        }
    for (i in 1:m) Sm[[i]] <- Sm[[i]]/eigen(Sm[[i]], symmetric = TRUE, 
        only.values = TRUE)$values[1]
    max.rank <- prod(d)
    r <- max.rank * r/d
    X <- tensor.prod.model.matrix(Xm)
    if (object$mp) {
        S <- tensor.prod.penalties(Sm)
        for (i in m:1) if (object$fx[i]) {
            S[[i]] <- NULL
            r <- r[-i]
        }
    }
    else {
        warning("single penalty tensor product smooths are deprecated and likely to be removed soon")
        S <- Sm[[1]]
        r <- object$margin[[i]]$rank
        if (m > 1) 
            for (i in 2:m) {
                S <- S %x% Sm[[i]]
                r <- r * object$margin[[i]]$rank
            }
        if (sum(object$fx) == m) {
            S <- list()
            object$fixed = TRUE
        }
        else {
            S <- list(S)
            object$fixed = FALSE
        }
        nr <- max.rank - r
        object$bs.dim <- max.rank
    }
    object$X <- X
    object$S <- S
    if (inter) 
        object$C <- matrix(0, 0, 0)
    else object$C <- C
    object$df <- ncol(X)
    object$null.space.dim <- prod(nr)
    object$rank <- r
    object$XP <- XP
    class(object) <- "tensor.smooth"
    object
}

