kmeansDisc <- function(E) {
    minmax <- function(x, min=0, max=1, quantiles=FALSE) {
        if (quantiles) {
            min = quantile(x, min, na.rm=TRUE)
            max = quantile(x, max, na.rm=TRUE)
        }
        x[x < min] <- min
        x[x > max] <- max
        x
    }

    qminmax <- function(x, min=0.01, max=0.99) {
        minmax(x, min, max, TRUE)
    }

    n <- nrow(E)
    out <- matrix(as.integer(0), nrow = n, ncol = ncol(E))
    for (i in 1:n) {
        kms <- kmeans(
            qminmax(E[i, ], 0.002, 0.998),
            centers = quantile(E[i, ], c(0.002, 0.5, 0.998)),
            iter.max = 1000
        )
        if (kms$ifault == 4) {
            kms <- kmeans(
                qminmax(E[i, ], 0.002, 0.998),
                centers = quantile(E[i, ], c(0.002, 0.5, 0.998)),
                iter.max = 10000,
                algorithm = "Lloyd"
            )
        }
        if( !(length(unique(kms$cluster)) == 3) ) stop(paste0("Gene ", i, " was not discretized into 3 groups"))
        out[i, ] <- kms$cluster
    }
    storage.mode(out) <- "integer"
    return(out)
}