require(ape)
require(mvtnorm)
require(motmot)
require(phangorn) # for Ancestors

# DONE have it report lambda values of trees/traits
essim <- function(phy, trait, nsim = 1000, es, log=FALSE, model=c("lambda", "BM"), return.es=FALSE) {
    model <- match.arg(model)
    # calculate ES speciation rates if not suppied as `es` argument
    if (missing(es)) {
        es <- compute_es(phy)
    }
    # check if all tips have traits. will work whether newly computed or passed in
    if (length(trait) < length(phy$tip.label)) {
        idx <- which(!phy$tip.label %in% names(trait));
        phy <- drop.tip(phy, idx);
    }
    if (log) {
        es <- log(es[phy$tip.label]) # log transform
    }
	trait <- trait[phy$tip.label]
	# Pearson's correlation between log inverse equal splits statistic and trait
	res <- cor.test(es, trait, method="pearson")
	# Brownian simulations
	if (model=="BM") {
        # Fit Brownian motion model to get diffusion rate and root state estimates
        vv <- vcv.phylo(as.phylo(phy))
        onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
        root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
        sig2 <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
        sims <- t(rmvnorm(nsim, sigma=sig2*vv))
        rownames(sims) <- rownames(vv)
        lambda <- 1
    }
    # Lambda simulations
    if (model=="lambda") {
        fit <- transformPhylo.ML(y=as.matrix(trait), phy=phy, model="lambda")
        lambda <- fit$Lambda[,"MLLambda"]
        sig2 <- fit$brownianVariance[1,1]
        anc <- fit$root.state
        sims <- rTrait(n=nsim, phy=phy, model="lambda", parameters=list(lambda=lambda, sigma2=sig2, ancestral.state=anc))
    }
	# Pearson's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(es[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
    names(sim.r) <- NULL

	# Calculate the two-tailed p value
	corr <- res$estimate
	# upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
	# lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
	# pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed
    pval <- mean(abs(sim.r) >= abs(corr))  # see REvell's ratebystate

    if (!return.es) { # output just rho and p value
		result <- list(rho=corr, p=pval, lambda=lambda, sig2=sig2, rhosims=sim.r)
		# names(result) <- c("rho", "P Value", "lambda", "sig2")
		return(result)
	} else { # output rho, p value, and list of es values
		result <- as.vector(c(corr, pval, list(es), lambda, sig2))
		names(result) <- c("rho", "P Value", "es", "lambda", "sig2")
		return(result)		
	}
}

compute_es <- function (phy, node=FALSE) {
    Ntip <- length(phy$tip.label);
    nnode <- Nnode(phy)
    rootnd <- Ntip + 1L;
    if (node) {
        ids <- c(1:Ntip, (Ntip+1) : (Ntip + nnode))    
    } else {
        ids <- c(1:Ntip)
    }
    es <- numeric(length(ids));
    for (k in 1:length(ids)) {
        # k = 11
        # get lineage
        lin <- c(ids[k], Ancestors(phy, ids[k]));
        # drop root node
        lin <- lin[-length(lin)];
        if (length(lin) == 0) {
            es[k] <- 0
            next
        }
        # set indices
        inds <- seq(1, length(lin), 1);
        # get els
        els <- phy$edge.length[match(lin, phy$edge[,2])];
        # compute in one
        es[k] <- 1/sum(els * (1 / 2^(inds-1)));
    }
    nms <- c(phy$tip.label, (Ntip(phy)+1) : (Ntip(phy) + Nnode(phy)) )
    names(es) <- nms[ids]
    return(es)
}

# test:
# phy <- rcoal(10)
# plot(phy, cex=.5)
# nodelabels(cex=.5)
# dr <- compute_es(phy, node=T)
# plot(phy, edge.width=dr[phy$edge[,2]])
# nodelabels(pch=1, cex=dr, col='red')


# here we want to see if rates of trait evolution covary with diversification rates
# doesn't make sense (to me) to treat species-specific rates as an evolving trait
# rather, we want to simulate several PC traits and get multivariate (RRphylo) rates 
# iteratively and then look at correlation between trait rates and speciation

# trait <- brain_pca$x[trees[[1]]$tip, 1:44]
# image(var(traits))
# phy=trees[[1]]
# trait=traits[,1:41]
# nsim=10
# es=speciation
# model="BM"
# trait = traits[,1]

essim_rates <- function(phy, trait, nsim = 1000, es, log=FALSE, model=c("lambda", "BM"), return.es=FALSE) {
    
    require(pbapply)

    model <- match.arg(model)

    # calculate ES speciation rates if not suppied as `es` argument
    if (missing(es)) {
        es <- compute_es(phy)
    }
    
    trait <- as.matrix(trait)
    
    # check if all tips have traits. will work whether newly computed or passed in
    if (nrow(trait) < length(phy$tip.label)) {
        idx <- which(!phy$tip.label %in% names(trait))
        phy <- drop.tip(phy, idx)
    }

    if (log) {
        es <- log(es) # log transform
    }

    es <- es[phy$tip.label]

    trait <- trait[phy$tip.label, , drop=F]

    rr <- RRphylo(phy, trait)
    id <- c((Ntip(phy)+1):(Ntip(phy)+Nnode(phy)), 1:Ntip(phy))
    tiprates <- abs(rr$rates[,1][match(phy$edge[,2], id)][phy$tip])

    tiprates <- log(tiprates)
# model="BM"
    # Pearson's correlation between log inverse equal splits statistic and trait
    obs <- cor.test(es, tiprates, method="pearson")

# plot(es~tiprates)

    # Brownian simulations 
    fit <- transformPhylo.ML(y=as.matrix(trait), phy=phy, model=model)
    sig2 <- diag(fit$brownianVariance)
    anc <- fit$root.state

    # sims <- lapply(1:nsim, rmvnorm, n=nrow(trait), mean=anc, sigma=sig2)
    sims <- lapply(1:nsim, function(i) {
        sapply(1:ncol(trait), function(x) {
            # x=44
            if (model=="BM") {
                res=rTrait(n=1, phy=phy, model="BM", parameters=list(sigma2=sig2[x], ancestral.state=anc[x]))
            }
            if (model=="lambda") {
                lambda <- fit$Lambda[,"MLLambda"]
                res=rTrait(n=1, phy=phy, model="lambda", parameters=list(lambda=lambda, sigma2=sig2[x], ancestral.state=anc[x]))
            }
            res
        })
    })

    # sims <- rTrait(n=nsim, phy=phy, model="lambda", parameters=list(lambda=lambda, sigma2=sig2, ancestral.state=anc))
    # phenogram(phy, trait[,44], ftype="off")
    # phenogram(phy, sims[[3]][,44], ftype="off", add=T, col='red')

    simrates <- pbsapply(1:length(sims), function(x) {
        rr <- RRphylo(phy, sims[[x]])
        abs(rr$rates[,1][match(phy$edge[,2], id)][phy$tip])
    })

    # ln-transform rates
    simrates <- log(simrates)

    # Pearson's correlations of simulated datasets
    sim.r <- sapply(1:nsim, function(x) cor(es, simrates[,x], method="pearson"))
    
    # Calculate the two-tailed p value
    corr <- obs$estimate
# hist(sim.r)
# abline(v=corr)

# phenogram(phy, tiprates, fsize=.5)
# phenogram(phy, simrates[, 1], ftype="off")

    # upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
    # lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
    # pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed
    pval <- mean(abs(sim.r) >= abs(corr))  # see REvell's ratebystate

    if (missing(return.es)) { # output just rho and p value
        result <- as.vector(c(corr, pval, lambda, sig2))
        names(result) <- c("rho", "P Value", "lambda", "sig2")
        return(result)
    } else { # output rho, p value, and list of es values
        result <- as.vector(c(corr, pval, list(es), lambda, sig2))
        names(result) <- c("rho", "P Value", "es", "lambda", "sig2")
        return(result)      
    }
}

essim_rates_multi <- function(phy, trait, nsim, es, log=FALSE) {
  # V <- phyl.vcv(X=trait, C=vcv(phy), lambda=lambda)$R
  fit <- transformPhylo.ML(y=trait, phy=kingtree, model="lambda")
  lambda <- fit$Lambda[1]
  V <- fit$brownianVariance
  # image(V)
  # image(cov(trait))
  # image(cov(simtraits[[2]]))
  a <- fit$root.state
  # a <- apply(trait, 2, fastAnc, tree=phy)
  # image(cov(transformPhylo.sim(phy=phy, n=7, model="lambda", lambda=lambda, rate=diag(V))))
  # image(cov(trait))
  phy2 <- geiger::rescale(phy, "lambda", lambda=lambda)
  # Simulate trait evolution under fit model (lambda~0)
  simtraits <- mvSIM(tree = phy2, model = "BM1", param = list(theta=a, sigma = V), nsim = nsim)
  # image(cov(simtraits[[3]]))
  # Recalculate trait evolution rates using ridge regression
  simrates <- pbsapply(seq_along(simtraits), function(x) {
    rr <- RRphylo(phy, simtraits[[x]])
    id <- c((Ntip(phy)+1):(Ntip(phy)+Nnode(phy)), 1:Ntip(phy))
    erates <- rr$rates[, 1][match(phy$edge[,2], id)]
    erates[phy$tip]
  })
  if (log) {
    es <- log(es)
  }
  # Compute null correlation between trait and speciation rates
  cors <- sapply(1:ncol(simrates), function(x) {
    cor(es, log(simrates[, x]))
  })
  list(straits=simtraits, srates=simrates, scors=cors)
}
