# simulation based test for a correlation between the state of x & the rate of y
# written by Liam J. Revell 2013
# modified by Chad M. Eliason 2020 for multivariate data (in y)
ratebystate_multi <- function(tree, x, y, nsim = 100, corr = c("pearson","spearman"), sim.method = "sim.corrs", model = c("BM", "lambda"), ...){
# tree <- kingtrees[[1]]
# plot(tree)
# head(kingdata)
# x=as.matrix(kingdata)[,"range"]
# y=kingdata[,1:6]
# library(motmot)
# y=brain_pca$x[tree$tip,1:44]
# transformPhylo.ML(tree, y=as.matrix(y), model="lambda")$Lambda[1]
# transformPhylo.ML(tree, y=as.matrix(x), model="lambda")$Lambda[1]
	model <- match.arg(model)
	if (model=="lambda") {
		tmp <- transformPhylo.ML(tree, y=y, model="lambda")
		# tmp$Lambda[1]
		Vy <- tmp$brownianVariance
		# tree <- geiger::rescale(tree, "lambda", lambda=lambday)
		# Vy <- phyl.vcv(y,vcv(tree),lambda=lambday)$R
		# image(Vy)
		# diag(Vy)
		tmp <- transformPhylo.ML(tree, y=as.matrix(x), model="lambda")#$Lambda[1]
		# tree <- geiger::rescale(tree, "lambda", lambda=lambdax)
		Vx <- tmp$brownianVariance
		# Vx <- phyl.vcv(as.matrix(x),vcv(tree),lambda=lambdax)$R
	}
	if (model=="BM") {
		V <- phyl.vcv(cbind(x,y),vcv(tree),lambda=1)$R
	}
	# corr<-corr[1]
	corr <- match.arg(corr)
	if(hasArg(sim.method)) sim.method<-list(...)$sim.method
	else sim.method<-"sim.corrs"
	if(hasArg(method)) method<-list(...)$method
	else method<-"by.node"
	if(hasArg(message)) message<-list(...)$message
	else message<-TRUE
	if(hasArg(logarithm)) logarithm<-list(...)$logarithm
	else logarithm<-FALSE
	if(!is.binary.tree(tree)) tree<-multi2di(tree)


# image(var(V[-1,-1]))
# image(var(y))

	# V <- phyl.vcv(y, vcv(tree), lambda=1)$R
	if (method=="by.branch"){
		aa<-c(x[tree$tip.label],fastAnc(tree,x))
		names(aa)[1:length(tree$tip)]<-1:length(tree$tip)
		aa<-rowMeans(matrix(aa[tree$edge],nrow(tree$edge),2))
		a<-vector()
		for(i in 1:tree$Nnode+length(tree$tip)){
			j<-which(tree$edge[,1]==i)
			a[i-length(tree$tip)]<-sum(aa[j]*tree$edge.length[j])/sum(tree$edge.length[j])
		}
		names(a)<-1:tree$Nnode+length(tree$tip)
	}	
	else a <- fastAnc(tree,x)
	if (logarithm) a <- exp(a)
	if (is.matrix(y)) {
		b <- sapply(1:ncol(y), function(i) {pic(phy = tree, x = y[, i])^2})
        b <- apply(b, 1, sum)
        b <- b[names(a)]
    } else {
        b<-pic(y,tree)[names(a)]^2
    }
	r.obs <- cor(a,b,method=corr)
	beta <- setNames(lm(b~a)$coefficients[2],NULL)
	foo <- function(tree,V){
		# if(sim.method=="fastBM") {
		# 	XY<-fastBM(tree,nsim=2)%*%sqrt(diag(diag(V)))
		# }
		if (sim.method=="sim.corrs") {
			XY <- sim.corrs(tree,V)
# plot(XY[,2:3])
# points(y[,1:2], pch=16)
		}
		a <- fastAnc(tree,XY[,1])
		if (is.matrix(y)) {
			b <- sapply(2:ncol(XY), function(i) {pic(phy = tree, x = XY[, i])^2})
        	b <- apply(b, 1, sum) # TODO check this
        	b <- b[names(a)]
		} else {
			b <- pic(XY[,2],tree)[names(a)]^2  # nodewise rates
		}
		# b <- picmulti(tree=tree, X=XY[,-1], scaled=T)$rates
		cor(a, b, method=corr)
	}
	# nsim=10
	# foo(tree,V)
	r.null <- c(r.obs, replicate(nsim-1, foo(tree,V)))
	P <- mean(abs(r.null)>=abs(r.obs))
	# if (message) {
	return(list(beta=beta, r=r.obs, r.null=r.null, P=P, corr=corr, method=method))
	# } else {
		# return(list(beta=beta,r=r.obs,P=P))
	# }
}

# function simulates rate by state evolution for x & y
# written by Liam J. Revell 2013
sim.ratebystate<-function(tree,sig2x=1,sig2y=1,beta=c(0,1),...){
	if(hasArg(method)) method<-list(...)$method
	else method<-"by.node"
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(logarithm)) logarithm<-list(...)$logarithm
	else logarithm<-FALSE
	x<-fastBM(tree,a=if(logarithm) beta[1] else 0,sig2=sig2x,internal=TRUE)
	if(method=="by.node") ss<-x[1:tree$Nnode+length(tree$tip.label)]
	else if(method=="by.branch") ss<-rowMeans(matrix(x[tree$edge],nrow(tree$edge),2))
	zz<-tree
	if(!logarithm) zz$edge.length<-beta[2]*zz$edge.length*(beta[1]+ss-min(ss))
	else zz$edge.length<-beta[2]*zz$edge.length*exp(ss)
	y<-fastBM(zz,sig2=sig2y)
	if(plot) phenogram(zz,x,type="b",color="blue",ftype="off",
		xlab="expected variance",ylab="independent variable (x)")
	x<-x[tree$tip.label]
	return(cbind(x,y))
}