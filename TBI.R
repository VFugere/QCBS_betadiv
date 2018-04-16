TBI <- function(mat1,mat2,method="%difference", pa.tr=FALSE, nperm=99, permute.sp=1, BCD=TRUE, replace=FALSE, test.BC=TRUE, test.t.perm=FALSE, save.BC=FALSE, seed.=NULL, clock=FALSE)
# Temporal beta diversity analysis. 
# TBI: temporal beta diversity index.
#
# Author:: Pierre Legendre
{
### Internal functions --

dissim <- function(mat1, mat2, n, method, tr=TRUE, BCD, ref)
# tr =TRUE : The species data have been transformed by decostand in function transform()
# BCD=TRUE : Method is {"ruzicka", "%difference"} and output table BCD was requested
# ref=TRUE : The function is called to compute the reference values of the TBI dissimil.
{
vecD = vector(mode="numeric",length=n)     # to receive the values D=(B+C)/den
if(ref & BCD) { 
	vecB = vector(mode="numeric",length=n) # to receive the values B/den
	vecC = vector(mode="numeric",length=n) # to receive the values C/den
	v.B  = vector(mode="numeric",length=n) # to receive the values B
	v.C  = vector(mode="numeric",length=n) # to receive the values C
	} else { vecB=NA; vecC=NA; v.B=NA; v.C=NA }
#
# Compute the dissimilarity between T1 and T2 for each object (site)
#
# 1. If method is {"euclidean", "chord", "hellinger"}
#    compute the Euclidean distance
if(any(method == c("euclidean", "chord", "hellinger")))  
	for(i in 1:n) vecD[i] = dist(rbind(mat1[i,], mat2[i,])) 
#
# 2. Compute the Ruzicka or %difference dissimilarity 
# "ruzicka"          # Quantitative form of Jaccard
# "%difference"      # Quantitative form of Sørensen
if(any(method == c("ruzicka", "%difference"))) { 
	for(i in 1:n) {
		tmp = RuzickaD(mat1[i,], mat2[i,], method=method, BCD=BCD, ref=ref) 
		if(ref & BCD) {
			vecB[i] <- tmp$B.den
			vecC[i] <- tmp$C.den
			v.B[i]  <- tmp$B
			v.C[i]  <- tmp$C    }
    	vecD[i] <- tmp$D
		}
	}
list(vecB=vecB, vecC=vecC, vecD=vecD, v.B=v.B, v.C=v.C)
}
###
transform <- function(mat, method)
{
if(method=="chord")     mat <- decostand(mat, "norm")
if(method=="hellinger") mat <- decostand(mat, "hellinger")
mat
}
###
permute.comm <- function(mat1, mat2, permute.sp=1)
# Function to permute community composition data, species by species independently.
# permute.sp=1 : permute a species in the same way in mat1 and mat2.
# permute.sp=2 : permute a species independently in mat1 and mat2.
# This permutation algorithm is 2x faster than   mat1.perm <- apply(mat1,2,sample)
{
n = nrow(mat1)
p = ncol(mat1)
mat1.perm = matrix(NA,n,p)
mat2.perm = matrix(NA,n,p)
for(j in 1:p) {
	order = sample(n)
	mat1.perm[,j] <- mat1[order,j]
	if(permute.sp==2) order = sample(n)
	mat2.perm[,j] <- mat2[order,j]
	}
out = list(mat1.perm=mat1.perm, mat2.perm=mat2.perm)
}
### End internal functions

###
A <- system.time({
#
# Set "seed." to a specified integer, e.g. 1234, in function call to repeat a calculation
if(!is.null(seed.)) set.seed(seed.)   
#
epsilon <- sqrt(.Machine$double.eps)
method <- match.arg(method, c("%difference", "ruzicka", "chord", "hellinger", "log.chord", "jaccard", "sorensen", "ochiai", "euclidean")) 
n = nrow(mat1)
p = ncol(mat1)
if((nrow(mat2)!=n) | (ncol(mat2)!=p)) stop("The matrices are not of the same size.")
#
if(method=="log.chord") { 
	mat1   <- log1p(mat1)   # log1p() transformation done only once, before permutations
	mat2   <- log1p(mat2)   # log1p() transformation done only once, before permutations
	method <- "chord" }
if(method=="jaccard") { pa.tr=TRUE; method="ruzicka" }
if(method=="sorensen") { pa.tr=TRUE; method="%difference" }
if(method=="ochiai") { pa.tr=TRUE; method="chord" }
#
if(pa.tr) {
	mat1 <- ifelse(mat1>0, 1, 0)
	mat2 <- ifelse(mat2>0, 1, 0) }
if(any(method == c("chord", "hellinger"))) {
	tr <- TRUE
	require(vegan)
	} else { tr <- FALSE }
test.B.C <- NA 
if( (any(method == c("ruzicka", "%difference"))) & BCD) { 
	BCD.mat <- matrix(0,n,3)
	if(method=="%difference") colnames(BCD.mat) <- 
			c("B/(2A+B+C)","C/(2A+B+C)","D=(B+C)/(2A+B+C)")
	if(method=="ruzicka")    colnames(BCD.mat) <- 
			c("B/(A+B+C)","C/(A+B+C)","D=(B+C)/(A+B+C)")
	rownames(BCD.mat) <- paste("Site",1:n,sep=".")
	Change = vector(mode="character",length=n)
	} else {
	BCD <- FALSE 
	BCD.mat <- NA 
	BCD.summ <- NA 
	}
###
# 1. Compute the reference D for each object from corresponding vectors in the 2 matrices.
if(tr) { 
	tmp <- dissim(transform(mat1,method),transform(mat2,method),n,method,tr,BCD,ref=TRUE)
	} else { tmp <- dissim(mat1, mat2, n, method, tr, BCD, ref=TRUE) 
	}
vecD.ref <- tmp$vecD
BC <- NA
if(BCD) { 
	BCD.mat[,1]<-tmp$vecB ; BCD.mat[,2]<-tmp$vecC ; BCD.mat[,3]<-tmp$vecD 
	for(i in 1:n) {
		if(tmp$vecB[i]>tmp$vecC[i]) Change[i]="–  " else Change[i]="+  " }
	BCD.summ = matrix(NA,1,6)
	colnames(BCD.summ) = c("mean(B/den)","mean(C/den)","mean(D)","B/(B+C)","C/(B+C)", 
		"Change")
	BCD.means = apply(BCD.mat,2,mean, na.rm=TRUE)  # Exclude the sites with value = NA
	BCD.summ[1,1:3] = BCD.means
	BCD.summ[1,4:5] = BCD.means[1:2]/BCD.means[3]
	BCD.summ = as.data.frame(BCD.summ)
	if(BCD.summ[1,1]>BCD.summ[1,2]) BCD.summ[1,6]="–  " else BCD.summ[1,6]="+  "
	rownames(BCD.summ) = ""
	#
	BCD.mat <- as.data.frame(BCD.mat)
	BCD.mat = cbind(BCD.mat,Change)
	#
	if((n>4) & test.BC) {   # Tests of significance of difference between B/den and C/den
		test.B.C = matrix(NA,1,4)
		rownames(test.B.C) = "Paired t.test"
		#
		# Paired t-test between the vectors of B and C values
		if(test.t.perm) {
			t.res <- t.paired.perm(tmp$vecB,tmp$vecC,nperm=nperm,alternative="two.sided")
			test.B.C[1,] <- c(t.res$estim, t.res$t.ref, t.res$p.param, t.res$p.perm)
			p.value <- t.res$p.param
			} else {
			t.res = t.test(tmp$vecB, tmp$vecC, paired=TRUE, alternative = "two.sided")
			test.B.C[1,] = c(t.res$estimate, t.res$statistic, t.res$p.value, NA)
			p.value <- t.res$p.value
			}
		signif. = ifelse(p.value>0.05, " ","*")
		test.B.C = as.data.frame(test.B.C)
		test.B.C = cbind(test.B.C, signif.)
		colnames(test.B.C) <- c("  mean(B-C)","Stat","p.param","p.perm","  p<=0.05")
		} else { test.B.C <- NA }
	# Matrix containing the observed values of B and C, in case they are needed later
	if(save.BC) {
		BC <- cbind(tmp$v.B, tmp$v.C)
		colnames(BC) <- c("B", "C")
		rownames(BC) <- paste("Site",1:n,sep=".")
		} 
	}
###
if(permute.sp!=3) {   # Permute the data separately in each column.
# 2. Permutation methods 1 and 2 --
# Permute *the raw data* by columns. Permute the two matrices in the same way, 
# saving the seed before the two sets of permutations through sample(). 
# Permutation test for each distance in vector D.
# seed: seed for random number generator, used by the permutation function 
#       sample(). It is reset to that same value before permuting the values in the  
#       columns of the second matrix. 
	if(nperm>0) {
		BCD <- FALSE
		nGE.D = rep(1,n)
		for(iperm in 1:nperm) {
			if(permute.sp==1) {    # Permutation method 1
				seed <- ceiling(runif(1,max=100000))
				# cat("seed =",seed,'\n')
				set.seed(seed)
				mat1.perm <- apply(mat1,2,sample,replace=replace)
				set.seed(seed)
				mat2.perm <- apply(mat2,2,sample,replace=replace)
			} else {  # Permutation method 2 - Do not force the permutations 
 					  # to start at the same point in the two matrices.
				mat1.perm <- apply(mat1,2,sample,replace=replace)
				mat2.perm <- apply(mat2,2,sample,replace=replace)
				}
# 3. Recompute transformations of the matrices and the D values of the paired vectors.
			if(tr) { tmp <- dissim(transform(mat1.perm,method), 
							transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
			} else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
			vecD.perm <- tmp$vecD
			ge <- which(vecD.perm+epsilon >= vecD.ref)
			if(length(ge)>0) nGE.D[ge] <- nGE.D[ge] + 1
			}
# 4. Compute the p-value associated with each distance (i.e. site).
		p.dist <- nGE.D/(nperm+1)
		} else { p.dist <- NA }   # if nperm=0

} else if(permute.sp==3) {   
# 2.bis  Permutation method 3 -- 
# Permute entire rows in each matrix separately.
	if(nperm>0) {
		BCD <- FALSE
		seed <- ceiling(runif(1,max=100000))
		set.seed(seed)
		nGE.D = rep(1,n)
		for(iperm in 1:nperm) {
			mat1.perm <- mat1[sample(n,replace=replace),]
			mat2.perm <- mat2[sample(n,replace=replace),]
			#
# 3.bis Recompute the D values of the paired vectors.
			if(tr) { tmp <- dissim(transform(mat1.perm,method), 
							transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
			} else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
			vecD.perm <- tmp$vecD
			ge <- which(vecD.perm+epsilon >= vecD.ref)
			if(length(ge)>0) nGE.D[ge] <- nGE.D[ge] + 1
			}
# 4.bis Compute the p-value associated with each distance.
		p.dist <- nGE.D/(nperm+1)
		} else { p.dist <- NA }   # if nperm=0
}
p.adj <- p.adjust(p.dist,"holm")
})
A[3] <- sprintf("%2f",A[3])
if(clock) cat("Computation time =",A[3]," sec",'\n')
#
list(TBI=vecD.ref, p.TBI=p.dist, p.adj=p.adj, BCD.mat=BCD.mat, BCD.summary=BCD.summ, t.test_B.C=test.B.C, BC=BC)
}

RuzickaD <- function(vec1, vec2, method="ruzicka", BCD=FALSE, ref=TRUE)
#
# Compute the Ruzicka dissimilarity (quantitative form of the Jaccard dissimilarity)
# or the percentage difference (quantitative form of the Sørensen dissimilarity).
# A single dissimilarity is computed because there are only two data vectors.
#
# Arguments --
# vec1, vec2 : data vectors (species abundance or presence-absence data)
# method == c("ruzicka", "%difference")
# BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
#             For the %difference, they are B/(2A+B+C), C/(2A+B+C), D/(2A+B+C).
#             For the Ruzicka D, they are B/(A+B+C), C/(A+B+C), D/(A+B+C).
# BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
# ref=TRUE  : Compute the reference values of D, B and C
#    =FALSE : Under permutation, compute only the value of D. Use separate code (shorter).
#
# License: GPL-2 
# Author:: Pierre Legendre, April 2015
{
# This algorithm is applicable to matrices Y containing two data vectors only
#
A <- sum(pmin(vec1, vec2))          # A = sum of minima from comparison of the 2 vectors
sum.Y <- sum(vec1, vec2)            # Sum of all values in the two vectors, (2A+B+C)
#
if(ref) {    # Compute the reference values of statistics D, B and C
	tmp = vec1 - vec2
	B = sum(tmp[tmp>0])                 # Sum of the species losses between T1 and T2
	C = -sum(tmp[tmp<0])                # Sum of the species gains between T1 and T2
	D = B+C                             # Dissimilarity

	# Under permutation, compute only the value of D. - Shorter computation time.
	} else { 
	D <- sum.Y-2*A                      # (B+C)
	}
# Compute the denominator (den) of the Ruzicka or %difference index
if(method == "ruzicka") { den <-(sum.Y-A)  # den = (A+B+C)
	} else { den <- sum.Y }                # den = (2A+B+C)
if(!BCD) { B <- NA ; C <- NA }
list(B.den=B/den, C.den=C/den, D=D/den, B=B, C=C)
}

t.paired.perm <- function(vec1, vec2, nperm=999, alternative="two.sided", silent=TRUE)
#
# This function computes a permutation test of comparison of the means
# of two paired vectors (related samples).
# For each object, permutations are restricted to the two related observations.
#
# Arguments of the function
#
#    vec1, vec2: the two vectors to be compared
#    nperm = number of permutations (default value: 999)
#    alternative = c("two.sided", "less", "greater"). Dafault value: "two.sided"
#    silent = FALSE: calculation results are printed to the R console.
#           = TRUE : calculation results are not printed to the R console.
#
# Values returned
# 
#    t.ref : reference value of the t-statistic
#    p.param : parametric p-value
#    p.perm : permutational p-value
#    nperm : number of permutations
#
# Example: Deer leg length data from Zar (1999, p. 162)
#
# deer = matrix( c(142,140,144,144,142,146,149,150,142,148,138,136,147,139,143,141,143, 145,136,146),10,2)
# rownames(deer) = c('Deer.1','Deer.2','Deer.3','Deer.4','Deer.5','Deer.6','Deer.7', 'Deer.8','Deer.9','Deer.10')
# colnames(deer) = c('Hind.leg', 'Fore.leg')
#
# res = t.paired.perm(deer[,1], deer[,2])   # Two-tailed test by default
#
# Compare results to:  res2 = t.test(deer[,1], deer[,2], paired=TRUE) 
#
#                          Pierre Legendre, November 2009
#                          Guillaume Blanchet, September 2015 (faster permutation code)
{
n1 <- length(vec1)
n2 <- length(vec2)
if(n1 != n2) stop("The two vectors have different lengths. They cannot be paired.")

tail <- match.arg(alternative, c("two.sided", "less", "greater"))

res = t.test(vec1, vec2, paired=TRUE, alternative=tail)
t.ref =  res$statistic

# Print these first results
if(!silent) { 
	cat('\nt-test comparing the means of two related samples','\n','\n')
	cat('Number of objects:',n1,'\n')
	cat('Mean of the differences:',res$estimate,'\n')
	cat('t statistic (paired observations):',t.ref,'\n')
	cat('95 percent confidence interval of t:',res$conf.int,'\n')
	cat('Degrees of freedom:',res$parameter,'\n')
	cat('Alternative hypothesis:',tail,'\n')
	cat('Prob (parametric):',res$p.value,'\n')
	}

# Perform the permutation test
# Permutations are restricted to the two related observations for each object.
nPGE <- 1

for(i in 1:nperm)
	{
## New permutation code, GB
	mat <- cbind(vec1,vec2)
	topermute <- rbinom(n1,1,0.5)
	mat[topermute==1,] <- mat[topermute==1,2:1]
## End new code, GB
	
	res.perm = t.test(mat[,1], mat[,2], paired=TRUE, alternative=tail)
	t.perm = res.perm$statistic
	
	if(tail == "two.sided") if( abs(t.perm) >= abs(t.ref) ) nPGE <- nPGE+1
	if(tail == "less")      if(t.perm <= t.ref) nPGE <- nPGE+1
	if(tail == "greater")   if(t.perm >= t.ref) nPGE <- nPGE+1
	}

# Compute and print the permutational p-value
P <- nPGE/(nperm+1)
if(!silent) cat('Prob (',nperm,'permutations):', formatC(P,digits=5,width=7,format="f"),'\n')
#
return(list(estim=res$estimate, t.ref=t.ref, p.param=res$p.value, p.perm=P, nperm=nperm))
}


# Examples -- 
# data(mite)
# res1 = TBI(mite[1:10,],mite[61:70,],method="%diff",nperm=999,permute.sp=1)
# 
# Example using permute.sp=3. This method is not recommended (low power).
# res2 = TBI(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=3)
