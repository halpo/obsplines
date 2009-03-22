setClass("SplineBasis",representation(knots="numeric", order="integer", Matrices="array", weights="numeric"))
setClass("OrthogonalSplineBasis", representation("SplineBasis", transformation="matrix"))

SplineBasis<-function(knots, order=4, keep.duplicates=FALSE) { 
	#degree<-order-1
	n<-length(knots)
	if(any(table(knots[order:(n-order+1)])>1)&&!keep.duplicates){
		warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
		knots<-unique(knots[order:(n-order+1)])
		knots<-knots[c(1,1,1,seq(length(knots)),length(knots),length(knots),length(knots))]
	}
	q<-n-order
	M<-FindSplineMatrices(knots,order)
	M2<-array(dim=c(order,q,n-2*order+1))
	for(i in order:(n-order)) {  #Identifying interior intervals
		M2[,,i-order+1]<-M[,,i]%*%I_Sel(i,order,n)
	}
	new("SplineBasis", knots=knots, order=as.integer(order), Matrices=M2, weights = rep(1,n-1))
}
OrthogonalizeBasis<-function(object,...){
	obase<-new("OrthogonalSplineBasis",object)
	M<-object@Matrices
	k<-object@order
	Delta<-1/Henkel(1:(2*k-1),k,k)
	
	di<-diff(object@knots[(k):(length(object@knots)-k+1)])
	d<-dim(M)
	s<-matrix(0,d[2],d[2])
	for(i in 1:d[3]) { 
		s<-s+di[i]*t(M[,,i])%*%Delta%*%M[,,i]
	}
	L<-solve(chol(s))
	N<-M
	for( i in 1:d[3]) { 
		N[,,i]<-M[,,i]%*%L
	}
	
	obase@Matrices<-N
	obase@transformation<-L

	return(obase)
}
OrthogonalSplineBasis<-function(knots,order=4)OrthogonalizeBasis(SplineBasis(knots,order))
dim.SplineBasis<-function(x)dim(x@Matrices)
EvaluateBasis<-function(object,x,...) { 
	stopifnot(is.numeric(x))
	dots<-list(...)
	if(length(x)>1){
		results<-t(sapply(x,EvaluateBasis, object=object))
	} else {
		M<-object@Matrices
		knots<-object@knots
		order<-object@order
		w<-object@weights
		if(x < knots[order] | x>knots[length(knots)-order+1])return(rep(NA,dim(object)[2]))
		ind<- x<knots
		if(all(ind)|all(!ind))  {
			if(x==knots[length(knots)-order+1]) 
				return(rep(1,order)%*%M[,,dim(M)[3]]) 
			else 
				return(rep(0,n+1))
		}
		i<-which(ind)[1]-1
		u<-(x-knots[i])/(knots[i+1]-knots[i])
		U<-u^(0:(order-1))
		return(w[i]*U%*%M[,,i-order+1])
	}
} 
deriv.SplineBasis<-function(object,l=1) {
	d<-dim(object)
	DM<-DerivativeMatrix(d[1])
	M<-object@Matrices
	k<-object@order
	knots<- object@knots
	
	N<-array(0,dim=d)
	for( i in 1:d[3]) { 
		N[,,i]<-t(MatrixPower(DM,l))%*%M[,,i]
	}
	
	w<-object@weights*diff(knots)^-l
	new("SplineBasis",knots=knots, order=object@order, Matrices=N, weights=w)
}


print.SplineBasis<-function(object) { 
	cat("Spline Basis\n")
	cat("Order: ",object@order,"\n",
		"Degree: ",object@order-1,"\n",
		"Knots: ",paste(object@knots,collapse=" "),"\n",sep="")
	invisible(object)
}
print.OrthogonalSplineBasis<-function(object) { 
	cat("Orthogonalized ")
	invisible(print.SplineBasis(object))
}
plot.SplineBasis<-function(x,y,xlab, ylab, main,...) { 
	basis<-x
	dots<-list(...)
	plotdata<-seq(basis@knots[basis@order],basis@knots[length(basis@knots)-basis@order+1],length=1000)
	xlabel<-if(!hasArg(xlab))"" else dots$xlab
	ylabel<-if(!hasArg(ylab))"Basis Functions" else dots$ylab
	ptitle<-if(!hasArg(main)) paste(substitute(x)) else dots$main
	matplot(plotdata,evaluate(basis,as(plotdata,"numeric")),type='l',xlab=xlabel,ylab=ylabel,main=ptitle,...)
}
OuterProdSecondDerivative<-function(basis){
	M<-basis@Matrices
	k<-basis@order
	d<-dim(M)
	Delta<-1/Henkel(1:(2*k-1),k,k)
	D2<-MatrixPower(DerivativeMatrix(k),2)
	OPSD<-0
	di<-diff(basis@knots[(k):(length(basis@knots)-k+1)])
	for(i in 1:d[3])OPSD<-OPSD+di[i]*t(M[,,i])%*%D2%*%Delta%*%t(D2)%*%M[,,i]
	OPSD
}

#  Basis Methods
setGeneric("orthogonalize",function(object,...)standardGeneric("orthogonalize"))
setMethod("orthogonalize",signature("SplineBasis"),OrthogonalizeBasis,valueClass="OrthogonalSplineBasis")
setMethod("dim","SplineBasis", dim.SplineBasis)
setMethod("show","SplineBasis",print.SplineBasis)
setMethod("show","OrthogonalSplineBasis",print.OrthogonalSplineBasis)
setGeneric("evaluate",function(object, x,...)standardGeneric("evaluate"))
setMethod("evaluate",signature("SplineBasis","numeric"),function(object, x, ...)EvaluateBasis(object=object, x=x, ...))
setMethod("plot",signature(x="SplineBasis",y="missing"),plot.SplineBasis)
setMethod("deriv", signature("SplineBasis"), function(expr,...)deriv.SplineBasis(object=expr,...))

#Other Helper Functions
DerivativeMatrix<-function(n){
	A<-rbind(0,diag(x=1,nrow=n-1,ncol=n))
	B<-diag(1:n,nrow=n,ncol=n)
	A%*%B
}
MatrixPower<-function(A,n){
	n<-as.integer(n)
	stopifnot(n>=1)
	stopifnot(ncol(A)==nrow(A))
	B<-diag(1,nrow(A))
	for(i in seq_len(n))B<-B%*%A
	return(B)
}
I_Sel<-function(i,k,n){
	cbind(
		matrix(0,nrow=k,ncol=i-(k)),
		diag(nrow=k),
		matrix(0,nrow=k,ncol=n-k-i)
		)
}
FindSplineMatrices<-function(knots,k){
n<-length(knots)-1 # n+1= number of knots

D0<-function(i,j,k){
	a<-knots[i]-knots[j]
	b<-knots[j+k-1]-knots[j]
	rtn<-a/b
	rtn[a==0]<-0
	rtn
}
D1<-function(i,j,k){
	a<-knots[i+1]-knots[i]
	b<-knots[j+k-1]-knots[j]
	rtn<-a/b
	rtn[a==0]<-0
	rtn
} 
M<-function(k,i) { 
	if(k==1) return(1);
	if(i<k|i>length(knots)-k) return(matrix(0,k,k));
	rbind(M(k-1,i),0)%*%(cbind(diag(x=1-D0(i,(i-k+2):i,k),k-1),0)+cbind(0,diag(x=D0(i,(i-k+2):i,k),k-1)))+
	rbind(0,M(k-1,i))%*%(cbind(diag(x= -D1(i,(i-k+2):i,k),k-1),0)+cbind(0,diag(x=D1(i,(i-k+2):i,k),k-1)))
} 
MFinal<-array(dim=c(k,k,length(knots)))
for(i in 1:length(knots)){
	MFinal[,,i]<-M(k,i)
}
class(MFinal)<-'SplineMatrices'

MFinal
}
Henkel<-function(x,nrow=length(x),ncol=length(x)){
	Z<-matrix(x[1:ncol],nrow=nrow,ncol=ncol,byrow=T)
	for(i in 1:(nrow-1)){
		Z[i+1,]<-x[-(1:i)][1:ncol]
	}
	Z
}

OBasis<-function(...)OrthogonalSplineBasis(...)
