#' @importFrom parallel mcmapply  
intervaler <- function(x,w){
	if(length(x) < 2) stop("size of vector inferior to 2")
	if(length(x) != length(w) ) stop("vectors x and w of different size")
	i <- 1:(length(x)-1)
	res <- lapply(i, function(i)c(x[i], x[i+1], w[i], w[i+1]))
	return(res)
}

borner_inf <- function(vect)
	vect[3] <= 0 & vect[4] >= 0
borner_sup <- function(vect)
	vect[3] >= 0 & vect[4] <= 0


binf_time <- function(uninterval ,fun ,nbdiv ,xmax ,fmax ,nbiter=0, p_res=NULL){
	print("**************************************************")
	print(nbiter)
	print("**************************************************")
	max_iter <- 20
	z <- seq(uninterval[1], uninterval[2], length.out=nbdiv)
	pas <- (max(z) - min(z)) / nbdiv
	if(!is.null(p_res)){
		z_p <- p_res$x
		w_p <- p_res$likel
		z <- z[!(z %in% z_p)]
	}
	w <- unlist(mclapply(z,fun))
	if(nbiter==0){
		iz <- which(z == xmax)
		print("checking")
		print(w[iz])
		w[iz] <- fmax
	}
	if(!is.null(p_res)){
		z=c(z,z_p)
		w <- c(w,w_p)
		w <- w[order(z)]
		z <- z[order(z)]
	}
	res <- list("x"=z,"likel"=w,"pas"=pas)
	if(nbiter >= max_iter)
		return(list("x"=z,"likel"=w,"pas"=pas))
	if(all(w>0)){
		print("all value are positive")
		int_r <- uninterval[2]-uninterval[1]	
		extended_interval <- c(uninterval[1]-int_r,uninterval[1])
		res <- binf_time(extended_interval, fun, nbdiv=nbdiv, xmax=xmax, nbiter=nbiter+1,p_res=res)
		return(res)
	}
	if(all(w<0)){
		stop("all values are negative")
	}
	intervalles <- intervaler(z,w)
	bornes_inf <- Filter(borner_inf,intervalles)
	if(length(bornes_inf) == 0)
		return(list("x"=z,"likel"=w))
	xinf <- min(unlist(lapply(bornes_inf,"[",1)))
	min_binf <- which(lapply(bornes_inf,"[",1)==xinf)
	min_binf <- bornes_inf[[min_binf]]
	ic_inf <- min_binf[1]
	condition <- (xmax-ic_inf)>20*pas
	if(!condition){
		res <- binf_time(min_binf[1:2], fun, nbdiv=2*nbdiv, xmax=xmax, nbiter=nbiter+1,p_res=res)
		return(res)
	}else{
		return(list("x"=z,"likel"=w,"pas"=pas))
	}
}

bsup_time <- function(uninterval,fun,nbdiv,xmax,fmax,nbiter=0,p_res=NULL){
	print("**************************************************")
	print(nbiter)
	print("**************************************************")
	max_iter <- 20
	z <- seq(uninterval[1],uninterval[2],length.out=nbdiv)
	pas=(max(z)-min(z))/nbdiv
	if(!is.null(p_res)){
		z_p <- p_res$x
		w_p <- p_res$likel
		z <- z[!(z%in%z_p)]
	}
	w <- unlist(mclapply(z,fun))
	if(nbiter==0){
		iz <- which(z == xmax)
		print("checking")
		print(w[iz])
		w[iz] <- fmax
	}
	if(!is.null(p_res)){
		z <- c(z,z_p)
		w <- c(w,w_p)
		w <- w[order(z)]
		z <- z[order(z)]
	}
	res <- list("x"=z,"likel"=w,"pas"=pas)
	if(nbiter >= max_iter)
		return(list("x"=z,"likel"=w,"pas"=pas))
	if(all(w>0)){
		print("all value are positive")
		int_r <- uninterval[2]-uninterval[1]	
		extended_interval <- c(uninterval[2],int_r+uninterval[2])
		res <- bsup_time(extended_interval, fun, nbdiv=nbdiv, xmax=xmax, nbiter=nbiter+1,p_res=res)
		return(res)
	}
	if(all(w<0)){
		stop("all values are negative")
	}
	intervalles <- intervaler(z,w)
	bornes_sup <- Filter(borner_sup,intervalles)
	if(length(bornes_sup) == 0)
		return(list("x"=z,"likel"=w))
	xsup <- max(unlist(lapply(bornes_sup,"[",1)))
	max_bsup <- which(lapply(bornes_sup,"[",1)==xsup)
	max_bsup <- bornes_sup[[max_bsup]]
	ic_sup <- max_bsup[2]
	condition <- (ic_sup-xmax)>20*pas
	if(!condition){
		res <- bsup_time(max_bsup[1:2], fun, nbdiv=2*nbdiv, xmax=xmax, nbiter=nbiter+1,p_res=res)
		return(res)
	}else{
		return(list("x"=z,"likel"=w,"pas"=pas))
	}
}



get_ic <- function(int_inf,int_sup){
	c("ic_inf"=int_inf[1], "ic_sup"=int_sup[2])
}
	

