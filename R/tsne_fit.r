


####################################################################
# Functions for using t-SNE.
####################################################################



# Internal function to validate control list
check_ctrl <- function(control_list){
	stopifnot(is.list(control_list))
	stopifnot(c('momentum_values', 'momentum_iter') %in% names(control_list))
	stopifnot(max(control_list$momentum_iter) >= control_list$max_iter)
	stopifnot(min(control_list$momentum_iter) > 1)
}



#' Dimension reduction with t-SNE
#'
#' t-Distributed Stochastic Neighbor Embedding (t-SNE).
#'
#' tsne_from_dist assumes that dist_ is a distance matrix interpretable as a
#' nearest neighbor distribution. If you only have a distance matrix, you can use
#' the \code{\link{dist_to_sne}} function to create a SNE dissimilarity matrix.
#'
#' Function \code{tsne} provides a convenient interface to use t-SNE using raw data.
#' Under the hood it calls the \code{\link{dist_sne}} and \code{\link{tsne_dist}} functions.
#'
#' @param xdata A matrix or data frame.
#' @param perplexity The perplexity. Usually a number between 5 and 50. Default is 25.
#' @param dist_ A matrix or \code{dist} object woth pairwise dissimilarities. All elements must be positive and sum to 1.
#' @param dimensions The number of dimensions of the t-SNE embedding. Default is 2, which is usually adequate.
#' @param rng_seed Provide a seed for generating intiial values for the coordiantes. The final t-SNE result is highly sensitive for starting values, so to get reproducible results you should provide a seed.
#' @param verbose The amount of output during the t-SNE estimation. Can be 0 (no output), 1 (some output, including a progress bar) and 2 (detailed output, mostly usefull for debugging purposes)
#' @param control A list of control parameters. See 'Details'.
#'
#'
#'
#' @section Control:
#' The control argument is a list that can supply any of the following components:
#'
#' \tabular{ll}{
#' \code{max_iter}: \tab The number of iterations. Default is 1000. \cr
#' \code{step_size}: \tab The gradient descent step size. Default is 100, but it can sometimes be usefull to set this to a lower value.\cr
#' \code{momentum_values}: \tab A vector of values that determines the momentum. The default is c(0.5, 0.8), following the original t-SNE paper.\cr
#' \code{momentum_iter}: \tab A vector of values that determines the iterations the momentum_values should be used. Default is c(250, 1000), following the original t-SNE paper. \cr
#' \code{early_exageration}: \tab Logical. Should the 'early exageration' trick be used? Default is TRUE, following the original t-SNE paper.\cr
#' }
#'
#'
#' @return
#'
#' A list of class 'tsne' with the following elements:
#'
#'\tabular{ll}{
#' \code{par}: \tab A matrix with the t-SNE embedding. \cr
#' \code{trace}: \tab A vector of the Kullback-Leibler divergences.
#' If you used early exaggeration (default) you will see a large jump in the
#' values at iteration 50. \cr
#'}
#'
#' @seealso
#'
#' The \code{\link{tsne_phyloseq}} function provides a nice interface for using t-SNE with microbiota data.
#'
#' @section References:
#'\itemize{
#' \item L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008.
#' \item L.J.P. van der Maaten. The t-SNE FAQ \url{https://lvdmaaten.github.io/tsne/}
#'}
#'
#' @examples
#' xx <- make_swiss_roll(150)
#' tsne_res <- tsne(xx[,1:3])
#' plot(tsne_res$par, col=xx[,4], pch=16)
#'
#' @export
tsne_dist <- function(dist_, dimensions=2, verbose=0, rng_seed=NULL, control=list()){

	# This already assumes that you have calculated a distance matrix from the data.

	if ('dist' %in% class(dist_)){
		# Distance matrix. Convert to regular matrix.
		pp_mat <- as.matrix(dist_)
	} else if ('matrix' %in% class(dist_)){
		pp_mat <- dist_
	}

	stopifnot(all(pp_mat >= 0))
	# make sure the dist_is a sne matrix.
	ok_pp <- sum(pp_mat) < 1.0001 & sum(pp_mat) > 0.999
	if (!ok_pp){
		stop('Input dist_ not a proper t-SNE dissimilarity matrix. sum(as.matrix(dist_)) should equal 1.')
	}

	stopifnot(dimensions < ncol(pp_mat))

	default_control <- list(max_iter=1000,
													step_size=100,
													momentum_values = c(0.5, 0.8),
													momentum_iter = c(250, 1000),
													early_exageration=TRUE)
	control <- utils::modifyList(default_control, control)
	check_ctrl(control)


	if (verbose > 0){
		cat('Gradient descent options:\n')
		cat(sprintf('Step size = %f\n', control$step_size))
	}

	nobs <- ncol(pp_mat)

	if (control$early_exageration){
		early_exag_factor <- 4
		early_exag_iter <- 50
	}

	if (verbose > 0){
		cat(sprintf('Momentum: %.2f (iter <= %d)\n', control$momentum_values, control$momentum_iter))
	}

	pp_mat_raw <- pp_mat

	KL_trace <- numeric(control$max_iter)

	# Initial parameter estimates.
	par_old <- rep(0, dimensions*nobs) # for Gradient descent momentum

	if (!is.null(rng_seed)){
		set.seed(rng_seed)
	}

	par_init <- stats::rnorm(dimensions*nobs, mean = 0, sd=0.001) # Current parameters.

	# Reset random seed.
	set.seed(seed=NULL)

	par_cur <- par_init

	cur_momentum <- 1

	if (verbose == 1){
		tpb <- utils::txtProgressBar(min=1, max=control$max_iter, initial=1)
	}

	for (kk in 1:control$max_iter){

		if (verbose == 1){
			utils::setTxtProgressBar(tpb, kk)
		}

		# Early exaggeration
		if (control$early_exageration){
			if (kk <= early_exag_iter){
				pp_mat <- pp_mat_raw * early_exag_factor
			} else if (kk == early_exag_iter+1){
				pp_mat <- pp_mat_raw
			}
		} else {
			pp_mat <- pp_mat_raw
		}


		# Momentum
		if (kk <= control$momentum_iter[cur_momentum]){
			momentum <- control$momentum_values[cur_momentum]
			if (kk == control$momentum_iter[cur_momentum]){
				# Update the momentum.
				cur_momentum <- cur_momentum + 1
			}
		}

		par_cur_mat <- matrix(par_cur, ncol=dimensions)

		## # calculate t SNE dissimilarities for yy
		QQ_num <- matrix(0, ncol=nobs, nrow=nobs)  		# Raw t dissimilarities
		QQ <- matrix(0, ncol=nobs, nrow=nobs) 				# Normalized t-dissimilarities
		for (ii in 1:nobs){
			numerators <- 1 / (1 + (colSums((par_cur_mat[ii,] - t(par_cur_mat))^2)))
			numerators[ii] <- 0 # set dissimialrities for i = j to 0.
			QQ_num[ii,] <- numerators
			QQ[ii,] <- numerators / sum(numerators)
		}

		# Symmetrize
		QQ_sym <- symmetrize_sne_matrix(QQ)
		diag(QQ_sym) <- 1

		# KL divergence
		KL_trace[kk] <- sum(pp_mat * log(pp_mat / QQ_sym), na.rm=TRUE)

		if (verbose >= 2){
			cat(sprintf('[%d / %d] Function value: %f\n', kk, control$max_iter, KL_trace[kk]))
		}

		# Calculate the gradient. Addapted from the tsne_p.m matlab script.
		LL <- (pp_mat - QQ_sym)* QQ_num
		gradient_ <- 4 * ((diag(colSums(LL)) - LL) %*% par_cur_mat)

		# Update parameters
		par_new <- par_cur - (control$step_size * gradient_) + (momentum * (par_cur - par_old))

		par_old <- par_cur
		par_cur <- par_new
	}

	out <- list(par=par_cur, trace=KL_trace)
	class(out) <- 'tsne'

	return(out)
}




#' @rdname tsne_dist
#' @export
tsne <- function(xdata, perplexity=25, dimensions=2, verbose=0, rng_seed=NULL, control=list()){

	if ('data.frame' %in% class(xdata)){
		all_numeric <- all(sapply(xdata, is.numeric))
		if (!all_numeric) {
			stop('Not all variables in xdata is numeric.')
		}
	} else if ('matrix' %in% class(xdata)){
		matrix_numeric <- is.numeric(xdata)
		if (!matrix_numeric) {
			stop('xdata is not numeric.')
		}
	}

	dist_ <- dist_sne(xdata, perplexity = perplexity)
	res <- tsne_dist(dist_=dist_, dimensions=dimensions,
													 control=control, rng_seed=rng_seed,
													 verbose=verbose)

	return(res)
}


#' Swiss roll data set
#'
#' Simulate swiss roll data. The swiss roll data set is used to test various dimension reduction algorithms.

#'
#' @param nn The number of data points to be generated.
#'
#'
#' @return A data frame with four variables. The three first are the swiss
#' roll data points, the fourth is vector fo colors useful for plotting.
#'
#' @section Reference:
#' Dinoj Surendran, Swiss Roll Dataset
#' \url{http://people.cs.uchicago.edu/~dinoj/manifold/swissroll.html}
#'
#' @examples
#' xx <- make_swiss_roll(150)
#' pairs(xx[,1:3], col=xx[,4], pch=16)
#'
#' @export
make_swiss_roll <- function(nn){

	means1 <- c(7.5, 7.5, 12.5, 12.5)
	means2 <- c(7.5, 12.5, 7.5, 12.5)

	xx1 <- stats::rnorm(nn, mean=means1)
	xx2 <- stats::rnorm(nn, mean=means2)

	xx <- data.frame(x1=xx1 * cos(xx1), x2=xx2, x3=xx1 * sin(xx1))
	# make colors according to x1.
	mycol <- grDevices::colorRamp(c('red', 'cadetblue', 'black'))
	xx1_tmp <- xx[,1] + abs(min(xx[,1]))
	xx1_tmp <- xx1_tmp / max(xx1_tmp)
	xcolor <- mycol(xx1_tmp)

	xx$plot_colors <- grDevices::rgb(xcolor[,1]/255, xcolor[,2]/255, xcolor[,3]/255)

	return(xx)
}


