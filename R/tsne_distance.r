
################################################
# This script contains functions for creating and
# handling distance matrices.
################################################

#
# The SNE symmetry operation for distance matrices
symmetrize_sne_matrix <- function(sne_matrix){
	nobs <- nrow(sne_matrix)
	sne_matrix_sym <- sne_matrix
	for (ii in 1:nobs){
		sne_matrix_sym[ii,] <- (sne_matrix[ii,] + sne_matrix[,ii]) / (2*nobs)
	}

	sym_matrix_sum <- sum(sne_matrix_sym)
	if (sym_matrix_sum > 1.001 | sym_matrix_sum < 0.999){
		stop()
	}

	return(sne_matrix_sym)
}


#' Calculate similarity matrices for a t-SNE analysis
#'
#' t-SNE uses a similarity measure that can be interpreted as the proability that
#' two points are nearest neighbors. These functions calculate these probabilities.
#'
#' \code{dist_sne} computes the SNE similarities using euclidean distances, as in the
#' originla t-SNE formulation.
#'
#' \code{dist_to_sne} converts an already computed dissimilarity matrix to a matrix
#' of SNE similarities.
#'
#' @param xdata A matrix or data frame.
#' @param dist_ A matrix of distances or \code{\link{dist}} object.
#' @param perplexity The perplexity. Usually a number between 5 and 50. Default is 25.
#'
#' @examples
#' data(iris)
#' iris_distances <- dist(iris[,1:4], method = 'manhattan')
#' iris_sne_distances <- dist_to_sne(iris_distances)
#'
#' @export
dist_sne <- function(xdata, perplexity=25){

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


	stopifnot(perplexity > 0)

	# Number of observations.
	nobs <- nrow(xdata)
	xdata <- as.matrix(xdata)

	# function for use with uniroot to find the sigma which yields the desired preplexity.
	perplex_error <- function(sigma, dist_vec, which_zero, perplex){
		# Calculate the nearest neighbor distribtuion.
		numerator_temp <- exp((-1 * dist_vec) / (2*sigma^2))
		numerator_temp[which_zero] <- 0
		pp_temp <- numerator_temp / sum(numerator_temp)

		cur_perpl <- 2^(sum(pp_temp * log2(pp_temp), na.rm=TRUE) * -1)
		perplexerr <- cur_perpl - perplex
		return(perplexerr)
	}

	# The SNE distance matrix
	PP <- matrix(0, ncol=nobs, nrow=nobs)

	for (jj in 1:nobs){

		eucl_dists <- colSums((xdata[jj,] - t(xdata))^2)
		uniroot_res <- tryCatch({stats::uniroot(perplex_error, interval=c(0.0001, 999),
																		 dist_vec=eucl_dists, which_zero=jj, perplex=perplexity)},
														error=function(e) return(0))
		if (is.numeric(uniroot_res)){
			browser()
			stop(sprintf('Failed to tune perplexity at observation %d. Try using a smaller perlexity.', jj))
		}

		sigmaa <- uniroot_res$root

		numerator_temp <- exp((-1 * eucl_dists) / (2*sigmaa^2))
		numerator_temp[jj] <- 0
		PP[jj,] <- numerator_temp / sum(numerator_temp)

	}

	PP <- symmetrize_sne_matrix(PP)

	return(PP)
}



#' @rdname dist_sne
#' @export
dist_to_sne <- function(dist_, perplexity=25){

	stopifnot(perplexity > 0)

	if ('dist' %in% class(dist_)){
		dist_ <- as.matrix(dist_)
	} else if ('matrix' %in% class(dist_)){
		# Should probably check for symmetry etc.
		dist_ <- dist_
	} else if ('data.frame' %in% class(dist_)){
		dist_ <- as.matrix(dist_)
	}

	# Number of observations.
	nobs <- nrow(dist_)

	# function for use with uniroot to find the sigma which yields the desired preplexity.
	perplex_error <- function(sigma, dist_vec, which_zero, perplex){
		# Calculate the nearest neighbor distribtuion.
		numerator_temp <- exp((-1 * dist_vec) / (2*sigma^2))
		numerator_temp[which_zero] <- 0
		pp_temp <- numerator_temp / sum(numerator_temp)

		cur_perpl <- 2^(sum(pp_temp * log2(pp_temp), na.rm=TRUE) * -1)
		perplexerr <- cur_perpl - perplex
		return(perplexerr)
	}

	# The SNE distance matrix
	PP <- matrix(0, ncol=nobs, nrow=nobs)

	for (jj in 1:nobs){

		uniroot_res <- tryCatch({stats::uniroot(perplex_error, interval=c(0.0001, 999),
													 dist_vec=dist_[jj,], which_zero=jj, perplex=perplexity)},
													 error=function(x) return(0))

		if (is.numeric(uniroot_res)){
			stop(sprintf('Failed to tune perplexity at observation %d. Try using a smaller perlexity.', jj))
		}

		sigmaa <- uniroot_res$root

		numerator_temp <- exp((-1 * dist_[jj,]) / (2*sigmaa^2))
		numerator_temp[jj] <- 0

		PP[jj,] <- numerator_temp / sum(numerator_temp)

	}

	PP <- symmetrize_sne_matrix(PP)

	return(PP)
}

