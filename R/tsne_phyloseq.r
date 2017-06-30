
###############################################################################
# This script contains functions for using t-SNE with
# microbiota data stored in phyloseq obejcts.
###############################################################################


#' Dimension reduction of microbiota data with t-SNE
#'
#' t-Distributed Stochastic Neighbor Embedding (t-SNE) for microbiota data.
#'
#' The \code{philr_options} argument can conatin two elements, part.weights and ilr.weights.
#' The default is 'enorm.x.gm.counts' and 'blw.sqrt', respectively, since these are
#' used in the philr vignette. See \code{\link[philr]{philr}} for more details.
#'
#' @param physeq  Phylogenetic sequencing data (\code{\link[phyloseq]{phyloseq-class}}).
#' @param distance The dissimilarity measure to use. Either one of \code{bray} (default) \code{unifrac}, \code{wunifrac}, \code{jsd}, \code{dpcoa} or \code{philr}.
#' @param perplexity The perplexity. Usually a number between 5 and 50. Default is 25.
#' @param dimensions The number of dimensions of the t-SNE embedding. Default is 2, which is usually adequate.
#' @param precomputed_distance Matrix or dist object of a precomputed dissimilarity matrix. If you use other ordination techniques or trying different perplexities it will be faster to compute the dissimilarities beforehand.
#' @param pseudocounts The number of psuedocounts, which is added to all elements in the otu table. Only applied when distance is philr.
#' @param verbose The amount of output during the t-SNE estimation. Can be 0 (no output), 1 (some output, including a progress bar) and 2 (detailed output, mostly usefull for debugging purposes)
#' @param rng_seed Provide a seed for generating intiial values for the coordiantes. The final t-SNE result is highly sensitive for starting values, so to get reproducible results you should provide a seed.
#' @param philr_options A list of options to use with philr distances. See 'Details'.
#' @param control A list of control parameters. See the 'Details' section for \code{\link{tsne}}.
#'
#'
#' @return
#'
#' A list of class tsne_phyloseq.
#'
#' \tabular{ll}{
#'
#' \code{distance_matrix} \tab The distance matrix computed using the method
#' specified in the \code{distance} argument. \cr
#'
#' \code{sne_distance_matrix} \tab The SNE dissimilarity matrix computed from the
#' distance_matrix with the perplexity given via the perplexity argument. \cr
#'
#' \code{names} \tab Samples names from the phyloseq object. \cr
#'
#' \code{tsne} \tab A list of class \code{tsne} as returned from the \code{\link{tsne}} function. Contains the t-SNE layout and some fit diagnostics,  \cr
#'
#' }
#'
#'
#' @section References:
#'\itemize{
#' \item L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605, 2008.
#' \item L.J.P. van der Maaten. The t-SNE FAQ \url{https://lvdmaaten.github.io/tsne/}
#' \item Paul J. McMurdie & Susan Holmes phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data \url{https://doi.org/10.1371/journal.pone.0061217}
#'}
#'
#' @seealso
#'
#' \code{\link{plot_tsne_phyloseq}}, \code{\link{add_tsne_to_sample_data}},
#' \code{\link{tsne}}.
#'
#' Other dimension reduction methods for microbiota data is available
#' via the \code{\link[phyloseq]{ordinate}} function in the phyloseq package.
#'
#' @examples
#' library(phyloseq)
#' data("GlobalPatterns")
#' tsne_res <- tsne_phyloseq(GlobalPatterns, distance='wunifrac',
#'    perplexity = 8, verbose=0, rng_seed = 3901)
#' plot_tsne_phyloseq(GlobalPatterns, tsne_res,
#'    color = 'SampleType', title='t-SNE (Weighted UniFrac)')
#'
#' @export
tsne_phyloseq <- function(physeq, distance='bray', perplexity=25,
													dimensions=2, precomputed_distance=NULL,
													pseudocounts=1, verbose=1,
													rng_seed=NULL,
													philr_options=list(), control=list()){

	stopifnot(dimensions >= 2)
	stopifnot(pseudocounts >= 0)
	stopifnot(perplexity >= 1)

	if (!requireNamespace("phyloseq", quietly = TRUE)) {
		stop("phyloseq needed for this function to work. Please install it.",
				 call. = FALSE)
	}

	if (is.null(physeq@phy_tree) & distance %in% c('unifrac', 'wunifrac', 'philr')){
		stop(sprintf("Phylogenetic tree not available in physeq obejct. distance = %s does not work.", distance),
				 call. = FALSE)
	}

	if (verbose >= 1){
		cat(sprintf('t-SNE for microbiota \n', distance))

		if (is.null(precomputed_distance)){
			cat(sprintf('Calculating distance matrix using %s \n', distance))
		}

	}

	if (is.null(precomputed_distance)){

		if (distance == 'philr'){
			if (!requireNamespace("philr", quietly = TRUE)) {
				stop("You need the philr package ot use the philr distance. Please install it.",
						 call. = FALSE)
			}

			# Philr options
			default_philr_options <- list(part.weights='enorm.x.gm.counts',
																		ilr.weights='blw.sqrt')
			philr_options <- utils::modifyList(default_philr_options, philr_options)
			check_philr_options(philr_options)

			# Some philr preprocessing, adapted from the philr vignette.
			phyloseq::taxa_names(physeq) <- paste('Otu', phyloseq::taxa_names(physeq), sep='')
			phyloseq::phy_tree(physeq) <- ape::makeNodeLabel(phyloseq::phy_tree(physeq), method="number", prefix='n')
			physeq_transformed <- phyloseq::transform_sample_counts(physeq, function(x) x+pseudocounts)

			otu_tab_tmp <- t(phyloseq::otu_table(physeq_transformed)@.Data)

			philr_res <- philr::philr(df=otu_tab_tmp, tree=phyloseq::phy_tree(physeq_transformed),
																part.weights=philr_options$part.weights,
																ilr.weights=philr_options$ilr.weights)

			distance_matrix <- stats::dist(philr_res, method='euclidean')
		} else if (distance %in% c('unifrac', 'wunifrac', 'bray', 'jsd', 'dpcoa')) {
			distance_matrix <- phyloseq::distance(physeq, method=distance)
		} else {
			#phyloseq::distanceMethodList$vegdist
			stop('Unsupported distance.')
		}

	} else {

		if ('dist' %in% class(precomputed_distance)){
			# Distance matrix. Convert to regular matrix.
			distance_matrix <- as.matrix(precomputed_distance)
		} else if ('matrix' %in% class(precomputed_distance)){
			distance_matrix <- precomputed_distance
		}

		stopifnot(ncol(distance_matrix) == nrow(distance_matrix))
		stopifnot(ncol(distance_matrix) == length(physeq@sam_data@row.names))

		cat(sprintf('Precomputed distance matrix used. Distance=%s will be ignored\n', distance))

	}

	# Make SNE conditional distance matrix
	sne_dist <- dist_to_sne(distance_matrix, perplexity = perplexity)


	sne_res <- tsne_dist(sne_dist, dimensions=dimensions,
											 control=control,
											 verbose=verbose, rng_seed=rng_seed)


	out <- list(distance_matrix=distance_matrix,
							sne_distance_matrix=sne_dist,
							names=physeq@sam_data@row.names,
							tsne=sne_res)

	class(out) <- 'tsne_phyloseq'
	return(out)
}


# Function that validates the philr control list.
check_philr_options <- function(philropt){
	names_in_list <- names(philropt)
	ok_names <-  c('part.weights', 'ilr.weights')
	names_valid <- names_in_list %in% ok_names
	if (any(!names_valid)){
		warning(sprintf('philr option %s is not applicable and is ignored.', names_in_list[!names_valid]))
	}
}



#' Plot microbiota t-SNE results using ggplot2
#'
#' This function provides an interface for combining a t-SNE ordination
#' with sample meta data from a phyloseq object.
#'
#' This function is fashioned after the \code{\link[phyloseq]{plot_ordination}}
#' function in the phyloseq package.
#'
#' @param physeq  Phylogenetic sequencing data (\code{\link[phyloseq]{phyloseq-class}}).
#' @param tsne_obj An list of class \code{tsne_phyloseq}.
#' @param axes The axes to plot. Must be a length 2 numeric.
#' @param color Character string indicating which variable in sample_data should be used as the ggplot2 color aesthetic
#' @param shape Character string indicating which variable in sample_data should be used as the ggplot2 shape aesthetic
#' @param title Plot title.
#' @param justDF If FALSE (default) the plot will be displayed and a ggplot2
#' object will be returned. If TRUE a data.frame useful will be returned instead.
#'
#' @seealso
#' \code{\link{tsne_phyloseq}}. \code{\link{add_tsne_to_sample_data}}.
#'
#' @examples
#' library(phyloseq)
#' data("GlobalPatterns")
#' tsne_res <- tsne_phyloseq(GlobalPatterns, distance='wunifrac',
#'   perplexity = 8, verbose=0, rng_seed = 3901)
#' plot_tsne_phyloseq(GlobalPatterns, tsne_res,
#'   color = 'SampleType', title='t-SNE (Weighted UniFrac)')
#'
#'
#' @export
plot_tsne_phyloseq <- function(physeq, tsne_obj, axes = 1:2, color = NULL,
															 shape = NULL, title = NULL,
															 justDF = FALSE){

	stopifnot('tsne_phyloseq' %in% class(tsne_obj))
	stopifnot('phyloseq' %in% class(physeq))
	stopifnot(is.numeric(axes))
	stopifnot(max(axes) <= ncol(tsne_obj$tsne$par))
	stopifnot(min(axes) >= 1)
	stopifnot(length(phyloseq::sample_names(physeq)) == nrow(tsne_obj$tsne$par))

	# Variable names.
	ok_names <- phyloseq::sample_data(physeq)@names

	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		stop("ggplot2 needed for this function to work. Please install it.",
				 call. = FALSE)
	}

	plot_df <- data.frame(tsne_obj$tsne$par)

	if (!is.null(color)){
		stopifnot(color %in% ok_names)
		plot_df[[color]] <- phyloseq::sample_data(physeq)[[color]]
	}

	if (!is.null(shape)){
		stopifnot(shape %in% ok_names)
		plot_df[[shape]] <- phyloseq::sample_data(physeq)[[shape]]
	}

	if (justDF){
		return(plot_df)
	} else {
		gg_obj <- ggplot2::ggplot(plot_df)
		gg_obj <- gg_obj + ggplot2::aes(x=plot_df[[axes[1]]], y=plot_df[[axes[2]]])
		gg_obj <- gg_obj + ggplot2::aes_string(color=color, shape=shape)
		gg_obj <- gg_obj + ggplot2::geom_point()
		gg_obj <- gg_obj + ggplot2::xlab(colnames(plot_df)[axes[1]])
		gg_obj <- gg_obj + ggplot2::ylab(colnames(plot_df)[axes[2]])
		gg_obj <- gg_obj + ggplot2::ggtitle(title)

		return(gg_obj)
	}

}




#' Add the results from t-SNE to the phyloseq sample data
#'
#' @param physeq  Phylogenetic sequencing data (\code{\link[phyloseq]{phyloseq-class}}).
#' @param tsne_obj An list of class \code{tsne_phyloseq}.
#' @param prefix Character string for making new column names.
#'
#' @return
#' A phyloseq object where the t-SNE coordiantes are added to the sample data frame.
#'
#' @export
add_tsne_to_sample_data <- function(physeq, tsne_obj, prefix='tsne'){

	stopifnot('tsne_phyloseq' %in% class(tsne_obj))
	stopifnot('phyloseq' %in% class(physeq))
	stopifnot(length(phyloseq::sample_names(physeq)) == nrow(tsne_obj$tsne$par))
	stopifnot(length(prefix) == 1)

	dimensions <- ncol(tsne_obj$tsne$par)
	column_names <- paste(prefix, 1:dimensions, sep='_')

	if (all(column_names %in% physeq@sam_data@names)){
		stop('Columnames alrdeady in sample_data. Try using a different prefix.')
	}

	if (all(tsne_obj$names != physeq@sam_data@names)){
		stop('The names in the tsne obejct does not match those in the physeq object.')
	}


	physeq@sam_data@names <- c(physeq@sam_data@names, column_names)

	ncol_sample_data <- length(physeq@sam_data@names)
	for (ii in 1:dimensions){
		physeq@sam_data@.Data[[ncol_sample_data+ii]] <- tsne_obj$tsne$par[,ii]
	}

	return(physeq)

}





