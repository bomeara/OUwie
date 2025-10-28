
	### Code by Priscilla Lau et al. they supplied to demonstrate the issue; slight tweaks by Brian O'Meara to handle zero alpha

	#library(pracma)
	#library(tidyverse)
	#library(tibble)
	

	parentNode <- function(tree, x) {
		m <- which(tree$edge[, 2] == x)
		return(tree$edge[m, 1])
	}

	# find nodes along a lineage towards root node by providing initial child (presumably tip) node
	nodesAlongLineage <- function(tree, old_node, young_node) {
		k <- young_node
		while (young_node != old_node) {
			k <- c(k, parentNode(tree, young_node))
			young_node <- tail(k, n = 1)
		}
		return(k)
	}

	# find subedges of a lineage
	lineage.constructor <- function(tree, root_node, e) {
		nodes <- nodesAlongLineage(tree, root_node, e)
		edges <- which(tree$edge[, 2] %in% nodes) # from root to tip
		subedge_lengths <- rev(unlist(lapply(edges, function(i) {
			tree$maps[[i]]
		}))) # tip to root

		state_changes <- names(subedge_lengths) # from tip to root

		return(tibble::tibble(state = state_changes, time_span = subedge_lengths))
	}

	weights.lineage <- function(tree, alpha, e) {
		print(paste("Calculating weights for edge:", e))
		root_node = length(tree$tip.label) + 1
		lineage <- lineage.constructor(tree, root_node, e)
		lineage[["alpha"]] = unlist(alpha[lineage[["state"]]])

		W = matrix(0, ncol = length(alpha), nrow = 1)
		colnames(W) = sort(names(alpha))

		if (length(lineage[[1]]) > 1) {
			lineage <- lineage |>
				dplyr::mutate(
					exp1 = -1 * expm1(-1 * alpha * time_span),
					sum2_temp = -1 * alpha * time_span
				)
			lineage$exp1[length(lineage$exp1)] = 1
			lineage$sum2 = 0

			for (i in 2:length(lineage[[1]])) {
				lineage$sum2[i] = lineage$sum2_temp[i - 1]
				lineage$sum2_temp[i] = lineage$sum2[i] + lineage$sum2_temp[i]
			}

			all_weights = lineage |>
				dplyr::mutate(exp_final = exp1 * exp(sum2)) |>
				dplyr::group_by(state) |>
				dplyr::summarise(weight = sum(exp_final))

			for (i in 1:nrow(all_weights)) {
				W[, all_weights$state[i]] = all_weights$weight[i]
			}
		} else {
			W[, lineage$state[1]] = 1
		}

		return(W)
	}

	# combine to form weight matrix
	weight.matrix.lau <- function(tree, alpha) {
		ntip = length(tree$tip.label)
		weight_matrix = matrix(0, nrow = ntip, ncol = length(alpha))
		rownames(weight_matrix) <- tree$tip.label
		colnames(weight_matrix) <- c(sort(names(alpha)))
		for (i in 1:ntip) {
			weight_matrix[i, ] <- weights.lineage(tree, alpha, i)
		}
		return(weight_matrix)
	}

	cov.accum <- function(tree, mrca_node, alpha, sigma2) {
		root_node = length(tree$tip.label) + 1
		if (mrca_node == root_node) {
			cov_accum = 0.0
		} else {
			nodes <- nodesAlongLineage(tree, root_node, mrca_node)
			edges <- which(tree$edge[, 2] %in% nodes) # from root to mcra_node
			subedge_lengths <- rev(unlist(lapply(edges, function(i) {
				tree$maps[[i]]
			}))) # from mcra_node to root

			subedge_lengths <- tibble::tibble(
				state = names(subedge_lengths),
				time_span = subedge_lengths,
				alpha = unlist(alpha[names(subedge_lengths)]),
				sigma2 = unlist(sigma2[names(subedge_lengths)])
			) |>
				dplyr::mutate(
					exp1 = -1 * expm1(-2 * alpha * time_span),
					sum2_temp = -2 * alpha * time_span
				)
			subedge_lengths$sum2 = 0

			if (length(subedge_lengths[[1]]) == 1) {
				subedge_lengths = subedge_lengths |>
					dplyr::mutate(cov = compute.piecewise.with.zero.possible(exp1, sigma2, alpha, time_span))
				cov_accum = subedge_lengths$cov[[1]]
			} else {
				for (i in 2:length(subedge_lengths[[1]])) {
					subedge_lengths$sum2[i] = subedge_lengths$sum2_temp[i - 1]
					subedge_lengths$sum2_temp[i] = subedge_lengths$sum2[i] +
						subedge_lengths$sum2_temp[i]
				}
				cov_accum = subedge_lengths |>
					dplyr::mutate(exp3 = exp1 * exp(sum2)) |>
					dplyr::group_by(state) |>
					dplyr::summarise(
						sum4 = sum(compute.piecewise.with.zero.possible(exp3, sigma2, alpha, time_span))
					) |>
					dplyr::reframe(sum_final = sum(sum4)) |>
					unlist() |>
					unname()
			}
		}
		return(cov_accum)
	}
	
	compute.piecewise.with.zero.possible <- function(exp_final, sigma2, alpha, time_span) {
		cov_value <- sigma2 / (2 * alpha) * exp_final
		baduns <- which(is.na(cov_value)) # occurs when alpha is zero, so we go to BM
		if(length(baduns) > 0){
			cov_value[baduns] <- sigma2[baduns] * time_span[baduns]
		}
		return(cov_value)
	}
		

	cov.loss <- function(tree, mrca_node, alpha, tip) {
		if (mrca_node == tip) {
			cov_loss_rate = 0
		} else {
			nodes <- nodesAlongLineage(tree, mrca_node, tip)
			nodes <- head(nodes, n = -1)
			edges <- which(tree$edge[, 2] %in% nodes) # from root to mcra_node
			subedge_lengths <- rev(unlist(lapply(edges, function(i) {
				tree$maps[[i]]
			}))) # from mcra_node to root
			subedge_lengths <- tibble::tibble(
				time_span = subedge_lengths,
				alpha = alpha[names(subedge_lengths)]
			)
			cov_loss_rate = subedge_lengths |>
				dplyr::mutate(sum1 = -1 * alpha * time_span) |>
				dplyr::reframe(sum_final = sum(sum1))
		}
		return(cov_loss_rate)
	}

	vcv.pairwise <- function(tree, alpha, sigma2, tip1, tip2) {
		mrca_node <- ape::mrca(tree)[tip1, tip2]

		cov_accum <- cov.accum(tree, mrca_node, alpha, sigma2)
		cov_loss1 <- cov.loss(tree, mrca_node, alpha, tip1)
		cov_loss2 <- cov.loss(tree, mrca_node, alpha, tip2)
		cov = cov_accum * exp(cov_loss1 + cov_loss2)
		return(unlist(unname(cov)))
	}

	vcv.matrix.lau <- function(tree, alpha, sigma2) {
		ntip <- length(tree$tip.label)
		V <- matrix(nrow = ntip, ncol = ntip)
		j = ntip
		while (j != 0) {
			for (i in 1:ntip) {
				V[i, j] <- vcv.pairwise(tree, alpha, sigma2, i, j)
				V[j, i] <- V[i, j]
			}
			j = j - 1
		}
		colnames(V) <- tree$tip.label
		rownames(V) <- tree$tip.label
		return(V)
	}

	sd_logL_vcv <- function(tree, continuousChar, alpha, sigma2, theta) {
		alpha = alpha[sort(names(alpha))]
		sigma2 = sigma2[sort(names(sigma2))]
		theta = theta[sort(names(theta))]
		states = names(theta)
		theta = matrix(theta, nrow = length(theta))
		rownames(theta) = states

		ntip <- length(tree$tip.label)
		V = vcv.matrix(tree, alpha, sigma2)

		W = weight.matrix(tree, alpha)

		C = chol(V) # upper triangular matrix
		L = t(C) # lower triangular matrix
		log_det_V = 0
		for (i in 1:ntip) {
			log_det_V = log_det_V + log(L[i, i])
		}
		log_det_V = log_det_V * 2.0

		y = NULL
		for (species in tree$tip.label) {
			y[species] = as.numeric(continuousChar[species])
		}

		# inverse of L
		r = solve(L) %*% y - solve(L) %*% W %*% theta

		res = 0.0
		res = res - (ntip / 2) * log(2 * pi)
		res = res - 0.5 * log_det_V
		res = res - 0.5 * pracma::dot(r, r)

		return(res)
	}

	# End of Lau et al. functions