#' Function to convert a tree with node labels and tip labels to a simmap-style tree
#' #' with mapped edges.
#' @param tree a phylogenetic tree, in ape phylo format and with internal nodes labeled denoting the ancestral selective regimes.
#' @param data a data.frame containing species information.
#' @param shift.point numeric value between 0 and 1 indicating the relative position along branches where regime shifts occur (default is 0.5, the midpoint).
#' @return A phylo object with mapped edges.
node_label_tree_to_simmap_tree <- function(tree, data, shift.point = 0.5) {
	# get the regimes from the data and node labels
	regimes <- c(data[, 2], tree$node.label)

	# create a simmap-style tree
	simmap_tree <- tree
	# create mapped edges
	simmap_tree$mapped.edge <- matrix(0, nrow = nrow(simmap_tree$edge), ncol = length(unique(regimes)))
	colnames(simmap_tree$mapped.edge) <- sort(unique(regimes))
	simmap_tree$maps <- vector("list", nrow(simmap_tree$edge))

	# fill in the mapped edges
	for (i in sequence(nrow(simmap_tree$edge))) {
		tipward_node <- simmap_tree$edge[i, 2]
		rootward_node <- simmap_tree$edge[i, 1]
		brlen <- simmap_tree$edge.length[i]
		tipward_state <- regimes[tipward_node]
		rootward_state <- regimes[rootward_node]
		if (tipward_state == rootward_state) {
			# no shift along this edge
			simmap_tree$mapped.edge[i, tipward_state] <- brlen
			simmap_tree$maps[[i]] <- c(tipward_state = brlen)
			names(simmap_tree$maps[[i]]) <- as.character(tipward_state)
		} else {
			# there is a shift along this edge
			shift_time <- brlen * shift.point
			simmap_tree$mapped.edge[i, tipward_state] <- brlen - shift_time
			simmap_tree$mapped.edge[i, rootward_state] <- shift_time
			simmap_tree$maps[[i]] <- c(
				rootward_state = shift_time,
				tipward_state = brlen - shift_time
			)
			names(simmap_tree$maps[[i]]) <- c(as.character(rootward_state), as.character(tipward_state))
		}
	}
	class(simmap_tree) <- c("simmap", "phylo")
	attr(simmap_tree, "map.order") <- "right-to-left"
	return(simmap_tree)
}