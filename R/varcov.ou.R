varcov.ou <- function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, assume.station=TRUE, shift.point=.5, corrected=TRUE, enhanced_tree=NULL) {
	if(corrected & simmap.tree) {
		if(is.null(enhanced_tree)) {
			enhanced_tree <- create_enhanced_tree_structure(phy)
		}
		return(varcov.ou.enhanced.tree (phy, enhanced_tree, Rate.mat, root.state, root.age, scaleHeight))
	} else {
		return(varcov.ou.original(phy, edges, Rate.mat, root.state, simmap.tree, root.age, scaleHeight, assume.station, shift.point))
	}
}



#OU variance-covariance matrix generator

#written by Jeremy M. Beaulieu



varcov.ou.original <- function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, assume.station=TRUE, shift.point=.5){
    
    if(assume.station == TRUE){
        alpha=Rate.mat[1,1]
        sigma=Rate.mat[2,1]
        vcv <- quickVCV(phy=phy, alpha=alpha, sigma.sq=sigma, scaleHeight=scaleHeight)
    }else{
        if(is.null(root.state)) {
            root.state <- which(edges[dim(edges)[1],]==1)-5
            edges <- edges[-1*dim(edges)[1],]
        }
        n=max(phy$edge[,1])
        ntips=length(phy$tip.label)
        if(simmap.tree==TRUE){
            k=length(colnames(phy$mapped.edge))
        }
        if(simmap.tree==FALSE){
            mm<-dim(edges)
            k<-length(6:mm[2])
        }
        
        pp <- prop.part(phy)
        oldregime=root.state
        nodevar1=rep(0,max(edges[,3]))
        nodevar2=rep(0,max(edges[,3]))
        alpha=Rate.mat[1,]
        sigma=Rate.mat[2,]
        n.cov1=matrix(rep(0,n), n, 1)
        n.cov2=matrix(rep(0,n), n, 1)
        
        if(simmap.tree==TRUE){
            regimeindex<-colnames(phy$mapped.edge)
            for(i in 1:length(edges[,1])){
                anc = edges[i, 2]
                desc = edges[i, 3]
                
                if(scaleHeight==TRUE){
                    currentmap<-phy$maps[[i]]/max(MakeAgeTable(phy, root.age=root.age))
                }
                else{
                    currentmap <- phy$maps[[i]]
                }
                oldtime=edges[i,4]
                for (regimeindex in 1:length(currentmap)){
                    regimeduration <- currentmap[regimeindex]
                    newtime <- oldtime+regimeduration
                    regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
                    nodevar1[i] <- nodevar1[i]+alpha[regimenumber]*(newtime-oldtime)
                    nodevar2[i] <- nodevar2[i]+sigma[regimenumber]*((exp(2*alpha[regimenumber]*newtime)-exp(2*alpha[regimenumber]*oldtime))/(2*alpha[regimenumber]))
                    oldtime <- newtime
                    newregime <- regimenumber
                }
                oldregime=newregime
                n.cov1[edges[i,3],]=nodevar1[i]
                n.cov2[edges[i,3],]=nodevar2[i]
            }
        }
        if(simmap.tree==FALSE){
            for(i in 1:length(edges[,1])){
                anc <- edges[i,2]
                oldtime <- edges[i,4]
                newtime <- edges[i,5]
                if(anc%in%edges[,3]){
                    start <- which(edges[,3]==anc)
                    oldregime <- which(edges[start,6:(k+5)]==1)
                }
                else{
                    #For the root:
                    oldregime=root.state
                }
                newregime=which(edges[i,6:(k+5)]==1)
                if(oldregime==newregime){
                    nodevar1[i] <- alpha[oldregime]*(newtime-oldtime)
                    nodevar2[i] <- sigma[oldregime]*((exp(2*alpha[oldregime]*newtime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
                }
                else{
                    shifttime <- newtime-((newtime-oldtime)*shift.point)
                    epoch1a <- alpha[oldregime]*(shifttime-oldtime)
                    epoch1b <- sigma[oldregime]*((exp(2*alpha[oldregime]*shifttime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
                    oldtime <- shifttime
                    newtime <- newtime
                    epoch2a <- alpha[newregime]*(newtime-oldtime)
                    epoch2b <- sigma[newregime]*((exp(2*alpha[newregime]*newtime)-exp(2*alpha[newregime]*oldtime))/(2*alpha[newregime]))
                    nodevar1[i] <- epoch1a+epoch2a
                    nodevar2[i] <- epoch1b+epoch2b
                }
                oldregime <- newregime
                n.cov1[edges[i,3],] <- nodevar1[i]
                n.cov2[edges[i,3],] <- nodevar2[i]
            }
        }
        vcv1 <- mat.gen(phy,n.cov1,pp)
        vcv2 <- mat.gen(phy,n.cov2,pp)
        if(any(abs(diff(alpha)) > 0)){
            species.variances <- diag(vcv1)
            species.total.variances <- matrix(0, dim(vcv1)[2], dim(vcv1)[2])
            count=0
            for(i in 1:dim(vcv1)[2]) {
                for(j in 1:dim(vcv1)[2]){
                    species.total.variances[i,j] <- exp(-(species.variances[i] + species.variances[j]))
                    count=count+1 # the count is always watching
                }
            }
            vcv <- species.total.variances * vcv2
        }else{
            if(is.null(root.age)){
                root.age <- max(branching.times(phy))
            }
            vcv <- exp(-2*alpha[1]*max(root.age)) * vcv2
        }
    }
    vcv
}

# to allow easier traversal and storage of tree structure
# assumes a simmap tree
create_enhanced_tree_structure <- function(phy) {
	new_structure <- data.frame(matrix(nrow=0, ncol=13))
	colnames(new_structure) <- c("edge_num", "rootward_treenode", "tipward_treenode", "tip_label", "rootward_mapnode", "tipward_mapnode", "segment_length", "parent_row", "regime", "tip_number", "rootward_height", "tipward_height", "segment_index")
	#phy <- ape::reorder.phylo(phy, "postorder") # which goes from root to tips if we start at the bottom and work up. But reordering in ape does not reorder the maps!!!
	heights <- phytools::nodeHeights(phy)
	for(i in rev(sequence(length(phy$edge[,1])))) {
		edge_num <- i
		rootward_treenode <- phy$edge[i, 1]
		tipward_treenode <- phy$edge[i, 2]
		rootward_height_treenode <- heights[i, 1]
		rootward_height_segment <- rootward_height_treenode
		tip_label <- ifelse(tipward_treenode <= length(phy$tip.label), phy$tip.label[tipward_treenode], NA)
		for (segment_index in sequence(length(phy$maps[[i]]))) {
			# for each regime segment on this edge, we need a row
			rootward_mapnode <- ifelse(segment_index == 1, 0, segment_index - 1)
			tipward_mapnode <- segment_index
			segment_length <- phy$maps[[i]][segment_index]
			tipward_height_segment <- rootward_height_segment + segment_length
			parent_row <- NA
			if (segment_index > 1) {
				parent_row <- nrow(new_structure)
			} else {
				# find the parent row
				# if (rootward_treenode %in% phy$edge[, 2]) {
				# 	parent_edge <- which(phy$edge[, 2] == rootward_treenode)
				# 	parent_row <- max(which(
				# 		new_structure$edge_num == parent_edge &
				# 			new_structure$tipward_mapnode ==
				# 				max(new_structure$tipward_mapnode[
				# 					new_structure$edge_num == parent_edge
				# 				])
				# 	))
				# } else {
				# 	parent_row <- NA
				# }
			}
			tip_label_local <- ifelse(segment_index == length(phy$maps[[i]]), tip_label, NA) # only the last segment gets the tip label
			tip_number <- ifelse(
				segment_index == length(phy$maps[[i]]) &
					tipward_treenode <= length(phy$tip.label),
				tipward_treenode,
				NA
			) 
			regime <- names(phy$maps[[i]])[segment_index]
			new_structure <- rbind(new_structure, data.frame(edge_num, rootward_treenode, tipward_treenode, tip_label, rootward_mapnode, tipward_mapnode, segment_length, parent_row, regime, tip_number, rootward_height_segment, tipward_height_segment, segment_index))
			rootward_height_segment <- tipward_height_segment
		}
	}
	for (i in sequence(nrow(new_structure))) {
		if (new_structure$segment_index[i] == 1) {
			if (new_structure$rootward_treenode[i] %in% phy$edge[, 2]) {
				parent_edge <- which(
					phy$edge[, 2] == new_structure$rootward_treenode[i]
				)
				parent_row <- max(which(
					new_structure$edge_num == parent_edge &
						new_structure$tipward_mapnode ==
							max(new_structure$tipward_mapnode[
								new_structure$edge_num == parent_edge
							])
				))
				new_structure$parent_row[i] <- parent_row
			}
		}
	}
	rownames(new_structure) <- NULL
	return(new_structure)

}

# Added by Brian O'Meara incorporating insights provided by Priscilla Lau. 

varcov.ou.enhanced.tree <- function(phy, enhanced_tree, Rate.mat, root.state, root.age=NULL, scaleHeight=FALSE){
	ntax <- max(enhanced_tree$tip_number, na.rm=TRUE)
	enhanced_tree$alpha <- Rate.mat[1,as.numeric(enhanced_tree$regime)]
	enhanced_tree$sigma_squared <- Rate.mat[2,as.numeric(enhanced_tree$regime)]
	enhanced_tree$alpha_t <- enhanced_tree$alpha * enhanced_tree$segment_length
	enhanced_tree$exp_2alpha_t <- exp(2 * enhanced_tree$alpha_t)
	vcv <- matrix(0, ntax, ntax)
	rownames(vcv) <- phy$tip.label
	colnames(vcv) <- phy$tip.label
	mrca_cache <- ape::mrca(phy)
	root_node_row <- which(is.na(enhanced_tree$parent_row))[1]
	for (row_index in sequence(ntax)) {
		for (col_index in sequence(ntax)) {
			if(col_index < row_index) {
				next
			}
			
			# along each descendant path
			
			mrca_node <- mrca_cache[row_index, col_index]
			
			
			left_sum <- traverse_down_for_vcv(enhanced_tree, row_index, mrca_node)
			right_sum <- traverse_down_for_vcv(enhanced_tree, col_index, mrca_node)
			exp_part <- exp(-(left_sum + right_sum))
			
			# from mrca to root
			stem_part <- traverse_to_root_from_mrca(enhanced_tree, mrca_node)
			vcv[row_index, col_index] <- stem_part * exp_part			
		}
	}
	
	# fill in lower triangle
	for (row_index in sequence(ntax)) {
		for (col_index in sequence(ntax)) {
			if(col_index < row_index) {
				vcv[row_index, col_index] <- vcv[col_index, row_index]
			}
		}
	}
	
	return(vcv)
}

traverse_down_for_vcv <- function(enhanced_tree, start_tip, end_node) {
	running_sum <- 0
	current_row <- max(which(enhanced_tree$tip_number == start_tip))
	while(!is.na(current_row)) {
		running_sum <- running_sum + enhanced_tree$alpha_t[current_row]
		current_row <- enhanced_tree$parent_row[current_row]
		if(is.na(current_row)) {
			break
		}
		if(enhanced_tree$tipward_treenode[current_row] == end_node) {
			break
		}
	}
	return(running_sum)
}

traverse_to_root_from_mrca <- function(enhanced_tree, end_node) {
	running_sum <- 0
	matching_rows <- which(enhanced_tree$tipward_treenode == end_node)
	if(length(matching_rows) == 0) {
		return(0)
	}
	current_row <- max(matching_rows)
	while (!is.na(current_row)) {
		running_sum <- running_sum +
			enhanced_tree$sigma_squared[current_row] *
				(exp(
					2 *
						enhanced_tree$alpha[current_row] *
						enhanced_tree$tipward_height_segment[current_row]
				) -
					exp(
						2 *
							enhanced_tree$alpha[current_row] *
							enhanced_tree$rootward_height_segment[current_row]
					)) /
				(2 * enhanced_tree$alpha[current_row])
		current_row <- enhanced_tree$parent_row[current_row]
		if (is.na(current_row)) {
			break
		}
	}
	return(running_sum)
}

## Quick VCV maker of OU1 and OUM -- since alpha and sigma.sq are constants, and since the regime does not matter, we can just do a simple plug and chug. It's a tad slower, but mostly for testing purposes.
quickVCV <- function(phy, alpha, sigma.sq, scaleHeight){
    phy$node.label <- NULL
    vcv <- matrix(0, Ntip(phy), Ntip(phy))
    split.times <- branching.times(phy)
    if(scaleHeight == TRUE){
        split.times <- split.times/max(split.times)
    }
    tot.time <- max(split.times)
    for(i in 1:Ntip(phy)){
        for(j in 1:Ntip(phy)){
            if(i == j){
                dij <- 0
                tij <- tot.time
                vcv[i,j] <- (sigma.sq /(2 * alpha)) * exp(-2 * alpha * dij)
            }else{
                split <- getMRCA(phy, tip=c(phy$tip.label[i], phy$tip.label[j]))
                dij <- split.times[which(names(split.times)==split)]
                tij <- tot.time - dij
                vcv[i,j] <- (sigma.sq /(2 * alpha)) * exp(-2 * alpha * dij)
            }
        }
    }
    return(vcv)
}


##Matrix generating function taken from vcv.phylo in ape:
mat.gen <- function(phy,piece.wise,pp){
    phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    ep <- piece.wise[,1]
    comp <- numeric(n + phy$Nnode)
    mat <- matrix(0, n, n)
    
    for (i in length(anc):1) {
        focal <- comp[anc[i]]
        comp[des[i]] <- focal + ep[des[i]]
        j <- i - 1L
        while (anc[j] == anc[i] && j > 0) {
            left <- if (des[j] > n) pp[[des[j] - n]] else des[j]
            right <- if (des[i] > n) pp[[des[i] - n]] else des[i]
            mat[left, right] <- mat[right, left] <- focal
            j <- j - 1L
        }
    }
    diag.elts <- 1 + 0:(n - 1)*(n + 1)
    mat[diag.elts] <- comp[1:n]
    
    mat
}

