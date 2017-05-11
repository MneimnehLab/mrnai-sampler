library(cluster)

# command lines arguments
# the first arg is the name of the input file contains the distance matrix
args <- commandArgs(trailingOnly = TRUE)
# print(args[1])


# read the input file as a table
inputTable <- read.table(args[1])


# convert the table into a distance matrix
distMatrix <- as.dist(inputTable)


# count the number of rows in the matrix
# this is essentially the number of structs
# why - 1 ???
rows <- dim(inputTable)[1] - 1


# iterations is the number of times we want to cut the tree
# it is 1 less than num of structs, since min # clusters = 2 and max = n 
# (therefore n - 2 + 1)
iters <- rows-1


# perform hierarchical clustering on the input data
# the result is kind of a dendogram

dendrogram <- hclust(distMatrix, method="complete")		# options: "single"/"complete"/?  but single looks better from some reason
# dendrogram <- hclust(distMatrix, method="single")		# options: "single"/"complete"/?  but single looks better from some reason
# x11()
 postscript("filename.eps", horizontal=F, width=60, height=12, paper="special", onefile=F)
 plot(dendrogram)



# groups <- cutree(dendrogram, h=0.95)
groups <- cutree(dendrogram, h=0.99)
# dput(groups) 
# for (i in groups)
# {
# 	cat(i, " ")
# }
# cat("\n")
write(groups, file = "clusters.out", append = FALSE, sep = " ")


# cutree cuts the dendogram in to k clusters
# since here k is a vector, the result of cutree will be a list (X) of list (Y) of clusters
# so X[i] is a list of i+1 clusters

######`# ---> comment from here to end
# groups <- cutree(dendrogram, k=2:rows) 


# # for each element in each cluster in X[i], we will compute the silhouette
# # the mean of each X[i] (i.e. each cut) is the mean of the sil. of each element in X[i]
# # so cutMeans is a list of means for each cut
# cutMeans <- rep(NA, iters)


# for (i in 1:iters ) 
# {
# 	# groups[,i] -> a cut consisting of i+1 clusters.
# 	#               a complex data structure but essentially a labelling of each point (elem) by its cluster num
	
# 	# comptute the silhouette of each element in this cut
# 	# sils contains a labelling of each point and the silhouette of each (individual!) point
# 	sils <- silhouette(groups[,i], distMatrix)
	
	
# 	# compute the mean of EACH CLUSTER
# 	# aggregate will first group points by their clusters using sils[,"cluster"], 
# 	# then compute the mean (FUN=mean) silhouette of each cluster
# 	# so agg -> list of i+1 means (for i+1 clusters in the i^th cut)
# 	# why [2] ???
# 	agg = aggregate(sils[,"sil_width"], by=list(sils[,"cluster"]), FUN=mean)[2]
	

# 	# we are interested in the silhouette of the cut, which is the mean of the
# 	# silhouettes of each cluster
# 	# cutMeans[i] -> silhouette of i^th cut (consisting of i+1 clusters)
# 	cutMeans[i] <- mean(agg[,1])
# 	# cat("i = ", i, ",  mean = ", cutMeans[i], "\n")
# }


# # plot a the curve of cut means
# # 
# # x11()
# # plot(cutMeans)
# # lines(cutMeans)

# # theoretically(?), the cut with the silhouette index closest to 1 has the best fit
# # since -1 <= s <= 1, pick the cut with the maximum silhouette index
# # wm -> index in array, so best num of clusters is wm+1
# wm <- which.max(cutMeans)
# cat("Optimal num of clusters: ", wm+1, "\n")
# # for (i in groups[,wm])
# # {
# # 	cat(i, " ")
# # }
# # cat("\n")


# # we need to the clustering details so that we can label actual structs
# # so output the labeling of each point/element in the input dataset
# # since the best cut is at index wm, get the labeling from groups[,wm]
# # write the output of file "clusters"
# write(groups[,wm], file = "clusters", append = FALSE, sep = " ")


# # If the script terminates, the plots disappear (if generated)
# # to keep the plots alive, the scripts needs to be kept alive
# # wait for input to script alive
# # 

# # message("Press Return To Continue")
# # invisible(readLines("stdin", n=1))
