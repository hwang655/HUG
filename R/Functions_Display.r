# functions for display
# image.matrix: display a matrix
# PlotGraphs: plot multiple undirected graphs in a single graph

# Display matrix in an image
image.matrix = function(mat, names.row = NULL, names.col = NULL, ...){
    mat = as.matrix(mat)
	n   = nrow(mat)
	p   = ncol(mat)
	if (is.null(names.row)) names.row = 1:n
	if (is.null(names.col)) names.col = 1:p
	image.plot(names.col, names.row, t(mat), 
		xlim = c(0.5, p+0.5), ylim = c(n + 0.5, 0.5), ...)
}

# generate ggplot2-like colors
gg_color_hue = function(n){
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

# icovs = Omegas; graph.layout = NULL; v.keep = FALSE; v.label = NULL; rmd = FALSE; name.icovs = NULL
# Plot Multiple Graphs in a Single Graph
PlotGraphs = function(icovs, graph4layout = NULL, cmpl.layout = NULL, v.keep = NULL, 
                      v.name = NULL, v.label = NULL, rmd = TRUE, name.icovs = NULL, 
                      legend = TRUE, layout.only = FALSE, conflict.disp = FALSE, ...) {
    if (!is.list(icovs)) {
        type.icovs = attr(class(icovs), "package")
        if (!is.null(type.icovs)) {
            if (type.icovs == "Matrix") {
                icovs = as.matrix(icovs)
            }
        }
        stopifnot(is.matrix(icovs))
        icovs == (icovs != 0) + 0
        if (!all(icovs == t(icovs))) {
            cat("Not a symmetric graph. Taking union...\n")
            icovs = icovs | t(icovs)
        }
        icovs = list(icovs)
    }
    K = length(icovs)
    p = nrow(icovs[[1]])
    lower.p = which(lower.tri(matrix(0, p, p)))
    name.icovs = if (is.null(name.icovs)) paste("GM", 1:K) else name.icovs
    graphs = HUG2::SummaryGraphs(icovs)
    cmpl.union = graphs$graph.union
    cmpl.common = graphs$graph.common
    cmpl.diff = graphs$graph.diff
    cmpl.each = graphs$graph.each
    any.diff = any(cmpl.diff[lower.p] > 0)
    colors.K = c("black", gg_color_hue(K + any.diff))
    if (!is.null(v.keep)) {
        connected = v.keep
    } else {
        diag(cmpl.union) = 0
        connected = which(Matrix::colSums(cmpl.union) > 0)
    }
    graph.union = cmpl.union[connected, connected]
    graph.common = cmpl.common[connected, connected]
    graph.diff = cmpl.diff[connected, connected]
    v.label = if (is.null(v.label)) NA else v.label[connected]
    v.name = if (is.null(v.name)) NA else v.name[connected]
    set.seed(41)
    if (is.null(cmpl.layout)) {
        graph4layout = if (is.null(graph4layout)) cmpl.union else graph4layout
        cmpl.layout = layout_nicely(graph.adjacency(as.matrix(graph4layout), 
                                                    mode = "undirected", diag = FALSE))
    }
    layout.grid = cmpl.layout[connected,]
    if (!layout.only) {
        if (!rmd) windows()
        for (k in 1:K) {
            graph.k = cmpl.each[[k]][connected, connected]
            g.k = graph.adjacency(as.matrix(graph.k), mode = "undirected", diag = FALSE)
            plot(g.k, layout = layout.grid, edge.color = colors.K[k+1+any.diff], 
                 edge.width = 1, edge.lty = 1, 
                 vertex.label = v.name, vertex.color = v.label, vertex.size = 3, 
                 vertex.frame.color = "grey")
            par(new = TRUE)
        }
        g.common = graph.adjacency(as.matrix(graph.common), mode = "undirected", diag = FALSE)
        plot(g.common, layout = layout.grid, edge.color = colors.K[1], 
             edge.width = 2, edge.lty = 1, vertex.label = v.name, 
             vertex.color = v.label, vertex.size = 3, vertex.frame.color = "grey")
        if (any.diff & conflict.disp) {
            par(new = TRUE)
            g.diff = graph.adjacency(as.matrix(graph.diff), mode = "undirected", diag = FALSE)
            plot(g.diff, layout = layout.grid, edge.color = colors.K[2], 
                 edge.width = 2, edge.lty = 1, vertex.label = v.name, 
                 vertex.color = v.label, vertex.size = 3, vertex.frame.color = "grey")
            if (legend) {
                legend("topleft", legend = c("Common", "Conflicts", name.icovs), 
                       col = colors.K, lty = c(1, 1, rep(1, K)), lwd = c(2, 2, rep(0.7, K)))
            }
        }
        else if ((K > 1) & legend) {
            legend("topleft", legend = c("Common", name.icovs), 
                   col = colors.K, lty = c(1, rep(1, K)), lwd = c(2, rep(1, K)))
        }
        par(new = FALSE)
    }
    return(invisible(cmpl.layout))
}