
# -------------------------------------
# Color palettes for plotting functions
# -------------------------------------

# rcolorbrewer "set" palette
set.colors = c(brewer.pal(9, 'Set1'), brewer.pal(7, 'Set2'), brewer.pal(12, 'Set3')[c(3,9,8)], 'violetred4')
set.colors[6] = 'khaki2'
set.colors[8] = 'lightskyblue2'
set.colors = rep(set.colors, 10)

# material.io "heat" palette
material.heat = function(n){
    colorRampPalette(
        c(
        "#283593", # indigo 800
        "#3F51B5", # indigo
        "#2196F3", # blue
        "#00BCD4", # cyan
        "#4CAF50", # green
        "#8BC34A", # light green
        "#CDDC39", # lime
        "#FFEB3B", # yellow
        "#FFC107", # amber
        "#FF9800", # orange
        "#FF5722"  # deep orange
        )
    )(n)
}

# NMF library "heatmap" palette
nmf.colors = c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")
