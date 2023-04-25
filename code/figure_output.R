# this is output for figures

pdf(file = "results/volcano.pdf", width = 4, height = 4)
print(volcano_plot)
dev.off()

pdf(file = "results/heatmap.pdf", width = 6.3, height = 5)
print(heat_plot)
dev.off()
