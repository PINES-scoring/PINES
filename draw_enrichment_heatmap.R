# Process data
library(ggplot2)
library(optparse)

option_list = list(make_option(c("--f"), type="character", default=NA,help="full path to the PINES output file copy_of_feature_enrichment_pvalues.RData", metavar="character"),make_option(c("--p"), type="character", default="annotations/",help="Path to annotation folder", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.na(opt$f))
{
	stop("Enter the full path to the PINES output file copy_of_feature_enrichment_pvalues.RData")
}

load(opt$f)
cell_types = read.table(paste(opt$p,"cell_types.csv",sep=""),header=T,sep=",")
histone_data = copy_of_feature_enrichment_pvalues[grep(x=names(copy_of_feature_enrichment_pvalues),pattern="H3K")]
histone_names = matrix(unlist(strsplit(names(histone_data),split="_")),byrow=T,ncol=3)
histone_group = as.vector(unlist(cell_types[match(histone_names[,1],cell_types[,1]),3]))
dhs_data = copy_of_feature_enrichment_pvalues[grep(x=names(copy_of_feature_enrichment_pvalues),pattern="DNase")]
dhs_names = matrix(unlist(strsplit(names(dhs_data),split="_")),byrow=T,ncol=2)
dhs_group = as.vector(unlist(cell_types[match(dhs_names[,1],cell_types[,1]),3]))
weight_matrix = rbind(cbind(histone_group,histone_names[,2]),cbind(dhs_group,dhs_names[,2]))
colnames(weight_matrix) = c("cell_type","annotation")
df_weight_matrix = as.data.frame(weight_matrix)
df_weight_matrix$weight = c(as.vector(unlist(histone_data)),as.vector(unlist(dhs_data)))

# Plot heatmap
base_size = 9
pdf(file="feature_enrichment_pvalues.pdf",width = 14, height = 7)
print(ggplot(data = df_weight_matrix, aes(y = annotation, x = cell_type)) + geom_tile(aes(fill = weight),colour="black") + scale_fill_gradient(low = "darkseagreen1",high = "blue4") + theme_grey(base_size = base_size) + labs(x = "", y = "", fill="Enrichment -log10(p)") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(legend.position = "bottom", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 90, hjust = 0, colour = "grey50")))
dev.off()

