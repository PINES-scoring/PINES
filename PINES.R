# Load libraries
library(parallel)
library(compiler)
library(optparse)

# Get user inputs
option_list = list(make_option(c("--f"), type="character", default=NA,help="path to the file containing the RS IDs of the variants to be scored", metavar="character"),make_option(c("--e"), type="character", default=NA,help="path to the file containing the RS IDs of the variants used for enrichment-based weighting of the epigenetic annotations", metavar="character"),make_option(c("--w"), type="double", default=1,help="manual weighting constant", metavar="character"),make_option(c("--t"), type="character", default=NA,help="tissue to upweight: select from neuro, liver, immune, or heart", metavar="character"),make_option(c("--c"), type="integer", default=1,help="number of cores to use for parallel computations", metavar="character"),make_option(c("--p"), type="character", default="annotations/",help="Path to annotation folder", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.na(opt$f))
{
        stop("Please provide a list of variants to score!")
} else {
	scoring_variants = unique(as.vector(unlist(read.table(opt$f,header=F))))
}
if(!is.na(opt$e))
{
	gwas_variants = unique(as.vector(unlist(read.table(opt$e,header=F))))
} else {
	gwas_variants = NULL
}
if((!is.na(opt$t)) & (!(opt$t %in% c("neuro","liver","immune","heart"))))
{
	stop("If provided by the user, the --t argument must take a value of neuro, liver, immune, or heart")
} else {
	to_tissue_type = opt$t
}
if(!is.integer(opt$c))
{
	stop("If provided by the user, the --c argument must be an integer")
} else {
	ncores = opt$c
}
if(!is.numeric(opt$w))
{
	stop("If provided by the user, the --w argument must be numeric")
} else {
	to_constant = opt$w
}

# Check inputs
if(length(scoring_variants) == 0)
{
	stop("Please provide a list of variants to score!")
}
if(length(gwas_variants) > 0)
{
	fine_mapping_region = unique(c(scoring_variants,gwas_variants))
} 
if(length(gwas_variants) == 0)
{
	fine_mapping_region = unique(scoring_variants)
}

# Set up datasets
load(paste(opt$p,"chrom_interact.RData",sep=""))
load(paste(opt$p,"w_hm3.snplist.RData",sep=""))
load(paste(opt$p,"honey_badger_dhs.RData",sep=""))
load(paste(opt$p,"promoter_indexes.RData",sep=""))
load(paste(opt$p,"fantom_enhancers.RData",sep=""))
enhancer_regions_tmp1 = as.vector(unlist(read.table(paste(opt$p,"enhancer_coordinates.txt",sep=""))))
enhancer_regions_tmp2 = matrix(unlist(strsplit(x=enhancer_regions_tmp1,split=":")),ncol=2,byrow=T)
enhancer_regions_tmp3 = matrix(unlist(strsplit(x=enhancer_regions_tmp2[,2],split="-")),ncol=2,byrow=T)
enhancer_regions_tmp4 = matrix(unlist(strsplit(x=enhancer_regions_tmp2[,1],split="chr")),ncol=2,byrow=T)[,2]
enhancer_regions = cbind(enhancer_regions_tmp4,enhancer_regions_tmp3)
class(enhancer_regions) = "numeric"
load(paste(opt$p,"human_hg19_CTCF_binding_sites.RData",sep=""))
CTCF_regions_liver = CTCF_regions[which(CTCF_regions[,4] == "Liver"),]
tss_data = read.table(paste(opt$p,"mart_export.txt",sep=""),header=T,sep="\t")

# Main annotation function
annotate_SNP = cmpfun(function(SNP_ID)
{
	# Get Haploreg location info
	initial_query = system(paste("grep -F -m 1 '",SNP_ID,";' ",opt$p,"endline_awk_short_unzipped_sorted_haploreg_v4.0_20151021.vcf",sep=""), intern=T)	
        if(length(initial_query) != 1)
        {
		return(0)
        }
       	else
        {
		parsed_initial_query = strsplit(strsplit(initial_query,split=";")[[1]],split=" ")[[1]]
		if(length(parsed_initial_query) != 3)
		{
			stop("HAPLOREG ERROR 2!!")
		}
		if(parsed_initial_query[3] != SNP_ID)
		{
			stop("HAPLOREG ERROR 2!")
		}
		SNP_chr = as.numeric(parsed_initial_query[1])
		SNP_pos = as.numeric(parsed_initial_query[2])
        }

	cell_types = read.table(paste(opt$p,"cell_types.csv",sep=""),header=T,sep=",")
	annotation_matrix = cbind(cell_types,0,0,0,0,0,0,0)
	colnames(annotation_matrix)[5:11] = c("Chrom_HMM_15","Chrom_HMM_25","H3K4me1_Enh","H3K4me3_Pro","H3K27ac_Enh","H3K9ac_Pro","DNase")

	# Extract SNP from HaploReg database
	db_query = paste(system(paste("tabix ",opt$p,"sorted_haploreg_v4.0_20151021.vcf.gz ",SNP_chr,":",SNP_pos-100,"-",SNP_pos+100," | grep \t",SNP_ID,"\t",sep=""), intern=T),collapse="")
	if(db_query == "")
	{
		return(0)
	}
        if(length(db_query) != 1)
        {
		stop("ERROR!")
        }
	if(length(as.vector(unlist(gregexpr(text=db_query,pattern="rs")[[1]]))) != 1)
	{
                return(0)
	}

	# Check intron/exon/intergenic
        locfun_query = grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="DBSNP=",value=T)
        if(length(locfun_query) > 0)
        {
                locfun_annot = strsplit(locfun_query,split="DBSNP=")[[1]][2]
        }
        else
        {
                locfun_annot = "intergenic"
        }
	if(!(locfun_annot %in% c("intergenic","INT")))
	{
                return(0)
	}

	# Get haploreg allele info
	haploreg_allele_A = strsplit(db_query,split="\t")[[1]][4]
        haploreg_allele_B = strsplit(db_query,split="\t")[[1]][5]
	haploreg_location = strsplit(db_query,split="\t")[[1]][2]
        haploreg_chr = strsplit(db_query,split="\t")[[1]][1]

	# Annotate with FANTOM 5 enhancer status
        is_fantom5_enhancer = 0
        chrom_index = which(fantom_enhancers[,1] == as.numeric(haploreg_chr))
        if(length(chrom_index) > 0)
        {
                if(length(which((fantom_enhancers[chrom_index,2] <= as.numeric(haploreg_location)) & (fantom_enhancers[chrom_index,3] >= as.numeric(haploreg_location)))) > 0)
                {
                        is_fantom5_enhancer = 1
                }
        }

	# Annotate with enhancer status
	is_enhancer = 0
	chrom_index = which(enhancer_regions[,1] == as.numeric(haploreg_chr))
	if(length(chrom_index) > 0)
	{
		if(length(which((enhancer_regions[chrom_index,2] <= as.numeric(haploreg_location)) & (enhancer_regions[chrom_index,3] >= as.numeric(haploreg_location)))) > 0)
		{
			is_enhancer = 1
		}
	}

	# Annotate with CTCF status
        is_CTCF = 0
        chrom_index = which(CTCF_regions_liver[,1] == paste("chr",haploreg_chr,sep=""))
        if(length(chrom_index) > 0)
        {
                if(length(which((CTCF_regions_liver[chrom_index,2] <= as.numeric(haploreg_location)) & (CTCF_regions_liver[chrom_index,3] >= as.numeric(haploreg_location)))) > 0)
                {
                        is_CTCF = 1
                }
        }

	# Annotate with ENCODE promoter status
        is_promoter = 0
        chrom_index = which(promoter_indexes[,1] == paste("chr",haploreg_chr,sep=""))
        if(length(chrom_index) > 0)
        {
                if(length(which((promoter_indexes[chrom_index,2] <= as.numeric(haploreg_location)) & (promoter_indexes[chrom_index,3] >= as.numeric(haploreg_location)))) > 0)
                {
                        is_promoter = 1
                }
        }

	# Annotate with distance to nearest TSS
	tss_min_distance = -1
        chrom_index = which(tss_data[,1] == as.numeric(haploreg_chr))
	if(length(chrom_index) > 0)
	{
		tss_min_distance = min(abs(tss_data[chrom_index,2] - as.numeric(haploreg_location)))
	}

	# Annotate HoneyBadger DHS
	honey_dhs_query = which((honey_badger_dhs[,1] == paste("chr",haploreg_chr,sep="")) & (honey_badger_dhs[,2] <= as.numeric(haploreg_location)) & (honey_badger_dhs[,3] >= as.numeric(haploreg_location)))
	if(length(honey_dhs_query) > 0)
	{
		is_honey_dhs = 1
	}
	else
	{
		is_honey_dhs = 0
	}

	# Annotate MAF
        maf_query = grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="AF=",value=T)
	if(length(maf_query) > 0)
	{
		maf_annot = as.numeric(strsplit(strsplit(maf_query,split="AF=")[[1]][2],split=",")[[1]])
	}
	else
	{
		maf_annot = c(0,0,0,0)
	}

        # Placeholders
        gwava_scores = c(1,1,1)
        Eigen_score = 1
        CADD_score = 1
        dann_score = 1
        linsight_score = 1
        fathmm_score = 1

	# Annotate conservation
	gerp_annot = length(grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="GERP",value=T))
        siphy_annot = length(grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="SIPHY",value=T))

	# Annotate 4C interaction region
	if(length(which((chrom_interact$InteractorAStart <= as.numeric(haploreg_location)) & (chrom_interact$InteractorAEnd >= as.numeric(haploreg_location)) & (chrom_interact$InteractorAChr == paste("chr",haploreg_chr,sep="")))) > 0)
	{
		c4_annot = 1
	}
	else
	{
		if(length(which((chrom_interact$InteractorBStart <= as.numeric(haploreg_location)) & (chrom_interact$InteractorBEnd >= as.numeric(haploreg_location)) & (chrom_interact$InteractorBChr == paste("chr",haploreg_chr,sep="")))) > 0)
		{
			c4_annot = 1
		}
		else
		{
			c4_annot = 0
		}
	}

	# Annotate histone modifications
	histone_query = grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="ROADMAP_HISTONE_MARK_PEAKS",value=T)
	if(length(histone_query) > 0)
	{
		histone_annot = matrix(unlist(strsplit(strsplit(strsplit(histone_query,split="=")[[1]][2],split=",")[[1]],split="\\|")),ncol=2,byrow=T)
		for(i in 1:nrow(histone_annot))
		{
			row_index = match(histone_annot[i,1],annotation_matrix[,1])
			col_index = match(histone_annot[i,2],colnames(annotation_matrix))

			if(is.na(row_index) | is.na(col_index))
			{
				stop("Histone ERROR!")
			}

			annotation_matrix[row_index,col_index] = 1
		}
	}

	# Annotate DNase regions
	dnase_query = grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="ROADMAP_DNASE",value=T)
	if(length(dnase_query) > 0)
	{
		dnase_annot = strsplit(strsplit(dnase_query,split="=")[[1]][2],split=",")[[1]]
		for(i in 1:length(dnase_annot))
		{
                        row_index = match(dnase_annot[i],annotation_matrix[,1])

			if(is.na(row_index))
			{
				stop("DNase ERROR!")
			}

			annotation_matrix[row_index,11] = 1
		}
	}

	# Annotate HMM 15
	hmm15_query = grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="ROADMAP_CHROMHMM_15STATE",value=T)
	if(length(hmm15_query) > 0)
	{
		hmm15_annot = matrix(unlist(strsplit(strsplit(strsplit(hmm15_query,split="=")[[1]][2],split=",")[[1]],split="\\|")),ncol=2,byrow=T)
       		for(i in 1:nrow(hmm15_annot))
        	{
			row_index = match(hmm15_annot[i,1],annotation_matrix[,1])

                        if(is.na(row_index))
                        {
                                stop("DNase ERROR!")
                        }

			if(strsplit(hmm15_annot[i,2],split="_")[[1]][1] %in% c(1,2,3,4,5,6,7,8))
			{
                		annotation_matrix[row_index,5] = 1
			}
			else
			{
				annotation_matrix[row_index,5] = -1
			}
        	}
	}

	# Annotate HMM 25
	hmm25_query = grep(strsplit(strsplit(db_query,split="\t")[[1]][8],split=";")[[1]],pattern="ROADMAP_CHROMHMM_25STATE",value=T)
	if(length(hmm25_query) > 0)
	{
		hmm25_annot = matrix(unlist(strsplit(strsplit(strsplit(hmm25_query,split="=")[[1]][2],split=",")[[1]],split="\\|")),ncol=2,byrow=T)
                for(i in 1:nrow(hmm25_annot))
                {
			row_index = match(hmm25_annot[i,1],annotation_matrix[,1])

                        if(is.na(row_index))
                        {
                                stop("DNase ERROR!")
                        }

                        if(strsplit(hmm25_annot[i,2],split="_")[[1]][1] %in% 1:20)
                        {
                                annotation_matrix[row_index,6] = 1
                        }
                        else
                        {
                                annotation_matrix[row_index,6] = -1
                        }
                }
	}

	# Condense annotation_matrix into a vector
	annotation_vector = NULL
	for(i in 1:nrow(annotation_matrix))
	{
		for(j in 5:ncol(annotation_matrix))
		{
			annotation_vector = c(annotation_vector,annotation_matrix[i,j])
			names(annotation_vector)[length(annotation_vector)] = paste(annotation_matrix[i,1],"_",colnames(annotation_matrix)[j],sep="")
		}
	}
	old_length = length(annotation_vector)
	annotation_vector = c(annotation_vector,maf_annot,gerp_annot,siphy_annot,c4_annot,is_honey_dhs,is_enhancer,is_CTCF,tss_min_distance,SNP_ID,locfun_annot,is_promoter,is_fantom5_enhancer,gwava_scores,Eigen_score,CADD_score,dann_score,linsight_score,fathmm_score)
	names(annotation_vector)[(old_length + 1):length(annotation_vector)] = c("maf_AFR","maf_AMR","maf_ASN","maf_EUR","gerp_conservation","siphy_conservation","interactions_4C","HoneyBadger_DHS","is_enhancer","is_CTCF","tss_min_distance","SNP_ID","locfun_annot","is_promoter","is_fantom5_enhancer","gwava_region","gwava_TSS","gwava_unmatched","Eigen_score","CADD_score","dann_score","linsight_score","fathmm_score")
	annotation_vector = c(SNP_chr,SNP_pos,annotation_vector)
	names(annotation_vector)[1:2] = c("chrom","location")

	if(length(annotation_vector) != 914)
        {
        	stop("Annotation length ERROR!")
        }

	write.table(t(annotation_vector),file=paste("annotation_vector_",SNP_ID,".txt",sep=""),sep=",",row.names=F,quote=F)

	return(0)		
})

# Genome-wide run
invisible(mclapply(X=fine_mapping_region,FUN=annotate_SNP,mc.cores=min(length(fine_mapping_region),ncores)))

# Collect results
system("ls annotation_vector_* | xargs -n 32 -P 1 cat >> combined_annotations.txt")

if(file.size("combined_annotations.txt") == 0)
{
	stop("The provided variants could not be annotated. Make sure that you are only attempting to score noncoding (intronic or intergenic) variants!")
}

new_annotations = read.table(file="combined_annotations.txt",header=T,sep=",")
rem_ind = which(new_annotations[,1] == "chrom")
if(length(rem_ind) > 0)
{
	combined_annotations = new_annotations[(-1)*rem_ind,]
} else {
	combined_annotations = new_annotations
}
save(combined_annotations,file="combined_annotations.RData")

# Clean workspace and free memory
rm(annotate_SNP,chrom_interact,CTCF_regions,CTCF_regions_liver,enhancer_regions,enhancer_regions_tmp1,enhancer_regions_tmp2,enhancer_regions_tmp3,enhancer_regions_tmp4,fantom_enhancers,final_SNP_matrix,honey_badger_dhs,promoter_indexes,tss_data)
invisible(gc())

# print(combined_annotations)

# Run scoring
data = as.matrix(local(get(load(paste(opt$p,"background_variants/combined_annotations.RData",sep="")))))
load(paste(opt$p,"eQTLs.RData",sep=""))
cell_types = read.table(paste(opt$p,"cell_types.csv",sep=""),header=T,sep=",")
tss_data = read.table(paste(opt$p,"mart_export.txt",sep=""),header=T,sep="\t")

# Incorporate additional annotations
if(!all(colnames(combined_annotations) == colnames(data)))
{
	stop("Annotation Error!")
}
if(length(which(!(combined_annotations[,903] %in% data[,903]))) > 0)
{
	to_add = as.matrix(combined_annotations[which(!(combined_annotations[,903] %in% data[,903])),])
	data = rbind(data,to_add)
	rm(to_add)
}

# Set indexes of variants to be returned to the user
test_variant_indexes = match(scoring_variants,data[,903])
index_of_na = which(is.na(test_variant_indexes))
if(length(index_of_na) > 0)
{
	test_variant_indexes = test_variant_indexes[(-1)*index_of_na]
}
if(length(test_variant_indexes) == 0)
{
	stop("The provided variants could not be annotated. Make sure that you are only attempting to score noncoding (intronic or intergenic) SNPs!")
}
rm(combined_annotations)
invisible(gc())

print("test_variant_indexes")
print(test_variant_indexes)

# Pick background indexes
pathogenic_variants = as.vector(unlist(read.table(paste(opt$p,"combined_regions_for_annotation.txt",sep=""))))
background_indexes = 1:nrow(data)
remove_pathogenic_indexes = which(data[background_indexes,903] %in% pathogenic_variants)
remove_eQTL_indexes = which(data[background_indexes,903] %in% eQTLs)
remove_enhancer_indexes = which(as.numeric(data[background_indexes,which(colnames(data) == "is_enhancer")]) == 1)
remove_promoter_indexes = which(as.numeric(data[background_indexes,which(colnames(data) == "is_promoter")]) == 1)
remove_CTCF_indexes = which(as.numeric(data[background_indexes,which(colnames(data) == "is_CTCF")]) == 1)
remove_fantom5_indexes = which(as.numeric(data[background_indexes,which(colnames(data) == "is_fantom5_enhancer")]) == 1)
remove_all_indexes = unique(c(remove_pathogenic_indexes,remove_eQTL_indexes,remove_enhancer_indexes,remove_promoter_indexes,remove_CTCF_indexes,remove_fantom5_indexes))
if(length(remove_all_indexes) > 0)
{
	background_indexes = background_indexes[(-1) * remove_all_indexes]
}
intronic_background_indexes = sample(background_indexes[which(data[background_indexes,904] == "INT")])
intergenic_background_indexes = sample(background_indexes[which(data[background_indexes,904] == "intergenic")])
maf_AFR = as.numeric(data[,892])
maf_AMR = as.numeric(data[,893])
maf_ASN = as.numeric(data[,894])
maf_EUR = as.numeric(data[,895])
common_variants = which((maf_AFR >= 0.05) & (maf_AMR >= 0.05) & (maf_ASN >= 0.05) & (maf_EUR >= 0.05) & (maf_AFR <= 0.95) & (maf_AMR <= 0.95) & (maf_ASN <= 0.95) & (maf_EUR <= 0.95))
intronic_background_indexes = intersect(intronic_background_indexes,common_variants)
intergenic_background_indexes = intersect(intergenic_background_indexes,common_variants)

print("Length intronic / intergenic background indexes")
print(length(intronic_background_indexes))
print(length(intergenic_background_indexes))

# Set TSS proximity indicator
TSS_dist = rep(0,times=nrow(data))
TSS_dist[which(as.numeric(data[,which(colnames(data) == "tss_min_distance")]) < 1000)] = 1

# Assemble shortened annotation matrix
annotation_matrix = cbind(data[,grep(x=colnames(data),pattern="H3K|DNase|conservation|HoneyBadger_DHS|Chrom_HMM|interactions_4C")],TSS_dist)
class(annotation_matrix) = "numeric"
annotation_matrix[which(annotation_matrix == -1)] = 0

# Remove all-zero columns in background data
zero_indexes = which(colMeans(annotation_matrix[c(intronic_background_indexes,intergenic_background_indexes),]) == 0)
if(length(zero_indexes) > 0)
{
        annotation_matrix = annotation_matrix[,(-1)*zero_indexes]
}
invisible(gc())

# Compute factor weights via GWAS enrichments (if available)
if(length(gwas_variants) > 0)
{
	num_draws_enrich = 200

	GWAS_hits = gwas_variants
	GWAS_variants = intersect(GWAS_hits,data[,903])
	index_of_GWAS_variants = match(GWAS_variants,data[,903])
	mafs_of_GWAS_variants = as.numeric(data[index_of_GWAS_variants,895])

	indexes_of_0_maf = which(mafs_of_GWAS_variants == 0)
	if(length(indexes_of_0_maf) > 0)
	{
		index_of_GWAS_variants = index_of_GWAS_variants[(-1) * indexes_of_0_maf]
		mafs_of_GWAS_variants = as.numeric(data[index_of_GWAS_variants,895])
	}

	matrix_of_matched_indexes = matrix(ncol=num_draws_enrich,nrow=length(index_of_GWAS_variants))
	for(i in 1:length(mafs_of_GWAS_variants))
	{
		diffs_in_maf = abs(as.numeric(data[,895]) - mafs_of_GWAS_variants[i])
		all_matched_variants = setdiff(x=which(diffs_in_maf < 0.05), y=index_of_GWAS_variants)
		matrix_of_matched_indexes[i,] = sample(x=all_matched_variants,size=num_draws_enrich,replace=F)
	}
	feature_enrichment_pvalues = vector(length=ncol(annotation_matrix))
	for(i in 1:ncol(annotation_matrix))
	{
		null_dist = vector(length=num_draws_enrich)
		for(k in 1:ncol(matrix_of_matched_indexes))
		{
			null_dist[k] = sum(annotation_matrix[matrix_of_matched_indexes[,k],i])
		}
		actual_overlap = sum(annotation_matrix[index_of_GWAS_variants,i])
		feature_enrichment_pvalues[i] = ppois(q=(actual_overlap-1),lambda=mean(null_dist),lower.tail=F)
	}

	feature_enrichment_pvalues = (-1) * log10(feature_enrichment_pvalues)

	copy_of_feature_enrichment_pvalues = feature_enrichment_pvalues
	names(copy_of_feature_enrichment_pvalues) = colnames(annotation_matrix)
	save(copy_of_feature_enrichment_pvalues,file="copy_of_feature_enrichment_pvalues.RData")

	invisible(gc())
}

# Set up cell type indexes
immune_cells_names = as.vector(unlist(cell_types[grep(x=cell_types[,4],pattern="Primary T|Primary B|Leukemia|Killer|Lymphoblastoid"),1]))
hdl_cells_names = as.vector(unlist(cell_types[grep(x=cell_types[,4],pattern="Liver|Hepato"),1]))
neuro_cell_names = as.vector(unlist(cell_types[grep(x=cell_types[,2],pattern="Brain"),1]))
heart_cell_names = as.vector(unlist(cell_types[grep(x=cell_types[,2],pattern="Muscle"),1]))

immune_cells = NULL
for(i in immune_cells_names)
{
        immune_cells = c(immune_cells,grep(pattern=i,x=colnames(annotation_matrix)))
}
hdl_cells = NULL
for(i in hdl_cells_names)
{
        hdl_cells = c(hdl_cells,grep(pattern=i,x=colnames(annotation_matrix)))
}
neuro_cells = NULL
for(i in neuro_cell_names)
{
	neuro_cells = c(neuro_cells,grep(pattern=i,x=colnames(annotation_matrix)))
}
heart_cells = NULL
for(i in heart_cell_names)
{
        heart_cells = c(heart_cells,grep(pattern=i,x=colnames(annotation_matrix)))
}

# Set up additional indexes
TSS_indexes = which(colnames(annotation_matrix) %in% c("TSS_dist","TSS_dist_pheno"))
interaction_4C_indexes = which(colnames(annotation_matrix) %in% c("interactions_4C"))
DHS_indexes = grep(x=colnames(annotation_matrix),pattern="HoneyBadger_DHS")
conservation_indexes = grep(x=colnames(annotation_matrix),pattern="conservation")

# Set spectral decomposition framework
intronic_background_means = colMeans(annotation_matrix[intronic_background_indexes,])
intergenic_background_means = colMeans(annotation_matrix[intergenic_background_indexes,])

cov_intronic = cov(annotation_matrix[intronic_background_indexes,])
cov_intergenic = cov(annotation_matrix[intergenic_background_indexes,])

spectral_decomp_intronic = eigen(cov_intronic,symmetric=T)
evalues_intronic = spectral_decomp_intronic$values
evectors_intronic = spectral_decomp_intronic$vectors

spectral_decomp_intergenic = eigen(cov_intergenic,symmetric=T)
evalues_intergenic = spectral_decomp_intergenic$values
evectors_intergenic = spectral_decomp_intergenic$vectors

smoothed_cov_intronic = 0
for(i in 1:30)
{
        smoothed_cov_intronic = smoothed_cov_intronic + evalues_intronic[i] * t(t(evectors_intronic[,i])) %*% t(evectors_intronic[,i])
}
inverse_smoothed_cov_intronic = 0
for(i in 1:30)
{
        inverse_smoothed_cov_intronic = inverse_smoothed_cov_intronic + (1 / evalues_intronic[i]) * t(t(evectors_intronic[,i])) %*% t(evectors_intronic[,i])
}

smoothed_cov_intergenic = 0
for(i in 1:30)
{
        smoothed_cov_intergenic = smoothed_cov_intergenic + evalues_intergenic[i] * t(t(evectors_intergenic[,i])) %*% t(evectors_intergenic[,i])
}
inverse_smoothed_cov_intergenic = 0
for(i in 1:30)
{
        inverse_smoothed_cov_intergenic = inverse_smoothed_cov_intergenic + (1 / evalues_intergenic[i]) * t(t(evectors_intergenic[,i])) %*% t(evectors_intergenic[,i])
}

invisible(gc())

# Get weight matrix
if(length(gwas_variants) > 0)
{
	w = diag(feature_enrichment_pvalues)
} 

if((length(gwas_variants) == 0) & (!is.na(to_constant)))
{
	w = diag(ncol(annotation_matrix))

	if(to_tissue_type == "neuro")
	{
		for(i in neuro_cells)
		{
       			w[i,i] = to_constant
		}
	}
        if(to_tissue_type == "liver")
        {
        	for(i in hdl_cells)
                {
                	w[i,i] = to_constant
                }
	}
        if(to_tissue_type == "immune")
        {
        	for(i in immune_cells)
                {
                	w[i,i] = to_constant
                }
	}
        if(to_tissue_type == "heart")
        {
        	for(i in heart_cells)
                {
                	w[i,i] = to_constant
                }
	}

	for(i in TSS_indexes)
	{
        	w[i,i] = 1
	}
	for(i in DHS_indexes)
	{
        	w[i,i] = 1
	}
	for(i in conservation_indexes)
	{
        	w[i,i] = 1
	}
	for(i in interaction_4C_indexes)
	{
        	w[i,i] = 1
	}
}

if((length(gwas_variants) == 0) & (is.na(to_constant)))
{
	w = diag(ncol(annotation_matrix))
}

# Get weighted Cholesky matrices
inverse_smoothed_cov_intronic_uw = inverse_smoothed_cov_intronic
inverse_smoothed_cov_intergenic_uw = inverse_smoothed_cov_intergenic
inverse_smoothed_cov_intronic_w = (w %*% inverse_smoothed_cov_intronic) %*% w
inverse_smoothed_cov_intergenic_w = (w %*% inverse_smoothed_cov_intergenic) %*% w

# Pre-compute length of all-1 vector for scoring
all_1_entry = rep(1,times=ncol(annotation_matrix))
len_all_1_entry_intronic_uw = sqrt(as.numeric(t(all_1_entry - intronic_background_means) %*% inverse_smoothed_cov_intronic_uw %*% (all_1_entry - intronic_background_means)))
len_all_1_entry_intergenic_uw = sqrt(as.numeric(t(all_1_entry - intergenic_background_means) %*% inverse_smoothed_cov_intergenic_uw %*% (all_1_entry - intergenic_background_means)))
len_all_1_entry_intronic_w = sqrt(as.numeric(t(all_1_entry - intronic_background_means) %*% inverse_smoothed_cov_intronic_w %*% (all_1_entry - intronic_background_means)))
len_all_1_entry_intergenic_w = sqrt(as.numeric(t(all_1_entry - intergenic_background_means) %*% inverse_smoothed_cov_intergenic_w %*% (all_1_entry - intergenic_background_means)))

# Save vectors for later
variant_type_data_vector = data[,904]
SNP_name_data_vector = data[,903]
variant_chr_data_vector = as.numeric(data[,1])
variant_pos_data_vector = as.numeric(data[,2])

# Test functions
get_mah_distance_w = cmpfun(function(index,type)
{
	if(type == "INT")
	{
		inv_sm_cov = inverse_smoothed_cov_intronic_w
		len_all_1 = len_all_1_entry_intronic_w
		back_means = intronic_background_means
	}else{
                inv_sm_cov = inverse_smoothed_cov_intergenic_w
                len_all_1 = len_all_1_entry_intergenic_w
                back_means = intergenic_background_means
	}
	x_annotations = annotation_matrix[index,]
	x_mu_sigma = t(x_annotations - back_means) %*% inv_sm_cov
	x_mu_sigma_y_mu = as.numeric(x_mu_sigma %*% (all_1_entry - back_means))
        len_x = sqrt(as.numeric(x_mu_sigma %*% (x_annotations - back_means)))
	if(len_x == 0)
	{
		return(1000)
	}
	theta = acos(x_mu_sigma_y_mu / (len_x * len_all_1)) / pi
	return(theta / (len_x * len_all_1))
})

get_mah_distance_uw = cmpfun(function(index,type)
{
        if(type == "INT")
        {
                inv_sm_cov = inverse_smoothed_cov_intronic_uw
                len_all_1 = len_all_1_entry_intronic_uw
                back_means = intronic_background_means
        }else{
                inv_sm_cov = inverse_smoothed_cov_intergenic_uw
                len_all_1 = len_all_1_entry_intergenic_uw
                back_means = intergenic_background_means
        }
        x_annotations = annotation_matrix[index,]
        x_mu_sigma = t(x_annotations - back_means) %*% inv_sm_cov
        x_mu_sigma_y_mu = as.numeric(x_mu_sigma %*% (all_1_entry - back_means))
        len_x = sqrt(as.numeric(x_mu_sigma %*% (x_annotations - back_means)))
        if(len_x == 0)
        {
                return(1000)
        }
        theta = acos(x_mu_sigma_y_mu / (len_x * len_all_1)) / pi
        return(theta / (len_x * len_all_1))
})

get_pvalues = cmpfun(function(i)
{
	variant_type = variant_type_data_vector[i]
	SNP_name = SNP_name_data_vector[i]
	variant_chr = toString(variant_chr_data_vector[i])
	variant_pos = toString(variant_pos_data_vector[i])

	mah_distance_w = get_mah_distance_w(index=i,type=variant_type)
	mah_distance_uw = get_mah_distance_uw(index=i,type=variant_type)

	if(variant_type == "INT")
        {
        	pvalue_method_weighted = length(which(background_distances_intronic_w <= mah_distance_w)) / length(background_distances_intronic_w)
                pvalue_method_unweighted = length(which(background_distances_intronic_uw <= mah_distance_uw)) / length(background_distances_intronic_uw)
	}else{
                pvalue_method_weighted = length(which(background_distances_intergenic_w <= mah_distance_w)) / length(background_distances_intergenic_w)
                pvalue_method_unweighted = length(which(background_distances_intergenic_uw <= mah_distance_uw)) / length(background_distances_intergenic_uw)
	}

	cat(SNP_name,variant_chr,variant_pos,pvalue_method_unweighted,pvalue_method_weighted,"\n",file=paste("pvalues_",SNP_name,".txt",sep=""))
        return(0)
})

rm(data)
invisible(gc())

# Get background score distributions
background_distances_intronic_uw = as.vector(unlist(mclapply(X=intronic_background_indexes,FUN=get_mah_distance_uw,type="INT",mc.cores=ncores)))
background_distances_intergenic_uw = as.vector(unlist(mclapply(X=intergenic_background_indexes,FUN=get_mah_distance_uw,type="intergenic",mc.cores=ncores)))
background_distances_intronic_w = as.vector(unlist(mclapply(X=intronic_background_indexes,FUN=get_mah_distance_w,type="INT",mc.cores=ncores)))
background_distances_intergenic_w = as.vector(unlist(mclapply(X=intergenic_background_indexes,FUN=get_mah_distance_w,type="intergenic",mc.cores=ncores)))
invisible(gc())

# Compute p-values
mclapply(X=test_variant_indexes,FUN=get_pvalues,mc.cores=min(length(test_variant_indexes),ncores))

# Create output file
cat("SNP_ID chr pos Unweighted_PINES Weighted_PINES\n",file="PINES_pvalues.txt")
system("ls pvalues_* | xargs -n 32 -P 1 cat >> PINES_pvalues.txt")
system("ls pvalues_* | xargs -n 32 -P 1 rm")
system("ls annotation_vector_* | xargs -n 32 -P 1 rm")
system("rm -f combined_annotations.txt")

print("PINES scoring complete!")

