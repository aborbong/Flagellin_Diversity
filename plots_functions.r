#Functions for flagellome diversity analyses 


#Setup for label sizes in plots
labels.x=theme(axis.title.x = element_text(size=14),
               axis.text.x = element_text(size=12))

labels.x.pca=theme(axis.title.x = element_text(size=14),
               axis.text.x = element_text(size=12)) 

labels.y=theme(axis.title.y = element_text(size=14),
               axis.text.y= element_text(size=12))

labels.legend=theme(legend.text = element_text (size=14),
                    legend.title = element_text(size=14))


#' Plot Color Continuous
#'
#' This function generates an ordination plot with continuous color mapping based on an input ordination matrix and color variable.
#'
#' @param physeq_obj A phyloseq object.
#' @param ordination.matrix Ordination matrix for plotting.
#' @param color_by Variable used for continuous color mapping.
#' @param axes Optional parameter specifying the axes to plot. Defaults to NULL, which plots the first two axes.
#'
#' @return An ordination plot with continuous color mapping.
#'
#' @import phyloseq
#' @import ggplot2
#' @import viridis
#'
#' @examples
#' plot_color_continuous(physeq_obj, ordination.matrix, color_by)
#'
#' @export
plot_color_continuous = function(physeq_obj, ordination.matrix, color_by, axes = NULL) {
  # Define axes to plot
  if (is.null(axes)) {
    axes = 1:2 
  }

  # Get the name of the metric for the plot title
  metric = deparse(substitute(ordination.matrix))
  title = strsplit(metric, ".", fixed = TRUE)[[1]][1]
  title = paste(toupper(substr(title, 1, 1)), substr(title, 2, nchar(title)), sep = "")

  # Create ordination plot
  ordination.plot <- plot_ordination(
    physeq = physeq_obj,
    ordination = ordination.matrix,
    type = "samples",
    color = color_by,
    axes = axes
  ) +
    theme_bw() +
    geom_point(size = 3) +
    theme(legend.position = "top") +
    scale_color_viridis() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labels.x + labels.y + labels.legend

  return(ordination.plot)
}
#' Plot Color Discrete
#'
#' This function generates an ordination plot with discrete color mapping based on an input ordination matrix and color variable.
#'
#' @param physeq_obj A phyloseq object.
#' @param ordination.matrix Ordination matrix for plotting.
#' @param color_by Variable used for discrete color mapping.
#' @param axes Optional parameter specifying the axes to plot. Defaults to NULL, which plots the first two axes.
#'
#' @return An ordination plot with discrete color mapping.
#'
#' @import phyloseq
#' @import ggplot2
#' @import cowplot
#'
#' @examples
#' plot_color_discrete(physeq_obj, ordination.matrix, color_by)
#'
#' @export
plot_color_discrete = function(physeq_obj, ordination.matrix, color_by, axes = NULL) {
  # Define axes to plot
  if (is.null(axes)) {
    axes = 1:2 
  }

  # Get the name of the metric for the plot title
  metric = deparse(substitute(ordination.matrix))
  title = strsplit(metric, ".", fixed = TRUE)[[1]][1]
  title = paste(toupper(substr(title, 1, 1)), substr(title, 2, nchar(title)), sep = "")

  # Create ordination plot
  ordination.plot <- plot_ordination(
    physeq = physeq_obj,
    ordination = ordination.matrix,
    type = "samples",
    color = color_by,
    axes = axes
  ) +
    theme_bw() +
    geom_point(size = 3,alpha=0.6) +
    theme(legend.position = "top") +
    scale_color_npg() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
         axis.title = element_text(size=14),
         axis.text = element_text(size=12)) 


  return(ordination.plot)
}

#' Plot Grid PCoA Test
#'
#' This function generates a grid of density plots for PCoA analysis based on input ordination vectors.
#'
#' @param input.ord Input ordination object.
#' @param ord.plot Ordination plot object.
#' @param mapping_value Variable name for mapping the fill color of density plots.
#'
#' @return A grid of density plots.
#'
#' @import ggplot2
#' @import cowplot
#' @import dplyr
#'
#' @examples
#' plot_grid_pcoa_test(input.ord, ord.plot, "mapping_value")
#'@export

plot_grid_pcoa <- function(input.ord, ord.plot, mapping_value, metadata_table, colors_df) {
  fill_by <- mapping_value
  
  # Ordination vectors
  ord.vector <- as.data.frame(input.ord$vectors) %>%
    rownames_to_column(var = "Sample") %>%
    inner_join(metadata_table)
  
  # PCoA
  pcoa.plot <- ord.plot +
      theme(legend.position = "none",
             axis.title = element_text(size=14),
             axis.text = element_text(size=12)) 
    
    
  # Density plot for x axis
  axis1 <- ggplot(data = ord.vector, aes(x = Axis.1, fill = .data[[fill_by]])) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = colors_df$color,
                     limits = colors_df$category) +  # Set colors using the provided data frame
    theme_void() +
    guides(fill = "none") +
    theme(
      axis.title.x = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Density plot for y axis
  axis2 <- ggplot(data = ord.vector, aes(x = Axis.2, fill = .data[[fill_by]])) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = colors_df$color,
                     limits = colors_df$category) +  # Set colors using the provided data frame
    theme_void() +
    guides(fill = "none") +
    coord_flip() +
    theme(
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Subpanels
  jaccard.axis1 <- axis1 + theme(plot.margin = unit(c(0.5, -0.5, 0, 1.7), "cm"))
  jaccard.cat <- pcoa.plot + theme(plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
  jaccard.axis2 <- axis2 + theme(plot.margin = unit(c(-0.6, 0.5, 2.7, 0), "cm"))
  
  # Grids for each panel
  axis1 <- plot_grid(axis1, pcoa.plot, ncol = 1, rel_heights = c(1, 3))
  axis2 <- plot_grid(NULL, axis2, ncol = 1, rel_heights = c(1, 3))
  
  pcoa.plot.grid <- plot_grid(axis1, axis2, ncol = 2, rel_widths = c(3, 1))
  return(pcoa.plot.grid)
}

# Function: getOrdinationVectors
# Description: Retrieves ordination vectors and merges them with a metadata table.
# Parameters:
#   - ordination: Ordination object containing vectors.
#   - metadata_table: Metadata table to be joined with ordination vectors.
# Returns:
#   - ord.vector: Data frame containing ordination vectors merged with metadata table.

getOrdinationVectors = function(ordination, metadata_table) {
  # Convert ordination vectors to a data frame, add rownames as "Sample" column,
  # and perform an inner join with the metadata_table
  ord.vector = as.data.frame(ordination$vectors) %>%
    rownames_to_column(var = "Sample") %>%
    inner_join(metadata_table)

  return(ord.vector)
}

# Function: wilcoxonOrdinationAxis
# Description: Performs pairwise Wilcoxon rank sum tests on ordination axis values against categories.
# Parameters:
#   - ord.vectors: Data frame containing ordination vectors.
# Returns:
#   - List of Wilcoxon rank sum test results for each ordination axis.

wilcoxonOrdinationAxis = function(ord.vectors, test_value) {
  # Perform pairwise Wilcoxon rank sum tests on each ordination axis against categories
  wilcoxon.axis1 = pairwise.wilcox.test(ord.vectors$Axis.1, ord.vectors[[test_value]], p.adjust.method = "fdr")
  wilcoxon.axis2 = pairwise.wilcox.test(ord.vectors$Axis.2, ord.vectors[[test_value]], p.adjust.method = "fdr")
  wilcoxon.axis3 = pairwise.wilcox.test(ord.vectors$Axis.3, ord.vectors[[test_value]], p.adjust.method = "fdr")
  wilcoxon.axis4 = pairwise.wilcox.test(ord.vectors$Axis.4, ord.vectors[[test_value]], p.adjust.method = "fdr")

  df.results = data.frame(
    axis1 = wilcoxon.axis1$p.value,
    axis2 = wilcoxon.axis2$p.value,
    axis3 = wilcoxon.axis3$p.value,
    axis4 = wilcoxon.axis4$p.value
  )    
   #return(list(axis1 = wilcoxon.axis1, axis2 = wilcoxon.axis2, axis3 = wilcoxon.axis3, axis4 = wilcoxon.axis4))  
  return(df.results)
}
   



#' Calculates group dissimilarity
#'
#' This function generates a list of dataframes of group dissimilarity measures from a list of distance methods
#'
#' @param physeq_obj Phyloseq object
#' @param dist_methods A vector containing the distance methods
#' @param mapping_value Variable name for group comparisons
#'
#' @return A list of data frames
#' @output[[1]]$data Data frame with distance measures
#' @output[[1]]$name Distance metric
#'
#' @import gridExtra
#' @import ggplot
#' @import MetagMisc
#'
#' @examples
#' dist_methods <- c("jaccard", "bray")
#' output <- compute_dissimilarity_groups_df(physeq, dist_methods, mapping_value = "mapping_variable")
#'
#' Access the dissimilarity data frames for each method
#' jaccard_data <- output[[1]]$data
#' jaccard_name <- output[[1]]$name
  
#' bray_data <- output[[2]]$data
#' bray_name <- output[[2]]$name
#'@export

compute_dissimilarity_groups_df <- function(physeq_obj, dist_methods, mapping_value) {
  dissimilarity_groups <- lapply(dist_methods, function(method) {
    diss.df <- phyloseq_group_dissimilarity(physeq_obj, 
                                            group = mapping_value, 
                                            method = method, 
                                            justDF = TRUE)
    list(name=method,data=diss.df)  # Assign method name as the name of the data frame
      
    return(list(name = method, data = diss.df))
  })
  
  return(dissimilarity_groups)
}


### Plot a grid with group dissimilarities based on a list of distance metrics
library(gridExtra)

compute_dissimilarity_groups_plots <- function(physeq_obj, dist_methods, mapping_value) {
  plots <- lapply(dist_methods, function(method) {
    diss.plot <- phyloseq_group_dissimilarity(physeq_obj, group = mapping_value, method = method, method_title = TRUE) +
      labels.x +
      labels.y +
      labels.legend +
      theme(legend.position = "bottom") +
      scale_fill_npg() +
      theme(axis.text.x = element_blank())

    return(diss.plot)
  })

  grid.plots = grid.arrange(grobs = plots, ncol = 2)  # Arrange plots in a grid with 2 columns

  return(grid.plots)  # Return NULL as the plots are displayed using grid.arrange
}

################################################################################
#' Convert phyloseq OTU count data into DGEList for edgeR package
#' 
#' Further details.
#' 
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#'  an \code{\link{otu_table-class}} object. 
#'  The latter is only appropriate if \code{group} argument is also a 
#'  vector or factor with length equal to \code{nsamples(physeq)}.
#'  
#' @param group (Required). A character vector or factor giving the experimental
#'  group/condition for each sample/library. Alternatively, you may provide
#'  the name of a sample variable. This name should be among the output of
#'  \code{sample_variables(physeq)}, in which case
#'  \code{get_variable(physeq, group)} would return either a character vector or factor.
#'  This is passed on to \code{\link[edgeR]{DGEList}},
#'  and you may find further details or examples in its documentation.
#'  
#' @param method (Optional). The label of the edgeR-implemented normalization to use.
#'  See \code{\link[edgeR]{calcNormFactors}} for supported options and details. 
#'  The default option is \code{"RLE"}, which is a scaling factor method 
#'  proposed by Anders and Huber (2010).
#'  At time of writing, the \link[edgeR]{edgeR} package supported 
#'  the following options to the \code{method} argument:
#'  
#'  \code{c("TMM", "RLE", "upperquartile", "none")}.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#' 
#' @examples
#' 
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}
################################################################################

### Perform differntially abundance analysis using EdgeR ####
### @param phyloseq_obj phyloseq_class object
### @param varianceThreshold a double vector, default=1e-5
### @param groupToCompare A character vector or factor giving the experimental group/condition for each sample/library

runEdgeR = function(phyloseq_obj,varianceThreshold=NULL,groupToCompare){

varianceThreshold = 1e-5
keepOTUs = names(which(apply(otu_table(phyloseq_obj), 1, var) > varianceThreshold))
filter_taxa = prune_taxa(keepOTUs, phyloseq_obj)

dge = phyloseq_to_edgeR(filter_taxa, group=groupToCompare)


# Perform binary test
et = exactTest(dge)
    
# Extract values from test results
tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
    

#Filter significant results based on FDR
    alpha = 0.001
    sigtab = res[(res$FDR < alpha), ]
    

#    print(tt$adjust.method)
#    print(tt$test)

    sigtab=filter(res,FDR<0.01)
    
# Check if sigtab is empty --> THIS PART IS NOT WORKING YET
#  if(nrow(sigtab) == 0) {
#    print("There are no differentially abundant features.")
#    stop("NULL") # Return NULL or you could return an empty list or similar depending on your needs
#  }

    group1 = levels(factor(tt$comparison))[1]
    group2 = levels(factor(tt$comparison))[2]
 
    # Add 'EnrichedIn' column
    sigtab = sigtab %>%
    mutate(EnrichedIn=ifelse(logFC>0,group2,group1))

    
#Sort results
sigtabgen = subset(sigtab, !is.na(Genus))

      # Sort results for each taxonomic level
  tax_levels <- c("Phylum", "Family", "Genus", "Species")
  for (level in tax_levels) {
    if (level %in% colnames(sigtab)) {
      x <- tapply(sigtab$logFC, sigtab[, level], function(x) max(x))
      x <- sort(x, decreasing = TRUE)
      sigtab[, level] <- factor(as.character(sigtab[, level]), levels = names(x))
    }
  }


    return(list(sigtab.final = sigtab,sigtabgen = sigtabgen))

    }

#Create tables with cols and rows annotations
annotationsForHeatmap <- function(ps_obj){

      annotationTable =data.frame(
                            sample=ps_obj@sam_data$Sample,
                            biome1=ps_obj@sam_data$biome_1,
                            biome2=ps_obj@sam_data$biome_2,
                            category=ps_obj@sam_data$category,
                            ncbi_order=ps_obj@sam_data$ncbi_order,
                            mammals=ps_obj@sam_data$mammals)

        rownames(annotationTable)=annotationTable$sample
    
    annotation.col = select(annotationTable,-sample)
    
    annotation.row=data.frame(ps_obj@tax_table)
    
    return(list(annotation.col = annotation.col,annotation.row = annotation.row))

}

#Create matrix for heatmap, applying log transformation to enhance visualization
createMatrixForHetmap = function(ps_obj){

matrix=as.matrix(ps_obj@otu_table[,2:ncol(ps_obj@otu_table)])
rownames(matrix)=rownames(ps_obj@otu_table)

matrixLog<-log((matrix*1000)+1)

    return(matrixLog)
}

#Function to plot specificing ordination matrix and mapping value
plot_color_by = function(physeq_obj,ordination.matrix,color_by){
    #Get the name of the metric for the plot title
    metric = deparse(substitute(ordination.matrix))
    title = strsplit(metric,".",fixed=TRUE)[[1]][1]
    title = paste(toupper(substr(title,1,1)),substr(title,2,nchar(title)),sep = "")
    #str_to_title(title) --> doesn't work, stringr package corrupt, meanwhile use another option
    #Create ordination plot
    ordination.plot  = plot_ordination(
                 physeq = physeq_obj,
                 ordination = ordination.matrix,
                 type = "samples",
                 color = color_by)+
                theme_bw()+
                geom_point(size = 3)+
                theme(legend.position = "top")+
                scale_color_npg() +
                ggtitle(title)+
                theme(plot.title = element_text(hjust = 0.5))+
                labels.x + labels.y + labels.legend
return(ordination.plot)

}

#Function to assign ncbi taxonomy from a metadata table containing the ncbi taxid

           assignTaxonomy = function(metadata.table){
               
library(taxonomizr)               
    accesionTaxa.db = "/ebio/abt3_projects2/databases_no-backup/NCBI_accession2taxid/accessionTaxa.sql"
    taxonomyTable = metadata.table %>%
                    mutate(ncbi_taxid=metadata.table$host_taxid)%>%
                  mutate(ncbi_phylum=getTaxonomy(metadata.table$host_taxid,accesionTaxa.db,desiredTaxa=c("phylum")))%>%
                  mutate(ncbi_class=getTaxonomy(metadata.table$host_taxid,accesionTaxa.db,desiredTaxa=c("class")))%>%
                  mutate(ncbi_order=getTaxonomy(metadata.table$host_taxid,accesionTaxa.db,desiredTaxa=c("order")))%>%
                  mutate(ncbi_family=getTaxonomy(metadata.table$host_taxid,accesionTaxa.db,desiredTaxa=c("family")))%>%
                  mutate(ncbi_genus=getTaxonomy(metadata.table$host_taxid,accesionTaxa.db,desiredTaxa=c("genus")))%>%
                  mutate(ncbi_species=getTaxonomy(metadata.table$host_taxid,accesionTaxa.db,desiredTaxa=c("species")))
    
    return(taxonomyTable)
}
           
# Function: getTaxSummary
# Description: Generates a summary table of taxonomic information from a given taxonomy data frame.
# Parameters:
#   - DA.taxonomy: The taxonomy data frame containing taxonomic information.
# Returns: 
#   - A tibble summarizing the taxonomic information.
           
getTaxSummary = function(DA.taxonomy){

        # Calculate the summary statistics
        taxSummary = tibble(EnrichedIn = DA.taxonomy$EnrichedIn[1],
                    totalEnriched = nrow(DA.taxonomy),
                    NoPhyla = length(unique(DA.taxonomy$Phylum)),
                    NoClass = length(unique(DA.taxonomy$Class)),
                    NoOrder = length(unique(DA.taxonomy$Order)),
                    NoFamily = length(unique(DA.taxonomy$Family)),
                    NoGenus = length(unique(DA.taxonomy$Genus)),
                    NoSpecies = length(unique(DA.taxonomy$Species)))

    # Return the taxonomic summary
    return(taxSummary)
}
           
# FilterSingleRepresentative: A function to filter rows with single representatives based on a specific column in an input table.

# Parameters:
#   input_table: A data frame representing the input table containing the data.
#   filter_column: A character vector specifying the column name to group and filter the data.

# Returns:
#   output_table: A data frame containing only the rows with single representatives for each unique value in the specified column.

FilterSingleRepresentative <- function(input_table, filter_column, subsample=FALSE, subsamplingDepth=NULL) {
  # Group the data by the specified column and filter rows with single representatives
  output_table <- input_table %>%
    group_by({{filter_column}}) %>%
    filter(n() == 1) %>%
    ungroup()
    
  # Check if subsampling is requested
  if (subsample) {
    # Check if the subsamplingDepth is provided and is valid
    if (is.null(subsamplingDepth) || subsamplingDepth <= 0 || subsamplingDepth >= nrow(output_table)) {
      stop("Invalid subsamplingDepth value. It should be a positive integer less than the number of rows in the filtered output_table.")
    }
    
    # Perform subsampling using dplyr's sample_n function
    output_table <- sample_n(output_table, subsamplingDepth)
  }
  
  # Print the number of rows in the output table
  print(nrow(output_table))
  
  # Return the filtered (and possibly subsampled) output table
  return(output_table)
}
           
           
## Export sequences to fasta
           
export_sequences_to_fasta <- function(sequence_fasta, seqList, category, output_file) {
  # Extract the relevant sequences based on the given category
  derep_category <- sequence_fasta[c(which(names(sequence_fasta) %in% seqList[[category]]$Accession))]
  
  # Write the sequences to a FASTA file
  write.fasta(sequences = derep_category, names = names(derep_category), nbchar = 80, file.out = output_file)
}
           
           
           
#Get sequence names from fasta file containing cds
#Input: fasta with cds sequences whose header has the cdsName_type_Accession in the first column of the header
#This function splits the first column of the header and create a mapping file with three columns: cds_name, type=cds, and protein accession
#Output: A data frame mapping cds_name to protein accessions
getSeqNames = function(fasta_file){
#     seq_names = getName(fasta_file)
    #Split the sequence names into three columns based on underscores
                seq_table = tibble(seqName=getName(fasta_file))%>%
                      separate(col=seqName,into=c("cds_name","type","Accession"),sep="_",remove="FALSE")
    return(seq_table)
    }

           

## Get counts matrix from a phyloseq object, based on the category value           
createCountsMatrix = function(ps_object){
    OTU1 = as(otu_table(ps_object),"matrix")
    #if(taxa_are_rows(ps_object)){OTU1 <- t(OTU1)}
    OTUdf = as.data.frame(OTU1)
    OTUdf = rownames_to_column(OTUdf,var="Accession")
    
    #Read metadata table from phyloseq object
    metadata = as_tibble(ps_object@sam_data)
    category.value = metadata$category[1]
    
    #Read taxonomy table from phyloseq object
    #Finish this section: Create the vector with corresponding taxonomy for each accession, to further add to the output.table
    #taxonomy = as_tibble(ps_object@tax_table) %>%
     #          rownames_to_column(var = "Accession")
#    species = select(taxonomy,Species)
    
    #Create table and assign category
    output.table = mutate(OTUdf,totalCount= rowSums(across(where(is.numeric))))%>%
                mutate(category = category.value) %>%
                mutate(NoSamples = nrow(metadata)) %>%
                mutate(ProportionInCategory = totalCount/NoSamples) %>%
                mutate(PresentInCategory = ifelse(ProportionInCategory>0,TRUE,FALSE))%>%
                select(Accession,category,totalCount,NoSamples,ProportionInCategory,PresentInCategory)%>%
                rename_with(~paste0("ProportionIn_",category.value),ProportionInCategory)%>%
                rename_with(~paste0("PresentIn_",category.value),PresentInCategory)%>%
                rename_with(~paste0("totalCount_",category.value),totalCount)

    #mutate(category.table,ExtractCategorypresent = ifelse(rowSums(select(.,-c(Sample,category)))>0,TRUE,FALSE))
    
    return(output.table)
}
           
#### Runs a PERMANOVA for a list of distance methods and defining a grouping variable
### This functions runs for default 999 permutations
          
           # Define the function with an additional argument for the grouping variable
run_permanova <- function(physeq_object, methods_list, grouping_variable) {
  # Initialize an empty list to store PERMANOVA results
  permanova_results <- list()
  
  # Loop over each method
  for (method in methods_list) {
    # Calculate distance matrix
    distance_matrix <- distance(physeq_object, method = method)
    
    # Run PERMANOVA using the specified grouping variable
    formula <- as.formula(paste("distance_matrix ~", grouping_variable))
    permanova_result <- adonis2(formula, 
                                data = as(sample_data(physeq_object), "data.frame"),
                                permutations = 999)
    
    # Add the results to the list
    permanova_results[[method]] <- permanova_result
  }
  
  # Return the list of PERMANOVA results
  return(permanova_results)
}
                  
                  
                  
                  
## This function creates a taxonomic counts summary per taxonomic hierarchies in a table with results of pariwise comparison in edgeR 
## Input: table sigtabgen produced with the function runEdgeR
## Input: table with the two levels of pairwise comparison and a column with 'EnrichedIn' category
                  
## TO FIX: levels for groups when the input table is empty (no diff. abundant found)
                  
getTaxSummary_pairwiseComparisons = function(DA.taxonomy){

        # Calculate the summary statistics
    group1 = levels(factor(DA.taxonomy$EnrichedIn))[1]
    group2 = levels(factor(DA.taxonomy$EnrichedIn))[2]
    tax_levels <- c("Phylum","Class","Order","Family", "Genus", "Species")
    
#    for (level in tax_levels){
        
        #Create a column with the count of unique values of each level, then append it to the table
 #       mutate(level = length(unique(level)))
        
  #  }
    
    tmp.group1 = filter(DA.taxonomy,EnrichedIn==paste(group1))
    taxSummary_group1 = tibble(EnrichedIn = tmp.group1$EnrichedIn[1],
                        totalEnriched = nrow(tmp.group1),
                        NoPhyla = length(unique(tmp.group1$Phylum)),
                        NoClass = length(unique(tmp.group1$Class)),
                        NoOrder = length(unique(tmp.group1$Order)),
                        NoFamily = length(unique(tmp.group1$Family)),
                        NoGenus = length(unique(tmp.group1$Genus)),
                        NoSpecies = length(unique(tmp.group1$Species))
                               )
    
    tmp.group2 = filter(DA.taxonomy,EnrichedIn==paste(group2))
    taxSummary_group2 = tibble(EnrichedIn = tmp.group2$EnrichedIn[1],
                        totalEnriched = nrow(tmp.group2),
                        NoPhyla = length(unique(tmp.group2$Phylum)),
                        NoClass = length(unique(tmp.group2$Class)),
                        NoOrder = length(unique(tmp.group2$Order)),
                        NoFamily = length(unique(tmp.group2$Family)),
                        NoGenus = length(unique(tmp.group2$Genus)),
                        NoSpecies = length(unique(tmp.group2$Species))
                               )
    
       
    taxSummary = bind_rows(taxSummary_group1,taxSummary_group2)%>%
                    mutate(comparison = paste(group1,"_",group2))
    # Return the taxonomic summary
    return(taxSummary)
}
    