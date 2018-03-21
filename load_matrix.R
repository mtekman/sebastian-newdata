suppressWarnings(library(SingleCellExperiment))
suppressWarnings(library(scater))

# Set the raw data and save/load it
umi <- NULL
if (!exists("generate_matrix")){
   generate_matrix = F
}

if (!exists("matrix_destination")){
   matrix_destination = "matrix.rds"
}


#genes_of_interest <- c(
#    "Eomes", "Brachyury", "Mesp1",                # meso
#    "Pou5f1", "nanog",                            # pluripotent
#    "Sox1", "Sox2", "Pou3f1", "zfp462", "slc7a3"  # neuroectoderm
#) 

if (is.null(umi) && generate_matrix){
    message("Regenerating matrix from source file ", matrix_destination ,"...")
    count_matrix <- read.table("source/E725.matrix.Seb_NewData_E725.3.quantif", 
                               sep="\t", 
                               comment.char = '*',   # important!
                               header = T, 
                               row.names = 1)

    umi <- as.matrix(count_matrix)
    sce <- SingleCellExperiment(assays = list(counts = umi))

    # Create SCE and annotate
    # Assign known/related groups
    #is.meso <- rownames(sce) %in% c("Eomes", "Brachyury", "Mesp1")
    #is.pluri <- rownames(sce) %in% c("Pou5f1", "nanog")
    #is.neuro <- rownames(sce) %in% c("Sox1", "Sox2", "Pou3f1", "zfp462", "slc7a3")
    #rowData(sce)$is_mesoderm <- is.meso
    #rowData(sce)$is_pluripotenz <- is.pluri
    #rowData(sce)$is_neuroectoderm <- is.neuro

    sce <- getBMFeatureAnnos(
        sce, 
        filters = "ensembl_gene_id",
        attributes = c(
            "ensembl_gene_id",              # Gene stable ID
            "external_gene_name",           # Casual name
            "external_transcript_name",     # Transcript-specific name
            "gene_biotype",                 # Gene biotype
            "transcript_biotype",           # Trans type
            "description",                  # Gene description
            "band",                         # Karyotype band
            "refseq_mrna",
            "go_id",                        # Go Term accession (cellular domains)
            "go_linkage_type",              # Go Term evidence code
            "name_1006",                    # Go Term name
            "definition_1006",              # Go Term definition
            "namespace_1003"                # Go domain             
        ),
        feature_symbol = "mgi_symbol",
        feature_id = "ensembl_gene_id",
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "mmusculus_gene_ensembl", 
        host = "www.ensembl.org"
    )

    saveRDS(sce, matrix_destination)

} else {
    message("Loading matrix from source file ", matrix_destination, "...")
    sce <- readRDS(matrix_destination)
}


message(dim(sce)[1], " genes x ", dim(sce)[2], " cells. (", 
        length(unique(colnames(sce))), ") unique barcodes.")
