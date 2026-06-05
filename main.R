library(optparse)

#---------------------------
# Argumentos
#---------------------------
option_list <- list(

  make_option("--variants",
              type    = "character",
              default = NULL,
              help    = "Ruta al CSV de variantes [requerido]"),

  make_option("--features",
              type    = "character",
              default = NULL,
              help    = "Ruta al CSV de features (requerido para --plot-type protein)"),

  make_option("--plot_type",
              type    = "character",
              default = "protein",
              help    = "Tipo de plot: protein | single_gene | multi_gene [default: %default]"),

  make_option("--gene_name",
              type    = "character",
              default = "GENE",
              help    = "Nombre del gen [default: %default]"),

  make_option("--protein_length",
              type    = "integer",
              default = NULL,
              help    = "Longitud de la proteina en aa (requerido para --plot-type protein)"),

  make_option("--transcript_id",
              type    = "character",
              default = NULL,
              help    = "RefSeq transcript ID (ej. NM_014727.3) para single_gene y multi_gene"),

make_option("--gene_list",
              type    = "character",
              default = "all",
              help    = "Genes a plotear: 'all', lista 'KMT2B,DNMT3A' o ruta a archivo .txt [default: %default]"), 

 make_option("--genome",
              type    = "character",
              default = "hg38",
              help    = "Ensamble genomico del CSV de variantes: hg38 | hg19 [default: %default]"),

  make_option("--grid",
              type    = "logical",
              default = FALSE,
              help    = "Dividir plot por variant_type: TRUE | FALSE [default: %default]"),

  make_option("--label_type",
              type    = "character",
              default = "protein",
              help    = "Tipo de etiqueta en lollipops: protein | cds [default: %default]"),

  make_option("--output",
              type    = "character",
              default = "output/lolliplot.png",
              help    = "Ruta del plot de salida [default: %default]"),
  
  make_option("--enrich",
              type    = "logical",
              default = FALSE,
              help    = "Enriquecer variantes desde ClinVar y NCBI [default: %default]"),
  make_option("--enrich_output",
              type    = "character",
              default = NULL,
              help    = "Ruta para guardar CSV enriquecido (ej. output/enriched.csv) [opcional]"),

  make_option("--exons",
              type    = "character",
              default = NULL,
              help    = "Exones a mostrar: numeros o rangos ej. '1,5-7,10' [default: todos]"),

  make_option("--count",
              type    = "logical",
              default = FALSE,
              help    = "Deduplicar variantes y escalar lollipops por frecuencia [default: %default]")

)

opt <- parse_args(OptionParser(option_list = option_list))

#---------------------------
# Validacion general
#---------------------------
if (is.null(opt$variants)) stop("--variants es requerido")
if (!file.exists(opt$variants)) stop("No se encuentra: ", opt$variants)

plot_type <- opt$plot_type

if (!plot_type %in% c("protein", "single_gene", "multi_gene")) {
  stop("--plot_type debe ser: protein | single_gene | multi_gene")
}

if (!opt$genome %in% c("hg38", "hg19")) {
  stop("--genome debe ser: hg38 | hg19")
}

#---------------------------
# Cargar funciones comunes
#---------------------------
source("R/parse_variants.R")
source("R/validate.R")

message("Leyendo variantes...")
variants_raw <- read.csv(opt$variants)

#---------------------------
# Enriquecimiento opcional
#---------------------------
if (opt$enrich) {
  message("Enriqueciendo variantes desde ClinVar y NCBI...")
  source("R/enrich_variants.R")
  variants <- enrich_variants(variants_raw, genome = opt$genome)
  variants$protein_pos <- as.numeric(
    sub("p\\.[A-Za-z]+([0-9]+).*", "\\1", variants$protein_change)
  )
} else {
  variants <- parse_variants(variants_raw)
}

# exportar CSV enriquecido si se especifica
  if (!is.null(opt$enrich_output)) {
    dir.create(dirname(opt$enrich_output),
               showWarnings = FALSE, recursive = TRUE)
    
    ext <- tools::file_ext(opt$enrich_output)
    sep <- switch(ext,
                  "tsv" = "\t",
                  "txt" = "\t",
                  ","
    )
    write.table(variants, opt$enrich_output,
                sep = sep, row.names = FALSE, quote = FALSE)
    message("CSV enriquecido guardado en: ", opt$enrich_output)
  }


#---------------------------
# Deduplicacion opcional
#---------------------------
if (opt$count) {
  message("Deduplicando variantes...")
  source("R/dedup_variants.R")
  variants <- dedup_variants(variants)
}

#---------------------------
# Plot tipo: protein
#---------------------------
if (plot_type == "protein") {

  if (is.null(opt$features)) stop("--features es requerido para --plot_type protein")
  if (!file.exists(opt$features)) stop("No se encuentra: ", opt$features)
  if (is.null(opt$protein_length)) stop("--protein_length es requerido para --plot_type protein")

  source("R/parse_features.R")
  source("R/plot_protein_lolliplot.R")

  message("Parseando features...")
  features_raw <- read.csv(opt$features)
  features     <- parse_features(features_raw)

  meta <- list(
    gene    = list(name   = opt$gene_name,
                   chr    = unique(variants$chr)[1],
                   start  = min(variants$pos, na.rm = TRUE),
                   end    = max(variants$pos, na.rm = TRUE),
                   strand = "+"),
    protein = list(id     = opt$gene_name,
                   length = opt$protein_length)
  )

  message("Validando variantes...")
  validate_variants(variants, meta)

  message("Generando plot proteico...")
  p <- plot_protein_lolliplot(variants, features, meta, grid = opt$grid)
}

#---------------------------
# Plot tipo: single_gene
#---------------------------
# en plots genomicos el default de label es cds
  if (opt$label_type == "protein") {
    message("Nota: usando label_type=cds para plot genomico")
    opt$label_type <- "cds"
  }

if (plot_type == "single_gene") {

  if (is.null(opt$transcript_id)) {
    stop("--transcript_id es requerido para --plot_type single_gene")
  }

  source("R/fetch_transcript.R")
  source("R/plot_gene_lolliplot.R")

  # liftover si datos en hg19
  if (opt$genome == "hg19") {
    message("Detectado --genome hg19, realizando liftover a hg38...")
    source("R/hgvs_to_genomic.R")

    library(rtracklayer)
    library(GenomicRanges)

    chain_file <- "/tmp/hg19ToHg38.over.chain"
    chain_gz   <- paste0(chain_file, ".gz")

    if (!file.exists(chain_file)) {
      message("Descargando chain file...")
      download.file(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        chain_gz, quiet = TRUE
      )
      R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
    }

    chain   <- import.chain(chain_file)
    valid   <- !is.na(variants$pos)
    gr_hg19 <- GRanges(
      seqnames = variants$chr[valid],
      ranges   = IRanges(start = variants$pos[valid],
                         end   = variants$pos[valid])
    )
    gr_hg38 <- liftOver(gr_hg19, chain)
    lifted  <- sapply(gr_hg38, function(x) if (length(x) == 1) start(x) else NA)
    variants$pos[valid] <- lifted
    message("Liftover completado")
  }

# obtener citobandas para ideograma
message("Obteniendo estructura del transcrito: ", opt$transcript_id)
  struct <- fetch_transcript_structure(opt$transcript_id)

  # filtrar exones si se especifica --exons
  source("R/parse_exon_selection.R")
  exon_sel    <- parse_exon_selection(opt$exons)
  exon_filter <- filter_exons(struct, exon_sel)
  struct      <- exon_filter$structure
  breaks      <- exon_filter$breaks

  # citobandas
  message("Obteniendo citobandas...")
  source("R/fetch_cytobands.R")
  source("R/plot_ideogram.R")
  chr   <- unique(struct$chr)[1]
  bands <- tryCatch(fetch_cytobands(chr), error=function(e) NULL)

  # filtrar variantes dentro del transcrito
  tx_start    <- min(struct$start)
  tx_end      <- max(struct$end)
  variants_in <- variants[!is.na(variants$pos) &
                           variants$pos >= tx_start &
                           variants$pos <= tx_end, ]

  n_out <- nrow(variants) - nrow(variants_in)
  if (n_out > 0) {
    message("Nota: ", n_out, " variantes fuera del transcrito fueron excluidas")
  }

  # excluir variantes fuera de los exones/intrones seleccionados
  if (!is.null(exon_sel)) {
    n_before <- nrow(variants_in)

    # obtener rangos validos (exones seleccionados + intrones entre consecutivos)
    valid_ranges <- exon_filter$structure %>%
      filter(type %in% c("exon", "intron")) %>%
      select(start, end)

    # filtrar variantes que caen en alguno de esos rangos
    variants_in <- variants_in[!is.na(variants_in$pos) &
      sapply(variants_in$pos, function(p) {
        any(p >= valid_ranges$start & p <= valid_ranges$end)
      }), ]

    n_excl <- n_before - nrow(variants_in)
    if (n_excl > 0) {
      message("Nota: ", n_excl,
              " variantes excluidas por caer fuera de los exones seleccionados")
    }
  }

message("Generando plot genomico...")
  p <- plot_gene_lolliplot(
    variants             = variants_in,
    transcript_structure = struct,
    gene_name            = opt$gene_name,
    label_type           = opt$label_type,
    grid                 = opt$grid,
    cytobands            = bands,
    breaks               = breaks
  )
}

#---------------------------
# Plot tipo: multi_gene (Fase 3)
#---------------------------

if (plot_type == "multi_gene") {

  source("R/plot_multi_gene.R")

  #---------------------------
  # parsear gene_list
  #---------------------------
  if (is.null(opt$gene_list) || opt$gene_list == "all") {
    gene_list <- "all"
  } else if (file.exists(opt$gene_list)) {
    message("Leyendo lista de genes desde: ", opt$gene_list)
    gene_list <- trimws(readLines(opt$gene_list))
    gene_list <- gene_list[nchar(gene_list) > 0]
    message("Genes en archivo: ", paste(gene_list, collapse=", "))
  } else {
    gene_list <- trimws(strsplit(opt$gene_list, ",")[[1]])
    message("Genes especificados: ", paste(gene_list, collapse=", "))
  }

  # liftover si datos en hg19
  if (opt$genome == "hg19") {
    message("Detectado --genome hg19, realizando liftover a hg38...")
    library(rtracklayer)
    library(GenomicRanges)

    chain_file <- "/tmp/hg19ToHg38.over.chain"
    chain_gz   <- paste0(chain_file, ".gz")

    if (!file.exists(chain_file)) {
      message("Descargando chain file...")
      download.file(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        chain_gz, quiet = TRUE
      )
      R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
    }

    chain   <- import.chain(chain_file)
    valid   <- !is.na(variants$pos)
    gr_hg19 <- GRanges(
      seqnames = variants$chr[valid],
      ranges   = IRanges(start = variants$pos[valid],
                         end   = variants$pos[valid])
    )
    gr_hg38 <- liftOver(gr_hg19, chain)
    lifted  <- sapply(gr_hg38, function(x) if (length(x) == 1) start(x) else NA)
    variants$pos[valid] <- lifted
    message("Liftover completado")
  }

  message("Generando plot multi-gen...")
  p <- plot_multi_gene(
    variants   = variants,
    gene_list  = gene_list,
    genome     = opt$genome,
    grid       = opt$grid,
    label_type = opt$label_type
  )
}

#---------------------------
# Guardar output
#---------------------------
message("Guardando output en: ", opt$output)
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

# alto dinamico segun plot_type y numero de genes
plot_height <- switch(plot_type,
  "protein"     = 6,
  "single_gene" = 7,
  "multi_gene"  = {
    n_genes <- if (identical(gene_list, "all")) {
      length(unique(variants$gene))
    } else {
      length(gene_list)
    }
    n_genes * 7
  }
)

ggplot2::ggsave(opt$output, p,
                width  = 14,
                height = plot_height,
                dpi    = 300,
                limitsize = FALSE)
