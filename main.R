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
              help    = "Ruta del plot de salida [default: %default]")
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
variants     <- parse_variants(variants_raw)

#---------------------------
# Plot tipo: protein
#---------------------------
if (plot_type == "protein") {

  if (is.null(opt$features)) stop("--features es requerido para --plot_type protein")
  if (!file.exists(opt$features)) stop("No se encuentra: ", opt$features)
  if (is.null(opt$protein_length)) stop("--protein_length es requerido para --plot_type protein")

  source("R/parse_features.R")
  source("R/plot_lolliplot.R")

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
  p <- plot_protein_lolliplot(variants, features, meta)
}

#---------------------------
# Plot tipo: single_gene
#---------------------------
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

  message("Obteniendo estructura del transcrito: ", opt$transcript_id)
  struct <- fetch_transcript_structure(opt$transcript_id)

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

  message("Generando plot genomico...")
  p <- plot_gene_lolliplot(
    variants             = variants_in,
    transcript_structure = struct,
    gene_name            = opt$gene_name,
    label_type           = opt$label_type,
    grid                 = opt$grid
  )
}

#---------------------------
# Plot tipo: multi_gene (Fase 3)
#---------------------------
if (plot_type == "multi_gene") {
  stop("--plot-type multi_gene estara disponible en la Fase 3")
}

#---------------------------
# Guardar output
#---------------------------
message("Guardando output en: ", opt$output)
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(opt$output, p, width = 14, height = 6, dpi = 150)
message("Listo!")
