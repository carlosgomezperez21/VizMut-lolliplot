library(optparse)

#---------------------------
# Argumentos
#---------------------------
option_list <- list(

  make_option("--variants",
              type = "character",
              default = NULL,
              help = "Ruta al CSV de variantes [requerido]"),

  make_option("--features",
              type = "character",
              default = NULL,
              help = "Ruta al CSV de features [requerido]"),

  make_option("--gene",
              type = "character",
              default = "KMT2B",
              help = "Nombre del gen [default: %default]"),

  make_option("--protein_length",
              type = "integer",
              default = 2715,
              help = "Longitud de la proteína en aa [default: %default]"),

  make_option("--output",
              type = "character",
              default = "output/lolliplot.png",
              help = "Ruta del plot de salida [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

#---------------------------
# Validación de argumentos
#---------------------------
if (is.null(opt$variants)) stop("--variants es requerido")
if (is.null(opt$features)) stop("--features es requerido")

if (!file.exists(opt$variants)) stop(paste("No se encuentra el archivo:", opt$variants))
if (!file.exists(opt$features)) stop(paste("No se encuentra el archivo:", opt$features))

#---------------------------
# Cargar funciones
#---------------------------
source("R/parse_variants.R")
source("R/parse_features.R")
source("R/validate.R")
source("R/plot_lolliplot.R")

#---------------------------
# Metadata
#---------------------------
meta <- list(
  gene = list(
    name   = opt$gene,
    chr    = "chr19",
    start  = 35718003,
    end    = 35738878,
    strand = "+"
  ),
  protein = list(
    id     = opt$gene,
    length = opt$protein_length
  )
)

#---------------------------
# Pipeline
#---------------------------
message("Leyendo datos...")
variants_raw <- read.csv(opt$variants)
features_raw <- read.csv(opt$features)

message("Parseando variantes...")
variants <- parse_variants(variants_raw)

message("Parseando features...")
features <- parse_features(features_raw)

message("Validando variantes...")
validate_variants(variants, meta)

message("Generando plot...")
p <- plot_protein_lolliplot(variants, features, meta)

message("Guardando output en: ", opt$output)
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
ggsave(opt$output, p, width = 14, height = 6, dpi = 150)

message("¡Listo!")
