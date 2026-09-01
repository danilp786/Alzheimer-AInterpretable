#!/usr/bin/env bash

###############################################################################
#PARÁMETROS QUE DEBE ADAPTAR EL USUARIO
#---------------------------------------
#
#--trim-length
#    Longitud final utilizada por Deblur.
#    Debe determinarse a partir de la distribución de longitudes
#    después del merge.
#
#--primer-error-rate
#    Tolerancia de error para la identificación de primers.
#    Ajustar según la calidad y variabilidad esperada de los primers.
#
#--threads
#    Número de hilos disponibles en el sistema.
#
#--classifier
#    Clasificador SILVA compatible con la región amplificada.
###############################################################################

set -euo pipefail


###############################################################################
# PIPELINE 16S V3-V4: SRA -> ASVs (Deblur) -> SILVA 138 taxonomy
#
# Flujo:
#   FASTQ paired-end
#       ↓
#   Importación QIIME 2
#       ↓
#   Eliminación de primers 341F / 805R
#       ↓
#   Eliminación de pares sin primers
#       ↓
#   Merge R1 + R2
#       ↓
#   Deblur
#       ↓
#   ASVs
#       ↓
#   Clasificación taxonómica SILVA
#
# Uso:
#   ./pipeline.sh \
#       --input /ruta/fastq \
#       --metadata /ruta/SraRunTable.csv \
#       --classifier /ruta/silva_classifier.qza \
#       --trim-length 250 \
#       --output /ruta/resultado
#
###############################################################################

THREADS=8

INPUT=""
METADATA=""
CLASSIFIER=""
TRIM_LENGTH=""
OUTPUT=""

###############################################################################
# ARGUMENTOS
###############################################################################

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            INPUT="$2"
            shift 2
            ;;

        --metadata)
            METADATA="$2"
            shift 2
            ;;

        --classifier)
            CLASSIFIER="$2"
            shift 2
            ;;

        --trim-length)
            TRIM_LENGTH="$2"
            shift 2
            ;;

        --output)
            OUTPUT="$2"
            shift 2
            ;;

        --threads)
            THREADS="$2"
            shift 2
            ;;

        *)
            echo "Argumento desconocido: $1"
            exit 1
            ;;
    esac
done

###############################################################################
# COMPROBACIÓN DE ARGUMENTOS
###############################################################################

if [[ -z "$INPUT" ||
      -z "$METADATA" ||
      -z "$CLASSIFIER" ||
      -z "$TRIM_LENGTH" ||
      -z "$OUTPUT" ]]; then

    echo ""
    echo "Uso:"
    echo "./pipeline.sh \\"
    echo "    --input /ruta/fastq \\"
    echo "    --metadata /ruta/SraRunTable.csv \\"
    echo "    --classifier /ruta/silva_classifier.qza \\"
    echo "    --trim-length 250 \\"
    echo "    --output /ruta/resultado \\"
    echo "    --threads 8"
    echo ""
    exit 1
fi

###############################################################################
# VALIDACIÓN DE PARÁMETROS
###############################################################################

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "ERROR: --threads debe ser un número entero mayor que 0."
    exit 1
fi

if ! [[ "$TRIM_LENGTH" =~ ^[0-9]+$ ]] || [[ "$TRIM_LENGTH" -lt 1 ]]; then
    echo "ERROR: --trim-length debe ser un número entero mayor que 0."
    exit 1
fi

###############################################################################
# DIRECTORIOS
###############################################################################

mkdir -p "$OUTPUT"

mkdir -p "$OUTPUT/01_import"
mkdir -p "$OUTPUT/02_quality"
mkdir -p "$OUTPUT/03_primers"
mkdir -p "$OUTPUT/04_merge"
mkdir -p "$OUTPUT/05_deblur"
mkdir -p "$OUTPUT/06_taxonomy"

###############################################################################
# COMPROBACIONES
###############################################################################

command -v qiime >/dev/null 2>&1 || {
    echo "ERROR: QIIME 2 no está disponible."
    exit 1
}

[[ -d "$INPUT" ]] || {
    echo "ERROR: No existe el directorio FASTQ: $INPUT"
    exit 1
}

[[ -f "$METADATA" ]] || {
    echo "ERROR: No existe el metadata: $METADATA"
    exit 1
}

[[ -f "$CLASSIFIER" ]] || {
    echo "ERROR: No existe el clasificador SILVA:"
    echo "$CLASSIFIER"
    exit 1
}

###############################################################################
# INFORMACIÓN DE LA EJECUCIÓN
###############################################################################

echo ""
echo "============================================================"
echo "CONFIGURACIÓN"
echo "============================================================"
echo "Input:          $INPUT"
echo "Metadata:       $METADATA"
echo "Classifier:     $CLASSIFIER"
echo "Trim length:    $TRIM_LENGTH"
echo "Threads:        $THREADS"
echo "Output:         $OUTPUT"
echo "============================================================"
echo ""

{
    echo "date=$(date -Iseconds)"
    echo "hostname=$(hostname)"
    echo "threads=$THREADS"
    echo "trim_length=$TRIM_LENGTH"
    echo "input=$INPUT"
    echo "metadata=$METADATA"
    echo "classifier=$CLASSIFIER"
    qiime info
} > "$OUTPUT/run-info.txt"

###############################################################################
# 1. IMPORTACIÓN DE FASTQ
###############################################################################

echo ""
echo "============================================================"
echo "1. IMPORTACIÓN DE FASTQ"
echo "============================================================"

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$INPUT" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "$OUTPUT/01_import/demux-paired.qza"

###############################################################################
# 2. CONTROL DE CALIDAD INICIAL
###############################################################################

echo ""
echo "============================================================"
echo "2. CONTROL DE CALIDAD INICIAL"
echo "============================================================"

qiime demux summarize \
    --i-data "$OUTPUT/01_import/demux-paired.qza" \
    --o-visualization "$OUTPUT/02_quality/demux-summary.qzv"

###############################################################################
# 3. ELIMINACIÓN DE PRIMERS
#
# 341F:
#   CCTACGGGNGGCWGCAG
#
# 805R:
#   GACTACHVGGGTATCTAATCC
#
# Se eliminan únicamente los pares donde se hayan encontrado
# los primers.
###############################################################################

echo ""
echo "============================================================"
echo "3. ELIMINACIÓN DE PRIMERS"
echo "============================================================"

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "$OUTPUT/01_import/demux-paired.qza" \
    --p-front-f CCTACGGGNGGCWGCAG \
    --p-front-r GACTACHVGGGTATCTAATCC \
    --p-error-rate 0.1 \
    --p-match-adapter-wildcards \
    --p-discard-untrimmed \
    --o-trimmed-sequences "$OUTPUT/03_primers/primer-trimmed.qza" \
    --verbose

###############################################################################
# 4. MERGE R1 + R2
###############################################################################

echo ""
echo "============================================================"
echo "4. MERGE DE READS FORWARD + REVERSE"
echo "============================================================"

qiime vsearch merge-pairs \
    --i-demultiplexed-seqs "$OUTPUT/03_primers/primer-trimmed.qza" \
    --o-merged-sequences "$OUTPUT/04_merge/merged.qza"

qiime demux summarize \
    --i-data "$OUTPUT/04_merge/merged.qza" \
    --o-visualization "$OUTPUT/04_merge/merged-summary.qzv"

###############################################################################
# 5. DEBLUR
#
# TRIM_LENGTH se proporciona mediante:
#
#   --trim-length 250
#
###############################################################################

echo ""
echo "============================================================"
echo "5. DEBLUR"
echo "============================================================"

qiime deblur denoise-16S \
    --i-demultiplexed-seqs "$OUTPUT/04_merge/merged.qza" \
    --p-trim-length "$TRIM_LENGTH" \
    --p-sample-stats \
    --p-jobs "$THREADS" \
    --o-representative-sequences "$OUTPUT/05_deblur/rep-seqs-deblur.qza" \
    --o-table "$OUTPUT/05_deblur/table-deblur.qza" \
    --o-stats "$OUTPUT/05_deblur/deblur-stats.qza"

###############################################################################
# 6. VISUALIZACIÓN DE LAS STATS DE DEBLUR
###############################################################################

echo ""
echo "============================================================"
echo "6. STATS DE DEBLUR"
echo "============================================================"

qiime deblur visualize-stats \
    --i-deblur-stats "$OUTPUT/05_deblur/deblur-stats.qza" \
    --o-visualization "$OUTPUT/05_deblur/deblur-stats.qzv"

###############################################################################
# 7. RESUMEN DE LA TABLA DE ASVs
###############################################################################

echo ""
echo "============================================================"
echo "7. RESUMEN DE TABLA DE ASVs"
echo "============================================================"

qiime feature-table summarize \
    --i-table "$OUTPUT/05_deblur/table-deblur.qza" \
    --m-sample-metadata-file "$METADATA" \
    --o-visualization "$OUTPUT/05_deblur/table-deblur.qzv"

###############################################################################
# 8. VISUALIZACIÓN DE LAS SECUENCIAS REPRESENTATIVAS
###############################################################################

echo ""
echo "============================================================"
echo "8. VISUALIZACIÓN DE LAS SECUENCIAS REPRESENTATIVAS"
echo "============================================================"

qiime feature-table tabulate-seqs \
    --i-data "$OUTPUT/05_deblur/rep-seqs-deblur.qza" \
    --o-visualization "$OUTPUT/05_deblur/rep-seqs-deblur.qzv"

###############################################################################
# 9. CLASIFICACIÓN TAXONÓMICA CON SILVA
#
# El clasificador se proporciona mediante:
#
#   --classifier /ruta/silva_classifier.qza
#
###############################################################################

echo ""
echo "============================================================"
echo "9. CLASIFICACIÓN TAXONÓMICA SILVA"
echo "============================================================"

qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$OUTPUT/05_deblur/rep-seqs-deblur.qza" \
    --o-classification "$OUTPUT/06_taxonomy/taxonomy.qza"

###############################################################################
# 10. VISUALIZACIÓN DE TAXONOMÍA
###############################################################################

echo ""
echo "============================================================"
echo "10. VISUALIZACIÓN DE TAXONOMÍA"
echo "============================================================"

qiime metadata tabulate \
    --m-input-file "$OUTPUT/06_taxonomy/taxonomy.qza" \
    --o-visualization "$OUTPUT/06_taxonomy/taxonomy.qzv"

###############################################################################
# FIN
###############################################################################

echo ""
echo "============================================================"
echo "PIPELINE FINALIZADO"
echo "============================================================"
echo ""

echo "ASV table:"
echo "$OUTPUT/05_deblur/table-deblur.qza"
echo ""

echo "ASV representative sequences:"
echo "$OUTPUT/05_deblur/rep-seqs-deblur.qza"
echo ""

echo "Deblur statistics:"
echo "$OUTPUT/05_deblur/deblur-stats.qzv"
echo ""

echo "Taxonomy:"
echo "$OUTPUT/06_taxonomy/taxonomy.qza"
echo ""

echo "Taxonomy visualization:"
echo "$OUTPUT/06_taxonomy/taxonomy.qzv"
echo ""
