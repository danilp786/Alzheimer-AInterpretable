#!/usr/bin/env bash

set -euo pipefail


###############################################################################
# PIPELINE 16S rRNA V3-V4
###############################################################################

FASTQ_DIR=""
METADATA=""
CLASSIFIER=""
OUTPUT="results"


###############################################################################
# 1. ARGUMENTOS
###############################################################################

while [[ $# -gt 0 ]]; do

    case "$1" in

        --fastq-dir)
            FASTQ_DIR="$2"
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

        --output)
            OUTPUT="$2"
            shift 2
            ;;

        *)
            echo "ERROR: argumento desconocido: $1"
            exit 1
            ;;

    esac

done


###############################################################################
# 2. COMPROBAR ARGUMENTOS
###############################################################################

if [[ -z "$FASTQ_DIR" ]]; then
    echo "ERROR: falta --fastq-dir"
    exit 1
fi

if [[ -z "$METADATA" ]]; then
    echo "ERROR: falta --metadata"
    exit 1
fi

if [[ -z "$CLASSIFIER" ]]; then
    echo "ERROR: falta --classifier"
    exit 1
fi


###############################################################################
# 3. COMPROBAR ARCHIVOS Y QIIME 2
###############################################################################

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: no existe el directorio FASTQ:"
    echo "$FASTQ_DIR"
    exit 1
fi

if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: no existe el metadata:"
    echo "$METADATA"
    exit 1
fi

if [[ ! -f "$CLASSIFIER" ]]; then
    echo "ERROR: no existe el clasificador:"
    echo "$CLASSIFIER"
    exit 1
fi

if ! command -v qiime >/dev/null 2>&1; then
    echo "ERROR: QIIME 2 no está disponible."
    echo "Ejecuta primero:"
    echo
    echo "conda activate qiime2-2021.4"
    exit 1
fi


###############################################################################
# 4. CREAR DIRECTORIOS
###############################################################################

mkdir -p "$OUTPUT"

mkdir -p "$OUTPUT/01_import"
mkdir -p "$OUTPUT/02_quality"
mkdir -p "$OUTPUT/03_primers"
mkdir -p "$OUTPUT/04_merged"
mkdir -p "$OUTPUT/05_deblur"
mkdir -p "$OUTPUT/06_taxonomy"
mkdir -p "$OUTPUT/07_metadata"


###############################################################################
# 5. CREAR METADATA V3-V4
###############################################################################

echo
echo "======================================================"
echo "5. FILTRANDO METADATA V3-V4"
echo "======================================================"


python3 <<PY

import csv

input_file = "$METADATA"
output_file = "$OUTPUT/07_metadata/metadata_V3-V4.tsv"

with open(input_file, "r", encoding="utf-8-sig", newline="") as f:

    reader = csv.DictReader(f)

    if "Run" not in reader.fieldnames:
        raise SystemExit("ERROR: no existe la columna Run")

    if "Library Name" not in reader.fieldnames:
        raise SystemExit("ERROR: no existe la columna Library Name")

    rows = []

    for row in reader:

        library = row["Library Name"]

        if "16S_V3/4" in library:
            rows.append(row)


if len(rows) == 0:
    raise SystemExit("ERROR: no se encontraron muestras V3-V4")


fieldnames = reader.fieldnames


with open(output_file, "w", encoding="utf-8", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow(["sample-id"] + fieldnames)

    for row in rows:

        writer.writerow(
            [row["Run"]] +
            [row[field] for field in fieldnames]
        )


print()
print("Muestras V3-V4 encontradas:", len(rows))
print("Metadata:", output_file)

if len(rows) != 101:

    print()
    print("ADVERTENCIA:")
    print("Se esperaban 101 muestras V3-V4.")
    print("Se encontraron:", len(rows))

PY


###############################################################################
# 6. CREAR MANIFEST PAIRED-END
###############################################################################

echo
echo "======================================================"
echo "6. CREANDO MANIFEST PAIRED-END"
echo "======================================================"


MANIFEST="$OUTPUT/01_import/manifest.tsv"


echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" \
    > "$MANIFEST"


python3 <<PY

import csv
import os

metadata_file = "$OUTPUT/07_metadata/metadata_V3-V4.tsv"
fastq_dir = "$FASTQ_DIR"
manifest_file = "$MANIFEST"


with open(metadata_file, "r", encoding="utf-8") as f:

    reader = csv.DictReader(f, delimiter="\t")

    samples = [row["sample-id"] for row in reader]


with open(manifest_file, "a", encoding="utf-8") as out:

    for sample in samples:

        forward_candidates = [
            os.path.join(fastq_dir, sample + "_1.fastq.gz"),
            os.path.join(fastq_dir, sample + "_R1.fastq.gz"),
            os.path.join(fastq_dir, sample + "_1.fastq"),
            os.path.join(fastq_dir, sample + "_R1.fastq"),
        ]

        reverse_candidates = [
            os.path.join(fastq_dir, sample + "_2.fastq.gz"),
            os.path.join(fastq_dir, sample + "_R2.fastq.gz"),
            os.path.join(fastq_dir, sample + "_2.fastq"),
            os.path.join(fastq_dir, sample + "_R2.fastq"),
        ]

        forward = None
        reverse = None

        for f in forward_candidates:

            if os.path.exists(f):
                forward = os.path.abspath(f)
                break

        for f in reverse_candidates:

            if os.path.exists(f):
                reverse = os.path.abspath(f)
                break

        if forward is None:
            raise SystemExit(
                f"ERROR: no se encontró R1 para {sample}"
            )

        if reverse is None:
            raise SystemExit(
                f"ERROR: no se encontró R2 para {sample}"
            )

        out.write(
            f"{sample}\t{forward}\t{reverse}\n"
        )


print("Manifest creado correctamente.")
print("Muestras:", len(samples))

PY


###############################################################################
# 7. IMPORTAR FASTQ PAIRED-END
###############################################################################

echo
echo "======================================================"
echo "7. IMPORTANDO FASTQ PAIRED-END"
echo "======================================================"


qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "$OUTPUT/01_import/demux-paired.qza"


###############################################################################
# 8. VISUALIZAR CALIDAD DE LOS READS ORIGINALES
###############################################################################

echo
echo "======================================================"
echo "8. CALIDAD DE READS ORIGINALES"
echo "======================================================"


qiime demux summarize \
    --i-data "$OUTPUT/01_import/demux-paired.qza" \
    --o-visualization "$OUTPUT/02_quality/demux-paired.qzv"


###############################################################################
# 9. ELIMINACIÓN DE PRIMERS
###############################################################################

echo
echo "======================================================"
echo "9. ELIMINANDO PRIMERS"
echo "======================================================"


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
# 10. MERGE FORWARD + REVERSE
###############################################################################

echo
echo "======================================================"
echo "10. MERGE FORWARD + REVERSE"
echo "======================================================"


qiime vsearch join-pairs \
    --i-demultiplexed-seqs "$OUTPUT/03_primers/primer-trimmed.qza" \
    --o-joined-sequences "$OUTPUT/04_merged/joined.qza"


###############################################################################
# 11. VISUALIZAR READS MERGEADOS
###############################################################################

echo
echo "======================================================"
echo "11. RESUMEN DE READS MERGEADOS"
echo "======================================================"


qiime demux summarize \
    --i-data "$OUTPUT/04_merged/joined.qza" \
    --o-visualization "$OUTPUT/04_merged/joined.qzv"


###############################################################################
# 12. FILTRADO DE CALIDAD PHRED
###############################################################################

echo
echo "======================================================"
echo "12. FILTRADO DE CALIDAD PHRED"
echo "======================================================"


qiime quality-filter q-score \
    --i-demux "$OUTPUT/04_merged/joined.qza" \
    --o-filtered-sequences "$OUTPUT/04_merged/filtered.qza" \
    --o-filter-stats "$OUTPUT/04_merged/filter-stats.qza"


###############################################################################
# 13. ESTADÍSTICAS DEL FILTRADO DE CALIDAD
###############################################################################

echo
echo "======================================================"
echo "13. ESTADÍSTICAS DEL FILTRADO DE CALIDAD"
echo "======================================================"


qiime metadata tabulate \
    --m-input-file "$OUTPUT/04_merged/filter-stats.qza" \
    --o-visualization "$OUTPUT/04_merged/filter-stats.qzv"


###############################################################################
# 14. DEBLUR
###############################################################################

echo
echo "======================================================"
echo "14. DEBLUR"
echo "======================================================"


TRIM_LENGTH=440


qiime deblur denoise-16S \
    --i-demultiplexed-seqs "$OUTPUT/04_merged/filtered.qza" \
    --p-trim-length "$TRIM_LENGTH" \
    --p-sample-stats \
    --o-representative-sequences \
        "$OUTPUT/05_deblur/rep-seqs-deblur.qza" \
    --o-table \
        "$OUTPUT/05_deblur/table-deblur.qza" \
    --o-stats \
        "$OUTPUT/05_deblur/deblur-stats.qza"


###############################################################################
# 15. ESTADÍSTICAS DE DEBLUR
###############################################################################

echo
echo "======================================================"
echo "15. DEBLUR STATS"
echo "======================================================"


qiime deblur visualize-stats \
    --i-deblur-stats "$OUTPUT/05_deblur/deblur-stats.qza" \
    --o-visualization "$OUTPUT/05_deblur/deblur-stats.qzv"


###############################################################################
# 16. RESUMEN DE LA TABLA ASV
###############################################################################

echo
echo "======================================================"
echo "16. RESUMEN TABLA ASV"
echo "======================================================"


qiime feature-table summarize \
    --i-table "$OUTPUT/05_deblur/table-deblur.qza" \
    --m-sample-metadata-file "$OUTPUT/07_metadata/metadata_V3-V4.tsv" \
    --o-visualization "$OUTPUT/05_deblur/table-deblur.qzv"


###############################################################################
# 17. CLASIFICACIÓN TAXONÓMICA CON SILVA 138
###############################################################################

echo
echo "======================================================"
echo "17. CLASIFICACIÓN TAXONÓMICA SILVA 138"
echo "======================================================"


qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$OUTPUT/05_deblur/rep-seqs-deblur.qza" \
    --o-classification "$OUTPUT/06_taxonomy/taxonomy.qza"


###############################################################################
# 18. VISUALIZAR TAXONOMÍA
###############################################################################

echo
echo "======================================================"
echo "18. VISUALIZANDO TAXONOMÍA"
echo "======================================================"


qiime metadata tabulate \
    --m-input-file "$OUTPUT/06_taxonomy/taxonomy.qza" \
    --o-visualization "$OUTPUT/06_taxonomy/taxonomy.qzv"


###############################################################################
# 19. FIN
###############################################################################

echo
echo
echo "======================================================"
echo "PIPELINE COMPLETADA"
echo "======================================================"

echo
echo "Metadata V3-V4:"
echo "  $OUTPUT/07_metadata/metadata_V3-V4.tsv"

echo
echo "Tabla ASV:"
echo "  $OUTPUT/05_deblur/table-deblur.qza"

echo
echo "Secuencias representativas:"
echo "  $OUTPUT/05_deblur/rep-seqs-deblur.qza"

echo
echo "Deblur stats:"
echo "  $OUTPUT/05_deblur/deblur-stats.qza"

echo
echo "Taxonomía:"
echo "  $OUTPUT/06_taxonomy/taxonomy.qza"

echo
echo "======================================================"
