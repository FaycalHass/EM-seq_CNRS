// this consumes fastq files, adapter trims, then aligns reads to the specified reference using bwa-meth
// mode can be set to tile_fastqs or run_fastqs depending on whether the system should map each tile's reads in distinct jobs then combine (tile_fastqs)
// or all reads for a library in a single job (run_fastqs)
nextflow.enable.dsl=1

flowcell = params.flowcell
genome = params.genome
params.tmp_dir = '/tmp'
outputPath = params.outdir //='output
fastq_mode = params.fastq_mode = 'run_fastqs'
println "Processing " + flowcell + "... => " + outputPath

fastq_glob = params.fastq_glob ?: '*.{1,2}*.fastq*'

println fastq_glob
Channel.fromFilePairs(fastq_glob)
    .map{ lib, read -> [flowcell: flowcell, 
                       library:lib,
                       insert_read1:read[0], 
                       insert_read2:read[1], 
                       barcode:'N', 
                       lane:'all', 
                       tile:'all' ]
    }.set{fq_set_channel}

process mapping {
    cpus fastq_mode == 'tile-fastq' ? 4 : 16
    // errorStrategy 'retry'
    tag { [flowcell, lib] }
    conda "bioconda::bwameth=0.2.2 bioconda::seqtk=1.3 bioconda::sambamba=0.7.0 bioconda::fastp=0.20.1 bioconda::mark-nonconverted-reads=1.1"

    input:
        // val fq_set from fq_set_channel
        set val(fcell), val(lib), file(read1), file(read2), val(barcode), val(lane), val(tile) from fq_set_channel

    output:
        set val(lib), file("*.aln.bam") into aligned_files
        set val(lib), file("*.nonconverted.tsv") into nonconverted_counts
        set val(lib), file("*_fastp.json") into fastp_log_files

    shell:
    '''
    inst_name=$(zcat -f '!{read1}' | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
    fastq_barcode=$(zcat -f '!{read1}' | head -n 1 | sed -r 's/.*://')

    if [[ "${inst_name:0:2}" == 'A0' || "${inst_name:0:2}" == 'NS' || \
          "${inst_name:0:2}" == 'NB' || "${inst_name:0:2}" == 'VH' || \
          "${inst_name: -2:2}" == 'NX' ]] ; then
       trim_polyg='--trim_poly_g'
       echo '2-color instrument: poly-g trim mode on'
    else
       trim_polyg=''
    fi

    seqtk mergepe <(zcat -f "!{read1}") <(zcat -f "!{read2}") > "!{lib}_${fastq_barcode}!{fcell}_!{lane}_!{tile}.fastq"

    fastp -l 2 \
        -Q ${trim_polyg} \
        --interleaved_in \
        --overrepresentation_analysis \
        -i "!{lib}_${fastq_barcode}!{fcell}_!{lane}_!{tile}.fastq" \
        -j "!{lib}_fastp.json" \
        -o "!{lib}_${fastq_barcode}!{fcell}_!{lane}_!{tile}_R1_cleaned.fastq" \
        --out2 "!{lib}_${fastq_barcode}!{fcell}_!{lane}_!{tile}_R2_cleaned.fastq"\
        2> fastp.stderr

    bwameth.py -p -t !{task.cpus} \
        --read-group "@RG\\tID:${fastq_barcode}\\tSM:!{lib}" \
        --reference !{genome} "!{lib}_${fastq_barcode}!{fcell}_!{lane}_!{tile}_R1_cleaned.fastq" "!{lib}_${fastq_barcode}!{fcell}_!{lane}_!{tile}_R2_cleaned.fastq" \
        > "!{lib}_${fastq_barcode}_!{fcell}_!{lane}_!{tile}.aln_bwameth.sam" \
        2> "!{lib}_${fastq_barcode}!{fcell}_!{lane}_!{tile}.log.bwamem"

    mark-nonconverted-reads.py \
        --bam "!{lib}_${fastq_barcode}_!{fcell}_!{lane}_!{tile}.aln_bwameth.sam" \
        --out "!{lib}_${fastq_barcode}_!{fcell}_!{lane}_!{tile}.aln_markNon.sam" \
        2> "!{lib}_${fastq_barcode}_!{fcell}_!{lane}_!{tile}.nonconverted.tsv"

    sambamba view -t 2 -S -f bam \
    -o "!{lib}_${fastq_barcode}_!{fcell}_!{lane}_!{tile}.aln.bam" \
    "!{lib}_${fastq_barcode}_!{fcell}_!{lane}_!{tile}.aln_markNon.sam" \
    2> sambamba.stderr
    '''
}

process mergeAndMarkDuplicates {
    cpus 8
    errorStrategy 'retry'
    tag { library }
    publishDir "${outputPath}", mode: 'copy', pattern: '*.{md.bam}*'
    conda "bioconda::samtools=1.9 bioconda::samblaster=0.1.24 bioconda::sambamba=0.7.0"

    input:
        set val(library), file(libraryBam) from aligned_files.groupTuple()

    output:
        set val(library), file('*.md.bam'), file('*.md.bam.bai') into md_bams
        file('*.samblaster') into samblaster_logs

    shell:
    '''
    samtools cat  -b <( find . -name '*.aln.bam' ) \
    | samtools view -h /dev/stdin \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t 2 -l 0 -S -f bam /dev/stdin \
    | sambamba sort --tmpdir=!{params.tmp_dir} -t !{task.cpus} -m 20GB -o !{library}.md.bam /dev/stdin

    '''
}

    // we need to send the same file into multiple channels for downstream 
    // analysis whether they came from aligned bams or fastqs
    md_bams.into {  md_files_for_mbias; md_files_for_extract; md_files_for_fastqc; 
                    md_files_for_samstats; md_files_for_picard; 
                    md_files_for_goleft; md_files_for_picard_gc; md_files_for_samflagstats; 
                    md_files_for_aggregate; md_files_for_human_reads;
                 }
                 

process methylDackel_mbias {
    cpus 8
    errorStrategy 'retry'
    tag {library}
    conda "bioconda::methyldackel=0.6.1 conda-forge::pigz=2.8"

    input:
        tuple library, file(md_file), file(md_bai) from md_files_for_mbias.groupTuple()

    output:
        file('*.svg') into mbias_output_svg
        file('*.tsv') into mbias_output_tsv
        tuple library, file('*.tsv') into mbias_for_aggregate

    shell:
        '''
        echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > !{library}_combined_mbias.tsv
        chrs=(`samtools view -H !{md_file} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

        for chr in ${chrs[*]}; do
            for context in CHH CHG CpG; do
                arg=''
                if [ $context = 'CHH' ]; then
                arg='--CHH --noCpG'
                elif [ $context = 'CHG' ]; then
                arg='--CHG --noCpG'
                fi
                # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
                # not sure why we need both --keepDupes and -F, probably a bug in mbias
                join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
                <( \
                    MethylDackel mbias --noSVG $arg -@ !{task.cpus} -r $chr !{genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                <( \
                    MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ !{task.cpus} -r $chr !{genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                | sed "s/^/${chr}\t${context}\t/" \
                >> !{library}_combined_mbias.tsv
            done
        done
        # makes the svg files for trimming checks
        MethylDackel mbias -@ !{task.cpus} --noCpG --CHH --CHG -r ${chrs[0]} !{genome} !{md_file} !{library}_chn
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

        MethylDackel mbias -@ !{task.cpus} -r ${chrs[0]} !{genome} !{md_file} !{library}_cpg
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;

        '''

}

    process methylDackel_extract {
        cpus 8
        tag {library}
        publishDir "${outputPath}", mode: 'copy'
        conda "bioconda::methyldackel=0.6.1 conda-forge::pigz=2.8"

        input:
            tuple library, file(md_file), file(md_bai) from md_files_for_extract.groupTuple()

        output:
            tuple library, file('*.methylKit.gz') into extract_output

        shell:
        '''
        MethylDackel extract --methylKit --nOT 0,0,0,5 --nOB 0,0,5,0 -@ !{task.cpus} --CHH --CHG -o !{library} !{genome} !{md_file}
        pigz -p !{task.cpus} *.methylKit
        '''

    }

    process select_human_reads {
        cpus 8 
        tag {library}
        conda "bioconda::sambamba=0.7.1 bioconda::bedtools=2.29.2"

        input:
            tuple library, file(md_file), file(md_bai) from md_files_for_human_reads.groupTuple()

        output:
            tuple library, file('*.human.bam') into human_bams_gc
            tuple library, file('*.human.bam') into human_bams_inserts

        shell:
        '''
        sambamba view -t !{task.cpus} -l 0 -f bam !{md_file} chr1 chr2 chr3 chr4 chr5 chr6 \
                                                  chr7 chr8 chr9 chr10 chr11 chr12 \
                                                  chr13 chr14 chr15 chr16 chr17 chr18 \
                                                  chr19 chr20 chr21 chr22 chrX chrY \
        > !{md_file}.human.bam
        '''

    }

    process runFastQC {
        cpus 1
        errorStrategy 'retry'
        tag { library }
        conda "bioconda::fastqc=0.11.8"

        input:
            tuple library, file(md_file), file(md_bai) from md_files_for_fastqc.groupTuple()

        output:
            file('*_fastqc.zip') into fastqc_results
            tuple library, file('*_fastqc.zip') into fastqc_results_for_aggregate

        shell:
        '''
        fastqc -f bam !{md_file}
        '''

    }

    process sum_nonconverted_reads {	

        input:	
            tuple library, file(count_files) from nonconverted_counts.groupTuple()	

        output:	
            file('*-nonconverted-counts.tsv') into cat_nonconversions	
            tuple library, file('*-nonconverted-counts.tsv') into nonconverted_counts_for_aggregate

        shell:	
        '''	
        files=(*.tsv)	
        paste *.tsv | awk -v numFiles=${#files[@]} -v OFS='\t' '	
        {	
        row = sep = ""	
        for(i=1; i < NF/numFiles; ++i) { row = row sep $i; sep = OFS }	
        sum = $(NF/numFiles) # last header col. / (1st) data col. to sum	
        for(i=2; i<=numFiles; ++i) sum += $(NF/numFiles * i) # add other cols.	
        printf "%s%s%s\\n", row, OFS, sum	
        }' > tmp-counts.tsv	
        awk '{print "!{library}\t" $0}' tmp-counts.tsv > !{library}-nonconverted-counts.tsv	
        '''    	
    }	

    process combine_nonconversion {	
        publishDir "${outputPath}", mode: 'copy'	

        input:	
            file ('*') from cat_nonconversions

        output:	
            file ("combined-nonconverted.tsv")	

        shell:	
        '''	
        cat *.tsv > combined-nonconverted.tsv	
        '''	

    }

    process samtools_flagstats {
        cpus 2
        errorStrategy 'retry'
        tag { library }
        conda "bioconda::samtools=1.9"

        input:
            tuple library, file(md_file), file(md_bai) from md_files_for_samflagstats.groupTuple(by: 0)

        output:
            file('*.flagstat') into flagstats
            file('*.idxstat') into idxstats
            tuple library, file('*.flagstat') into flagstats_for_aggregate
            tuple library, file('*.idxstat') into idxstats_for_aggregate
            

        shell:
        '''
        samtools flagstat -@!{task.cpus} !{md_file} > !{md_file}.flagstat
        samtools idxstats !{md_file} > !{md_file}.idxstat
        '''
    }

    process samtools_stats {
        cpus 2
        errorStrategy 'retry'
        tag { library }
        conda "bioconda::samtools=1.9"

        input:
            tuple library, file(md_file),file(md_bai) from md_files_for_samstats.groupTuple()

        output:

            file('*.samstat') into samstats

        shell:
        '''
        samtools stats -@!{task.cpus} !{md_file} > !{md_file}.samstat
        '''
    }

    process picard_gc_bias {
        cpus 1
        errorStrategy 'retry'
        tag { library }
        conda "bioconda::picard=2.20.7"

        input:
            tuple library, file(md_file), file(md_bai) from md_files_for_picard_gc.groupTuple()

        output:
            file('*gc_metrics') into picard_gc_stats

        shell:
        '''
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=LENIENT I=!{md_file} O=!{md_file}.gc_metrics S=!{md_file}.gc_summary_metrics CHART=!{md_file}.gc.pdf R=!{genome}
        '''
    }

    process picard_stats {

        cpus 4
        errorStrategy 'retry'
        tag { library }
        conda "bioconda::picard=2.20.7"

        input:
            tuple library, file(md_file), file(md_bai) from md_files_for_picard.groupTuple()

        output:
            file('*_metrics') into picard_stats

        shell:
        '''
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=!{md_file} O=!{md_file}.insertsize_metrics MINIMUM_PCT=0 HISTOGRAM_FILE=/dev/null
        '''
    }

    process human_gc_bias {
        cpus 1
        errorStrategy 'retry'
        tag { library }
        conda "bioconda::picard=2.20.7"

        input:
            tuple library, file(md_file) from human_bams_gc.groupTuple()

        output:
            tuple library, file('*gc_metrics') into human_gc_stats_for_aggregate

        shell:
        '''
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=LENIENT I=!{md_file} O=!{md_file}.gc_metrics S=!{md_file}.gc_summary_metrics CHART=!{md_file}.gc.pdf R=!{genome}
        '''
    }

    process human_insert_size {

        cpus 4
        errorStrategy 'retry'
        tag { library }
        conda "bioconda::picard=2.20.7"

        input:
            tuple library, file(md_file) from human_bams_inserts.groupTuple()

        output:
            tuple library, file('*insertsize_metrics') into human_stats_for_aggregate

        shell:
        '''
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=!{md_file} O=!{md_file}.insertsize_metrics MINIMUM_PCT=0.0001 HISTOGRAM_FILE=/dev/null
        '''
    }

    process goleft {
        cpus 1
        conda "bioconda::goleft=0.2.0"

        input:
            tuple library, file(md_file), file(md_bai) from md_files_for_goleft.groupTuple()

        output:
            file("${library}/*-indexcov.ped") into goleft_ped
            file("${library}/*-indexcov.roc") into goleft_roc

        shell:
        '''
            goleft indexcov --directory !{library} *.bam
        '''
    }

    process multiqc {
        cpus 1
        publishDir "${outputPath}", mode: 'copy'
        conda "bioconda::multiqc=1.14.0"

        input:
            file('*') from fastqc_results.flatten().toList()
            file('*') from flagstats.flatten().toList()
            file('*') from idxstats.flatten().toList()
            file('*') from samstats.flatten().toList()
            file('*') from picard_stats.flatten().toList()
            file('*') from picard_gc_stats.flatten().toList()
            file('*') from goleft_ped.flatten().toList()
            file('*') from goleft_roc.flatten().toList()
            file('*') from samblaster_logs.flatten().toList()
            file('*') from fastp_log_files.flatten().toList()

        output:
            file "*report.html"

        shell:
        '''
        for file in $(cat input.* | sed -e 's/\\[//g' | sed -e 's/, \\|\\]/\\n/g'); do ln -s ${file} ./; done
        cat <<CONFIG > multiqc_config.yaml 
    title: Bwameth Alignment Summary - !{flowcell}
    extra_fn_clean_exts:
        - '.md'
        - '_fastp'
    custom_plot_config:
        picard_insert_size:
            xmax: 1000
    table_columns_placement:
        Samtools Stats:
            raw_total_sequences: 10
            reads_mapped_percent: 20
            reads_properly_paired_percent: 30
            reads_MQ0_percent: 35
        Samblaster:
            pct_dups: 40
        Picard:
            summed_median: 50
    table_columns_visible:
        Picard:
            PCT_PF_READS_ALIGNED: False
            summed_mean: False
        Samtools Stats:
            reads_mapped: False
            mapped_passed: False
            non-primary_alignments: False
            reads_MQ0_percent: True
        Samtools Flagstat:
            mapped_passed: False
        samtools_idxstats_always:
            - plasmid_puc19c
            - phage_lambda
        FastQC:
            percent_duplicates: False
            total_sequences: False
            avg_sequence_length: False
            percent_fails: False
            total_sequences: False
    CONFIG

        multiqc -ip  .
        '''
    }

    process combine_mbias_tsv {
        publishDir "${outputPath}", mode: 'copy', pattern: 'combined*'

        input:
            file(tsv) from mbias_output_tsv

        output:
            file "combined-mbias.tsv"

        shell:
        '''
            echo -ne 'flowcell\tlibrary\t' > combined-mbias.tsv
            ls *_mbias.tsv | head -n 1 | xargs head -n 1 >> combined-mbias.tsv
            for f in *_mbias.tsv; do 
                filebase=`basename "${f}" _combined_mbias.tsv`
                paste <( yes "!{flowcell}	${filebase}" | head -n `nl "$f" | tail -n 1 | cut -f 1` ) "$f" | tail -n +2  >> combined-mbias.tsv
            done
        '''
    }

    process combine_mbias_svg {
        publishDir "${outputPath}", mode: 'copy', pattern: 'combined*'
        conda 'conda-forge::cairosvg=2.4.2 conda-forge::ghostscript=9.22'

        input:
            file(svg) from mbias_output_svg.groupTuple()

        output:
            file "combined-mbias.pdf"

        shell:
        '''
        for f in *.svg; do
            cairosvg <(sed s/-nan/-1/ $f) -o $f.pdf
        done
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=combined-mbias.pdf *.pdf
        '''
    }









































# Charger les bibliothèques nécessaires
library(methylKit)
library(ggplot2)

# Fonction principale
analyzeMethylation <- function(input_file, output_dir, sample_id = "Sample1", assembly = "hg38",
                               treatment = 0, context = "CpG", lo.count = 10, hi.perc = 99.9) {
  
  # Créer le dossier de sortie si nécessaire
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Gestion d'erreur pour la lecture du fichier
  cat("Lecture du fichier de méthylation...\n")
  meth <- tryCatch({
    methRead(input_file, sample.id = sample_id, assembly = assembly,
             treatment = treatment, context = context)
  }, error = function(e) {
    stop(paste("Erreur lors de la lecture du fichier :", e$message))
  })
  
  # Sauvegarde du graphique de statistiques brutes
  cat("Génération du graphique de statistiques brutes...\n")
  png(file.path(output_dir, paste0(sample_id, "_raw_methylation_stats.png")), width = 800, height = 600)
  getMethylationStats(meth, plot = TRUE)
  dev.off()
  
  # Sauvegarde de l'histogramme de couverture
  cat("Génération de l'histogramme de couverture...\n")
  png(file.path(output_dir, paste0(sample_id, "_coverage_histogram.png")), width = 800, height = 600)
  hist(getData(meth)$coverage, breaks = 100, main = "Distribution de la couverture", xlab = "Couverture")
  dev.off()
  
  # Résumé statistique
  cat("Résumé statistique de la couverture :\n")
  print(summary(getData(meth)$coverage))
  
  # Filtrage des données
  cat("Filtrage des données par couverture...\n")
  meth.filt <- filterByCoverage(meth, lo.count = lo.count, hi.perc = hi.perc)
  
  # Sauvegarde du graphique après filtrage
  png(file.path(output_dir, paste0(sample_id, "_filtered_methylation_stats.png")), width = 800, height = 600)
  getMethylationStats(meth.filt, plot = TRUE, both.strands = FALSE)
  dev.off()
  
  # Extraire les données
  meth.data <- getData(meth.filt)
  
  # Sauvegarde des résultats
  cat("Sauvegarde des résultats...\n")
  write.csv(meth.data, file = file.path(output_dir, paste0(sample_id, "_meth_data_filtered.csv")), row.names = FALSE)
  saveRDS(meth.filt, file = file.path(output_dir, paste0(sample_id, "_meth_filtered.rds")))
  
  cat("Analyse terminée avec succès pour : ", sample_id, "\n")
}

# Exemple d'utilisation
analyzeMethylation(
  input_file = "/Isiprod1/project/SB_fhassani/EM-seq/GSM_kit/GSM2772524_NA12878_CpG_methylkit_original.txt.gz",
  output_dir = "/Isiprod1/project/SB_fhassani/EM-seq/test_CpG/figure_GSM2772524_NA12878_CpG",
  sample_id = "GSM2772524_NA12878_CpG",
  treatment = 0,
  assembly = "hg19",
  context = "CpG"
) 























































































# ------------------------------------------------------------------------
#  EM-seq Methylation Pipeline – tout-en-un (PNG version, options forcées)
#  (c) 2025  –  Adapté pour : /Isiprod1/project/SB_fhassani/EM-seq
# ------------------------------------------------------------------------
#  ► Toutes les étapes optionnelles (annotation & enrichissement GO) sont
#    désormais **forcées** – si les packages manquent, le script tentera
#    de les installer via BiocManager.
#  ► Tous les graphiques sont exportés en PNG (300 dpi).
# ------------------------------------------------------------------------
#  ► À lancer :
#     Rscript  EM_seq_pipeline_full_png.R
# ------------------------------------------------------------------------

## ───────────────────────── 0. PARAMETERS ──────────────────────────────
projDir   <- "/Isiprod1/project/SB_fhassani/EM-seq"   # racine projet
assembly  <- "hg38"
sample.id <- c("CTRL", "TRAITE")
file.list <- c(
  "output_test_SRR13953195/SRR13953195_CpG.methylKit.gz",
  "GSM_kit/GSM2772524_NA12878_CpG_percent.txt.gz"
)
file.list <- as.list(file.path(projDir, file.list))

treatment <- c(0, 1)

# -------- NOUVEAU : dossier de sortie centralisé --------
outDir <- "/Isiprod1/project/SB_fhassani/EM-seq/test_CpG/test_result_GSM"
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

# Seuils de différenciation
site.diff <- 25; site.qval <- 0.01
win.diff  <- 10; win.qval <- 0.05

cores <- max(1, parallel::detectCores() - 1)
opt   <- options(mc.cores = cores)

## ───────────────────────── 1. PACKAGES ────────────────────────────────
message("[1] Loading packages …")
suppressPackageStartupMessages({
  library(methylKit); library(ggplot2); library(data.table)
})

# ---- Forcer l'installation des packages optionnels si absents --------
optional_pkgs <- c("genomation", "GenomicFeatures", "clusterProfiler", "org.Hs.eg.db", "enrichplot")
for (p in optional_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(sprintf("[1] Installing missing package: %s", p))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}
# (Re)charger après installation
suppressPackageStartupMessages({
  lapply(optional_pkgs, function(pkg) {
    suppressWarnings(require(pkg, character.only = TRUE))
  })
})

## helper for high‑res PNG
gsave_png <- function(plot, filename, width = 7, height = 7, dpi = 300) {
  ggsave(filename = filename, plot = plot, width = width, height = height,
         dpi = dpi, units = "in")
}

## ───────────────────────── 2. INPUT QC ────────────────────────────────
message("[2] Reading methylation calls …")
stopifnot(length(file.list) == length(sample.id),
          length(sample.id) == length(treatment))

myobj <- methRead(location = file.list,
                  sample.id = as.list(sample.id),
                  assembly  = assembly,
                  treatment = treatment,
                  context   = "CpG")

# Couverture avant filtre — un PNG par échantillon
for (i in seq_along(myobj)) {
  png(file.path(outDir, sprintf("QC_coverage_before_filter_%s.png", sample.id[i])),
      width = 7, height = 7, units = "in", res = 300)
  getCoverageStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
}

## ───────────────────── 3. FILTER & NORMALISE ──────────────────────────
message("[3] Filtering & normalising …")
myobj.filt <- filterByCoverage(myobj, lo.count = 10, hi.perc = 99.9)
myobj.norm <- normalizeCoverage(myobj.filt)

# Couverture après filtre
for (i in seq_along(myobj.norm)) {
  png(file.path(outDir, sprintf("QC_coverage_after_filter_%s.png", sample.id[i])),
      width = 7, height = 7, units = "in", res = 300)
  getCoverageStats(myobj.norm[[i]], plot = TRUE, both.strands = FALSE)
  dev.off()
}

## ───────────────────── 4. UNION & DIFFERENTIAL ────────────────────────
message("[4] Unite CpG positions …")

meth.united <- unite(myobj.norm, destrand = TRUE)

diffMeth <- calculateDiffMeth(meth.united,
                              test = "fast.fisher",
                              overdispersion = "MN",
                              mc.cores = cores)

dmr.hyper <- getMethylDiff(diffMeth, difference = site.diff,
                           qvalue = site.qval, type = "hyper")
dmr.hypo  <- getMethylDiff(diffMeth, difference = site.diff,
                           qvalue = site.qval, type = "hypo")

## ───────────────────── 5. WINDOW (DMR) ANALYSIS ───────────────────────
message("[5] Differential methylation (1 kb windows) …")
tiles        <- tileMethylCounts(myobj.norm, win.size = 1000, step.size = 1000,
                                 mc.cores = cores)
tiles.united <- unite(tiles, destrand = TRUE)
tiles.diff   <- calculateDiffMeth(tiles.united, test = "fast.fisher",
                                  mc.cores = cores)
tiles.sig    <- getMethylDiff(tiles.diff, difference = win.diff,
                              qvalue = win.qval)

## ───────────────────── 6. SAVE RESULTS ────────────────────────────────
message("[6] Writing outputs …")
write.csv(dmr.hyper, file.path(outDir, "CpGs_hypermethylated.csv"), row.names = FALSE)
write.csv(dmr.hypo , file.path(outDir, "CpGs_hypomethylated.csv"), row.names = FALSE)

saveRDS(diffMeth ,  file.path(outDir, "diffMeth_perSite.rds"))
saveRDS(tiles.diff, file.path(outDir, "diffMeth_tiles.rds"))

if (nrow(dmr.hyper) > 0) {
  bed <- data.table(chr = dmr.hyper$chr, start = dmr.hyper$start, end = dmr.hyper$end)
  fwrite(bed[1:min(1000,.N)], file.path(outDir, "top_hyper_CpG.bed"),
         sep = "\t", col.names = FALSE)
}
if (nrow(dmr.hypo) > 0) {
  bed <- data.table(chr = dmr.hypo$chr, start = dmr.hypo$start, end = dmr.hypo$end)
  fwrite(bed[1:min(1000,.N)], file.path(outDir, "top_hypo_CpG.bed"),
         sep = "\t", col.names = FALSE)
}

## ───────────────────── 7. PLOTS INTERPRETATION ────────────────────────
message("[7] Creating interpretation plots …")
if (nrow(dmr.hyper) > 0) {
  p.volcano <- ggplot(dmr.hyper, aes(meth.diff, -log10(qvalue))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_vline(xintercept = c(site.diff, -site.diff), lty = 2) +
    geom_hline(yintercept = -log10(site.qval), lty = 2) +
    labs(title = "CpG hyperméthylés", x = "∆ méthylation (%)", y = "-log10(qvalue)") +
    theme_minimal()
  gsave_png(p.volcano, file.path(outDir, "plot_volcano_hyper.png"))
}

df.all <- as.data.frame(diffMeth)
if (nrow(df.all) > 0) {
  df.all$chr <- factor(df.all$chr, levels = paste0("chr", c(1:22, "X", "Y")))
  df.all$pos <- (df.all$start + df.all$end)/2
  p.manhattan <- ggplot(df.all, aes(pos/1e6, -log10(qvalue), colour = chr)) +
    geom_point(alpha = 0.4, size = 0.3) +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    theme_minimal(base_size = 8) +
    theme(legend.position = "none") +
    labs(x = "Position (Mb)", y = "-log10(qvalue)")
  gsave_png(p.manhattan, file.path(outDir, "plot_manhattan.png"), width = 12, height = 4)
}

## ───────────────────── 8. ANNOTATION (FORCÉ) ──────────────────────────
if (nrow(tiles.sig) > 0) {
  message("[8] Annotating DMRs …")
  suppressPackageStartupMessages({library(genomation); library(GenomicFeatures)})
  txdb <- makeTxDbFromUCSC(genome = assembly, tablename = "refGene")
  ann.tiles <- annotateWithGeneParts(as(tiles.sig, "GRanges"), getGeneParts(txdb))
  write.table(as.data.frame(ann.tiles), file.path(outDir, "DMR_tiles_annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  png(file.path(outDir, "plot_annotation_bar.png"), width = 7, height = 7, units = "in", res = 300)
  plotTargetAnnotation(ann.tiles, col = "cornflowerblue", border = NA)
  dev.off()
} else {
  message("[8] No tiles for annotation – skipped.")
}

## ───────────────────── 9. ENRICHISSEMENT GO (FORCÉ) ───────────────────
if (exists("ann.tiles") && length(ann.tiles) > 0) {
  message("[9] GO enrichment …")
  suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)})
  genes <- unique(ann.tiles$geneId[ann.tiles$feature == "promoter"])
  if (length(genes) > 0) {
    ego   <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
    p.go   <- dotplot(ego, showCategory = 20) + ggtitle("GO Biological Process (promoters DMR)")
    gsave_png(p.go, file.path(outDir, "GO_BP_enrichment.png"), width = 8, height = 6)
    saveRDS(ego, file.path(outDir, "GO_enrich.rds"))
  } else {
    message("⏩ GO enrichment ignoré (aucun gène/promoteur détecté).")
  }
} else {
  message("⏩ GO enrichment ignoré (packages ou annotations manquants).")
}

## ───────────────────── 10. SESSION INFO ──────────────────────────────
writeLines(capture.output(sessionInfo()), file.path(outDir, "sessionInfo.txt"))
options(opt)
message("
✔ Pipeline terminé : tous les résultats sont dans → ", outDir)
