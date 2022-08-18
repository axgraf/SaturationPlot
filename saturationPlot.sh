#!/bin/sh

OUTPUT_PDF="saturation_plot.pdf"
COV_THRESHOLD=10
FILES_SAMPLE="sample1.htseq.count.tsv,sample1 sample2.htseq.count.tsv,sample2"


R --no-restore --no-save --slave --quiet --args $OUTPUT_PDF $COV_THRESHOLD $FILES_SAMPLE < saturationPlot.r

