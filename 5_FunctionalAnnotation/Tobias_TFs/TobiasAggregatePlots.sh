#!/bin/bash
source /home/rna2/PROGRAMS/miniconda3/etc/profile.d/conda.sh
conda activate tobias
cd "/mnt/Data0/PROJECTS/CROPSeq/EnhancerPredictionModels/Results/Tobias/"
OUTPUT="Plots/AggregatePlots/"
OUTPUT_TXT="Footprint/Aggregate_txt/"
SIGNALS="Footprint/NHA_ATAC_S3.filtered_corrected.bw"
HITS=../../Data/PeakData/Hit_Enhancers.bed
NON_HITS=../../Data/PeakData/NonHit_Enhancers.bed
PEAKS="/mnt/Data0/PROJECTS/CROPSeq/PreprocessedData/Bulk_ATAC_RNAseq/GOK8505-GOK8837/GOK8837A1/Mapping/MACS2/NHA_ATAC_S3.filtered.BAM_peaks.narrowPeak"

#Compare Hits, Non_hits & All peaks  - unbound TFs
UNBOUND_TFs="Footprint/BINDetect/unbound_combined.bed"
TOBIAS PlotAggregate --TFBS $UNBOUND_TFs \
--signals $SIGNALS \
--TFBS-labels "Unbound_TFs" \
--regions $HITS $NON_HITS \
--region-labels "Hits" "Non-Hits" \
--output  $OUTPUT/Aggregate_Hits_vs_Non_Hits_Unbound_Only.pdf \
--output-txt $OUTPUT_TXT/Aggregate_Hits_vs_Non_Hits_Unbound_Only.txt \
--share_y both --plot_boundaries

#Compare Hits, Non_hits & All peaks  - bound TFs
BOUND_TFs="Footprint/BINDetect/combined.bed"
TOBIAS PlotAggregate --TFBS $BOUND_TFs \
--signals $SIGNALS \
--TFBS-labels "Bound_TFs" \
--regions $HITS $NON_HITS \
--region-labels "Hits" "Non-Hits" \
--output  $OUTPUT/Aggregate_Hits_vs_Non_Hits_Bound_Only.pdf \
--output-txt $OUTPUT_TXT/Aggregate_Hits_vs_Non_Hits_Bound_Only.txt \
--share_y both --plot_boundaries


ALL_TFs="Footprint/BINDetect/combined_all.bed"
#Compare Hits, Non_hits & All peaks - all TFs
TOBIAS PlotAggregate --TFBS $ALL_TFs \
--signals $SIGNALS \
--TFBS-labels "All_TFs" \
--regions $HITS $NON_HITS $PEAKS \
--region-labels "Hits" "Non-Hits" "All_Peaks" \
--output  $OUTPUT/Aggregate_Hits_vs_Non_Hits_vs_All.pdf \
--output-txt $OUTPUT_TXT/Aggregate_Hits_vs_Non_Hits_vs_All.txt \
--share_y both --plot_boundaries

