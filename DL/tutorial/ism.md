## In-silico saturation mutagenesis
# Input to this function is super similar to the exon_test function

### Example for human
python -u exon_ism.py  \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0/human_ISM \
--rbp 0 --seq 1 --splice 1 \
--head 1 --multihead True \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0
