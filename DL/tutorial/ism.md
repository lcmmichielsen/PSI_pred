## Making predictions

### Example making predictions with model trained on glia, sequence only
### Two mandatory inputs: model directory and data directory
python -u exon_test.py --save \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0/run0 \
--rbp 0 --seq 1 --splice 1 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0/run0 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0

### Example making predictions with model trained on glia, sequence + RBP
python -u exon_test.py --save \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqRBP_var01/glia/fold0/run0 \
--rbp 1 --seq 1 --splice 1 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqRBP_var01/glia/fold0/run0 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqRBP_var01/glia/fold0

### Example making predictions with model trained on glia, human + mouse
### Head 0 --> first input when training, mouse in this case
### Head 1 --> second input, human in this case

### Predictions mouse
python -u exon_test.py --save \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0/mouse \
--rbp 0 --seq 1 --splice 1 \
--head 0 --multihead True \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_mouse_seqonly_var01/glia/fold0

### Predictions human
python -u exon_test.py --save \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0/human \
--rbp 0 --seq 1 --splice 1 \
--head 1 --multihead True \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0
