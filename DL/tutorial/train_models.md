## Training models

### Example training model on glia, sequence only
### Two mandatory inputs: params file and data directory
python -u exon_train.py \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0/run0 \
--rbp 0 --seq 1 --splice 1 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0/params.json \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0

### Example training model on glia, sequence + RBP
### Two mandatory inputs: params file and data directory
python -u exon_train.py \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqRBP_var01/glia/fold0/run0 \
--rbp 1 --seq 1 --splice 1 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqRBP_var01/glia/fold0/params.json \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqRBP_var01/glia/fold0

### Example training on glia, human + mouse
python -u exon_train.py \
-o /exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
--rbp 0 --seq 1 --splice 1 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0/params_multitask.json \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_human_seqonly_var01/glia/fold0 \
/exports/humgen/lmichielsen/PSI_pred/test_output/tfr_records_mouse_seqonly_var01/glia/fold0 

