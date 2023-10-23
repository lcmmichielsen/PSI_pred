# Training models
Here, we show how to train the DL models using the previously created TFRecords. There are two mandatory inputs for the `exon_train.py` function: the params.json file and the path to the directory with the TFRecords.

## Human data, sequence only
```
python -u exon_train.py \
-o test_output/tfr_records_human_seqonly_var01/glia/fold0/run0 \
--rbp 0 --seq 1 --splice 1 \
test_output/tfr_records_human_seqonly_var01/glia/fold0/params.json \
test_output/tfr_records_human_seqonly_var01/glia/fold0
```

## Human data, sequence + RBP
```
python -u exon_train.py \
-o test_output/tfr_records_human_seqRBP_var01/glia/fold0/run0 \
--rbp 1 --seq 1 --splice 1 \
test_output/tfr_records_human_seqRBP_var01/glia/fold0/params.json \
test_output/tfr_records_human_seqRBP_var01/glia/fold0
```
## Multihead model, trained on human and mouse
```
python -u exon_train.py \
-o test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
--rbp 0 --seq 1 --splice 1 \
test_output/tfr_records_human_seqonly_var01/glia/fold0/params_multitask.json \
test_output/tfr_records_human_seqonly_var01/glia/fold0 \
test_output/tfr_records_mouse_seqonly_var01/glia/fold0
```

