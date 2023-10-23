# Making predictions
Here, we show how to make predictions on the test fold using our previously trained models. To the `exon_test.py` script, there are 2 mandatory inputs: the path to the trained model and the path to the TFRecords.

## Human data, sequence only
```
python -u exon_test.py --save \
-o test_output/tfr_records_human_seqonly_var01/glia/fold0/run0 \
--rbp 0 --seq 1 --splice 1 \
test_output/tfr_records_human_seqonly_var01/glia/fold0/run0 \
test_output/tfr_records_human_seqonly_var01/glia/fold0
```

## Human data, sequence + RBP
```
python -u exon_test.py --save \
-o test_output/tfr_records_human_seqRBP_var01/glia/fold0/run0 \
--rbp 1 --seq 1 --splice 1 \
test_output/tfr_records_human_seqRBP_var01/glia/fold0/run0 \
test_output/tfr_records_human_seqRBP_var01/glia/fold0
```
## Multihead model
Here, we have to specify which head we want to make predictions with. Head 0 corresponds to the first input while training the model (mouse data) and head 1 to the second input.

### Predictions mouse
```
python -u exon_test.py --save \
-o test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0/mouse \
--rbp 0 --seq 1 --splice 1 \
--head 0 --multihead True \
test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
test_output/tfr_records_mouse_seqonly_var01/glia/fold0
```

### Predictions human
```
python -u exon_test.py --save \
-o test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0/human \
--rbp 0 --seq 1 --splice 1 \
--head 1 --multihead True \
test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
test_output/tfr_records_human_seqonly_var01/glia/fold0
```
