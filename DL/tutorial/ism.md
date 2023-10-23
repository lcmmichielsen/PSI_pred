# In-silico saturation mutagenesis
The input to this script is very similar to the `exon_test.py` function.

## Multihead model, ISM for human sequences
```
python -u exon_ism.py  \
-o test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0/human_ISM \
--rbp 0 --seq 1 --splice 1 \
--head 1 --multihead True \
test_output/tfr_records_human_mouse_seqonly_var01/glia/fold0/run0 \
test_output/tfr_records_human_seqonly_var01/glia/fold0
```
