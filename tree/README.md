# Single Cell Tree Construction

## Build Tree by Maximum Likelihood Method

[IQ-TREE](https://github.com/Cibiv/IQ-TREE)

```bash
iqtree2 -nt 10 -s input.fa -m GTR2+FO+I+R10 -alrt 1000 -bb 1000 -asr --mlrate -wslmr -wspmr
```

Refer [run_iqtree.sh](run_iqtree.sh) script for more details.

## Fix Tree

Refer [set_outgroup.py](set_outgroup.py) script for more details.

Refer [iqtree_internal_state_table_to_sequence.py](iqtree_internal_state_table_to_sequence.py) script for more details.

## Visualize Tree


