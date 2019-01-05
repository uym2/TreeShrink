## `pasta_treeshrink_pipeline.py`

This script creates three scripts that can be used to run the following pipeline:

1. Run PASTA (default: 1 iteration)
2. Run TreShrink (default: `-b 33 -q 0.05`
3. Run PASTA on TreeShrink filtered sequences

Prerequisites:
1. Have TreeShrink 1.2.1 or later
2. Have PASTA 1.8.2 or later

You need to have a directory with a structure like this:

~~~
avian/1000/:
input.faa

avian/10000/:
input.faa

avian/10001/:
input.faa

avian/10002/:
input.faa
~~~

Then, you can run:

* For DNA:

~~~
python pasta_treeshrink_pipeline.py -i avian -a input.faa
~~~

* for AA sequences:

~~~
python pasta_treeshrink_pipeline.py -i avian -a input.faa -p
~~~

This will create three files:

1. `firststep.sh`: all commands needed to run the initial PASTA jobs
2. `secondstep.sh`: the command needed to run TreeShrink
3. `thirdstep.sh`: all commands needed to run the second PASTA job. 

The lines in the first and third line should be run independently. You can do this using `bash thirdstep.sh` to run the jobs in serial, using `parallel < thirdstep.sh` to run them in parallel, or if you have a cluster, you can use whatever system your cluser provides. 

At the end, under each gene directory, you will have something like the following. The most important files are:

- `filtered.tre`: the final filtered gene tree
- `filtered.aligned.fasta`: the final filtered alignment 
- `firstround_shrunk_RS_0.05.txt`: the list of filtered species

~~~
avian/10052:
filtered.aligned.fasta                                            filtered_temp_iteration_initialsearch_seq_unmasked_alignment.gz   firstround_shrunk_RS_0.05.txt
filtered.err.txt                                                  filtered_temp_iteration_initialsearch_tree.tre                    firstround_temp_iteration_0_seq_alignment.txt
filtered.out.txt                                                  filtered_temp_name_translation.txt                                firstround_temp_iteration_0_seq_unmasked_alignment.gz
filtered.score.txt                                                filtered_temp_pasta_config.txt                                    firstround_temp_iteration_0_tree.tre
filtered.tre                                                     firstround.err.txt                                                firstround_temp_iteration_initialsearch_seq_alignment.txt
filtered_temp_iteration_0_seq_alignment.txt                       firstround.marker001.input.faa.aln                                firstround_temp_iteration_initialsearch_seq_unmasked_alignment.gz
filtered_temp_iteration_0_seq_unmasked_alignment.gz               firstround.marker001.input.faa.idx                                firstround_temp_iteration_initialsearch_tree.tre
filtered_temp_iteration_0_tree.tre                                firstround.marker001.input.faa_shrunk0.05.aln                    firstround_temp_name_translation.txt
filtered_temp_iteration_1_seq_alignment.txt                       firstround.out.txt                                                firstround_temp_pasta_config.txt
filtered_temp_iteration_1_seq_unmasked_alignment.gz               firstround.score.txt                                              input.faa
filtered_temp_iteration_1_tree.tre                                firstround.tre
filtered_temp_iteration_initialsearch_seq_alignment.txt           firstround_shrunk_0.05.tre
~~~

Here is the help of the script:

~~~bash 
> python related_scripts/pasta_treeshrink_pipeline.py -h
usage: pasta_treeshrink_pipeline.py [-h] -i INDIR [-a INPUT] [-1 PTARGS1]
                                    [-2 TSARGS] [-3 PTARGS2] [-n NS] [-p] [-c]

optional arguments:
  -h, --help            show this help message and exit
  -i INDIR, --indir INDIR
                        The parent input directory where the trees (and
                        alignments) can be found
  -a INPUT, --input INPUT
                        The name of the input sequence files. Each
                        subdirectory under it must contain a sequence file
                        with this name. Default: input.fasta
  -1 PTARGS1, --ptargs1 PTARGS1
                        A set of arguments passed to first PASTA run; e.g., '
                        --iter-limit 1 --mask-gappy-sites=20 -d Protein'
  -2 TSARGS, --tsargs TSARGS
                        A set of arguments passed to tree shrink; e.g., '-q
                        0.01 -b 20'
  -3 PTARGS2, --ptargs2 PTARGS2
                        A set of arguments passed to first PASTA run; e.g., '
                        --iter-limit 3 --mask-gappy-sites=3 -d Protein'
  -n NS, --ns NS        Number of species'
  -p, --protein         Are inputs proteins
  -c, --cleanup         Cleanup the directories
~~~