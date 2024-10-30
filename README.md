# ProFiT (PRObe FInding Tool)
Code to create probe set that hits target sequences without hitting negative sequences for the paper [A novel pipeline for exploring plasmid-bacterial hosts dynamics at the community and the single cell level in gut microbiomes]()

Usage example: `python main.py -t <targets fasta> -n <negatives fasta> -o <output directory> -l <length of probe> -d <max # degenerate bases in probe> -c <# probes to cover target> -f <max # times a negative sequence may be hit>`

Command line arguments
```
Output possible probe sets

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET, --target TARGET
                        fasta file of the target sequences to find probes for
                        (default: None)
  -n NEGATIVE, --negative NEGATIVE
                        fasta file of the negtive sequence set that probes
                        cannot hit (default: None)
  -o OUTDIR, --outdir OUTDIR
                        output directory (default: None)
  -l PROBE_LEN, --probe_len PROBE_LEN
                        length of probes (default: 20)
  -d MAX_DEGENERATE, --max_degenerate MAX_DEGENERATE
                        maximum number of degenerate bases allowed in a probe
                        (default: 1)
  -c MIN_COVERAGE, --min_coverage MIN_COVERAGE
                        minimum number of probes required to consider a target
                        as covered (default: 15)
  -f MAX_FALSE, --max_false MAX_FALSE
                        maximum number of times a negative sequence may be hit
                        (default: 3)
```
Example output files can be found in the files `probes_list.tsv` and `probes_hits.tsv`. `probes_list.tsv` is a **tab-separated** file listing in order the sequence of each probe, the sequence names of any new target sequences hit by adding that probe to the list and the sequnece names of any new negative sequence that would be hit by adding that probe to the set. `probes_list.tsv` is a **tab-separated** file listing the sequence of each probe and the sequence name and location of each positive and negative hit of the probe in the target and negative seqeuences.
