# Primer Finder
A small python script to find out which primers match your sequence.

This uses fuzzy regex matching to determin alignments, allowing only substitution errors. Only primers matching a single time at the best score will be reported.

## Dependencies
Python 3 (tested with 3.9) and the following non-standard packages (available via Conda):
- [Regex](https://pypi.org/project/regex/), an alternative package to perform regular expression matching. Tested using verion 2022.9.13.
- [Plotly](https://github.com/plotly/plotly.py), for generating the figure. Tested using version 5.15.0.
- [Biopython](https://biopython.org/) for reverse-complementing because I am lazy. Tested using version 1.79.

## Usage
The only absolute requirement is a primer data file in the same folder as the script. The primer file is tab-separated with double-quote-enclosed fields. The primer name is field 2, the primer sequence is field 8, the presence of a dye is in field 17 (1-indexed). Run using:

```python <path to primer_finder.py> ```

This will look for a primer data file in the same folder as the script, read it in and then prompt you to paste in the sequence to analyse. Warning: the max amount of characters for the input prompt may be limited by your terminal (ca 1000 for macOS).

### Additional parameters
In addition to the minimal run shown above, you can tweak a few aspects using parameters. A brief overview and description of parameters can be obtained by running:

```python <path to primer_finder.py> -h```

#### Input parameters:
```-p <path>```	path to the primer data file, use this to point the script to a text file containing the data exported from the primer database.

```-s <path>```	ideally this is a path to a text file containing your sequence in fasta format. You can also just paste in your sequence instead, though be wary of size limits.

#### Analysis parameters:
```-e <N>```	this is the edge parameter, defining how many nucleotides must be a perfect match at the 3’ end. By default, this is 5.

```-mm <N>```	this is the mismatch parameter, defining how many mismatches are allowed in the primer alignment. By default, this is 2.

```-d <N>```	this is the degenerate parameter, defining how many degenerate positions are allowed in the primer. By default, this is 3.

#### Output parameters:
```-o <path>```	this is the folder you want the results to be saved to. If you don’t specify this parameter, you will only get the temporary graph displayed in your browser.

```-go```		only in conjunction with -o, save only the graph output.

```-to```		only in conjunction with -o, save only the text output.

