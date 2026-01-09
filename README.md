# G4Boost
This repository contains documentation code for G4 quadruplex prediction tool [G4Boost](https://link.springer.com/article/10.1186/s12859-022-04782-z).

## Overview
G4Boost is a machine learning–based framework for predicting the formation and stability of G-quadruplex structures from nucleotide sequence data.

The framework was originally developed and validated for biologically relevant G4 prediction tasks and has since been applied and cited in high-impact studies investigating both RNA and DNA G-quadruplex biology.

## Why G4Boost?
G4Boost was designed to address key limitations of motif-based G-quadruplex predictors:
- Predicts G-quadruplex formation and relative structural stability beyond simple motif screening
- Supports both classification and regression-style prediction tasks
- Enables flexible feature engineering and model retraining
- Has been independently validated and cited by subsequent high-impact studies
- Suitable for integration into bioinformatics pipelines and downstream analyses

## Supported G-Quadruplex Types
G4Boost supports prediction tasks involving both RNA and DNA G-quadruplex–forming sequences.
- RNA G-quadruplexes: Commonly analyzed in transcriptomic and post-transcriptional regulatory contexts, including translation, splicing, and RNA stability.
- DNA G-quadruplexes: Frequently studied in genomic regulatory regions such as promoters, enhancers, and replication-associated loci to investigate transcriptional control, genome organization, and replication dynamics.

As with any computational prediction method, G4Boost outputs are best interpreted in the context of complementary experimental data and biological knowledge.


--------------------
## Installation
G4Boost requires Python 3 and the dependencies listed in requirements.txt.

### Clone the repository and obtain the source code
```
git clone https://github.com/hbusra/G4Boost.git
cd G4Boost
```
Alternatively, users may download the repository as a ZIP file from GitHub (Click Code → Download ZIP) and extract it locally instead of using git clone.


### Create a dedicated environment and install the required dependencies

#### Option 1: Conda (MacOS)

Using conda is the most reliable installation method across platforms and is strongly recommended for macOS users due to OpenMP requirements of XGBoost.

```bash
conda env create -f environment.yml
conda activate g4boost_env
```

#### Option 2: pip (Linux)

```bash
pip install -r requirements.txt
```

Packaging and simplified installation options (pip and conda) are planned for future releases.


## Quick Start

Prepare an input file (FASTA):
```
head test_sequence.fa
>test_seq_1
GGGCAGAAGGGAGGGCTGGGG
```

To run G4Boost:
```
python3 G4Boost.py -f test_sequence.fa
```


The outputs include two files

1) a gff file for the position of the putative G4s
2) a csv file for the secondary structure prediction scores

```
head test_sequence.fa.gff
test_seq_1	0	21	test_seq_1_0_21	21	+	GGGCAGAAGGGAGGGCTGGGG

head test_sequence.fa.g4scores.csv
seq	seq_length	g4motif	length	loops	G-quartet	maxlbase	minlbase	G	C	GG	CC	g4_pred	g4_prob	mfe_pred
test_seq_1	21	GGGcagaaGGGaGGGctgGGG	21	3	3	5	1	66	9	42	0	1	0.99999976	-13.413148
```
 
To change parameters:
```
python3 G4Boost.py -f test_sequence.fa --maxloop 7 --minloop 1 --maxG 4 --minG 3 --loops 3 --noreverse --quiet
```

--------------------
## Documentation
Prediction models was constructed using the sklearn module. Prebuilt regression and classification models are available through GitHub.
Locate regression and classification models (G4Boost_regressor.json and G4Boost_classifier.json) into the same directory with G4Boost.py



	required python packages:
		pandas==1.4.3
		xgboost==1.4.2


	Options:
	--version             show program's version number and exit  
	-h, --help            show this help message and exit
	-f FASTA_FILE, --fasta FASTA_FILE
			      Query sequences in fasta format
	--classifier JSON_FILE
			      Classification model in .json format [optional]
	--regressor JSON_FILE
			      Regression model in .json format [optional]
	-N MAX_LOOP_LENGTH, --maxloop MAX_LOOP_LENGTH
			      The maximum number of bases allowed within a loop region [default=12]
	-n MIN_LOOP_LENGTH, --minloop MIN_LOOP_LENGTH
			      The minimum number of bases allowed within a loop region [default=1]
	-G MAX_G_BASES, --maxG MAX_G_BASES
			      The maximum number of consecutive guanine bases allowed within a G-stem [default=7]
	-g MIN_G_BASES, --minG MIN_G_BASES
			      The minimum number of consecutive guanine bases allowed within a G-stem [default=1]
	-l NUMBER_OF_LOOPS, --loops NUMBER_OF_LOOPS
			      The maximum number of flexible loops separating the G-stems [default=11]
	--noreverse
			      Do not search the reverse complement of the input fasta.
	-q, --quiet
			      Do not print progress report (i.e. sequence names as they are scanned).


## Citation
If you use G4Boost in your research, please cite:

Cagirici H.B. et al.  
G4Boost: a machine learning framework for G-quadruplex prediction  
BMC Bioinformatics (2022)


Additional studies citing G4Boost include publications in [Nature Communications (2024)](https://www.nature.com/articles/s41467-024-54958-9) and [Cell (2024)](https://www.cell.com/cell/fulltext/S0092-8674(24)01134-6).

