# G4Boost
G-quadruplex identification and secondary structure stability prediction tool 

G4Boost requires python3.

Prediction models was constructed using the sklearn module. Prebuilt regression and classification models are available through GitHub.\n


	required python packages:
		pandas
		xgboost

Program Usage

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
			      The maximum number of fleaxible loops separating the G-stems [default=11]
	--noreverse
			      Do not search the reverse complement of the input fasta.
	-q, --quiet
			      Do not print progress report (i.e. sequence names as they are scanned).



--------------------
Prepare an input file (FASTA):
	head test_sequence.fa
	>test_seq_1
  GGGCAGAAGGGAGGGCTGGGG

Locate regression and classification models into the same directory with G4Boost.py


--------------------
To run G4Boost:

	python3 G4Boost.py -f test_sequence.fa


--------------------
The outputs include two files
1) a gff file for the position of the putative G4s
2) a csv file for the secondary structure predcition scores

	head test_sequence.fa.gff
	test_seq_1	0	21	test_seq_1_0_21	21	+	GGGCAGAAGGGAGGGCTGGGG
  
  head test_sequence.fa.g4scores.csv
  seq	seq_length	g4motif	length	loops	G-quartet	maxlbase	minlbase	G	C	GG	CC	g4_pred	g4_prob	mfe_pred
  test_seq_1	21	GGGcagaaGGGaGGGctgGGG	21	3	3	5	1	66	9	42	0	1	0.99999976	-13.413148
 
--------------------
To change parameters:

	python3 G4Boost.py -f test_sequence.fa --maxloop 7 --minloop 1 --maxG 4 --minG 3 --loops 3 --noreverse --quiet


--------------------
