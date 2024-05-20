# CSE185_Project_Group32

## Description

The Gibbs Sampler Motif Finding Tool is a command-line Python application designed to identify motifs within a given set of DNA sequences. It implements the Gibbs sampling algorithm, which iteratively selects k-mers from DNA sequences to find the most significant motifs. This tool is modeled after the [MEME](https://meme-suite.org/meme/doc/meme.html?man_type=web) (Multiple EM for Motif Elicitation) tool but specifically uses a randomized approach to motif discovery. 

## Installation

### Requirements
- Python 3.8 or higher
- NumPy (optional)

### Setup
Clone this repository to your local machine using: 


## Usage
To run the Gibbs Sampler, you need to provide a file containing DNA sequences. You will also need to specify the motif length (k), the number of sequences (t), and the number of iterations (n).

### Command-Line Arguments
- `-f, --file`: Path to the file containing DNA sequences.
- `-k, --motif_length`: Length of the motif to find.
- `-t, --num_sequences`: Number of sequences to use.
- `-n, --iterations`: Number of iterations for the Gibbs sampler.


