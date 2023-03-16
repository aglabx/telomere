# Telomere: Eukaryotic Telomeres Annotation Tool

**Telomere** is a powerful and user-friendly tool designed for annotating telomeres in eukaryotic genome assemblies. The tool employs fast algorithm to identify and characterize telomeric sequences, providing valuable insights into the structure and function of telomeres in various eukaryotic organisms.

## Table of Contents

- Overview
- Installation
- Usage
- Output
- Contributing
- Citation
- License

## Overview

Telomeres are essential for the stability and integrity of linear eukaryotic chromosomes. They are typically composed of tandem repeats of short DNA sequences that protect the ends of chromosomes from degradation and fusion. Accurate annotation and visualization of telomeres is crucial for understanding their roles in various biological processes, such as aging and cancer development. Telomere Annotator is designed to facilitate this process by providing a comprehensive and accurate telomere annotation solution for eukaryotic genome assemblies.

Key features of *Telomere* include:

- Identification of telomeric sequences in both assembled and unassembled eukaryotic genomes
- Characterization of telomere structure, including repeat composition and length
- Estimation of telomere length from sequencing data
- Visualization tools for easy interpretation of results
- A new quality control metric for T2T assemblies

# Installation

To install *Telomere*, please follow the instructions below:

Clone the repository:

```bash
git clone https://github.com/aglabx/telomere.git
```

Change to the repository directory:

```bash
cd telomere
```

Install the required dependencies using the provided requirements file:

```bash
pip install -r requirements.txt
```

Run the setup script to install Telomere Annotator:

```bash
python setup.py install
```

## Usage

After installing Telomere Annotator, you can use it with the following command:

```bash
telomere-annotator [options] -i <input_file> -o <output_directory>
```

### Options

-i, --input: Path to the input genome assembly file (in FASTA format)
-o, --output: Path to the output directory where results will be saved
--min-repeat-length: Minimum length of telomeric repeats (default: 5)
--max-repeat-length: Maximum length of telomeric repeats (default: 8)
--min-telomere-length: Minimum length of telomeric sequences (default: 1000)
--monomer: Telomere monomer sequence (TTAGGG)
-v, --verbose: Display verbose output


## Output

Telomere Annotator generates the following output files:

- telomere_summary.tsv: A tab-separated file containing a summary of identified telomeres, including their location, repeat composition, and length.
- telomere_sequences.fasta: A FASTA file containing the extracted telomeric sequences.
- telomere_annotation.gff3: A GFF3 file containing the genomic coordinates of identified telomeres.
- telomere_report.html: An interactive HTML report containing a summary of results, visualizations, and comparative analysis
- telomere_karyotype.svg: A visualization of telomere localization on chromosomes

## Contributing

We welcome contributions to improve and expand Telomere. To contribute, please follow these steps:

- Fork the repository.
- Create a new branch with a descriptive name.
- Make your changes and commit them to the branch.
- Open a pull request with a clear description of your changes.


## License

Telomere is released under the GPLv3 License.

## Citation

If you use Telomere in your research, please cite the following publication:

Daniil Listas, Aleksey Komissarov (2023). Telomere: A Robust Annotation and Quality Control Toolkit for Telomeres in Chromosome-level Eukaryotic Assemblies. xxx, vol. xx, no. x, pp. xxx-xxx.

## Contact

If you have any questions or suggestions, please feel free to open an issue on GitHub or contact us at ad3002@gmail.com.

