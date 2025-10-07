# BioPhi CDR Grafting CLI Command - Implementation Guide

## Overview

This document describes the implementation of a new `biophi cdrgraft` command-line interface for antibody humanization using CDR grafting methodology. CDR grafting is an alternative humanization approach to Sapiens that transplants CDR loops from input antibodies onto human germline frameworks.

## Implementation Date
October 7, 2025

## Files Created/Modified

### 1. New CLI Command Implementation

**File Locations:**
- Conda CLI: `/opt/miniforge3/envs/biophi_env/lib/python3.12/site-packages/biophi/humanization/cli/cdrgraft.py`
- Docker Source: `/home/dwang/Desktop/biophi/BioPhi/biophi/humanization/cli/cdrgraft.py`

**Description:** Complete CLI implementation for CDR grafting with support for:
- FASTA-only output (fast mode)
- Full report generation with OASis humanness evaluation
- Interactive mode for single antibody humanization
- Automatic or manual germline selection
- Optional Vernier zone backmutations
- Optional Sapiens refinement iterations

### 2. CLI Main Entry Point Updates

**File Locations:**
- Conda CLI: `/opt/miniforge3/envs/biophi_env/lib/python3.12/site-packages/biophi/common/cli/main.py`
- Docker Source: `/home/dwang/Desktop/biophi/BioPhi/biophi/common/cli/main.py`

**Changes Made:**
```python
# Added import
from biophi.humanization.cli.cdrgraft import cdrgraft

# Registered new command
main.add_command(cdrgraft)
```


## Installation Instructions

### Patching an Existing Conda BioPhi Installation

If you have BioPhi installed via conda and want to add the CDR grafting functionality, follow these steps:

#### Step 1: Clone or Download This Repository

```bash
# Clone the repository
git clone https://github.com/qdev1915/BioPhi.git
cd BioPhi
```

#### Step 2: Identify Your Conda Environment

First, find where your BioPhi conda environment is located:

```bash
# Activate your BioPhi environment
conda activate biophi_env  # or your environment name

# Find the Python site-packages location
python -c "import site; print(site.getsitepackages()[0])"
```

This will output something like:
```
/opt/miniforge3/envs/biophi_env/lib/python3.12/site-packages
```

#### Step 3: Copy the CDR Grafting CLI File

Copy the new `cdrgraft.py` command to your conda environment:

```bash
# Set your site-packages path (adjust based on Step 2 output)
SITE_PACKAGES="/opt/miniforge3/envs/biophi_env/lib/python3.12/site-packages"

# Copy the CDR grafting CLI file
sudo cp biophi/humanization/cli/cdrgraft.py \
    ${SITE_PACKAGES}/biophi/humanization/cli/cdrgraft.py

# Set proper permissions
sudo chmod 644 ${SITE_PACKAGES}/biophi/humanization/cli/cdrgraft.py
```

#### Step 4: Update the Main CLI Entry Point

Update the `main.py` file to register the new command:

```bash
# Backup the original file
sudo cp ${SITE_PACKAGES}/biophi/common/cli/main.py \
    ${SITE_PACKAGES}/biophi/common/cli/main.py.backup

# Copy the updated main.py
sudo cp biophi/common/cli/main.py \
    ${SITE_PACKAGES}/biophi/common/cli/main.py

# Set proper permissions
sudo chmod 644 ${SITE_PACKAGES}/biophi/common/cli/main.py
```

**Alternatively**, you can manually edit the file:

```bash
sudo nano ${SITE_PACKAGES}/biophi/common/cli/main.py
```

Add the import line after the other imports:
```python
from biophi.humanization.cli.cdrgraft import cdrgraft
```

Add the command registration after `main.add_command(sapiens)`:
```python
main.add_command(cdrgraft)
```

#### Step 5: Verify Installation

Test that the command is available:

```bash
biophi cdrgraft --help
```

You should see the help message for the CDR grafting command.

#### Step 6: Test with Sample Data

```bash
# Create a test directory
mkdir -p ~/biophi_test
cd ~/biophi_test

# Create a sample FASTA file (or use your own)
cat > test_antibody.fa << 'FASTA'
>test_VH
QVQLQESGPGLVKPSETLSLTCTVSGGSISSGGYYWSWIRQPPGKGLEWIGYIYYSGSTYGNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCARDTFWSGYYYGMDVWGQGTTVTVSS
>test_VL
DIQMTQSPSSLSASVGDRVTITCRASQGISSALAWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQFNSYPLTFGGGTKVEIK
FASTA

# Run CDR grafting
biophi cdrgraft test_antibody.fa --fasta-only --output cdrgraft_output.fa

# View results
cat cdrgraft_output.fa
```

#### Troubleshooting Installation

**Permission Denied Errors:**
- Use `sudo` when copying files to system directories
- Ensure you have write permissions to the conda environment

**Command Not Found:**
- Verify your conda environment is activated: `conda activate biophi_env`
- Check file was copied correctly: `ls -l ${SITE_PACKAGES}/biophi/humanization/cli/cdrgraft.py`
- Restart your terminal or re-activate the conda environment

**Import Errors:**
- Clear Python cache: `find ${SITE_PACKAGES} -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null`
- Verify `main.py` was updated correctly
- Check for syntax errors in the copied files

### Installing from Source (Development Mode)

If you want to install BioPhi with CDR grafting support from source:

```bash
# Clone this repository
git clone https://github.com/qdev1915/BioPhi.git
cd BioPhi

# Create a new conda environment (optional but recommended)
conda create -n biophi_dev python=3.12
conda activate biophi_dev

# Install in development mode
pip install -e .

# Install dependencies if needed
pip install -r requirements.txt

# Verify installation
biophi cdrgraft --help
```

### Docker Installation

For Docker deployments, rebuild the image from this repository:

```bash
# Clone this repository
git clone https://github.com/qdev1915/BioPhi.git
cd BioPhi

# Build Docker image
docker build -t biophi:cdrgraft .

# Run container
docker run -it biophi:cdrgraft biophi cdrgraft --help
```


## Command Syntax

```bash
biophi cdrgraft [OPTIONS] [INPUTS]...
```

### Arguments

- `INPUTS` - Input FASTA file path(s). If not provided, creates an interactive session.

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--output` | TEXT | stdout | Output directory path. With --fasta-only, output FASTA file path. |
| `--fasta-only` | FLAG | False | Output only a FASTA file with humanized sequences (speeds up processing) |
| `--oasis-db` | TEXT | None | OAS peptide database connection string (required to run OASis) |
| `--scheme` | TEXT | kabat | Numbering scheme: one of imgt, aho, chothia, kabat |
| `--cdr-definition` | TEXT | kabat | CDR definition: one of imgt, chothia, kabat, north |
| `--heavy-v-germline` | TEXT | auto | Heavy chain V germline gene (auto for automatic selection) |
| `--light-v-germline` | TEXT | auto | Light chain V germline gene (auto for automatic selection) |
| `--backmutate-vernier` | FLAG | True | Backmutate Vernier zone residues to parental |
| `--sapiens-iterations` | INTEGER | 0 | Additional Sapiens iterations after CDR grafting |
| `--limit` | INTEGER | None | Process only first N records |

## Usage Examples

### 1. Basic CDR Grafting (FASTA output only)

```bash
biophi cdrgraft input.fa --fasta-only --output humanized.fa
```

**Output:** Single FASTA file with humanized sequences

### 2. CDR Grafting with Full Report

```bash
biophi cdrgraft input.fa --output ./report/ \
  --oasis-db /opt/biophi/data/OASis_9mers_v1.db
```

**Output:** 
- `humanized.fa` - Humanized sequences
- `CDRGraft_report.xlsx` - Excel report with OASis scores and mutation details

### 3. CDR Grafting with Sapiens Refinement

```bash
biophi cdrgraft input.fa \
  --sapiens-iterations 2 \
  --output ./report/ \
  --oasis-db /opt/biophi/data/OASis_9mers_v1.db
```

**Output:** CDR grafted sequences followed by 2 iterations of Sapiens refinement

### 4. Custom Germline Selection

```bash
biophi cdrgraft input.fa \
  --heavy-v-germline IGHV3-23*01 \
  --light-v-germline IGKV1-39*01 \
  --output ./report/
```

**Output:** CDRs grafted onto specified human germlines

### 5. CDR Grafting with Different Numbering Schemes

```bash
# Using IMGT numbering and CDR definition
biophi cdrgraft input.fa \
  --scheme imgt \
  --cdr-definition imgt \
  --fasta-only \
  --output humanized_imgt.fa

# Using Chothia numbering and CDR definition
biophi cdrgraft input.fa \
  --scheme chothia \
  --cdr-definition chothia \
  --fasta-only \
  --output humanized_chothia.fa
```

### 6. Straight Graft (No Vernier Backmutations)

```bash
biophi cdrgraft input.fa \
  --no-backmutate-vernier \
  --fasta-only \
  --output humanized_straight.fa
```

**Note:** By default, `--backmutate-vernier` is enabled. Use `--no-backmutate-vernier` to disable.

### 7. Interactive Mode

```bash
biophi cdrgraft
```

Prompts for antibody sequences interactively.

## CDR Grafting Methodology

### How It Works

1. **Input Parsing**: Reads VH and/or VL sequences from FASTA files
2. **Germline Selection**: Automatically selects the most similar human germline or uses specified germline
3. **CDR Grafting**: Transplants CDR loops from input sequence onto human germline framework
4. **Vernier Zone Backmutation** (optional): Restores key framework residues that support CDR structure
5. **Sapiens Refinement** (optional): Further optimizes humanization using deep learning
6. **OASis Evaluation** (optional): Assesses humanness using 9-mer peptide database

### Comparison with Sapiens

| Feature | Sapiens | CDR Grafting |
|---------|---------|-------------|
| Method | Deep learning mutations | Framework replacement |
| CDR Preservation | Optional | Always preserved |
| Framework Source | Optimized by AI | Human germline |
| Vernier Zone | Not specifically addressed | Optionally backmutated |
| Speed | Fast (1-3 iterations typical) | Fast (single operation) |
| Customization | Iteration count | Germline selection, Vernier, Sapiens refinement |

## Output Files

### FASTA Only Mode (`--fasta-only`)

**File:** Single FASTA file with humanized sequences

**Example:**
```
>test1_VH VH (Humanized test1_VH CDR_Grafted_kabat_Vernier_None BioPhi)
QVQLQESGPGLVKPSDTLSLTCTVSGFSLSSYGVSWIRQPPGKGLEWIGIIWDDGATDYA
SWAKSRSTISRDTSKNQVSLKLSSVTAADTAVYYCEKGGSAYIWGQGTMVTVSS
>test1_VL VL (Humanized test1_VL CDR_Grafted_kabat_Vernier_None BioPhi)
AQQLTQSPSSLSASVGDRVTITCQASQNVYSNNRLAWYQQKPGKAPKQLIYGASNLASGV
PSRFSGSGSGADFTLTISSLQPEDFATYYCLGEFNCYAGDCFAWGPGTKVDIK
```

### Full Report Mode (default)

**Files:**
1. `humanized.fa` - Humanized sequences in FASTA format
2. `CDRGraft_report.xlsx` - Excel workbook with:
   - **Overview Sheet**: Antibody names, OASis scores, germlines, mutation counts

## Implementation Details

### Core Functions

#### `cdrgraft()`
Main CLI command function that handles argument parsing and delegates to appropriate processing function.

#### `cdrgraft_fasta_only()`
Processes input files and outputs only humanized FASTA sequences without OASis analysis or Excel reports.

#### `cdrgraft_full_report()`
Processes input files, generates humanized sequences, runs OASis analysis (if database provided), and creates Excel reports.

#### `cdrgraft_interactive()`
Interactive mode for single antibody humanization with real-time feedback.

### Core Humanization Logic

The CLI uses existing BioPhi core functions:

```python
from biophi.humanization.methods.humanization import (
    humanize_antibody, 
    CDRGraftingHumanizationParams
)

# Create parameters
params = CDRGraftingHumanizationParams(
    scheme='kabat',
    cdr_definition='kabat',
    heavy_v_germline='auto',
    light_v_germline='auto',
    backmutate_vernier=True,
    sapiens_iterations=0
)

# Humanize
result = humanize_antibody(vh=vh_chain, vl=vl_chain, params=params)
```

## Testing

### Test Command

```bash
cd /home/dwang/test
conda run -n biophi_env biophi cdrgraft input/test1.fa \
  --fasta-only \
  --output output/test1_cdrgraft.fa
```

### Expected Output

```
Settings:
- Humanization method: CDR Grafting
- Numbering scheme: kabat
- CDR definition: kabat
- Heavy V germline: auto
- Light V germline: auto
- Backmutate Vernier zone: Yes

Reading input files...
Processing 2 sequences...
Humanizing: 100%|██████████| 2/2 [00:02<00:00,  1.05s/it]
Writing output to output/test1_cdrgraft.fa...
Completed: 2 sequences humanized
```

## Installation Instructions

### For Conda Environment (CLI)

Files are already installed in:
```
/opt/miniforge3/envs/biophi_env/lib/python3.12/site-packages/biophi/
```

No additional installation needed. Command is immediately available:
```bash
conda activate biophi_env
biophi cdrgraft --help
```

### For Docker Image

Files are copied to Docker source location:
```
/home/dwang/Desktop/biophi/BioPhi/biophi/
```

To include in Docker image:
1. Rebuild Docker image from source
2. The new command will be available in the container

```bash
cd /home/dwang/Desktop/biophi/BioPhi
docker build -t biophi:latest .
```

## Troubleshooting

### Command Not Found

**Issue:** `biophi cdrgraft` command not recognized

**Solution:** 
1. Verify conda environment is activated: `conda activate biophi_env`
2. Check if file exists: `ls /opt/miniforge3/envs/biophi_env/lib/python3.12/site-packages/biophi/humanization/cli/cdrgraft.py`
3. Verify main.py includes the import and registration

### Module Import Errors

**Issue:** `ModuleNotFoundError: No module named 'biophi.humanization.cli.cdrgraft'`

**Solution:**
1. Ensure `cdrgraft.py` has correct permissions: `chmod 644 cdrgraft.py`
2. Clear Python cache: `find /opt/miniforge3/envs/biophi_env -type d -name __pycache__ -exec rm -rf {} +`
3. Restart terminal or re-activate conda environment

### Logo Display Error

**Issue:** `AssertionError: Logos need to have same amount of lines`

**Solution:** Ensure the ASCII logo has exactly 6 lines (including final empty line) to match BioPhi's logo format requirements.

## Future Enhancements

Potential improvements for future versions:

1. **Batch Processing**: Add parallel processing for large FASTA files
2. **Progress Tracking**: Enhanced progress bars with time estimates
3. **Germline Database**: Interactive germline selection with descriptions
4. **Comparison Mode**: Side-by-side comparison of Sapiens vs CDR Grafting results
5. **Web UI Integration**: Add CDR grafting to existing BioPhi web interface
6. **Custom Mutations**: Allow manual specification of additional backmutations
7. **Multi-format Output**: Support for additional output formats (JSON, CSV)

## References

- BioPhi Documentation: https://github.com/Merck/BioPhi
- Antibody Numbering: http://www.bioinf.org.uk/abs/
- OASis Database: https://academic.oup.com/bioinformatics/article/35/14/i576/5529249

## Support

For questions or issues:
1. Check BioPhi GitHub issues: https://github.com/Merck/BioPhi/issues
2. Review BioPhi documentation
3. Contact BioPhi maintainers

## Changelog

### Version 1.0 (October 7, 2025)
- Initial implementation of `biophi cdrgraft` command
- Support for FASTA-only and full report modes
- Integration with OASis humanness evaluation
- Automatic and manual germline selection
- Optional Vernier zone backmutations
- Optional Sapiens refinement after grafting
- Interactive mode support
