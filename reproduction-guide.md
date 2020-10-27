# Reproduce the results of "The origin, distribution, and genetic interactions of *KRAS* alleles across cancer types."

Follow the instructions below to reproduce all of the results of the paper.
Please open an [issue](https://github.com/jhrcook/comutation/issues) if you run into any problems.

### 1. Download the source code

Download or clone the repostory [jhrcook/comutation](https://github.com/jhrcook/comutation/tree/resubmission) from GitHub.

```bash
# Download the git repository.
git clone https://github.com/jhrcook/comutation.git
# Download as a zip.
wget https://github.com/jhrcook/comutation/archive/master.zip
unzip -q master.zip

```

### 2. Acquire the raw data

Unfortunately, due to data limits on GitHub, we are unable to upload the data with the source code.
Please contact Kevin M. Haigis (corresponding author) at `kevin_haigis at dfci dot harvard dot edu` for the raw data.


### 3. Install software

Ensure that all of the required software is installed.

The R packages and conda virtual environment can be installed using the following command.
Make sure to first set the project's directory as the current working directory.

```bash
./config/create-R-python-libs.sh
```

### 4. Configuration

There are multiple analyses that must be specifically configured to run on your system.

#### Snakemake configuration file for the Row-Column Test (RC-test) for reduced comutation

The RC-test is run in parallel using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow tool.
In the files ["src/20_21_run-rc-test-snakemake.sh"](src/20_21_run-rc-test-snakemake.sh) and ["src/20_63_nonallelespec_run-rc-test-snakemake.sh"](src/20_63_nonallelespec_run-rc-test-snakemake.sh), adjust the `snakemake` commands in accordance to your computing platform.
In particular, the `--drmaa` flag should be adjusted to match your platform's job-requesting command.
A full list of the command line options to `snakemake` can be found [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

#### GSEA

In the file ["src/10_35_gsea-depmap.sh"](src/10_35_gsea-depmap.sh), change the path in the variable `GSEA_PATH` to reflect the location of the '.jar' file on your computer.
For this analysis, GSEA was run as a batch array on a SLURM scheduler.
If that is not the system used on your computing platform, you will need to adjust the files ["src/10_34_submit-gsea-depmap.sh"](src/10_34_submit-gsea-depmap.sh) and ["src/10_35_gsea-depmap.sh"](src/10_35_gsea-depmap.sh) to align your platform.

### 5. Running analyses

The complete analysis is now ready to run using the following command.

```bash
./run-all-analyses.R
```
