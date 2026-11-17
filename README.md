# TCRAFT
TCRAFT is a a scalable and cost-effective platform for pooled assembly of synthetic TCR libraries. This software package converts a table of user-inputted TCR sequences into CDR3 oligonucleotide pools to order for the TCRAFT assembly process. 

The TCRAFT package requires prior installation of [Python](https://www.python.org/downloads/) 3.7 or higher as well as an updated version of pip.

## Quick start: the TCRAFT application (no Python needed)

If you do not want to install anything, download the ready-to-run application from the
[Releases page](https://github.com/birnbaumlab/TCRAFT/releases) and double-click it.
It contains everything TCRAFT needs.

| System | File to download | First launch |
|---|---|---|
| Windows | `TCRAFT.exe` | Windows SmartScreen may warn that the publisher is unknown: click **More info**, then **Run anyway**. |
| macOS | `TCRAFT-macOS.zip` | Unzip it, then **right-click** `TCRAFT.app` and choose **Open**, then **Open** again. Double-clicking the first time shows a Gatekeeper warning instead. |
| Linux | `TCRAFT-Linux` | In a terminal: `chmod +x TCRAFT-Linux`, then run it. A desktop Linux already has what it needs; on a bare server install, add the Qt libraries first: `sudo apt install libegl1 libgl1 libxkbcommon-x11-0`. |

These warnings appear because the application is not code-signed; signing requires a paid
developer certificate. They only show up the first time.

Using the application:
1. **Generate oligos** tab: choose your TCR list (CSV). TCRAFT immediately shows how many TCRs
   it read, whether the required columns are present, which column holds your clonotype names,
   and how many duplicate TCRs it found.
2. Adjust the options if you wish, then press **Generate CDR3 oligos**. Progress and messages
   appear live.
3. When it finishes, a summary lists the oligos per pool; **Open results folder** takes you to
   the files to order.
4. The **Validate assembly** tab checks a generated oligo file by simulating Golden Gate assembly.


## Installation

### Option 1: Install from GitHub
This option requires you to have [Git](https://git-scm.com/downloads) installed.
1. Open Terminal or Command Prompt on your machine.
2. Enter the following command: `pip install git+https://github.com/birnbaumlab/TCRAFT.git`

### Option 2: Build the TCRAFT package locally
1. Download or clone the TCRAFT git repository locally.
2. In Terminal or Command Prompt, navigate to the base directory of the package (`/TCRAFT`) and execute the following command: `pip install .`
  
### Note: Using TCRAFT command line interface
In order to use the command line interface to run TCRAFT, you will need to either install TCRAFT inside of a [Conda](https://www.anaconda.com/download/success) environment or make sure that the folder where your Python packages are installed is added to your `$PATH` environment variable.

If you are unfamiliar with these steps or are not interested in this, TCRAFT can also be imported and run from within a Python script. Both the Python and Command Line versions are identical, and instructions are provided for both options below.

### Note: Install Dependencies Manually
By default, pip will automatically find and install all Python packages that are necessary for TCRAFT to run. In the event that dependency errors arise during installation, you can also manually install the requisite packages in your Python environment by running the following line in Terminal or Command Prompt:

`pip install numpy pandas tqdm biopython`

To also get the graphical interface when installing from source, use `pip install ".[gui]"` and start it with the `TCRAFT-gui` command.

TCRAFT uses the following package versions:
- numpy>=1.26.4
- pandas>=2.2.3
- tqdm>=4.67.0
- biopython>=1.84
- PySide6>=6.8 (graphical interface only)

## How to use the TCRAFT software package
Once installed, TCRAFT can be executed from the either the terminal or a Python script using two simple functions: TCRAFT-generate and TCRAFT-validate. TCRAFT-generate takes in a CSV-formatted list of TCR sequences to assemble and generates lists of oligonucleotide sequences to order in pooled format for the TCRAFT assembly process. TCRAFT-validate takes in the list of oligonucleotides created by TCRAFT-generate and simulates TCRAFT Golden Gate assembly of these oligos into complete TCR sequences. TCRAFT-validate then checks all of these simulated assemblies to ensure that the desired TCR sequence is correctly assembled and does not contain any errors.

Below are detailed instructions on how to use this package to prepare CDR3 oligo pools for your TCRAFT experiment:


1. **Prepare a CSV table of TCR sequences to assemble.** 
The CSV table *must* contain the following columns or else TCRAFT-generate will throw an error:
    - `V_alpha` (TCR $\alpha$ chain V region IMGT ID)
    - `V_beta`  (TCR $\beta$ chain V region IMGT ID)
    - `J_alpha` (TCR $\alpha$ chain J region IMGT ID)
    - `J_beta` (TCR $\beta$ chain J region IMGT ID)
    - `CDR3_alpha` (TCR $\alpha$ chain CDR3 amino acid sequence)
    - `CDR3_beta` (TCR $\beta$ chain CDR3 amino acid sequence)

    A sample TCR list is provided in this repository as an example of how to properly format and prepare the TCR list and can be used to test the TCRAFT pipeline.

2. **Generate CDR3 oligo pools by executing one of the following:**

    a. Command Line option:
        
        TCRAFT-generate /path/to/TCR/list.csv --output_dir /path/to/output/directory

    b. Python script option:
    
        from TCRAFT.Generate import generate_cdr3_oligos
        generate_cdr3_oligos(input_csv_path="/path/to/TCR/list.csv", output_dir="/path/to/output/directory")

    Note: The `output_dir` argument is optional and allows you to specify the output file directory. If not specified the output files will be stored in `./TCRAFT-generate_<current date>`. 

    TCRAFT-generate creates the following output files in the output directory:
    - `All_<Number of TCRs>_CDR3_oligos.csv`: Contains all successfully generated oligos
    - `Under350_All_CDR3_oligos.csv`: Contains only oligos with length <= 350bp (which can be ordered from suppliers such as Twist Biosciences).
    - `Over350_All_CDR3_oligos.csv`: Contains only oligos with length > 350bp (which can be ordered from alternative suppliers such as IDT).
    - `Under350_Pool<A or B>_CDR3_oligos.csv`: Split up the under 350bp oligos into pools that can be ordered together from an oligo pool supplier.
    - `Over350_Pool<A or B>_CDR3_oligos.csv`: Split up the oversize >350bp oligos into pools that can be ordered together from an oligo pool supplier.
    - `Oligo_Pool_Metadata.csv`: Table of oligo pool statistics which are helpful to plan out TCRAFT experimental protocols (ex: number of oligos per pool, average oligo size, etc.). This table can be directly pasted into the TCRAFT experiment planning workbook template to automatically compute resuspension volumes, mixing ratios, etc.
    - `Invalid_CDR3_Table.csv`: Contains TCRs that failed TCRAFT oligo generation. The most common reason for failure is due to specification of a V/J allele ID which is not supported by TCRAFT (TCRAFT only uses functional, coding alleles as annotated by IMGT).
    - `Under350_Pool<A or B>_CDR3_oligos_padded.csv`: The same pool, with every oligo padded at its 3' end up to the longest oligo of that pool, so that all oligos of a pool have the same length. The padding is taken from a shared random sequence and is checked never to create a Type-IIS restriction site or a homopolymer, including at the junction. The `Padding_added` and `Padding_status` columns record what happened to each oligo.
    - `Padding_report.csv`: Lists only the oligos that needed a padding window other than the shared one, or that could not be padded safely. Normally empty.
    - `All_TCR_deduplicated.csv`: The list of unique TCRs actually used, with the clonotype names as they appear in the oligo tables.
    - `All_TCR_removed.csv`: Audit trail of every duplicate TCR that was removed, together with `kept_clonotype`, the name of the clonotype that was kept in its place.

    **Clonotype names.** Every oligo carries a clonotype name in every output table. TCRAFT takes it from a column of your CSV named `clonotype`, `clonotype_id`, `clone_id`, `tcr_id`, `id`, `name` and similar (the match ignores case); use `--clonotype_id_column` to name the column yourself. If your file has no such column, names `TCR_0001`, `TCR_0002`, ... are generated. If two *different* TCRs share one name, the repeats are suffixed (`name`, `name_2`, ...) and a warning is printed.

    **Duplicate TCRs** are detected on the six required columns, ignoring case and surrounding or inner whitespace, so that `TRAV6 ` and `trav6` count as the same TCR. The first occurrence is kept.

    **Reproducibility.** By default TCRAFT uses a fixed random seed, so the same input file always yields exactly the same oligos. Pass `--seed <number>` to try another one if a few TCRs fail to generate, or `--no_seed` for a different result on every run. The seed used is recorded in `Oligo_Pool_Metadata.csv`.

3. **Validate the CDR3 oligo pools by executing one of the following:**

    a. Command Line option:
    
        TCRAFT-validate /path/to/TCRAFT/oligo/file.csv --output_dir /path/to/output/directory

    b. Python script option:
        
        from TCRAFT.Validate import validate_cdr3_oligos
        validate_cdr3_oligos(input_csv_path="path/to/TCRAFT/oligo/file.csv", output_dir="/path/to/output/directory")
        

    Note: Similar to above, the `output_dir` argument is optional. If not specified the output files will be stored in `./TCRAFT-validate_<current date>`. 

    TCRAFT-validate will simulate Golden Gate assembly of all inputted CDR3 oligos. If an oligo fails assembly validation, either due to presence of extraneous restriction sites or mismatch between the simulated and desired TCR protein sequence, the program will print the oligo/TCR that needs fixing. At the end, the program will output the percentage of inputted TCRs that passed assembly validation.

    TCRAFT-validate creates the following output files in the output directory:
    - `Reference_Assembled_TCR_Sequences.csv`: These are the full successfully validated TCR sequences obtained after `TCRAFT-validate` simulated TCRAFT assembly of the CDR3 oligos into full TCR sequences. This reference list is sometimes useful to cross-check downstream experimental sequencing data against the ground-truth list of assembled sequences.
    - `Failed_Assemblies.csv`: Contains the TCRs that failed assembly validation along with info that is helpful for debugging the failures.
    - `Assembly_Metadata.csv`: Similar to the `Oligo_Pool_Metadata.csv` file, this table contains useful statistics for experimental planning and can also be directly pasted into the TCRAFT experiment planning workbook template.

## Version
1.2.1

## Author
[Rachit Mukkamala](mailto:rsmukk@mit.edu)
