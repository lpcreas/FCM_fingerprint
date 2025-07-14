# Python implementation of the flow cytometric fingerprinting method proposed by Props et al. (2016)

Props, R., Monsieurs, P., Mysara, M., Clement, L. and Boon, N. (2016), 'Measuring the biodiversity of microbial communities by flow cytometry', _Methods in Ecology and Evolution_ __7__(11), 1376--1385.

This script processes flow cytometry data to generate phenotypic fingerprints and calculate diversity indices
using Hill numbers. The script is structured into several processing steps, including data transformation,
subsampling, normalization, kernel density estimation (KDE), and diversity index calculation.

## Modules

- `numpy`: For numerical operations.
- `matplotlib.pyplot`: For plotting.
- `FlowCal`: For flow cytometry data handling and gating.
- `scipy.stats`: For statistical operations, including KDE.

## Constants

- `datapath`: Path to the input `.fcs` file containing flow cytometry data.
- `figpath`: Path to save generated figures.
- `CHANNEL_ID`: Dictionary mapping channel identifiers to their respective names.
- `MIN_SAMPLE_AMOUNT`: Minimum number of samples for subsampling.
- `THRESH_AOBNOB`: Threshold for gating AOB and NOB populations.

## Functions

- `get_fingerprinting_channels(elem)`: Extracts relevant channels (FSC, SSC, C3) for fingerprinting.
- `subsample_fingerprint_channels(elem)`: Subsamples the data to a minimum sample amount.
- `log_minmax_normalize_fingerprint_channels(elem)`: Applies log transformation and min-max normalization.
- `calculate_2D_kde(m1, m2)`: Computes 2D kernel density estimation for two variables.
- `minmax_normaliser(Z)`: Normalizes a 2D KDE to the `[0, 1]` range.
- `get_kde_fingerprint(elem, plot=False)`: Generates a phenotypic fingerprint by concatenating normalized KDEs.
- `hill_number(signal, q=0)`: Calculates the Hill number (diversity index) for a given signal.
- `calculate_hill_numbers(elem, orders=[0, 1, 2])`: Computes Hill numbers for specified orders.

## Processing Steps

1. Data transformation to relative fluorescence intensity (RFI).
2. Gating to filter data based on aspect ratio and fluorescence thresholds.
3. Extraction of fingerprinting channels.
4. Subsampling to ensure uniform sample size.
5. Log transformation and normalization of data.
6. KDE computation for phenotypic fingerprinting.
7. Calculation of Hill numbers to quantify diversity.

## Outputs

- Phenotypic fingerprints for AOB and NOB populations.
- Hill numbers (diversity indices) for AOB and NOB populations.

## Usage

1. Ensure the required `.fcs` file is available at the specified `datapath`.
2. Install the required libraries using `pip` (Python 3.9) and the `requirements.txt` file:

    ```bash
    pip install -r requirements.txt
    ```

3. Adjust constants (e.g., `MIN_SAMPLE_AMOUNT`, `THRESH_AOBNOB`) as needed.
4. Run the script to process the data and generate outputs.
