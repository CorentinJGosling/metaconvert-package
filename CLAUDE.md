# metaConvert R Package Documentation

This document provides a comprehensive overview of the `metaConvert` R package, detailing its functionality, structure, and intended workflow. The analysis is based on an automated review of the source code.

## Core Purpose

The `metaConvert` package is a specialized toolkit for researchers conducting systematic reviews and meta-analyses. Its primary function is to automate the calculation and conversion of various effect size (ES) measures from a wide array of summary statistics that are commonly reported in primary research studies.

The package can estimate 11 different types of effect sizes (such as Cohen's d, Hedges' g, Odds Ratios, and Correlation Coefficients), making it a versatile tool for synthesizing evidence across studies that may report their findings in different formats.

## Project Structure

The package follows a standard R package layout, with a clear organization of its source code and documentation:

-   **`R/`**: Contains all the R source code. The files are logically structured:
    -   `main_*.R` files house the primary, user-facing functions like `convert_df`.
    -   A large number of `es_from_*.R` files each contain a specific formula to calculate an effect size from a particular set of inputs (e.g., `es_from_means_sd.R` calculates an ES from means and standard deviations).
    -   `internal_*.R` files contain helper functions used by other parts of the package.
-   **`man/`**: Holds the `.Rd` documentation files for every exported function, providing detailed help that is accessible via `?function_name` in R.
-   **`data/`**: Includes several example datasets (`df.haza`, `df.short`, etc.) that are used in the documentation and for testing purposes.
-   **`vignettes/`**: Contains the `Tutorial.Rmd` file, which serves as a long-form, narrative guide to using the package with practical examples.

## Data Structures

The package is designed to work with a single, specific data structure: a **wide-format `data.frame`**.

-   Each row in the dataframe must represent a single study or a single effect size calculation.
-   The columns must use specific, predefined names that the package recognizes (e.g., `n_exp`, `mean_exp`, `student_t`, `or`).

To facilitate data entry, the package provides a key utility function: `data_extraction_sheet()`. When called, this function generates a template `data.frame` that contains all possible column names that `metaConvert` can use. This template serves as the definitive guide to the expected input data format.

## Key Functions & Workflow

The intended workflow for using `metaConvert` is a straightforward, two-step process that efficiently handles complex calculations under the hood.

### Step 1: `convert_df()`

This is the main engine of the package. The user prepares their data in the wide format described above and passes the dataframe to this function.

-   `convert_df` iterates through each row of the input data.
-   For each row, it systematically attempts to apply *every possible `es_from_*` calculation function*.
-   The function returns a complex S3 object of class `metaConvert`. This object is a large list (containing over 70 dataframes), where each dataframe corresponds to the results from one specific calculation method. This list is intelligently ordered based on a default (or user-specified) hierarchy of evidence, which prioritizes more reliable data sources (e.g., calculations from means and SDs are ranked higher than calculations from a p-value).

### Step 2: `summary()`

The user then passes the `metaConvert` object generated in the previous step to the generic `summary()` function.

-   The package provides a specialized method, `summary.metaConvert`, to handle this object.
-   This method processes the hierarchical list of dataframes. For each study (row), it selects the **single best available effect size** based on the evidence hierarchy.
-   The final output is a clean, tidy dataframe containing one row per study with the most reliable effect size calculation, ready for subsequent meta-analysis.

### Ancillary Functions

In addition to this primary workflow, the package also exports:
-   All the individual **`es_from_*`** functions, allowing users to perform specific, one-off calculations if needed.
-   **`compare_df()`**: A utility to check for discrepancies between two data-entry files.
-   **`see_input_data()`**: A helper function that shows which input data columns are used by which calculation functions.
