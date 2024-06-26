# dea_edgeR

## Overview

This repository contains code and data for differential expression analysis using the edgeR package. The project aims to identify differentially expressed genes in Psoriatic Arthritis (PsA) samples.

## Table of Contents

- [Overview](#overview)
- [Table of Contents](#table-of-contents)
- [Installation](#installation)
- [Usage](#usage)


## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/montahac/dea_edgeR.git
    cd dea_edgeR
    ```

2. Ensure you have R and the necessary packages installed:
    ```r
    install.packages("edgeR")
    install.packages("readr")
    install.packages("dplyr")
    ```

## Usage

1. **Reading raw expression data**:
    ```r
    expr0 <- read.table(file = 'GSE205748_read_counts_PsA.csv', 
                        header = TRUE,
                        row.names = 1)
    ```

2. **Reading metadata (sample information)**:
    ```r
    meta1 <- read.table('GSE205748_series_matrix_edit.txt',
                        header = TRUE)
    ```

3. **Selecting rows and columns**:
    ```r
    # Select the first row
    meta1[1,]
    # Select the first column
    meta1[,1]
    ```

4. **Data Preparation**:
   - checking if column order of expr table is the same as row order of meta1
   ```r
    identical(colnames(expr0), meta1$Sample_code)
    ```

   - The purpose of this code is to ensure that the column names of the expression data (expr0) match exactly with the sample codes listed in the metadata (meta1). This is important in bioinformatics and other data analysis contexts because it ensures that the expression data is correctly aligned with the sample metadata. Misalignment could lead to incorrect analysis results.
  
   - Model matrix and expression check
  

   


