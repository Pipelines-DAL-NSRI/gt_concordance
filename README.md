# gt_concordance

## Description
This package calculates the concordance of data from two separate files. It accepts two input files, either xlsx or csv, and identifies overlapping data (based on the first column). It calculates the number of matching rows and creates a csv file and a barplot.

## Function
```
concordance(file1, file2, haplotype = FALSE)
```

## Parameters
@param *file1* and *file2* Should be excel or csv files. The first column of the input file should contain the samples while the markers are in the succeeding columns. See samplefile1.csv or samplefile2.csv as an example. 
@param *haplotypes* is a required parameter asking if the data are haplotypes. If set to FALSE (default), the function will homogenize the allele orders.

## Usage
concordance("samplefile1.xslx", "samplefile2.xlsx")
