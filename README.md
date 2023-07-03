# CNV Analysis From NGS Depth of Coverage

This project contains a Perl script configured for Copy Number Variation (CNV) analysis from Next-Generation Sequencing (NGS) depth of coverage data on tumor-normal pairs.

## Algorithm
The CNV analysis Perl script performs the following steps:

1. Input Parameters: The script expects command line parameters to be provided in pairs. Each pair consists of the normal sample and the tumor sample depth of coverage interval summary files, respectively.

2. Loading Annotations: The script loads annotations from the provided EXOME.interval_list.annotated file. This file contains information about genomic regions and associated genes.

3. Normalization: For each pair of samples, the script calculates the average coverage for both the normal and tumor samples. It then determines a normalizer value based on the ratio of normal average coverage to tumor average coverage.

4. Calculating Ratios: The script calculates the coverage ratio between the tumor and normal samples for each interval. If the normal sample's coverage for a specific interval is zero, the ratio is set to 1.

5. Statistical Analysis: The script calculates the average and standard deviation of the coverage ratios across all intervals. It uses these statistics to calculate the Z-score for each interval, indicating the deviation from the average.

6. Output Files: The script generates two output files for each pair of samples. The first file, ending in '.totals' contains detailed information about each interval, including sample names, genomic location, associated gene, normal depth, tumor-normalized depth, ratio, and Z-score. The second file, ending in '.cnv', contains filtered information based on the specified specificity threshold, marking regions as amplifications or deletions.

## Usage
1. Clone the repository: $ git clone https://github.com/zpollock/cnv_from_depth_of_coverage.git
2. Ensure Perl is installed on your system.
3. Download an Exome Interval List from GATK.  Place the file in the root directory and name it EXOME.interval_list.annotated. 
4. Provide the depth of coverage interval summary files as command line parameters in pairs, where the first entry is the normal sample file and the second entry is the tumor sample file.
```bash
$ perl z_cnv.pl patient_1_normal.interval_summary patient_1_tumor.interval_summary
```

5. The script will generate two output files for each pair, detailed above.

## Requirements
- Perl (version 5.X.X)
- Depth of Coverage Interval Summary Files (GATK format)
- EXOME.interval_list.annotated file (annotations for genomic regions, can be downloaded from gatk.broadinstitute.org)

## License
This project is licensed under the MIT License.

## Contributing
Contributions are welcome! Please feel free to contribute, report issues, or provide suggestions.
