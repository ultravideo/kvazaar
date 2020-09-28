Examples
========
Examples of external files for use with Kvazaar.

## Region of interest (roi) files
A simple text file can be used with the `--roi` switch to setup regions of interest for encoding.
Header row of the file will tell how many regions the encoded frames are divided (columns, rows).
The header must be followed by a data row with number entries equal to columns * rows.
The data row will tell the encoder which delta QP value will be assigned to each region.
The included example file will split frames into four regions with the top regions having a delta QP of +5
```
2 2
5 5 0 0
``` 