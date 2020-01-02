# lib-indexer
EDBV project for indexing books with unique labels in library shelves.

## Prerequisites

The program was tested on Matlab 9.3 R2017b and should work with all newer versions. 

Additionally a dataset or image of choice must be provided. The following criteria must be met for the program to work properly:

+ The image must display a section of a bookshelf
+ The bookshelf must contain books, which are labeled in the same matter as the TU Wien library books. 
+ The shelf must not be tilted more than +-30° in the image
+ The books must not be tilted more than +-5° (vertically) in the image
+ The resolution of the image must be high 

## Usage

To use the program, change the image path in "main.m" to match the image you want to be processed and run the program. 

## Output 

The program will output a .json file containing the results of the optical character recognition applied on each label.
The labels are indexed with **1 to n**, whereas n is the number of labels/books displayed in the image.

## To-Do

1. Complete OCR 
2. Change the program structure to release
3. Documentation
4. Prepare and provide a dataset
5. Test program on dataset and validate results and output

## Authors

Anand Eichner, Laurenz Fiala, Aleksandar Vucenovic
