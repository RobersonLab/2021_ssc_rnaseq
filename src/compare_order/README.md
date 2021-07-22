compare_order
=============

## Purpose
A simple C++ program for comparing the read identifiers of two sequencing files (Illumina formatted single-line or FASTQ format). It reads the files line by line together and compares their read identifiers. The goal is to confirm that reads are ordered the same in each file. The identifier must be separated from index/pair by a '#', a single space character, or a trailing '/', e.g. HWIFAKE:1:1:1:1#0/1 or HWIFAKE:1:1:1:1 0/1 or HWIFAKE:1:1:1:1/1.

## Compiling

### Linux
    make
	make install

### Windows
    g++ -o compare_order -O3 compare_order.cpp
	
### Running
Requires specifying three parameters:
--format scarf|fastq
--file1 filename
--file2 filename

```bash
compare_order --format scarf --file1 mate_1.txt --file2 mate_2.txt
compare_order --format fastq --file1 mate_1.fq --file2 mate_2.fq
```

### What about compressed FASTQ files?
Most FASTQ files in long-term storage will be compressed, but this simple program doesn't automagically handle FASTQ files. How can we get around that? Simple on Linux-like systems that support mkfifo. Assume your files are called read1.fastq.gz and read2.fastq.gz. Using mkfifo with gunzip -cd will use an intermediate file-in-file-out to stream the decompressed file to the program, but leave the compressed file intact. The '-cd' option for gunzip, however, is **critical** because otherwise the original file will be consumed in the process and disappear.

```bash
mkfifo r1.fastq r2.fastq
nohup gunzip -cd read1.fastq.gz > r1.fastq &
nohup gunzip -cd read2.fastq.gz > r2.fastq &
compare_order --format fastq --file1 r1.fastq --file2 r2.fastq > order.log
rm r1.fastq r2.fastq
```

