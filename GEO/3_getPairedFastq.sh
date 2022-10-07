#!/bin/sh

ls *_R[12]*.fastq.gz| sort | tr "\n" "\t" | sed 's/\t/\n/2; P; D'
