#!/bin/bash

for filename in *Tn5.bed_peaks.narrowPeak *Tn5.bed_peaks.xls *Tn5.bed_summits.bed; do
    if [[ "$filename" == *Tn5.bed_peaks.narrowPeak ]]; then
        new_filename="${filename%%.Tn5.bed_peaks.narrowPeak}_peaks.narrowPeak"
    elif [[ "$filename" == *Tn5.bed_peaks.xls ]]; then
        new_filename="${filename%%.Tn5.bed_peaks.xls}_peaks.xls"
    elif [[ "$filename" == *Tn5.bed_summits.bed ]]; then
        new_filename="${filename%%.Tn5.bed_summits.bed}_summits.bed"
    else
        new_filename="$filename"  
    fi
    
    mv "$filename" "$new_filename"
done
