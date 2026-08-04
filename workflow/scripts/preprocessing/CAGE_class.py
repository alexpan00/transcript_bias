#!/usr/bin/env python
import os
from bx.intervals import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree
from collections import defaultdict

class CAGEPeak:
    """
    A class to represent and query CAGE (Cap Analysis of Gene Expression) peaks from a BED file.

    Attributes
    ----------
    cage_bed_filename : str
        The filename of the BED file containing CAGE peak data.
    cage_peaks : defaultdict
        A dictionary where keys are tuples of (chromosome, strand) and values are IntervalTree objects containing intervals of peaks.

    Methods
    -------
    __init__(cage_bed_filename):
        Initializes the CAGEPeak object with the given BED filename and reads the BED file to populate the peaks.
    read_bed():
        Reads the BED file and populates the cage_peaks attribute with intervals of peaks.
    find(chrom, strand, query, search_window=10000):
        Queries the CAGE peaks to determine if a given position falls within a peak and calculates the distance to the nearest TSS.
    """    
    def __init__(self, cage_bed_filename):
        self._validate_input(cage_bed_filename)
        self.cage_bed_filename = cage_bed_filename
        self.cage_peaks = defaultdict(lambda: IntervalTree()) # (chrom,strand) --> intervals of peaks

        self.read_bed()

    def _validate_input(self, cage_bed_filename):
        if not cage_bed_filename.endswith('.bed'):
            raise ValueError("CAGE peak file must be in BED format.")
        if not os.path.exists(cage_bed_filename):
            raise FileNotFoundError(f"CAGE peak {cage_bed_filename} does not exist.")
        
    def read_bed(self):
        with open(self.cage_bed_filename, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                raw = line.strip().split('\t')
                if len(raw) < 6:
                    continue
                chrom = raw[0]
                start0 = int(raw[1])
                end1 = int(raw[2])
                peak_id = raw[3]
                strand = raw[5]
                tss0 = int((start0 + end1) / 2)
                # Store (tss0, peak_id) tuple as interval payload
                self.cage_peaks[(chrom, strand)].insert(start0, end1, (tss0, peak_id))

    def find(self, chrom, strand, query, search_window=10000):
        """
        :param chrom: Chromosome to query
        :param strand: Strand of the query ('+' or '-')
        :param query: Position to query
        :param search_window: Window around the query position to search for peaks
        :return: <True/False falls within a cage peak>, <nearest distance to TSS>
        If the distance is negative, the query is upstream of the TSS.
        If the query is outside of the peak upstream of it, the distance is NA
        """
        within_peak, dist_peak = 'FALSE', float('inf')
        peaks = self.cage_peaks[(chrom, strand)].find(query - search_window, query + search_window)

        for interval in peaks:
            start0 = interval.start
            end1 = interval.end
            val = interval.value
            if isinstance(val, tuple):
                tss0, peak_id = val
            elif isinstance(val, int):
                tss0, peak_id = val, ""
            else:
                tss0, peak_id = int((start0 + end1) / 2), str(val)

            # Checks if the TSS is upstream of a peak
            if (strand == '+' and start0 > query and end1 > query) or \
               (strand == '-' and start0 < query and end1 < query):
                continue
            # Checks if the query is within the peak and the distance to the TSS
            within_out = start0 <= query < end1 if strand == '+' else start0 < query <= end1
            d = (tss0 - query) * (-1 if strand == '-' else 1)
            w = 'TRUE' if within_out else 'FALSE'
            
            if within_peak != 'TRUE':
                within_peak = w
                if w == 'TRUE' or abs(d) < abs(dist_peak):
                    dist_peak = d
            else:
                if abs(d) < abs(dist_peak) and not (w == 'FALSE' and within_peak == 'TRUE'):
                    within_peak, dist_peak = w, d 
        if dist_peak == float('inf'):
            dist_peak = 'NA'
        return within_peak, dist_peak