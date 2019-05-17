#! /usr/bin/env python3
import numpy as np
import vcf
import vcf.utils
import pysam

__author__ = 'Alexander Benkö'


class Assignment2:
    
    def __init__(self, file_name):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        self.file_name = file_name
        self.vcf_file = vcf.Reader(open(self.file_name, 'r'))
        self.file_metadata = self.vcf_file.metadata

        self.avg_qual = []
        self.total_var = 0
        self.chrom = None
        self.seq_meth = None
        self.indel_cnt = 0
        self.snv_cnt = 0
        self.het_cnt = 0


        for record in self.vcf_file:
            if self.chrom == None:  ## Chromosome Nr.
                self.chrom = record.CHROM
            self.avg_qual.append(record.QUAL) ## Avg. Quality
            self.total_var += 1 ## Total Number of Variants
            if record.is_indel: ## Number of Indels
                self.indel_cnt += 1
            if record.is_snp: ## Number of SNVs
                self.snv_cnt += 1
            self.het_cnt += record.num_het ## Number of Heterozygous Variants
            self.var_caller = record.INFO['callsetnames'] ## Variant Caller

        self.avg_qual = np.mean(self.avg_qual)

        with open (self.file_name) as fh: ## Sequencing Method
            for line in fh:
                if "Illumina" in line and self.seq_meth == None:
                    self.seq_meth = "Illumina"


    def merge_chrs_into_one_vcf(self, sec_file, file_name_combo):
        '''
        Creates one VCF containing all variants of chr21 and chr22
        :return:
        '''
        vcf_22 = vcf.Reader(open(self.file_name), "r")
        vcf_21 = vcf.Reader(open(sec_file), "r")
        vcf_files = [vcf_21,vcf_22]
        vcf_Writer = vcf.Writer(open(f"{file_name_combo}.vcf", "w"), vcf_21)

        for i in vcf_files:
            for record in i:
                vcf_Writer.write_record(record)

        vcf_temp = vcf.Reader(open(f"{file_name_combo}.vcf", "r"))
        self.comb_var_total = 0
        for i in vcf_temp:
            self.comb_var_total += 1

        print(f"\nCombined Total Number of Variants: {self.comb_var_total}")



    def print_summary(self):
        print(f"\nFileformat: {self.file_metadata['fileformat']}")
        print(f"\nReference Sequence: {self.chrom}")
        print(f"\nSequencing Method: {self.seq_meth}")
        print(f"\nVariant Caller: {self.var_caller[1]}")
        print(f"\nAverage PHRED Quality: {self.avg_qual}")
        print(f"\nTotal Number of Variants: { self.total_var}")
        print(f"\nIndel Count: {self.indel_cnt}")
        print(f"\nSNV Count: {self.snv_cnt}")
        print(f"\nHeterozygous Variant Count: {self.het_cnt}")

    
def main():
    print("Assignment 2:\n")
    print(f"\nAuthor: {__author__}")
    assignment2 = Assignment2("chr22_new.vcf")
    assignment2.print_summary()
    assignment2.merge_chrs_into_one_vcf("chr21_new.vcf", "chr_21_22")
    print("\nDone with assignment 2")
        
        
if __name__ == '__main__':
    main()
   
    



