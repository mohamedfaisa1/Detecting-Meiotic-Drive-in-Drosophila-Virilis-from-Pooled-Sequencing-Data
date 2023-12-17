import argparse
import pysam



def calculate_af(cargs):
    # Chromosomes to calculate allele frequencies in
    contigs = ["Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_6", "Chr_X"]

    # Initialize VCF file information
    vcf_gz = cargs.vcf
    vcf_index = cargs.vcf_index
    vcf = pysam.VariantFile(filename=vcf_gz, index_filename=vcf_index, mode='rb')
    # Samples & MetaData
    v48 = "vir48"
    v47 = "vir47"
    GT = "GT"

    # Initialize Sam file information
    bam_file = cargs.bam
    bai_file = cargs.bai
    sam_file = pysam.AlignmentFile(filename=bam_file, index_filename=bai_file, mode='rb')

    # Initialize output file
    output_file = cargs.out_file
    # For each chromosome, calculate the allele frequencies
    for contig in contigs:
        with open(contig+"_"+output_file, "w") as out_file:
            out_file.write(f"chromosome,position,{v48}_snp,{v48}_snp_count,{v47}_snp,{v47}_snp_count\n")
            vcf_contig = vcf.fetch(contig)
            for vc in vcf_contig:
                alleles = vc.alleles
                ref_allele = vc.ref
                alt_alleles = vc.alts
                v48_gt = vc.samples[v48][GT]
                v47_gt = vc.samples[v47][GT]

                # get only SNPs at bi-allelic sites

                v48_snp = ''
                v47_snp = ''
                if (v48_gt != v47_gt) and (v48_gt[0] == v48_gt[1]) and (v47_gt[0] == v47_gt[1]) \
                        and (len(alleles) == 2) and (len(ref_allele) == 1) and (len(alt_alleles) == 1) \
                        and (len(ref_allele[0]) == 1) and (len(alt_alleles[0]) == 1) \
                        and v48_gt[0] != '.' and v47_gt[0] != '.':

                    alt_allele = alt_alleles[0]
                    if v48_gt[0] == 0:
                        v48_snp = ref_allele
                        v47_snp = alt_allele
                    elif v48_gt[0] == 1:
                        v48_snp = alt_allele
                        v47_snp = ref_allele
                    else:
                        continue

                    # to get position of bi-allelic site & initialize counts
                    snp_position = vc.pos
                    v47_count = 0
                    v48_count = 0

                    # compute parental allele counts at the bi-allelic site
                    pileups = sam_file.pileup(contig=contig, start=snp_position-1, stop=snp_position, truncate=True)
                    for pileup_col in pileups:
                        for pileup_read in pileup_col.pileups:
                            if not pileup_read.is_del and not pileup_read.is_refskip:
                                snp_index = pileup_read.query_position
                                SNP = pileup_read.alignment.query_sequence[snp_index]
                                if SNP == v47_snp:
                                    v47_count += 1
                                elif SNP == v48_snp:
                                    v48_count += 1
                                else:
                                    continue
                    # do not write if there are no counts at the bi-allelic site
                    if v48_count == 0 and v47_count == 0:
                        continue
                    else:
                        out_file.write(f"{contig},{snp_position},{v48_snp},{v48_count},{v47_snp},{v47_count}\n")

    return


def main():
    parser = argparse.ArgumentParser(
        description="Compute allele frequencies at bi-allelic sites in the Drosophila virilis genome of strains Vir47 and Vir48.")
    parser.add_argument("-vcf", action="store", dest="vcf", type=str,
                        default="Dvir.GATK_calls_8-24.vir48_vir47.vcf.gz")
    parser.add_argument("-vcf_index_file", action="store", dest="vcf_index", type=str,
                        default="Dvir.GATK_calls_8-24.vir48_vir47.vcf.gz.csi")
    parser.add_argument("-bam", action="store", dest="bam", type=str,
                        default="12925_11745_151591_HTCKHBGXK_V1A_CAGTGAAA_AGACCAGT_sorted.bam")
    parser.add_argument("-bai", action="store", dest="bai", type=str,
                        default="12925_11745_151591_HTCKHBGXK_V1A_CAGTGAAA_AGACCAGT_sorted.bam.bai")
    parser.add_argument("-o", action="store", dest="out_file", type=str,
                        default="allele_frequencies.csv")

    args = parser.parse_args()
    calculate_af(args)
    return


if __name__ == '__main__':
    main()