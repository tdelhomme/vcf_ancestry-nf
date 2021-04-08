#!/usr/bin/env nextflow

// Copyright (C) 2021 IRB Barcelona

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null
params.input_vcf = null
params.output_vcf = "FALSE"
params.KG_folder = null
params.ancestry = "EUR"
params.output_folder = "vcf_ancestry_output"

log.info ""
log.info "-----------------------------------------------------------------------"
log.info "VCF ancestry - nextflow pipeline"
log.info "-----------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run main.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input_vcf                FILE                 Input VCF file to be filtered on ancestry"
    log.info "--KG_folder               FOLDER               1000 genomes folder for reference"
    log.info "--fasta_ref                FILE                 Reference genome in fasta format"
    log.info ""
    log.info "Optional arguments:"
    log.info "--output_vcf               FILE                 Output VCF (default=renamed with ancestry)"
    log.info "--output_folder            FOLDER               Output folder (default=vcf_ancestry_output)"
    log.info "--ancestry                 STRING               Ancestry (from 1000 genomes) for filtering samples"
    log.info ""
    log.info "Flags:"
    log.info "--help                                          Display this message"
    log.info ""
    exit 1
}

assert (params.input_vcf != null) : "please provide the --input_vcf option"
assert (params.fasta_ref != null) : "please provide the --fasta_ref option"
assert (params.KG_folder != null) : "please provide the --KG_folder option"

vcf = file(params.input_vcf)
fasta_ref = file(params.fasta_ref)

process filter_VCF {

  publishDir params.output_folder, mode: 'copy', pattern: "*_arraysnps_merged.vcf.gz"

  input:
  file vcf
  file fasta_ref

  output:
  file "tmp.vcf.gz" into filt_vcf

  shell:
  '''
  # DR2 filtering
  bcftools norm -m - -Oz -f !{fasta_ref} !{vcf} | bcftools filter -i 'INFO/DR2>0.3' | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' > tmp0.vcf.gz

  # filtering on HWE and MAF
  plink --vcf tmp0.vcf.gz --maf 0.01 --hwe 1e-6 --make-bed --out filtered_vcf --recode vcf

  # LD pruning
  plink --vcf filtered_vcf.vcf --indep-pairwise 50 5 0.5 --out filtered_vcf_LD_prun
  plink --vcf filtered_vcf.vcf --extract filtered_vcf_LD_prun.prune.in --recode vcf --out filtered_vcf_prun
  bgzip -c filtered_vcf_prun.vcf > tmp.vcf.gz
  '''
}

process extract_common_SNPs {

  publishDir params.output_folder, mode: 'copy', pattern: "*_arraysnps_merged.vcf.gz"

  input:
  file vcf from filt_vcf

  output:
  file "common_snps" into common_snps

  shell:
  '''
  # get input VCF snp list
  zcat !{vcf} | grep -v "^#" | cut -f3 > vcf.snps

  # get 1KG snp list
  touch 1KG.snps
  for chr in {1..22}; do bcftools view !{params.KG_folder}/chr${chr}.1kg.phase3.v5a.vcf | grep -v "^#" | cut -f3 >> 1KG.snps ; done

  # Make intersection file
  grep -Fxf "vcf.snps" "1KG.snps" > intersection.snps

  # Extract the common scripts
  mkdir -p common_snps
  for chr in {1..22}; do plink --bfile !{params.KG_folder}/chr${chr}.1kg.phase3.v5a --extract intersection.snps --make-bed --out common_snps/1KG_intersection_chr${chr}; done
  plink --vcf !{vcf} --const-fid 0 --extract intersection.snps --make-bed --out common_snps/input_VCF_intersection
  '''
}

process PCA {

  publishDir params.output_folder, mode: 'copy', pattern: "*_arraysnps_merged.vcf.gz"

  input:
  file common_snps from common_snps

  output:

  shell:
  '''
  ls
  '''
}
