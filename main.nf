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
params.PED = null
params.ancestry = "EUR"
params.output_folder = "vcf_ancestry_output"
params.keep_chr = "FALSE"

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
    log.info "--KG_folder                FOLDER               1000 genomes folder for reference"
    log.info "--fasta_ref                FILE                 Reference genome in fasta format"
    log.info "--PED                      FILE                 Reference PED file for ancestry estimation"
    log.info ""
    log.info "Optional arguments:"
    log.info "--output_vcf               FILE                 Output VCF (default=renamed with ancestry)"
    log.info "--output_folder            FOLDER               Output folder (default=vcf_ancestry_output)"
    log.info "--ancestry                 STRING               Ancestry (from 1000 genomes) for filtering samples"
    log.info ""
    log.info "Flags:"
    log.info "--keep_chr                                           Chr with prefix"
    log.info "--help                                          Display this message"
    log.info ""
    exit 1
}

assert (params.input_vcf != null) : "please provide the --input_vcf option"
assert (params.fasta_ref != null) : "please provide the --fasta_ref option"
assert (params.KG_folder != null) : "please provide the --KG_folder option"
assert (params.PED != null) : "please provide the --PED option"

vcf = file(params.input_vcf)
fasta_ref = file(params.fasta_ref)
ped = file(params.PED)

process filter_VCF {

  input:
  file vcf
  file fasta_ref

  output:
  file "tmp.vcf.gz" into filt_vcf

  shell:
  '''
  # DR2 filtering
  # bcftools norm -m - -Oz -f !{fasta_ref} !{vcf} | bcftools filter -i 'INFO/DR2>0.3' | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' > tmp0.vcf.gz
  tabix -p vcf !{vcf}
  bcftools norm -m - -Oz -f !{fasta_ref} !{vcf} | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' > tmp0.vcf.gz

  # filtering on HWE and MAF
  plink --vcf tmp0.vcf.gz --maf 0.1 --hwe 1e-6 --make-bed --out filtered_vcf --recode vcf

  # LD pruning
  plink --vcf filtered_vcf.vcf --indep-pairwise 50 5 0.5 --out filtered_vcf_LD_prun
  plink --vcf filtered_vcf.vcf --extract filtered_vcf_LD_prun.prune.in --recode vcf --out filtered_vcf_prun
  bgzip -c filtered_vcf_prun.vcf > tmp.vcf.gz
  '''
}

process extract_common_SNPs {

  input:
  file vcf from filt_vcf

  output:
  file "common_snps" into common_snps

  shell:
  if(params.keep_chr!="FALSE"){ keep_chr="yes" } else {keep_chr="no"} 
  '''
  # get input VCF snp list
  zcat !{vcf} | grep -v "^#" | cut -f3 > vcf.snps
  sed -i 's/chr//' vcf.snps

  # get 1KG snp list
  touch 1KG.snps
  for chr in {1..22}; do bcftools view !{params.KG_folder}/chr${chr}.1kg.phase3.v5a.bcf | grep -v "^#" | cut -f3 >> 1KG.snps ; done
  sed -i 's/chr//' 1KG.snps

  # Make intersection file
  grep -Fxf "vcf.snps" "1KG.snps" > intersection.snps
  cp intersection.snps intersection.snps1
  if [ "!{keep_chr}" = "yes" ]; then cat intersection.snps | awk '{print "chr" $0}' > intersection.snps1 ; fi

  # Extract the common scripts
  mkdir -p common_snps
  for chr in {1..22}; do plink --bfile !{params.KG_folder}/chr${chr}.1kg.phase3.v5a --extract intersection.snps --make-bed --out common_snps/1KG_intersection_chr${chr}; done
  plink --vcf !{vcf} --const-fid 0 --extract intersection.snps1 --make-bed --out common_snps/input_VCF_intersection
  '''
}

process PCA {

  input:
  file common_snps from common_snps

  output:
  file "plink_eigenvec" into eigenvec

  shell:
  '''
  mkdir -p merge
  echo "common_snps/1KG_intersection_chr1" > merge/list
  for chr in {2..22}; do echo "common_snps/1KG_intersection_chr${chr}" >> merge/list; done
  echo "common_snps/input_VCF_intersection" >> merge/list

  plink --memory 4000 --merge-list merge/list --out merge/1KG_with_input_VCF

  mkdir -p pca && cd pca
  plink --bfile ../merge/1KG_with_input_VCF --pca
  cd .. && cp pca/plink.eigenvec plink_eigenvec
  '''
}


process PCA_analysis {

  publishDir params.output_folder, mode: 'copy', pattern: "*pdf"
  publishDir params.output_folder, mode: 'copy', pattern: "table_3PCs.txt"

  input:
  file eigenvec from eigenvec
  file ped

  output:
  file "*pdf" into plots
  file "table_3PCs.txt" into res

  shell:
  '''
  Rscript !{baseDir}/bin/analysis_PCA.R --eigenvec_file=!{eigenvec} --PED=!{ped}
  '''
}
