#!/bin/bash
#
INSTALL.sh
#
# BUILD_DB.sh
#
GENERATE_ADAPTERS.sh
GET_PHIX.sh
GET_MTDNA.sh
#
gunzip VDB.fa.gz
#
declare -a READS=("healthy-skin:MT18021B_S1_L001_R1_001.fastq.gz:MT18021B_S1_L001_R2_001.fastq.gz:MT18021B_S1_L002_R1_001.fastq.gz:MT18021B_S1_L002_R2_001.fastq.gz" "skin:MT18021C_S2_L001_R1_001.fastq.gz:MT18021C_S2_L001_R2_001.fastq.gz:MT18021C_S2_L002_R1_001.fastq.gz:MT18021C_S2_L002_R2_001.fastq.gz" "bone:MT18021D_S3_L001_R1_001.fastq.gz:MT18021D_S3_L001_R2_001.fastq.gz:MT18021D_S3_L002_R1_001.fastq.gz:MT18021D_S3_L002_R2_001.fastq.gz" "colon:MT18021E_S4_L001_R1_001.fastq.gz:MT18021E_S4_L001_R2_001.fastq.gz:MT18021E_S4_L002_R1_001.fastq.gz:MT18021E_S4_L002_R2_001.fastq.gz" "heart:MT18021F_S5_L001_R1_001.fastq.gz:MT18021F_S5_L001_R2_001.fastq.gz:MT18021F_S5_L002_R1_001.fastq.gz:MT18021F_S5_L002_R2_001.fastq.gz" "liver:MT18021G_S6_L001_R1_001.fastq.gz:MT18021G_S6_L001_R2_001.fastq.gz:MT18021G_S6_L002_R1_001.fastq.gz:MT18021G_S6_L002_R2_001.fastq.gz" "spleen:MT18022B_S7_L001_R1_001.fastq.gz:MT18022B_S7_L001_R2_001.fastq.gz:MT18022B_S7_L002_R1_001.fastq.gz:MT18022B_S7_L002_R2_001.fastq.gz" "kidney:MT18022C_S8_L001_R1_001.fastq.gz:MT18022C_S8_L001_R2_001.fastq.gz:MT18022C_S8_L002_R1_001.fastq.gz:MT18022C_S8_L002_R2_001.fastq.gz" "lung:MT18022D_S9_L001_R1_001.fastq.gz:MT18022D_S9_L001_R2_001.fastq.gz:MT18022D_S9_L002_R1_001.fastq.gz:MT18022D_S9_L002_R2_001.fastq.gz" "plasma:MT18022E_S10_L001_R1_001.fastq.gz:MT18022E_S10_L001_R2_001.fastq.gz:MT18022E_S10_L002_R1_001.fastq.gz:MT18022E_S10_L002_R2_001.fastq.gz" "blood:MT18022F_S11_L001_R1_001.fastq.gz:MT18022F_S11_L001_R2_001.fastq.gz:MT18022F_S11_L002_R1_001.fastq.gz:MT18022F_S11_L002_R2_001.fastq.gz" "bone-marrow:MT18022G_S12_L001_R1_001.fastq.gz:MT18022G_S12_L001_R2_001.fastq.gz:MT18022G_S12_L002_R1_001.fastq.gz:MT18022G_S12_L002_R2_001.fastq.gz" "teeth:MT18023B_S13_L001_R1_001.fastq.gz:MT18023B_S13_L001_R2_001.fastq.gz:MT18023B_S13_L002_R1_001.fastq.gz:MT18023B_S13_L002_R2_001.fastq.gz" "brain:MT18023C_S14_L001_R1_001.fastq.gz:MT18023C_S14_L001_R2_001.fastq.gz:MT18023C_S14_L002_R1_001.fastq.gz:MT18023C_S14_L002_R2_001.fastq.gz")

################################################################################
## LOOP TRHOUGH THE ARRAY READS AND SPLITS THE DATA INTO 5 FIELDS
for read in "${READS[@]}"
  do
  ORGAN_T=`echo $read | tr ':' '\t' | awk '{ print $1 }'`;
  SPL_R1A=`echo $read | tr ':' '\t' | awk '{ print $2 }'`;
  SPL_R2A=`echo $read | tr ':' '\t' | awk '{ print $3 }'`;
  SPL_R1B=`echo $read | tr ':' '\t' | awk '{ print $4 }'`;
  SPL_R2B=`echo $read | tr ':' '\t' | awk '{ print $5 }'`;
  echo "RUN ORGAN=$ORGAN_T F1=$SPL_R1A R1=$SPL_R2A F2=$SPL_R1B R2=$SPL_R2B";

  ## MERGE FILES
  echo "MERGGING FILES ...";
  rm -f FW_READS.fq.gz RV_READS.fq.gz
  zcat $SPL_R1A $SPL_R1B | gzip > FW_READS.fq.gz
  zcat $SPL_R2A $SPL_R2B | gzip > RV_READS.fq.gz

  #
  ./TRIM_FILTER_READS.sh
  #
  ./REMOVE_PHIX.sh
  #
  ./METAGENOMICS.sh $ORGAN_T
  #
  # profiles...
  ./EXTRACT_MTDNA.sh
  #
  ./ASSEMBLE_MT.sh $ORGAN_T
  #
  done
#

################################################################################
