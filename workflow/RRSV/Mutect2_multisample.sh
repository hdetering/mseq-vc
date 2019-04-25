#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user lauratomaslopezslurm@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 00:05:00
#SBATCH --mem 50G
#SBATCH -p thinnodes,thin-shared,cola-corta




# Reading configuration file

source ReadConfig.sh $1


# Loading modules
module load jdk/8u181
GATK=/mnt/netapp1/posadalab/APPS/gatk-4.1.1.0/gatk


# Selecting inputs and outputs

HEALTHY=ChineseSon.H
LIST=ChineseFinalReplicates.SampleList

REPLICATE=R$((((${SLURM_ARRAY_TASK_ID} - 1)%10)+1))
SCENARIO=S$((((${SLURM_ARRAY_TASK_ID} - 1)/10)+1))
echo $REPLICATE $SCENARIO

echo "I am running in the node"
hostname
echo "Healthy is $HEALTHY"
echo "Tumor is ${SCENARIO}.${REPLICATE}"
echo "Running Mutect2 for $HEALTHY and ${SCENARIO}.${REPLICATE}"


# Command

GERMRES=${RESDIR}/af-only-gnomad.raw.sites.b37.vcf.gz
PON=/mnt/netapp1/posadalab/Tama_Laura/PON/workdir/PON.TruSeq.vcf


time(

$GATK Mutect2\
  -R ${RESDIR}/${REF}.fa \
  -I ${WORKDIR}/ChineseSon.T1.${SCENARIO}.${REPLICATE}.bam \
  -I ${WORKDIR}/ChineseSon.T2.${SCENARIO}.${REPLICATE}.bam \
  -I ${WORKDIR}/ChineseSon.T3.${SCENARIO}.${REPLICATE}.bam \
  -I ${WORKDIR}/ChineseSon.T4.${SCENARIO}.${REPLICATE}.bam \
  -I ${WORKDIR}/ChineseSon.T5.${SCENARIO}.${REPLICATE}.bam \
  -I ${WORKDIR}/${HEALTHY}.chr21.bam \
  -normal ${HEALTHY} \
  --germline-resource ${GERMRES} \
  --af-of-alleles-not-in-resource 0.00003125 \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  -pon ${PON} \
  -L 21 \
  --output ${WORKDIR}/${SCENARIO}.${REPLICATE}.Mutect2_multisample.vcf

)


# Filtering
$GATK FilterMutectCalls \
        -V ${WORKDIR}/${SCENARIO}.${REPLICATE}.Mutect2_multisample.vcf \
        -O ${WORKDIR}/${SCENARIO}.${REPLICATE}.Mutect2_multisample.PASS.vcf \
        -R ${RESDIR}/${REF}.fa \
	--stats ${WORKDIR}/${SCENARIO}.${REPLICATE}.Mutect2_multisample.vcf.stats
