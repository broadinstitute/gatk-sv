CONTAINER_ID="0de4667c8078"

echo "localizing scripts"
podman cp src/bash_workflows/collect_counts.sh "${CONTAINER_ID}":/
podman cp src/bash_workflows/gather_sample_evidence.sh "${CONTAINER_ID}":/
podman cp src/bash_workflows/run_manta.sh "${CONTAINER_ID}":/
podman cp src/bash_workflows/collect_sv_evidence.sh "${CONTAINER_ID}":/
podman cp src/bash_workflows/scramble.sh "${CONTAINER_ID}":/
podman cp src/bash_workflows/run_whamg.sh "${CONTAINER_ID}":/
podman cp src/bash_workflows/standardize_vcf.sh "${CONTAINER_ID}":/
podman cp src/bash_workflows/realign_soft_clipped_reads.sh "${CONTAINER_ID}":/


echo "localizing cram files"
#podman cp NA12878.final.cram "${CONTAINER_ID}":/
#podman cp NA12878.final.cram.crai "${CONTAINER_ID}":/
podman cp downsampled_HG00096.final.cram "${CONTAINER_ID}":/
podman cp downsampled_HG00096.final.cram.crai "${CONTAINER_ID}":/

echo "localizing reference assembly files"
podman cp Homo_sapiens_assembly38.fasta "${CONTAINER_ID}":/
podman cp Homo_sapiens_assembly38.fasta.fai "${CONTAINER_ID}":/
podman cp Homo_sapiens_assembly38.dict "${CONTAINER_ID}":/

echo "localizing misc files"
#podman cp primary_contigs.list "${CONTAINER_ID}":/
#podman cp contig.fai "${CONTAINER_ID}":/
podman cp downsampled_primary_contigs.list "${CONTAINER_ID}":/
podman cp downsampled_contig.fai "${CONTAINER_ID}":/

#podman cp preprocessed_intervals.interval_list "${CONTAINER_ID}":/
podman cp downsampled_preprocessed_intervals.interval_list "${CONTAINER_ID}":/

#podman cp primary_contigs_plus_mito.bed.gz "${CONTAINER_ID}":/
#podman cp primary_contigs_plus_mito.bed.gz.tbi "${CONTAINER_ID}":/
podman cp downsampled_primary_contigs_plus_mito.bed.gz "${CONTAINER_ID}":/
podman cp downsampled_primary_contigs_plus_mito.bed.gz.tbi "${CONTAINER_ID}":/

#podman cp Homo_sapiens_assembly38.dbsnp138.vcf "${CONTAINER_ID}":/
podman cp downsampled_Homo_sapiens_assembly38.dbsnp138.vcf "${CONTAINER_ID}":/

podman cp hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz "${CONTAINER_ID}":/

#podman cp wham_whitelist.bed "${CONTAINER_ID}":/
podman cp downsampled_wham_whitelist.bed "${CONTAINER_ID}":/

podman cp Homo_sapiens_assembly38.fasta.64.alt "${CONTAINER_ID}":/
podman cp Homo_sapiens_assembly38.fasta.64.amb "${CONTAINER_ID}":/
podman cp Homo_sapiens_assembly38.fasta.64.ann "${CONTAINER_ID}":/
podman cp Homo_sapiens_assembly38.fasta.64.bwt "${CONTAINER_ID}":/
podman cp Homo_sapiens_assembly38.fasta.64.pac "${CONTAINER_ID}":/
podman cp Homo_sapiens_assembly38.fasta.64.sa "${CONTAINER_ID}":/

