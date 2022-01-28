bsub -q big -J test -o test.log -sla miket_sc "gsutil cp gs://talkowski-sv-gnomad-cromwell-execution/CleanVcfChromosomeLastTwoSteps/645b0267-5cef-4127-ba20-f2adf8aaa616/call-StitchFragmentedCnvs/gnomad-sv-v3.chr21.stitch_fragmented_cnvs.vcf.gz  ./stitch_fix/"
bsub -q big -J test -o test.log -sla miket_sc "gsutil cp gs://talkowski-sv-gnomad-cromwell-execution/CleanVcfChromosomeLastTwoSteps/3a9edace-6db7-4e91-b7f3-e0a32c0dbe1a/call-StitchFragmentedCnvs/gnomad-sv-v3.chr22.stitch_fragmented_cnvs.vcf.gz  ./stitch_fix/"
bsub -q big -J test -o test.log -sla miket_sc "gsutil cp gs://talkowski-sv-gnomad-cromwell-execution/CleanVcfChromosomeLastTwoSteps/58725c2c-92b4-418f-bf56-2eb20206a6e6/call-StitchFragmentedCnvs/gnomad-sv-v3.chrY.stitch_fragmented_cnvs.vcf.gz  ./stitch_fix/"

