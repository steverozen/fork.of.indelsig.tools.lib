vcf = ICAMS::ReadVCFs(
  "inst/Strelka.ID.GRCh37.s1.vcf",
  variant.caller = "strelka"
)[[1]]

colnames(vcf)[1] = "chr"
colnames(vcf)[2] = "position"
vcf$chr = paste0("chr", vcf$chr)

devtools::load_all()
debug(indel_classifier89)

indel_classifier89(vcf, "hg19")
