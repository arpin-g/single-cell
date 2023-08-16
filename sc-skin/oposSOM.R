# library(devtools)
# install_github("hloefflerwirth/oposSOM")

library(oposSOM)

sm <- fic_metacells

# seurat <- sm[[1]]
metacell.data <- sm[[2]]

env <- opossom.new(list(dataset.name = "2023_08_Reynolds_fic_som_1",
                        database.host = "https://apr2022.archive.ensembl.org/",
                        dim.1stLvlSom = 60,
                        standard.spot.modules = "overexpression",
                        spot.coresize.modules = 2,
                        spot.threshold.modules = 0.92 ) )

# Load input data
env$indata <- metacell.data

extract_part <- function(value) {
  parts <- unlist(strsplit(value, "_adult|_fetal"))
  return(parts[1])
}

# Define sample groups
env$group.labels <- sapply(colnames(env$indata), extract_part)
o <- order( env$group.labels)
env$group.labels <- env$group.labels[o]
env$indata <- env$indata[,o]

# execute
env <- opossom.run(env)

interest <- c("cell development",
              "cell differentiation",
              "cell migration",
              "cellular response to fibroblast growth factor stimulus",
              "collagen binding",
              "collagen catabolic process",
              "collagen fibril organization",
              "collagen trimer",
              "collagen-containing extracellular matrix",
              "developmental growth",
              "developmental process",
              "extracellular matrix",
              "extracellular matrix binding",
              "extracellular matrix disassembly",
              "extracellular matrix organization",
              "extracellular matrix structural constituent",
              "extracellular region",
              "fibroblast growth factor receptor signaling pathway",
              "immune response",
              "positive regulation of angiogenesis",
              "positive regulation of cell cycle",
              "positive regulation of cell differentiation",
              "positive regulation of cell migration",
              "positive regulation of cell population proliferation",
              "positive regulation of DNA replication",
              "positive regulation of interleukin-6 production",
              "positive regulation of Wnt signaling pathway",
              "regulation of cell cycle",
              "regulation of cell differentiation",
              "regulation of cell migration",
              "Wnt signaling pathway",
              "wound healing",
              "PROTEINATLAS_skin",
              "ESTIMATE_stromal_signature",
              "Bagaev_Cancer_associated_fibroblasts",
              "Bagaev_Matrix",
              "Bagaev_Matrix_remodeling",
              "1_TssP_Fibroblasts",
              "10_ReprPC_Fibroblasts",
              "11_K9K27me3_Fibroblasts",
              "12_Het_Fibroblasts",
              "13_HetRpts_Fibroblasts",
              "14_ZNF_Fibroblasts",
              "15_Quies_Fibroblasts",
              "2_TssA_Fibroblasts",
              "3_TssF_Fibroblasts",
              "4_TxTrans_Fibroblasts",
              "5_Tx_Fibroblasts",
              "6_EnhG_Fibroblasts",
              "7_Enh_Fibroblasts",
              "8_EnhP_Fibroblasts",
              "9_ReprPCWk_Fibroblasts",
              "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
              "HALLMARK_TGF_BETA_SIGNALING",
              "JAZAG_TGFB1_SIGNALING_UP",
              "JAZAG_TGFB1_SIGNALING_DN",
              "JAZAG_TGFB1_SIGNALING_VIA_SMAD4_UP",
              "JAZAG_TGFB1_SIGNALING_VIA_SMAD4_DN",
              "WILLERT_WNT_SIGNALING",
              "PLASARI_TGFB1_SIGNALING_VIA_NFIC_1HR_UP",
              "PLASARI_TGFB1_SIGNALING_VIA_NFIC_1HR_DN",
              "PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_UP",
              "PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_DN",
              "MEBARKI_HCC_PROGENITOR_WNT_UP_CTNNB1_DEPENDENT",
              "MEBARKI_HCC_PROGENITOR_WNT_DN_CTNNB1_DEPENDENT",
              "MEBARKI_HCC_PROGENITOR_WNT_UP_CTNNB1_INDEPENDENT",
              "MEBARKI_HCC_PROGENITOR_WNT_DN_CTNNB1_INDEPENDENT",
              "MEBARKI_HCC_PROGENITOR_WNT_UP_BLOCKED_BY_FZD8CRD",
              "MEBARKI_HCC_PROGENITOR_WNT_DN_BLOCKED_BY_FZD8CRD",
              "NABA_COLLAGENS",
              "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
              "REACTOME_COLLAGEN_FORMATION",
              "REACTOME_DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING",
              "BIOCARTA_TGFB_PATHWAY",
              "KEGG_WNT_SIGNALING_PATHWAY",
              "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
              "KEGG_TGF_BETA_SIGNALING_PATHWAY",
              "KEGG_VEGF_SIGNALING_PATHWAY",
              "PID_PDGFRB_PATHWAY",
              "PID_PDGFRA_PATHWAY",
              "PID_SMAD2_3PATHWAY",
              "WP_TGFBETA_SIGNALING_PATHWAY",
              "WP_IL6_SIGNALING_PATHWAY")

env$gs.def.list <- env$gs.def.list[names(env$gs.def.list) %in% interest]

env <- pipeline.genesetStatisticSamples(env)
env <- pipeline.genesetStatisticModules(env)

pipeline.genesetOverviews(env)

pipeline.summarySheetsSamples(env)
pipeline.summarySheetsModules(env)
