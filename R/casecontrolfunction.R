casecontrolfunction <- function(data_count){
  data_count <- data_count[, data_count[1, ] != "scaled_estimate"][-1, ]
  sampleID_count <- colnames(data_count)[-1]
  # leave sample 11B out
  index <- substr(sampleID_count, 14, 16)
  data_control <- data_count[, -1][, index == "11A"]
  data_tumor <- data_count[, -1][, index == "01A"]
  gene_id <- data_count[, 1]
  data_control <- cbind(gene_id, data_control)
  data_tumor <- cbind(gene_id, data_tumor)
  # leave gene ? out
  df <- data_control[, 1]
  n <- grep("\\?", df)
  data_control <- data_control[-n, ]
  dg <- data_tumor[, 1]
  m <- grep("\\?", dg)
  data_tumor <- data_tumor[-m, ]
  # sample matching
  sampleID_control <- colnames(data_control)[-1]
  sampleID_tumor <- colnames(data_tumor)[-1]
  ID_gene_normal <- substr(sampleID_control, 1, 12)
  ID_gene_tumor <- substr(sampleID_tumor, 1, 12)
  commonID <- ID_gene_tumor[ID_gene_tumor %in% ID_gene_normal]
  ID_gene_normal0 <- c("ID_gene_normal", ID_gene_normal)
  ID_gene_tumor0 <- c("ID_gene_tumor", ID_gene_tumor)
  data_tumor <- rbind(ID_gene_tumor0, data_tumor)
  data_control <- rbind(ID_gene_normal0, data_control)
  gene_id <- data_tumor[, 1]
  data_case0 <- data_tumor[, which(data_tumor[1, ] %in% commonID)]
  data_case1 <- cbind(gene_id, data_case0)
  gene_id <- data_control[, 1]
  data_control0 <- data_control[, which(data_control[1, ] %in% commonID)]
  data_control1 <- cbind(gene_id, data_control0)
  # merge the case-control table
  case_control_data <- merge(data_case1[-1, ], data_control1[-1, ], by = "gene_id")
  case_control_data[, 2:ncol(case_control_data)] <- lapply(case_control_data[, 2:ncol(case_control_data)], as.numeric)
  return(case_control_data)
  }

