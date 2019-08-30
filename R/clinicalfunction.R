clinicalfunction <- function(clinical,commonID22){
  clinical$gender <- ifelse(clinical$gender == "MALE", 1, 2)
  clinical$vital_status <- ifelse(clinical$vital_status == "Alive", 1, 2)
  clinical$tobacco_smoking_pack_years_smoked <- ifelse(clinical$tobacco_smoking_pack_years_smoked == "[Not Available]", NA, clinical$tobacco_smoking_pack_years_smoked)
  clinical$tobacco_smoking_pack_years_smoked <- as.numeric(clinical$tobacco_smoking_pack_years_smoked)
  for (i in 1:length(clinical$death_days_to)) {
    clinical$death_days_to[i] <- ifelse(clinical$death_days_to[i] == "[Not Available]", clinical$last_contact_days_to[i], clinical$death_days_to[i])
    clinical$death_days_to[i] <- ifelse(clinical$death_days_to[i] == "[Not Applicable]", clinical$last_contact_days_to[i], clinical$death_days_to[i])
  }
  clinical$death_days_to <- ifelse(clinical$death_days_to == "[Not Available]", NA, clinical$death_days_to)
  clinical$death_days_to <- ifelse(clinical$death_days_to == "[Discrepancy]", NA, clinical$death_days_to)
  clinical$death_days_to <- as.numeric(clinical$death_days_to)
  clinical <- clinical[!is.na(clinical$death_days_to), ]
  clinical$bcr_patient_barcode <- paste(clinical$bcr_patient_barcode, "01A", sep = "-")
  clinical0 <- clinical[clinical$bcr_patient_barcode %in% commonID22, ]
  commonID222 <- clinical0[, 3]
  clinical0$tobacco_smoking_pack_years_smoked <- ifelse(clinical0$tobacco_smoking_pack_years_smoked == "[Not Available]", NA, clinical0$tobacco_smoking_pack_years_smoked)
  clinical1 <- clinical0[!is.na(clinical0$tobacco_smoking_pack_years_smoked), ]
  clinical1 <- clinical1[clinical1$death_days_to != 0, ]
  clinical1$age_at_initial_pathologic_diagnosis <- as.numeric(clinical1$age_at_initial_pathologic_diagnosis)
  clinical1$followupage <- (-as.numeric(as.character(clinical1$birth_days_to)))/365
  clinical1$followupage <- floor(clinical1$followupage)
  clinical1$age <- ifelse(clinical1$followupage != clinical1$age_at_initial_pathologic_diagnosis, clinical1$followupage, clinical1$age_at_initial_pathologic_diagnosis)
  clinical1$age_at_initial_pathologic_diagnosis <- ifelse(is.na(clinical1$age), clinical1$age_at_initial_pathologic_diagnosis,
                                                          clinical1$age)
  clinical2 <- clinical1[!is.na(clinical1$age_at_initial_pathologic_diagnosis), ]
  return(clinical2)
}