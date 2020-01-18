
############## Read and prepare the input-data ##############

prepare_input <- function() {
  # read the modelpoints
  DT.input_raw <- fread("EB/DMWV__.txt", header = TRUE, skip = "!")
  
  # pre-select those columns which could be important for further use
  DT.input_selCol <- DT.input_raw[, .(
    DGO_SPCODE_I,
    AGE_AT_ENTRY_I,
    B2_PRAEMIE_I,
    BIL_RES_I,
    DURATIONIF_M_I,
    INIT_ACC_BON_I,
    INIT_POLSTAT_I,
    POL_TERM_Y_I,
    PREM_FREQ_I,
    PREM_PAYBL_Y_I,
    RABATT_KZ_I,
    SEX_I,
    STAMM_SEG_I,
    SUM_ASSURED_I,
    VERTRIEB_I,
    ZUGANG_ART_I,
    ZV_PRAEMIE1_I
  )]
  
  # Some columns need to be one hot encoded because they have categorical data (sex,...)
  # For this analysis the columns SEX_I and STAMM_SEG_I are one hot encoded.
  # See the commented code outside the function for the one hot encoding of 
  # some other variables.
  
  
  ########## one hot encoding for SEX_I   ########## 
  DT.oH_SEX_I_tmp <- DT.input_selCol[, .(DGO_SPCODE_I, SEX_I)]
  DT.oH_SEX_I_tmp[, SEX_I := as.character(SEX_I)]
  DT.oH_SEX_I_tmp[SEX_I == "0", SEX_I := "male"]
  DT.oH_SEX_I_tmp[SEX_I == "1", SEX_I := "female"]
  DT.oH_SEX_I <- dcast(data = DT.oH_SEX_I_tmp,
                             DGO_SPCODE_I ~ SEX_I,
                             fun = length)
  
  ########## one hot encoding for STAMM_SEG_I   ##########
  DT.oH_STAMM_SEG_I_tmp <- DT.input_selCol[, .(DGO_SPCODE_I, STAMM_SEG_I)]
  DT.oH_STAMM_SEG_I_tmp[, STAMM_SEG_I := as.character(STAMM_SEG_I)]
  DT.oH_STAMM_SEG_I_tmp[STAMM_SEG_I == "0", STAMM_SEG_I := "NoStammSeg"]
  DT.oH_STAMM_SEG_I_tmp[STAMM_SEG_I == "1", STAMM_SEG_I := "StammSeg"]
  DT.oH_STAMM_SEG_I <- dcast(data = DT.oH_STAMM_SEG_I_tmp,
                             DGO_SPCODE_I ~ STAMM_SEG_I,
                             fun = length)
   
  # add the one hot encoded columns 
  DT.input_tmp <- cbind(DT.input_selCol,
                        DT.oH_SEX_I[, .(female, male)],
                        DT.oH_STAMM_SEG_I[, .(NoStammSeg, StammSeg)])
  
  
  # take a sample which forms a homogeneous group 
  DT.input_tmp <- DT.input_tmp[INIT_POLSTAT_I == 1               # premium paying phase
                        & PREM_FREQ_I    == 12             # monthly payments
                        & RABATT_KZ_I    == "N"            # no discounts
                        & VERTRIEB_I     == "MAKLER"       # sold via broker
                        & ZUGANG_ART_I %in% c("N", "0"),]  # new policy or top up 
  
  # take only those columns which provide the key indicators of those policies
  DT.input <- DT.input_tmp[, .(
    DGO_SPCODE_I,   # Unique policy-ID for every contract 
    AGE_AT_ENTRY_I, 
    B2_PRAEMIE_I,   # Premium 
    BIL_RES_I,      # Actuarial reserve
    DURATIONIF_M_I, # Months since contract start
    INIT_ACC_BON_I, # Surpluses
    POL_TERM_Y_I,   # Duration of contract in years
    PREM_PAYBL_Y_I, # Premium payment period in years
    male,           # one hot from SEX_I
    female,         # one hot from SEX_I
    StammSeg,       # one hot from STAMM_SEG_I
    NoStammSeg,     # one hot from STAMM_SEG_I
    SUM_ASSURED_I,  # Sum assured
    ZV_PRAEMIE1_I   # Premium for riders
  )]
  
  return(DT.input)
}




# ########## one hot encoding for INIT_POLSTAT_I   ########## 
# oH_INIT_POLSTAT_I_tmp <- DT.input_selCol[, .(DGO_SPCODE_I, INIT_POLSTAT_I)]
# oH_INIT_POLSTAT_I_tmp <- oH_INIT_POLSTAT_I_tmp[, INIT_POLSTAT_I := as.character(INIT_POLSTAT_I)]
# oH_INIT_POLSTAT_I_tmp[INIT_POLSTAT_I == "1",  INIT_POLSTAT_I := "Polpämpfl"]
# oH_INIT_POLSTAT_I_tmp[INIT_POLSTAT_I == "3",  INIT_POLSTAT_I := "Polpämfr"]
# oH_INIT_POLSTAT_I_tmp[INIT_POLSTAT_I == "4",  INIT_POLSTAT_I := "Polpämfr2"]
# oH_INIT_POLSTAT_I <- dcast(data = oH_INIT_POLSTAT_I_tmp,
#                            DGO_SPCODE_I ~ INIT_POLSTAT_I,
#                            fun = length)
# 
# ########## one hot encoding for PREM_FREQ_I   ########## 
# oH_PREM_FREQ_I_tmp <- DT.input_selCol[, .(DGO_SPCODE_I, PREM_FREQ_I)]
# oH_PREM_FREQ_I_tmp <- oH_PREM_FREQ_I_tmp[, PREM_FREQ_I := as.character(PREM_FREQ_I)]
# oH_PREM_FREQ_I_tmp[PREM_FREQ_I == "1",  PREM_FREQ_I := "jährlich"]
# oH_PREM_FREQ_I_tmp[PREM_FREQ_I == "12", PREM_FREQ_I := "mtl"]
# oH_PREM_FREQ_I_tmp[PREM_FREQ_I == "2",  PREM_FREQ_I := "halbjäh"]
# oH_PREM_FREQ_I_tmp[PREM_FREQ_I == "4",  PREM_FREQ_I := "quartä"]
# oH_PREM_FREQ_I <- dcast(data = oH_PREM_FREQ_I_tmp,
#                         DGO_SPCODE_I ~ PREM_FREQ_I,
#                         fun = length)


########## one hot encoding for VERTRIEB_I   ########## 
# oH_VERTRIEB_I_tmp <- DT.input_selCol[, .(DGO_SPCODE_I, VERTRIEB_I)]
 # oH_VERTRIEB_I <- dcast(data = oH_VERTRIEB_I_tmp,
 #                       DGO_SPCODE_I ~ VERTRIEB_I,
 #                       fun = length)




# one hot encoding for RABATT_KZ_I --> nur Fälle mit "N" nehmen! --> kein Encoding
# DT.input_raw[, .N, by = .(RABATT_KZ_I)]



# one hot encoding for ZUGANG_ART_I --> in Prophet alles gleich behandelt solange
#                                       K oder P

# DT.input_raw[, .N, by = .(ZUGANG_ART_I)]

#" 0" = Neuzugang
#" 1" = Anpassungsbrief (S) / Dynamisierung (D)
#" 2" = Ersatz (mit Rückdatierung)
#" 4" = Reduktion
#" 5" = Reaktivierung
#" 8" = Indexklauselerhöhung bei Großleben Anpassung bei Kollektiv (S)
#" 9" = Interner Zugang (infolge Ablauf der Prämienzahlung)
#" K" = Konvertierung (S) / Ersatz bei prämienfreier Summe (D)
#" N" = Nachtrag, Aufstockung, Indexbrief
#" P" = Prämienfreie Einrechnung (S)
#" S" = Prolongation (nur Großleben) (S) / Verlängerung (D)
#" T" = Versorgertod
#" W" = Switch ohne Tarifänderung ( FLV )
#" Y" = Switch mit Tarifänderung ( FLV )
#" Z" = Summenabfall wegen Zwischenauszahlung
#"  " bei Änderungsarten, die keinen Zugang inkludieren


#input_selCol[, .N,
#      by = .(INIT_POLSTAT_I, PREM_FREQ_I, RABATT_KZ_I, VERTRIEB_I, ZUGANG_ART_I)] %>% fwrite("analyse.csv") 
