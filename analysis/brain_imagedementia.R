
library(dplyr)
#bids <- "P:/ostpre_bids"
#bids <- "/research/groups/ostpre/ostpre_bids"

qs::qload(file=file.path(bids,"clinical/mem_dates.qs"))
qs::qload(file=file.path(bids,"nifti/tags.qs"))

graph_titles <- list(
  'average_thickness'="Average thickness",
  'entorhinal_thickness'="Entorhinal thickness",
  'GM_volume'="GM volume",
  'hippocampus'="Hippocampus",
  'Jack_signature_CT'="Jack Signature CT",
  'quality_IQR'="Quality IQR",
  'TIV'="TIV",
  'ventricle'="Ventricle"
)

set_labels <- function(variable,value){
  return(graph_titles[value])
}


tags |> distinct(lomno1) |> count()   # Henkilöitä NIfTI konversion jälkeen
tags |> distinct(lomno1,ses) |> count()   # Sessioita NIfTI konversion jälkeen
tags |> count() # Erillisiä kuvia/"sekvenssejä"


# Irrotetaan metadatoista analyyseihin tarvittavat tagit
reltags <- tags |>
  mutate(
    StudyID=as.character(StudyInstanceUID),
    SeriesID=as.character(SeriesInstanceUID),
    MF=as.character(Manufacturer),
    SN=as.character(DeviceSerialNumber),
    MagStre=as.character(MagneticFieldStrength)
  ) |>
  select(lomno1,ses,StudyID,SeriesID,MagneticFieldStrength,Manufacturer,ManufacturersModelName,DeviceSerialNumber,StationName,MF,SN,MagStre) |>
  distinct() |>
  mutate(
    SN=case_when(
      SN=="NULL" ~ "25345",     # Manuaalinen korjaus eli jos SerialNumber (jota ei kylläkään käytetä varsinaisissa analyyseissa) puuttuu, niin näyttävät olevan hyvin samanlaisia kuin tämän laitteen (muut) tulokset
      TRUE ~ SN
    )
  )

imgdata_in <- readxl::read_excel(file.path(bids,"bids/derivatives/summary_measures/ostpre_CAT12_summary_measures.xlsx")) |>
  mutate(
    lomno1=as.numeric(gsub("sub-","",subj)),
    ses=as.numeric(gsub("ses-","",session))
  )

outliers <- readxl::read_excel(file.path(bids,"bids/derivatives/summary_measures/ostpre_outliers.xlsx"))

outliers |> count()
outliers |> distinct(subj) |> count()

outlierdata <- imgdata_in |>  
  inner_join(outliers, by=c("subj","session")) |>
  mutate(qual=case_when(
    quality_percent>68 ~ ">68",
    TRUE ~ "<=68"
  ))

outlierdata |> count(qual)

imgdata <- imgdata_in |>  
  anti_join(outliers, by=c("subj","session")) # Poistetaan poikkeavat kuvat (todettu myös manuaalisesti kelvottomiksi analyysiin)

sessions <- readr::read_delim(file.path(bids,"bids/sessions.csv"),delim=";",col_names=FALSE) |>
  mutate(
    lomno1=as.numeric(substr(X6,4,8)),
    ses=as.numeric(gsub("_","",substr(X6,14,15))),
    sesdate=as.Date(X1)
  ) |>
  arrange(lomno1,ses) |>
  select(lomno1,ses,sesdate,StudyID=X2,SeriesID=X3,SeriesName=X4) |>
  left_join(imgdata, by=c("lomno1","ses")) |>
  left_join(mem_dates |> rename(lomno1=LOMNO1),by=c("lomno1")) |>
  left_join(reltags, by=c("SeriesID")) |>
  mutate(
    memtype=case_when( # memtypen avulla ryhmittely
      dem-366 < sesdate ~ "Dementia",
      atc-366 < sesdate ~ "ATC-Dementia", 
      mci-366 < sesdate ~ "MCI",  # Mild Congnitive Impairment
      mem-366 < sesdate ~ "SMC",  # Subjective Memory Complaint
      dem-366 < mci ~ "Pre MCI-Dementia",
      dem-366 < mem ~ "Pre SMC-Dementia",
      is.na(dem) & !is.na(mci) ~ "Pre MCI",
      is.na(dem) & !is.na(mem) ~ "Pre SMC",
      !is.na(dem) ~ "Pre Dementia",
      is.na(dem) & is.na(mem) & is.na(atc) ~ "None",
      is.na(dem) & is.na(mem) & !is.na(atc) ~ "Pre ATC-Dementia"
    ),
    img_memtype=as.factor(memtype),
    age=floor(lubridate::interval(spvm,sesdate)/lubridate::dyears(1)),
    Age=lubridate::interval(spvm,sesdate)/lubridate::dyears(1),
    mt=factor(case_when(
      memtype %in% c("Dementia","ATC-Dementia") ~ "Dementia",
      grepl("^Pre|^None",memtype) ~ "None",   # None incl pre
      TRUE ~ memtype
    ), levels=c("Dementia","MCI","SMC","None")),
    ORIGPROT=case_when(
      MF %in% c("GE","Toshiba") ~ "Other",
      MF %in% c("Philips") ~ "Philips",
      TRUE ~ paste0(MF,SN)
      ),
    FS=case_when(
      MagStre==3 ~ "3Tesla",
      TRUE ~ "1.5Tesla"
    )
  )

sessions |> count()
sessions |> distinct(lomno1.x) |> count()

imgdata |> count()
imgdata |> distinct(lomno1) |> count()

imgdata_in |> count()
imgdata_in |> distinct(lomno1) |> count()


sessions |> count(MF,SN)
sessions |> filter(quality_percent>68) |> count(MF,SN) # Vain yli 68 laatuiset valitaan
sessions |> count(FS)

memdist <- sessions |>
  filter(quality_percent>10) |>
  count(age,mt) |>
  tidyr::pivot_wider(
    id_cols=age,
    names_from="mt",
    values_from="n",
    values_fill = 0
  ) |>
  select(age,2,3,5,4)

#openxlsx::write.xlsx(memdist,file="memagedist.xlsx")

agedist <- sessions |>
  count(age)


library(ggplot2)
ls <- sessions |>
  filter(quality_percent>68) |>
  tidyr::pivot_longer(
    cols=quality_percent:ventricle
  )


ls |>
  ggplot(
    aes(x='',y=value)
  ) +
  geom_boxplot(aes(fill=mt),position=position_dodge(.9), notch=TRUE) +
  facet_wrap(~name,scale="free", labeller=set_labels) +
  stat_summary(fun="mean",aes(group=mt),position=position_dodge(.9), geom="point", shape=4) +
  labs(x='', y='')
  
adni <- readxl::read_excel(file.path(bids, "bids/derivatives/summary_measures/ADNI_CAT12_ostprematched_summary_measures_v2.xlsx")) |>
  filter(DX_bl != "NA") |>
  mutate(
    img_memtype=factor(case_when(
      DX_bl %in% c("AD") ~ "Dementia",
      DX_bl %in% c("EMCI","LMCI") ~ "MCI",
      DX_bl %in% c("SMC") ~ "SMC",
      TRUE ~ "None"
    )),
    mt=factor(case_when(
      DX_bl %in% c("AD") ~ "Dementia - ADNI",
      DX_bl %in% c("EMCI","LMCI") ~ "MCI - ADNI",
      DX_bl %in% c("SMC") ~ "SMC - ADNI",
      TRUE ~ "None - ADNI"
    )),
    MF=case_when(
      MANUFACTURER==1 ~ "Siemens",
      MANUFACTURER==2 ~ "Philips",
      MANUFACTURER==3 ~ "GE",
    ),
    FS=case_when(
      grepl("^3",FLDSTRENG) ~ "3Tesla",
      TRUE ~ "1.5Tesla"
    )
  ) |>
  rename(adni_subj=subj)
  
  
  
library(ggplot2)
ls2 <- adni |>
  filter(quality_percent>68) |>
  tidyr::pivot_longer(
    cols=quality_percent:ventricle
  ) 

ls2 |>
  ggplot(
    aes(x='',y=value)
  ) +
  geom_boxplot(aes(fill=mt),position=position_dodge(.9), notch=TRUE) +
  facet_wrap(~name,scale="free", labeller=set_labels) +
  stat_summary(fun="mean",aes(group=mt),position=position_dodge(.9), geom="point", shape=4) +
  labs(x='', y='')


  
yhd <- adni |>
  bind_rows(sessions) |>
  filter(quality_percent>68)

yhd |> count(MF,FS)

library(ggplot2)

yhd |> 
  ggplot(
    aes(x=quality_percent, fill=MF)
  ) +
  geom_density(alpha=0.6)

yhd |> 
  ggplot(
    aes(x=quality_percent, fill=FS)
  ) +
  geom_density(alpha=0.6)


yls <- yhd |>
  tidyr::pivot_longer(
    cols=quality_percent:ventricle
  )

yls |>
  ggplot(
    aes(x='',y=value)
  ) +
  geom_violin(aes(fill=mt),position=position_dodge(.9), notch=TRUE) +
  facet_wrap(~name,scale="free", labeller=set_labels) +
  labs(x='', y='')



yhd_ana <- yhd |>
  mutate(
    subj=factor(case_when(
      is.na(subj) ~ sprintf("ADNI-%05d",adni_subj),
      TRUE ~ sub("sub","OSTPRE",subj)
    )),
    gr=factor(case_when(
      grepl("^ADNI",subj) ~ "ADNI",
      TRUE ~ "OSTPRE"
    )),
    mtyp=case_when(
      gr=="OSTPRE" ~ mt,
      TRUE ~ img_memtype
    ),
    mpf=factor(case_when(
      grepl("^Dementia",mtyp) ~ "3-Dementia",
      grepl("^MCI",mtyp) ~ "2-MCI",
      grepl("^SMC",mtyp) ~ "1-SMC",
      grepl("^None",mtyp) ~ "0-None"
    ), labels=c("None","SMC","MCI","Dementia")),
    cage=Age-75,
    MF=factor(MF),
    FS=factor(FS)
  ) |>
  select(gr,subj,session,Age,cage,mpf,MF,FS,TIV:ventricle)


yhd_ana |> count(gr)
yhd_ana |> filter(gr=="OSTPRE") |> distinct(subj) |> count()
yhd |> filter(!is.na(subj)) |> count()
yhd |> filter(!is.na(subj)) |> distinct(subj) |> count()

# Table 1 tietoja
yhd_ana |>
  group_by(gr) |>
  summarise(
    n=n(),
    nsub=n_distinct(subj),
    age=mean(Age,na.rm=TRUE),
    agesd=sd(Age,na.rm=TRUE),
  )

yhd_ana |>
  count(gr,mpf) |>
  group_by(gr) |>
  mutate(np=n/sum(n)*100)

yhd_ana |>
  count(gr,FS) |>
  group_by(gr) |>
  mutate(np=n/sum(n)*100)

yhd_ana |>
  count(gr,MF) |>
  group_by(gr) |>
  mutate(np=n/sum(n)*100)


md <- yhd_ana |>
#  select(gr,subj,Age,cage,mpf,MF,FS,TIV,meas=hippocampus) |>
#  select(gr,subj,Age,cage,mpf,MF,FS,TIV,meas=ventricle) |>
  select(gr,subj,Age,cage,mpf,MF,FS,TIV,meas=GM_volume) |>
#  select(gr,subj,Age,cage,mpf,MF,FS,TIV,meas=average_thickness) |>
#  select(gr,subj,Age,cage,mpf,MF,FS,TIV,meas=Jack_signature_CT) |>
#  select(gr,subj,Age,cage,mpf,MF,FS,TIV,meas=entorhinal_thickness) |>
  filter(!is.na(meas) & !is.na(Age)) |>
  mutate(
    smeas=scale(meas),
    nmeas=meas-mean(meas),
    sTIV=scale(TIV)
    )

md |> count(gr,mpf)


m1 <- lm(smeas ~ splines::bs(Age) + splines::bs(sTIV) + gr*mpf + MF + FS, data=md)
# m1 <- nlme::lme(smeas ~ splines::bs(Age) + splines::bs(sTIV) + gr*mpf + MF + FS, random= ~ 1 | subj, data=md)


mdp1 <- ggeffects::predict_response(m1, terms=c("Age [65:90]"), vcov_fun="vcovCR", vcov_type="CR0", vcov_args=list(cluster=md$subj))
# mdp1 <- ggeffects::predict_response(m1, terms=c("Age [65:90]","gr"), vcov_fun="vcovCR", vcov_type="CR0", vcov_args=list(cluster=md$subj))
plot(mdp1)

mdp2 <- ggeffects::predict_response(m1, terms=c("sTIV [-3:3]"), vcov_fun="vcovCR", vcov_type="CR0", vcov_args=list(cluster=md$subj))
plot(mdp2)

mdp3 <- ggeffects::predict_response(m1, terms=c("Age [65:85]","mpf", "gr"), vcov_fun="vcovCR", vcov_type="CR0", vcov_args=list(cluster=md$subj))
plot(mdp3)

mdp <- ggeffects::predict_response(m1, terms=c("mpf","gr"), vcov_fun="vcovCR", vcov_type="CR0", vcov_args=list(cluster=md$subj))
plot(mdp)

tdp <- ggeffects::test_predictions(mdp)
tdp # Täältä näkyvät kiinnostavat kontrastit (OSTPRE-OSTPRE ja ADNI-ADNI), mutta on ylimääräisiä mukana

dd_tpd <- ggeffects::test_predictions(
  mdp, 
  test=c(
    "(b1 - b4) = (b5 - b8)",   # None vs Dementia difference between ADNI and OSTPRE
    "(b1 - b3) = (b5 - b7)",   # None vs MCI difference between ADNI and OSTPRE
    "(b1 - b2) = (b5 - b6)",   # None vs SMC difference between ADNI and OSTPRE
    "(b3 - b4) = (b7 - b8)"    # MCI vs Dementia difference between ADNI and OSTPRE
    )
  )
dd_tpd # Tämä lisäksi OSTPRE-ADNI vertailu kiinnostavien asioiden osalta

