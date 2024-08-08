#R Code of the proyect: Visceral adipose tissue is associated with major adverse cardiovascular events in patients with premature coronary artery disease: a prospective sub-analysis of the GEA study.
## Data Analysis: Neftali Eduardo Antonio-Villa (neftalivilla@comunidad.unam.mx)
## Latest version of Analysis: 8-March-2024
## Any question regarding analysis, please contact Neftali Eduardo Antonio-Villa

#####Library#####

library(tidyverse)
library(readxl)
library(dplyr)
library(mice)
library(haven)
library(epiR)
library(ggpubr)
library(rstatix)
library(ggthemes)
library(patchwork)
library(gtsummary)
library(data.table)
library(rms)
library(patchwork)
library(ggpmisc)
library(forestplot)
library(pROC)
library(rmda)
library(ggplotify)
library(ggbreak)
library(survminer)
library(Epi)
library(lubridate)
library(MatchIt)


#####Dataset: Read Data and Setting Labels#####

setwd("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/VAT - MACE")

#Base Basal
base.original <- read_sav("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/Base GEA/BASE_BASAL para compartir_2.sav", encoding="latin1")
base.original<-janitor::clean_names(base.original)
base.original<-base.original%>%
  mutate(id=paste0(muestra,numero))

base.original<-labelled::remove_labels(base.original)%>%dplyr::filter(muestra=="GP")

#Base Prematuros (Casos)
base.prem <-read_sav("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/MCP-1/pacientes con EAC prematura.sav", encoding="latin1")
base.prem <- janitor::clean_names(base.prem)
base.prem <- base.prem%>%mutate(id=paste0(muestra,numero))%>%dplyr::select(-c(muestra,numero))

#Base MACE (Casos)
base.mace <-read_sav("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/MCP-1/SegGEAP_MACEyMortalidad_estatinas.sav", encoding="latin1")
base.mace <- janitor::clean_names(base.mace)
base.mace <- base.mace%>%dplyr::select(gp,mortalidad_total,mace_no_fatal,mace_fatal,seguimiento_anos,edad_eac,fecha_muerte_o_ultima_visita,muerte_string)
base.mace$muestra<-c("GP")
base.mace <- base.mace%>%mutate(id=paste0(muestra,gp))%>%dplyr::select(-c(muestra))


#Base de Fecha de ECA (NEAV)
base.NEAV.ECA<-read_excel("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/Base GEA/BASE_EDAD_NEAV.xlsx")
base.NEAV.ECA<-janitor::clean_names(base.NEAV.ECA)%>%
  mutate(muestra=c("GP"))%>%
  mutate(id=paste0(muestra,gp))%>%
  dplyr::select(id,fecha_primer_eca_neav)

#Base de Pacientes Nuevos
base.AIDA.nuevos<-read_excel("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/VAT - MACE/TIPO_MACE_V2.xlsx",sheet = "AÑADIR CON SEGUIMIENTO")
base.AIDA.nuevos<-janitor::clean_names(base.AIDA.nuevos)%>%
  mutate(muestra=c("GP"))%>%
  mutate(id=paste0(muestra,gp))%>%
  left_join(base.original%>%dplyr::select(id,fechade_entrevista),by = "id")%>%
  mutate(seguimiento_anos=case_when(!is.na(tipo_evento)~ as.numeric(as.Date(fecha_evento, format="%d/%m/%Y")-as.Date(fechade_entrevista, format="%d/%m/%Y"))/365,
                                    is.na(tipo_evento)~ as.numeric(as.Date(fecha_vis_segui, format="%d/%m/%Y")-as.Date(fechade_entrevista, format="%d/%m/%Y"))/365))%>%
  dplyr::select(id,tipo_evento,seguimiento_anos,no_ea_cp,descart)


##Tipos de MACE
base.tipo.mace<-read_excel("TIPO_MACE_V2.xlsx", sheet = "Tipo MACE")
base.tipo.mace <- janitor::clean_names(base.tipo.mace)
base.tipo.mace<-base.tipo.mace%>%
  mutate(muestra=c("GP"))%>%
  mutate(id=paste0(muestra,gp))%>%
  dplyr::select(id,tipo_mace)

#Base Final
base.original<-base.original%>%
  left_join(base.prem,by = "id")%>% #Casos Prematures
  left_join(base.mace,by = "id")%>% #Casos MACE
  left_join(base.AIDA.nuevos,by = "id") %>% #Nuevos Casos
  left_join(base.tipo.mace,by = "id")


base.original<-base.original %>% 
  mutate(seguimiento_anos = purrr::pmap_dbl(select(., c(seguimiento_anos.x, seguimiento_anos.y)), pmax, na.rm=TRUE))%>%
  mutate(no_ea_cp = purrr::pmap_dbl(select(., c(no_ea_cp.x, no_ea_cp.y)), pmin, na.rm=TRUE))%>%
  mutate(descartado = purrr::pmap_dbl(select(., c(descartado.x, descartado.y)), pmin, na.rm=TRUE))

#Depurar 
base<-base.original%>%
  dplyr::filter(no_ea_cp==0)%>%
  dplyr::filter(descartado==0)%>%
  dplyr::filter(!is.na(seguimiento_anos))%>%
  dplyr::filter(!is.na(gav))

nrow(base)


#####Dataset: Multiple Imputation Analysis#####

base2<-base%>%dplyr::select(id,circunfer_cintura,gav,edad_eac,edad,cldl,hldl,pcr,promsisto,promdiasto,
                              gas,rel_lhd_bazo,homa,adiponectina,rita,creat,apo_b)
base2_imp<-mice::mice(base2, m=5, maxit=5,seed = 123)
base2_imp_2<-complete(base2_imp,1)
base<-base%>%dplyr::select(-c(circunfer_cintura,gav,edad_eac,edad,cldl,hldl,pcr,promsisto,promdiasto,
                              gas,rel_lhd_bazo,homa,adiponectina,rita,creat,apo_b))%>%
  left_join(base2_imp_2,by="id")


#####Dataset: Recoding Variables####

#Age Categories
base$age_group <- base$edad %>%
  cut(c(-Inf, 45,55,65, Inf), labels = 1:4, right = F)

#MACE Typos
base$mace_no_fatal[is.na(base$mace_no_fatal)]<-0
base$mace_no_fatal[!is.na(base$tipo_evento)]<-1
base$mace_fatal[is.na(base$mace_fatal)]<-0

#MACE OVERALL
base$morbimortalidad_mace<-NULL
base$morbimortalidad_mace[base$mace_fatal==1]<-1
base$morbimortalidad_mace[base$mace_no_fatal==1]<-1
base$morbimortalidad_mace<-na.tools::na.replace(base$morbimortalidad_mace,0)

#MACE CATEGORIAS
base$MACE_CAT<-NULL
base$MACE_CAT[base$mace_no_fatal==1]<-1
base$MACE_CAT[base$mace_fatal==1]<-2
base$MACE_CAT[base$morbimortalidad_mace==0]<-0

#MACE OVERALL 
base$morbimortalidad_mace_2<-base$morbimortalidad_mace


#P>75
quantile(base$gav,probs = c(0.75))
base$gav_p75<-NULL
base$gav_p75[base$gav>=211]<-1
base$gav_p75[base$gav<211]<-0

#VAT Terciles

quantile(base$gav,probs = c(0.33,0.66))
base$gav_terciles_neav<-NULL
base$gav_terciles_neav[base$gav<142.00]<-1
base$gav_terciles_neav[base$gav>=142.00 & base$gav<194.92]<-2
base$gav_terciles_neav[base$gav>=194.92]<-3
base$gav_terciles_neav_2<-base$gav_terciles_neav

#Age ≥65 years
base$edad_65<-NULL
base$edad_65[base$edad>=65]<-1
base$edad_65[base$edad<65]<-0

#Non HDL ≥100 mg/dl
base$NO_HDL_100<-NULL
base$NO_HDL_100[base$cno_hdl>=100]<-1
base$NO_HDL_100[base$cno_hdl<100]<-0

#Glucose ≥100 mg/dl
base$gluc_100<-NULL
base$gluc_100[base$glucosa>=100]<-1
base$gluc_100[base$glucosa<100]<-0

#HDL-C
base$chdl_cat<-NULL
base$chdl_cat[base$sexo==1 & base$chdl<=50]<-1
base$chdl_cat[base$sexo==0 & base$chdl<=40]<-1
base$chdl_cat<-na.tools::na.replace(base$chdl_cat,0)

#Estatinas
base$Estatinas_REC<-NULL
base$Estatinas_REC[base$estatinas==1]<-1
base$Estatinas_REC[is.na(base$Estatinas_REC)]<-0

#AntiHas
base$uso_antihertensivo[is.na(base$uso_antihertensivo)]<-0

#FECHA DE PRIMER ECA
base$edad.primer.eca.neav<-base$edad_eac

#Dislipidemia
base$dislipidemia_neav
base$dislipidemia_neav[base$dislipidemia1==1]<-1
base$dislipidemia_neav[base$ct>=200]<-1
base$dislipidemia_neav[base$tg>=150]<-1
base$dislipidemia_neav[base$sexo==1 & base$hldl>=50]<-1
base$dislipidemia_neav[base$sexo==0 & base$hldl>=40]<-1
base$dislipidemia_neav<-na.tools::na.replace(base$dislipidemia_neav,0)

#Año de Entrevista
base$year<-as.numeric(format(as.Date(base$fechade_entrevista, format="%d/%m/%Y"),"%Y"))
base$year[is.na(base$year)]<-2008

#Typos
base$seguimiento_anos<-as.numeric(base$seguimiento_anos)

#CKD
base$TFG_CKD_EPI<-nephro::CKDEpi_RF.creat(base$creat, base$sexo, base$edad)

#Estimar Tiempo Seguimiento
# base$tiempo_seguimiento_ECA_neav<-as.numeric(as.Date(base$fecha_muerte_o_ultima_visita, format="%d/%m/%Y")-as.Date(base$fecha_primer_eca_neav, format="%d/%m/%Y"))/365
# base$tiempo_seguimiento_ECA_neav<-na.tools::na.replace(base$tiempo_seguimiento_ECA_neav,base$seguimiento_anos)

#Fecha de Reclutamiento Lexit
base$fechade_entrevista[is.na(base$fechade_entrevista)]<-"2008-07-03"
base$date_recruited<-decimal_date(as.Date(base$fechade_entrevista, format="%d/%m/%Y"))



base$AGE_CVD_CAT<-NULL
base$AGE_CVD_CAT[base$edad_eac<35]<-0
base$AGE_CVD_CAT[base$edad_eac>=35 & base$edad_eac<=65]<-1
base$AGE_CVD_CAT<-na.tools::na.replace(base$AGE_CVD_CAT,1)

base$BMI_CAT<-NULL
base$BMI_CAT[base$imc<25]<-1
base$BMI_CAT[base$imc>=25 & base$imc<=30]<-2
base$BMI_CAT[base$imc>=30]<-3

#LDL <70 mg/dl
base$cldl_cat<-NULL
base$cldl_cat[base$cldl<70]<-0
base$cldl_cat[base$cldl>=70]<-1

#Trigliceridos
base$tg_cat<-NULL
base$tg_cat[base$tg<150]<-0
base$tg_cat[base$tg>=150]<-1

#No-HDLC
base$nhdl_cat<-NULL
base$nhdl_cat[base$cno_hdl<130]<-0
base$nhdl_cat[base$cno_hdl>=130]<-1

#APO-B
quantile(base$apo_b,probs = c(0.75),na.rm = T)
base$apob_p75_neav<-NULL
base$apob_p75_neav[base$apo_b>=102]<-1
base$apob_p75_neav[base$apo_b<102]<-0

#HOMA IR
base$homa_p75_neav<-NULL
base$homa_p75_neav[base$homa>=7.145616]<-1
base$homa_p75_neav[base$homa<7.145616]<-0

#ADIPO-IT
quantile(base$rita,probs = c(0.75))
base$rita_p25_neav<-NULL
base$rita_p25_neav[base$adiponectina>=14.75652]<-1
base$rita_p25_neav[base$adiponectina<14.75652]<-0

#PCR
base$pcr_cat<-NULL
base$pcr_cat[base$pcr>=3]<-1
base$pcr_cat[base$pcr<3]<-0

#Adiponectina
quantile(base$adiponectina,probs = c(0.25))
base$adipo_p25_neav<-NULL
base$adipo_p25_neav[base$adiponectina>=3.05]<-1
base$adipo_p25_neav[base$adiponectina<3.05]<-0

#Total adipossity
quantile(base$gat,probs = c(0.75))
base$GAT_p75_neav<-NULL
base$GAT_p75_neav[base$gat>=530.5]<-1
base$GAT_p75_neav[base$gat<530.5]<-0

#Subcutanea
quantile(base$gas,probs = c(0.75))
base$GAS_p75_neav<-NULL
base$GAS_p75_neav[base$gas>=323]<-1
base$GAS_p75_neav[base$gas<323]<-0

#####Dataset: Label Factors#####

base$morbimortalidad_mace<-factor(base$morbimortalidad_mace,labels = c("Without MACE","Incident MACE"))
base$tabaquismo_cat<-factor(base$tabaquismo_cat,labels = c("Never","Previous","Current"))
base$gav_terciles_neav<-factor(base$gav_terciles_neav,labels = c("Lower-Tertile\n(<140 cm\u00b2)",
                                                                 "Middle-Tertile\n(140-194 cm\u00b2)",
                                                                 "Upper-Tertile\n(>194 cm\u00b2)"))

#####Dataset: Lexit Dataset#####

#Overall MACE
base_lexis <- Lexis(entry = list("period" = date_recruited, 
                                 "age" = edad, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + seguimiento_anos), 
                    data = base, exit.status = morbimortalidad_mace_2) 

base_fin2<-splitLexis(base_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
base_fin2$age.cat <- timeBand(base_fin2, "age", type = "factor")
base_fin2$premature<-ifelse(base_fin2$lex.Xst==T,1,0); base_fin2$premature[base_fin2$premature==1 & base_fin2$ageout>=75]<-0 

base_gea_lexit <- base_fin2 %>%
  rename("id.lex" = lex.id, 
         "fu_period" = lex.dur, 
         "entry_status" = lex.Cst, 
         "exit_status" = lex.Xst) %>% 
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)

#Time Lexit
base_gea_lexit$time_lexit<-base_gea_lexit$time_at_exit-base_gea_lexit$time_at_entry

##NON-Fatal MACE
base_lexis_NF_MACE <- Lexis(entry = list("period" = date_recruited, 
                                 "age" = edad, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + seguimiento_anos), 
                    data = base%>%dplyr::filter(mace_fatal==0), exit.status = mace_no_fatal) 

base_fin2_NF_MACE<-splitLexis(base_lexis_NF_MACE, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
base_fin2_NF_MACE$age.cat <- timeBand(base_fin2_NF_MACE, "age", type = "factor")
base_fin2_NF_MACE$premature<-ifelse(base_fin2_NF_MACE$lex.Xst==T,1,0)
base_fin2_NF_MACE$premature[base_fin2_NF_MACE$premature==1 & base_fin2_NF_MACE$ageout>=75]<-0 

base_gea_lexit_NF_MACE <- base_fin2_NF_MACE %>%
  rename("id.lex" = lex.id, 
         "fu_period" = lex.dur, 
         "entry_status" = lex.Cst, 
         "exit_status" = lex.Xst) %>% 
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)

#Time Lexit
base_gea_lexit_NF_MACE$time_lexit<-base_gea_lexit_NF_MACE$time_at_exit-base_gea_lexit_NF_MACE$time_at_entry


#Fatal MACE
base_lexis_F_MACE <- Lexis(entry = list("period" = date_recruited, 
                                         "age" = edad, "time_at_entry" = 0),
                            exit = list("period" = date_recruited + seguimiento_anos), 
                            data = base%>%dplyr::filter(mace_no_fatal==0), exit.status = mace_fatal) 

base_fin2_F_MACE<-splitLexis(base_lexis_F_MACE, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
base_fin2_F_MACE$age.cat <- timeBand(base_fin2_F_MACE, "age", type = "factor")
base_fin2_F_MACE$premature<-ifelse(base_fin2_F_MACE$lex.Xst==T,1,0)
base_fin2_F_MACE$premature[base_fin2_F_MACE$premature==1 & base_fin2_F_MACE$ageout>=75]<-0 

base_gea_lexit_F_MACE <- base_fin2_F_MACE %>%
  rename("id.lex" = lex.id, 
         "fu_period" = lex.dur, 
         "entry_status" = lex.Cst, 
         "exit_status" = lex.Xst) %>% 
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)

#Time Lexit
base_gea_lexit_F_MACE$time_lexit<-base_gea_lexit_F_MACE$time_at_exit-base_gea_lexit_F_MACE$time_at_entry


#####Dataset: Labels#####

setattr(base$edad, "label", "Age, (Years)")
setattr(base$sexo, "label", "Male Sex, (%)")
setattr(base$seguimiento_anos, "label", "Follow-up, (Years)")
setattr(base$dm, "label", "Diabetes, (%)")
setattr(base$hta, "label", "Hypertension, (%)")
setattr(base$dislipidemia_neav, "label", "Dyslipidemia, (%)")
setattr(base$edad.primer.eca.neav, "label", "Age at First CVD, (%)")
setattr(base$tabaquismo_cat, "label", "Smoking, (%)")
setattr(base$promsisto, "label", "Systolic BP, (mmHg)")
setattr(base$promdiasto, "label", "Diastolic BP, (mmHg)")
setattr(base$uso_antihertensivo, "label", "Antihypertensive use, (%)")
setattr(base$Estatinas_REC, "label", "Statin use, (%)")
setattr(base$glucosa, "label", "Glucose, (mg/dl)")
setattr(base$tg, "label", "Tryglycerides, (mg/dl)")
setattr(base$ct, "label", "Total Cholesterol, (mg/dl)")
setattr(base$cldl, "label", "LDL-C, (mg/dl)")
setattr(base$chdl, "label", "HDL-C, (mg/dl)")
setattr(base$cno_hdl, "label", "Non-HDL-C, (mg/dl)")

setattr(base$TFG_CKD_EPI, "label", "eGFR-CKD-EPI, (ml/min/1.73 m²)")
setattr(base$homa, "label", "HOMA-IR")
setattr(base$rita, "label", "ADIPO-IR")
setattr(base$adiponectina, "label", "Adiponectin, (μg/mL)")
setattr(base$pcr, "label", "C-reactive protein, (μg/mL)")
setattr(base$imc, "label", "BMI, (Kg/m²)")
setattr(base$circunfer_cintura, "label", "WC, (cm)")
setattr(base$gav, "label", "VAT, (cm²)")
setattr(base$gas, "label", "SAT, (cm²)")
setattr(base$gat, "label", "TAT, (cm²)")
setattr(base$apo_b, "label", "Apo-B, (mg/dl)")

setattr(base$agatston, "label", "Agaston Score, (HU)")
setattr(base$morbimortalidad_mace, "label", "MACE, (%)")

setattr(base$mace_no_fatal, "label", "Non-Fatal MACE, (%)")
setattr(base$mace_fatal, "label", "Fatal MACE, (%)")

setattr(base$aspirina, "label", "Aspirin use, (%)")
setattr(base$anticoagulantes_orales, "label", "Anticoagulants use, (%)")

setattr(base$gluc_100, "label", "Glucose ≥100 mg/dl, (%)")
setattr(base$tg_cat, "label", "Tryglicerides ≥150 mg/dl, (%)")
setattr(base$cldl_cat, "label", "LDL-C ≥70 mg/dl, (%)")
setattr(base$nhdl_cat, "label", "Non-HDL-C ≥130 mg/dl, (%)")
setattr(base$chdl_cat, "label", "Low-HDL-C, (%)")
setattr(base$apob_p75_neav, "label", "Apo-B ≥p75, (%)")
setattr(base$pcr_cat, "label", "C-reactive protein ≥3.0 mg/dl, (%)")
setattr(base$homa_p75_neav, "label", "HOMA-IR ≥p75, (%)")
setattr(base$rita_p25_neav, "label", "ADIPO-IR ≥p75, (%)")
setattr(base$adipo_p25_neav, "label", "ADIPO-IR ≤p25, (%)")
setattr(base$gasp75, "label", "SAT ≥p75, (%)")
setattr(base$gatp75, "label", "TAT ≥p75, (%)")
setattr(base$BMI_CAT, "label", "BMI Categories, (%)")


#####Analysis: Descriptive characteristics of the study sample (Table 1)######

base %>% 
  dplyr::select(gav_terciles_neav,
                edad,sexo,seguimiento_anos,edad_eac,
                imc,BMI_CAT,
                tabaquismo_cat,dm,hta,
                Estatinas_REC,anticoagulantes_orales,aspirina,
                glucosa,gluc_100,
                tg,tg_cat,
                ct,
                cldl,cldl_cat,
                chdl,chdl_cat,
                cno_hdl,nhdl_cat,
                apo_b,apob_p75_neav,
                pcr,pcr_cat,
                homa,homa_p75_neav,
                rita,rita_p25_neav,
                adiponectina,adipo_p25_neav,
                gas,gasp75,
                gat,gatp75,
                gav,
                morbimortalidad_mace,mace_no_fatal,mace_fatal)%>%
  tbl_summary(by = gav_terciles_neav,
              missing = "ifany")%>%
  bold_labels()%>%
  add_overall()%>%
  add_p()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Table_1_TERCILES.docx")


#####Analysis: Incidence of MACE Events by VAT Tertiles (Figure 2A)####

#Function to estimate mortality rate
get.death.rate <- function(x){
  a <- x %>% summarise(cases=n(), time=sum(seguimiento_anos,na.rm = T))
  time <- sum(a$time); cases <- a$cases[2]
  lambda<-(cases/time)*1000; lambda.se<-(sqrt(cases/time^2))*1000
  LI<-lambda-(qnorm(.975)*lambda.se); LS<-lambda+(qnorm(.975)*lambda.se)
  cases.L<-LI*time/1000; cases.U<-LS*time/1000
  data.frame(lambda, cases, cases.L, cases.U, time)} #Crude death rates
un.ag.st <- function(x){
  a <- data.frame(); for(i in 1:length(table(base$age_group))){
    b <- x %>% dplyr::filter(age_group==i) %>% get.death.rate()
    a <- rbind(a, b)}; c <- a %>% cbind(
      "age"=base$age_group %>% levels); c} #Uniformly age and sex standardized
un.ag.st.CI <- function(x){
  x %>% na.omit %>% mutate("w"=1/length(age)) %>%
    mutate("var"=(((w)^2*lambda)/time)) %>% summarise(
      lambda=mean(lambda), var=sum(var), lambda.se=sqrt(sum(var)),
      cases=sum(cases), cases.L=sum(cases.L), cases.U=sum(cases.U)) %>% 
    mutate("Lower"=lambda+(cases.L-cases)*sqrt(var/cases),
           "Upper"=lambda+(cases.U-cases)*sqrt(var/cases))} #Standardized confidence intervals


#Overall MACE
rates1.M <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group) %>% dplyr::filter(seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
rates1.F <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group) %>% dplyr::filter(seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
sup.tab.1.0 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

#Lower-Tertile
rates1.1.M <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group,gav_terciles_neav_2) %>% dplyr::filter(gav_terciles_neav_2==1,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
rates1.1.F <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group,gav_terciles_neav_2) %>% dplyr::filter(gav_terciles_neav_2==1, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
rates1A <- rbind(rates1.1.M, rates1.1.F) %>% un.ag.st.CI

#Middle-Tertile
rates1.2.M <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group,gav_terciles_neav_2) %>% dplyr::filter(gav_terciles_neav_2==2,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
rates1.2.F <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group,gav_terciles_neav_2) %>% dplyr::filter(gav_terciles_neav_2==2, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
rates1B <- rbind(rates1.2.M, rates1.2.F) %>% un.ag.st.CI

#Upper-Tertile
rates1.3.M <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group,gav_terciles_neav_2) %>% dplyr::filter(gav_terciles_neav_2==3,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
rates1.3.F <- base %>% select(seguimiento_anos,morbimortalidad_mace_2,sexo,age_group,gav_terciles_neav_2) %>% dplyr::filter(gav_terciles_neav_2==3, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(morbimortalidad_mace_2) %>% un.ag.st
rates1C <- rbind(rates1.3.M, rates1.3.F) %>% un.ag.st.CI


Figure1A.df<-round(rbind(rates1A,rates1B,rates1C),1)
Figure1A.df$group<-c("Lower\n(<140 cm²)",
                     "Middle\n(140-194 cm²)", 
                     "Upper\n(>194 cm²)")


Figure1A<-Figure1A.df %>%
  ggplot(aes(x = group, y = lambda, fill = group)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = .2, position = position_dodge(.9)) +
  theme_classic() +
  xlab("VAT-Tertiles") +
  ylab("Events per 1,000 Persons-Years\n(Age and Sex Standarized)") +
  scale_fill_manual(values = c("gray80", "gray50", "gray10")) +
  geom_text(
    aes(label = paste0(sprintf("%.1f", lambda), "\n", "(", sprintf("%.1f", Lower), "-", sprintf("%.1f", Upper), ")"), y = lambda + 8),
    position = position_dodge(0.9),
    vjust = .5,
    size = 5
  ) +
  scale_y_continuous(limits = c(-.1, 51)) +
  labs(title = "Incidence Rate of MACE", fill = "") +
  annotate("text", x = 2, y = 51, size = 5, fontface = 2, label = paste0("Overall Rate: ", sprintf("%.1f", sup.tab.1.0$lambda), " (95% CI: ", sprintf("%.1f", sup.tab.1.0$Lower), "-", sprintf("%.1f", sup.tab.1.0$Upper), ")"))+
  theme(legend.position = "none",
        plot.title = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


#####Analysis: Kaplan-Meier Survival Curves (Figure 2B)######

#Kaplan Meier Curve Analisis
mod1_km<-survfit(Surv(seguimiento_anos, morbimortalidad_mace_2) ~factor(gav_terciles_neav), data = base)
KM_fig1<-ggsurvplot(mod1_km, data = base, size = 1,palette =c("gray80","gray50","gray10"),conf.int = T,
                    risk.table = T,
                    fun = "event",
                    ggtheme = theme_classic(),
                    xlab="Time of Follow-Up, (Years)",
                    ylab="Cumulative Incidence, (%)",
                    legend.title = "VAT-Tertiles",
                    legend.labs = c("Lower\n(<140 cm\u00b2)",
                                    "Middle\n(140-194 cm\u00b2)",
                                    "Upper\n(>194 cm\u00b2)"),
                    ylim= c(0,1.0),
                    xlim= c(0,7),
                    break.y.by= c(0.1),
                    break.x.by= c(1),
                    pval = TRUE, 
                    pval.method = TRUE,
                    log.rank.weights = "1", 
                    pval.method.size = 3,
                    pval.coord = c(1, 0.10),
                    pval.method.coord = c(1, 0.14))

Figure1B <-KM_fig1 + theme_survminer(base_size = 12,
                                     base_family = "Arial",
                                     font.x = c(12, "plain" ), 
                                     font.y = c(12, "plain"),
                                     font.main = c(14, "plain"),
                                     font.caption = c(8, "plain"), 
                                     font.legend = c(8, "plain"),
                                     font.tickslab = c(8, "plain"))



Figure1A_TOP<-Figure1B$plot+
  theme_survminer(base_size = 9,base_family = "Arial")+
  ggbreak::scale_y_break(breaks = c(0.60,0.90),scales = 0.2,space = 0.50)+
  
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position = "top",
        axis.text.y.right =element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())+
  labs(title="Cumulative Incidence of MACE")

Figure1A_DOWN<-Figure1B$table

Figure1B<-Figure1A_TOP / Figure1A_DOWN + 
  plot_layout(heights =  c(2, 0.60))

Figure1A<-Figure1A+plot_layout(heights =  c(2))

Figure1<-ggarrange(Figure1A,Figure1B,ncol = 2,nrow = 1,labels = c("A","B"))

ggsave(Figure1,
       filename = "Figure2.png", 
       bg = "white",
       width = 35, 
       height = 18.5,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)

#####Analysis: Cox-Models to Evaluate the Association with VAT and MACE (Figure 3)#####

##Cox Proportional Hazard Regression Models
#Univariate Association
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo), data=base_gea_lexit);summary(m0)
cox.zph(m0); summary(m0); BIC(m0)

#Model 1
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit);summary(m1)
cox.zph(m1); summary(m1); BIC(m1); car::vif(m1)

#Model 2
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit);summary(m2)
cox.zph(m2); summary(m2); BIC(m2)

#Model 3
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit);summary(m3)
cox.zph(m3); summary(m3); BIC(m3); car::vif(m3)

#Forest Plot
#Standard Units
m0<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav, data=base_gea_lexit))
m1<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit))
m2<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
                    BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
                    gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit))
m3<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
                    BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
                    gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
                    homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit))

tid.gav.units<-tibble(subgroup=c("Model 0\n(Unadjusted)","Model 1\n(Clinical)","Model 2\n(Biochemical)","Model 3\n(Adipose Tissue Markers)"),
                     name=c("+1 cm\u00b2","+1 cm\u00b2","+1 cm\u00b2","+1 cm\u00b2"),
                     beta=c(m0$coefficients[1,1],
                            m1$coefficients[1,1],
                            m2$coefficients[1,1],
                            m3$coefficients[1,1]),
                     se=c(m0$coefficients[1,3],
                          m1$coefficients[1,3],
                          m2$coefficients[1,3],
                          m3$coefficients[1,3]),
                     pvalue=c(m0$coefficients[1,5],
                              m1$coefficients[1,5],
                              m2$coefficients[1,5],
                              m3$coefficients[1,5]))

#Scale Units
m0<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav), data=base_gea_lexit))
m1<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit))
m2<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
                    BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
                    gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit))
m3<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
                    BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
                    gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
                    homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit))

tid.scale.units<-tibble(subgroup=c("Model 0\n(Unadjusted)","Model 1\n(Clinical)","Model 2\n(Biochemical)","Model 3\n(Adipose Tissue Markers)"),
                      name=c("+1 Standard Deviation\n(66.5 cm\u00b2)",
                             "+1 Standard Deviation\n(66.5 cm\u00b2)",
                             "+1 Standard Deviation\n(66.5 cm\u00b2)",
                             "+1 Standard Deviation\n(66.5 cm\u00b2)"),
                      beta=c(m0$coefficients[1,1],
                             m1$coefficients[1,1],
                             m2$coefficients[1,1],
                             m3$coefficients[1,1]),
                      se=c(m0$coefficients[1,3],
                           m1$coefficients[1,3],
                           m2$coefficients[1,3],
                           m3$coefficients[1,3]),
                      pvalue=c(m0$coefficients[1,5],
                               m1$coefficients[1,5],
                               m2$coefficients[1,5],
                               m3$coefficients[1,5]))


#Tertiles
m0<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav), data=base_gea_lexit))
m1<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit))
m2<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
                    BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
                    gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit))
m3<-summary(coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
                    BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
                    gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
                    homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit))

tid.t2.units<-tibble(subgroup=c("Model 4\n(Adipose Tissue Markers)",
                                "Model 3\n(Biochemical)",
                                "Model 2\n(Clinical)",
                                "Model 1\n(Age-at-risk and Sex)"),
                        name=c("Middle-Tertile\n(140-194 cm\u00b2)",
                               "Middle-Tertile\n(140-194 cm\u00b2)",
                               "Middle-Tertile\n(140-194 cm\u00b2)",
                               "Middle-Tertile\n(140-194 cm\u00b2)"),
                        beta=c(m3$coefficients[1,1],
                               m2$coefficients[1,1],
                               m1$coefficients[1,1],
                               m0$coefficients[1,1]),
                        se=c(m3$coefficients[1,3],
                             m2$coefficients[1,3],
                             m1$coefficients[1,3],
                             m0$coefficients[1,3]),
                        pvalue=c(m3$coefficients[1,5],
                                 m2$coefficients[1,5],
                                 m1$coefficients[1,5],
                                 m0$coefficients[1,5]))

tid.t3.units<-tibble(subgroup=c("Model 4\n(Adipose Tissue Markers)",
                                "Model 3\n(Biochemical)",
                                "Model 2\n(Clinical)",
                                "Model 1\n(Age-at-risk and Sex)"),
                     name=c("Upper-Tertile\n(>194 cm\u00b2)",
                            "Upper-Tertile\n(>194 cm\u00b2)",
                            "Upper-Tertile\n(>194 cm\u00b2)",
                            "Upper-Tertile\n(>194 cm\u00b2)"),
                     beta=c(m3$coefficients[2,1],
                            m2$coefficients[2,1],
                            m1$coefficients[2,1],
                            m0$coefficients[2,1]),
                     se=c(m3$coefficients[2,3],
                          m2$coefficients[2,3],
                          m1$coefficients[2,3],
                          m0$coefficients[2,3]),
                     pvalue=c(m3$coefficients[2,5],
                              m2$coefficients[2,5],
                              m1$coefficients[2,5],
                              m0$coefficients[2,5]))


tid.t1.units<-tibble(subgroup=c("Model 4\n(Adipose Tissue Markers)",
                                "Model 3\n(Biochemical)",
                                "Model 2\n(Clinical)",
                                "Model 1\n(Age-at-risk and Sex)"),
       name=c("Lower-Tertile\n(<194 cm\u00b2)",
              "Lower-Tertile\n(<194 cm\u00b2)",
              "Lower-Tertile\n(<194 cm\u00b2)",
              "Lower-Tertile\n(<194 cm\u00b2)"),
       beta=c(rep(0,4)),
       se=c(rep(0,4)),
       pvalue=c(rep(1,4)))


tid<-rbind(tid.t3.units,tid.t2.units,tid.t1.units)
tid<-tid %>%  mutate(subgroup=ordered(subgroup, 
                                 levels=c("Model 4\n(Adipose Tissue Markers)",
                                          "Model 3\n(Biochemical)",
                                          "Model 2\n(Clinical)",
                                          "Model 1\n(Age-at-risk and Sex)")))
tid$name<-factor(tid$name,labels = c("Lower\n(<140 cm\u00b2)\n[Reference]","Middle\n(140-194 cm\u00b2)","Upper\n(<194 cm\u00b2)"))

Figure2A<-ggforestplot::forestplot(
  df = tid,
  estimate = beta,
  shape = subgroup,
  logodds = TRUE,
  xlab = "Hazard Ratio\n(HR, 95% CI)",
  ylab = "VAT-Tertile")+
  labs(subtitle  = "Cox Proportional Hazard Regression Models",
       title = "B")+
  theme_classic() +
  theme(legend.position = "top")+
  labs(shape="")

#Plot Restricted Cubic Spline
base.fig1A<-base_gea_lexit%>%dplyr::select(time_at_entry, time_at_exit, exit_status,morbimortalidad_mace_2,gav)
base.fig1A<-as.data.frame(base.fig1A)
dd <- datadist(base.fig1A) ; options(datadist = "dd")
cph <- cph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data = base.fig1A, x = TRUE, y = TRUE)
pdata <- Predict(cph, gav, fun = exp)

Figure2B<-ggplot(data = pdata) +
  geom_hline(aes(yintercept = 1), linetype = 1) +
  geom_vline(aes(xintercept = 176), linetype = 2)+
  theme_classic() +
  scale_x_continuous(breaks = seq(from = 30, to = 500, by = 20))+
  scale_y_continuous(breaks = seq(from = 0, to = 8, by = 1))+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
    labs(x = "Visceral Adiposse Tissue (cm\u00b2)",
         y = "Hazard Ratio (95% CI)",
         title = "A",
         subtitle = "Simulation-Based Plot")+
  annotate(geom="text", x=100, y=6, label="Linear Association\np value <0.0001")

Figure3<-Figure2B+Figure2A

ggsave(Figure3,
       filename = "Figure3.png", 
       bg = "white",
       width = 40, 
       height = 15,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)

#####Analysis: Tressholds Evaluation (Text)######

m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav, data=base_gea_lexit);summary(m0)

#Muerte All Vascular
#Overall Dataset
res.cut.all.1 <- survminer::surv_cutpoint(base, time = "seguimiento_anos", event = "morbimortalidad_mace_2",variables = c("gav"))
summary(res.cut.all.1)

#Lexit Dataset
res.cut.all.1 <- survminer::surv_cutpoint(base_gea_lexit, time = "time_lexit", event = "exit_status",variables = c("gav"))
summary(res.cut.all.1)


#####Analysis: Interaction with VAT and Cardiometabolic Risk Factors (Figure 4)#####

#Sex 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(sexo==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(sexo==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*sexo, data=base_gea_lexit)
summ1a <- summary(c1); summ1b <- summary(c2);summ1c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(sexo) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.MALE <- deaths$Deaths[deaths$sexo==1]
d.FEMALE <- deaths$Deaths[deaths$sexo==0]

#Age 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(edad_65==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(edad_65==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*edad_65, data=base_gea_lexit)
summ2a <- summary(c1); summ2b <- summary(c2);summ2c <- summary(c3)

## Deaths by age definition ##
deaths <- base %>% group_by(edad_65) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.AGE.65 <- deaths$Deaths[deaths$edad_65==1]
d.AGE.non65 <- deaths$Deaths[deaths$edad_65==0]

#Age at CVD 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(AGE_CVD_CAT==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(AGE_CVD_CAT==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*AGE_CVD_CAT, data=base_gea_lexit)
summ3a <- summary(c1); summ3b <- summary(c2);summ3c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(AGE_CVD_CAT) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.CVD.EDAD.35<- deaths$Deaths[deaths$AGE_CVD_CAT==1]
d.CVD.EDAD.35.65 <- deaths$Deaths[deaths$AGE_CVD_CAT==0]

#BMI Categories
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(BMI_CAT==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(BMI_CAT==2))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(BMI_CAT==3))
c4<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*BMI_CAT, data=base_gea_lexit)
summ4a <- summary(c1); summ4b <- summary(c2);summ4c <- summary(c3);summ4d <- summary(c4)

## Deaths by sex definition ##
deaths <- base %>% group_by(BMI_CAT) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NORMO<- deaths$Deaths[deaths$BMI_CAT==1]
d.SOBRE <- deaths$Deaths[deaths$BMI_CAT==2]
d.OBES <- deaths$Deaths[deaths$BMI_CAT==3]

#Smoking Status
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(tabaquismo_cat=="Never"))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(tabaquismo_cat %in% c("Previous","Current")))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*(tabaquismo_cat=="Never"), data=base_gea_lexit)
summ5a <- summary(c1); summ5b <- summary(c2);summ5c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(tabaquismo_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NEVER <- deaths$Deaths[deaths$tabaquismo_cat=="Never"]
d.PRE.CURR <- sum(deaths$Deaths[deaths$tabaquismo_cat=="Previous" | deaths$tabaquismo_cat=="Current"])

#Diabetes
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(diabetes==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(diabetes==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*diabetes, data=base_gea_lexit)
summ6a <- summary(c1); summ6b <- summary(c2);summ6c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(diabetes) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.DIABETES <- deaths$Deaths[deaths$diabetes==1]
d.NON.DIABETES <- deaths$Deaths[deaths$diabetes==0]

#Hypertension Status
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(hta==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(hta==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*hta, data=base_gea_lexit)
summ7a <- summary(c1); summ7b <- summary(c2);summ7c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(hta) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.HTA <- deaths$Deaths[deaths$hta==1]
d.NON.HTA <- deaths$Deaths[deaths$hta==0]

#Statins Use
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(Estatinas_REC==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(Estatinas_REC==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*Estatinas_REC, data=base_gea_lexit)
summ8a <- summary(c1); summ8b <- summary(c2);summ8c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(Estatinas_REC) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ESTATINA <- deaths$Deaths[deaths$Estatinas_REC==1]
d.NON.ESTATINA <- deaths$Deaths[deaths$Estatinas_REC==0]

#Anticoagulation
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(anticoagulantes_orales==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(anticoagulantes_orales==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*anticoagulantes_orales, data=base_gea_lexit)
summ9a <- summary(c1); summ9b <- summary(c2);summ9c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(anticoagulantes_orales) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ANTICOAG <- deaths$Deaths[deaths$anticoagulantes_orales==1]
d.NON.ANTICOAG <- deaths$Deaths[deaths$anticoagulantes_orales==0]

#Aspirin
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(aspirina==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(aspirina==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*aspirina, data=base_gea_lexit)
summ10a <- summary(c1); summ10b <- summary(c2);summ10c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(aspirina) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ASA <- deaths$Deaths[deaths$aspirina==1]
d.NON.ASA <- deaths$Deaths[deaths$aspirina==0]

#Glucose
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(gluc_100==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(gluc_100==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*gluc_100, data=base_gea_lexit)
summ11a <- summary(c1); summ11b <- summary(c2);summ11c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(gluc_100) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GLUC <- deaths$Deaths[deaths$gluc_100==1]
d.NON.GLUC <- deaths$Deaths[deaths$gluc_100==0]

#Tryglycerides
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(tg_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(tg_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*tg_cat, data=base_gea_lexit)
summ12a <- summary(c1); summ12b <- summary(c2);summ12c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(tg_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.TG <- deaths$Deaths[deaths$tg_cat==1]
d.NON.TG <- deaths$Deaths[deaths$tg_cat==0]


#LDL-C
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(cldl_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(cldl_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*cldl_cat, data=base_gea_lexit)
summ13a <- summary(c1); summ13b <- summary(c2);summ13c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(cldl_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.LDL <- deaths$Deaths[deaths$cldl_cat==1]
d.NON.LDL <- deaths$Deaths[deaths$cldl_cat==0]


#NON-HDL
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(nhdl_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(nhdl_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*nhdl_cat, data=base_gea_lexit)
summ14a <- summary(c1); summ14b <- summary(c2);summ14c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(nhdl_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NHDL <- deaths$Deaths[deaths$nhdl_cat==1]
d.NON.NHDL <- deaths$Deaths[deaths$nhdl_cat==0]

#APO-B p75 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(apob_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(apob_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*apob_p75_neav, data=base_gea_lexit)
summ15a <- summary(c1); summ15b <- summary(c2); summ15c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(apob_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.APOB <- deaths$Deaths[deaths$apob_p75_neav==1]
d.NON.APOB <- deaths$Deaths[deaths$apob_p75_neav==0]


#CRP
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(pcr_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(pcr_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*pcr_cat, data=base_gea_lexit)
summ16a <- summary(c1); summ16b <- summary(c2);summ16c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(pcr_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.PCR <- deaths$Deaths[deaths$pcr_cat==1]
d.NON.PCR <- deaths$Deaths[deaths$pcr_cat==0]


#HOMA-2IR
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(homa_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(homa_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*homa_p75_neav, data=base_gea_lexit)
summ17a <- summary(c1); summ17b <- summary(c2);summ17c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(homa_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.HOMA <- deaths$Deaths[deaths$homa_p75_neav==1]
d.NON.HOMA <- deaths$Deaths[deaths$homa_p75_neav==0]

#RITA
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(rita_p25_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(rita_p25_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*rita_p25_neav, data=base_gea_lexit)
summ18a <- summary(c1); summ18b <- summary(c2); summ18c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(rita_p25_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.RITA <- deaths$Deaths[deaths$rita_p25_neav==1]
d.NON.RITA <- deaths$Deaths[deaths$rita_p25_neav==0]


#Adiponectin
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(adipo_p25_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(adipo_p25_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*adipo_p25_neav, data=base_gea_lexit)
summ19a <- summary(c1); summ19b <- summary(c2); summ19c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(adipo_p25_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ADIPONEC <- deaths$Deaths[deaths$adipo_p25_neav==1]
d.NON.ADIPONEC <- deaths$Deaths[deaths$adipo_p25_neav==0]


#GAT p>75
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(GAT_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(GAT_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*GAT_p75_neav, data=base_gea_lexit)
summ20a <- summary(c1); summ20b <- summary(c2); summ20c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(GAT_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GAT <- deaths$Deaths[deaths$GAT_p75_neav==1]
d.NON.GAT <- deaths$Deaths[deaths$GAT_p75_neav==0]

#GAS p>75
#Stratification

c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(GAS_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit%>%dplyr::filter(GAS_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*GAS_p75_neav, data=base_gea_lexit)
summ21a <- summary(c1); summ21b <- summary(c2); summ21c <- summary(c3)


## Deaths by sex definition ##
deaths <- base %>% group_by(GAS_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GAS <- deaths$Deaths[deaths$GAS_p75_neav==1]
d.NON.GAS <- deaths$Deaths[deaths$GAS_p75_neav==0]

#Dataset Merging
F.TAB_WWI <- rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4]),c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4]), c(summ1c$coefficients[3,2], summ1c$conf.int[3,3:4]),
  c(summ2a$coefficients[1,2], summ2a$conf.int[1,3:4]),c(summ2b$coefficients[1,2], summ2b$conf.int[1,3:4]),c(summ2c$coefficients[3,2], summ2c$conf.int[3,3:4]),
  c(summ3a$coefficients[1,2], summ3a$conf.int[1,3:4]),c(summ3b$coefficients[1,2], summ3b$conf.int[1,3:4]),c(summ3c$coefficients[3,2], summ3c$conf.int[3,3:4]),
  c(summ4a$coefficients[1,2], summ4a$conf.int[1,3:4]),c(summ4b$coefficients[1,2], summ4b$conf.int[1,3:4]),c(summ4c$coefficients[1,2], summ4c$conf.int[1,3:4]),c(summ4d$coefficients[3,2], summ4d$conf.int[3,3:4]),
  c(summ5a$coefficients[1,2], summ5a$conf.int[1,3:4]),c(summ5b$coefficients[1,2], summ5b$conf.int[1,3:4]),c(summ5c$coefficients[3,2], summ5c$conf.int[3,3:4]),
  c(summ6a$coefficients[1,2], summ6a$conf.int[1,3:4]),c(summ6b$coefficients[1,2], summ6b$conf.int[1,3:4]),c(summ6c$coefficients[3,2], summ6c$conf.int[3,3:4]),
  c(summ7a$coefficients[1,2], summ7a$conf.int[1,3:4]),c(summ7b$coefficients[1,2], summ7b$conf.int[1,3:4]),c(summ7c$coefficients[3,2], summ7c$conf.int[3,3:4]),
  c(summ8a$coefficients[1,2], summ8a$conf.int[1,3:4]),c(summ8b$coefficients[1,2], summ8b$conf.int[1,3:4]),c(summ8c$coefficients[3,2], summ8c$conf.int[3,3:4]),
  c(summ9a$coefficients[1,2], summ9a$conf.int[1,3:4]),c(summ9b$coefficients[1,2], summ9b$conf.int[1,3:4]),c(summ9c$coefficients[3,2], summ9c$conf.int[3,3:4]),
  c(summ10a$coefficients[1,2], summ10a$conf.int[1,3:4]),c(summ10b$coefficients[1,2], summ10b$conf.int[1,3:4]),c(summ10c$coefficients[3,2], summ10c$conf.int[3,3:4]),
  c(summ11a$coefficients[1,2], summ11a$conf.int[1,3:4]),c(summ11b$coefficients[1,2], summ11b$conf.int[1,3:4]),c(summ11c$coefficients[3,2], summ11c$conf.int[3,3:4]),
  c(summ12a$coefficients[1,2], summ12a$conf.int[1,3:4]),c(summ12b$coefficients[1,2], summ12b$conf.int[1,3:4]),c(summ12c$coefficients[3,2], summ12c$conf.int[3,3:4]),
  c(summ13a$coefficients[1,2], summ13a$conf.int[1,3:4]),c(summ13b$coefficients[1,2], summ13b$conf.int[1,3:4]),c(summ13c$coefficients[3,2], summ13c$conf.int[3,3:4]),
  c(summ14a$coefficients[1,2], summ14a$conf.int[1,3:4]),c(summ14b$coefficients[1,2], summ14b$conf.int[1,3:4]),c(summ14c$coefficients[3,2], summ14c$conf.int[3,3:4]),
  c(summ15a$coefficients[1,2], summ15a$conf.int[1,3:4]),c(summ15b$coefficients[1,2], summ15b$conf.int[1,3:4]),c(summ15c$coefficients[3,2], summ15c$conf.int[3,3:4]),
  c(summ16a$coefficients[1,2], summ16a$conf.int[1,3:4]),c(summ16b$coefficients[1,2], summ16b$conf.int[1,3:4]),c(summ16c$coefficients[3,2], summ16c$conf.int[3,3:4]),
  c(summ17a$coefficients[1,2], summ17a$conf.int[1,3:4]),c(summ17b$coefficients[1,2], summ17b$conf.int[1,3:4]),c(summ17c$coefficients[3,2], summ17c$conf.int[3,3:4]),
  c(summ18a$coefficients[1,2], summ18a$conf.int[1,3:4]),c(summ18b$coefficients[1,2], summ18b$conf.int[1,3:4]),c(summ18c$coefficients[3,2], summ18c$conf.int[3,3:4]),
  c(summ19a$coefficients[1,2], summ19a$conf.int[1,3:4]),c(summ19b$coefficients[1,2], summ19b$conf.int[1,3:4]),c(summ19c$coefficients[3,2], summ19c$conf.int[3,3:4]),
  c(summ20a$coefficients[1,2], summ20a$conf.int[1,3:4]),c(summ20b$coefficients[1,2], summ20b$conf.int[1,3:4]),c(summ20c$coefficients[3,2], summ20c$conf.int[3,3:4]),
  c(summ21a$coefficients[1,2], summ21a$conf.int[1,3:4]),c(summ21b$coefficients[1,2], summ21b$conf.int[1,3:4]),c(summ21c$coefficients[3,2], summ21c$conf.int[3,3:4])) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper")) %>% 
  mutate("Cause"=c(rep(1,3),
                   rep(2,3),
                   rep(3,3),
                   rep(4,4),
                   rep(5,3),
                   rep(6,3),
                   rep(7,3),
                   rep(8,3),
                   rep(9,3),
                   rep(10,3),
                   rep(11,3),
                   rep(12,3),
                   rep(13,3),
                   rep(14,3),
                   rep(15,3),
                   rep(16,3),
                   rep(17,3),
                   rep(18,3),
                   rep(19,3),
                   rep(20,3),
                   rep(21,3)), 
         "Strat"=rep(1:64,1)) %>% 
  mutate("Deaths"=c(d.MALE,d.FEMALE,"", #1
                    d.AGE.65,d.AGE.non65,"", #2
                    d.CVD.EDAD.35,d.CVD.EDAD.35.65,"", #3
                    d.NORMO,d.SOBRE,d.OBES,"", #4
                    d.NEVER,d.PRE.CURR,"", #5
                    d.DIABETES,d.NON.DIABETES,"", #6
                    d.HTA,d.NON.HTA,"", #7
                    d.ESTATINA,d.NON.ESTATINA,"", #8
                    d.ANTICOAG,d.NON.ANTICOAG,"", #9
                    d.ASA,d.NON.ASA,"", #10
                    d.LDL,d.NON.LDL,"", #11
                    d.GLUC,d.NON.GLUC,"", #12
                    d.TG,d.NON.TG,"", #13
                    d.NHDL,d.NON.NHDL,"", #14
                    d.APOB,d.NON.APOB,"", #15
                    d.PCR,d.NON.PCR,"", #16
                    d.HOMA,d.NON.HOMA,"", #17
                    d.RITA,d.NON.RITA,"", #18
                    d.ADIPONEC,d.NON.ADIPONEC,"", #19
                    d.GAT,d.NON.GAT,"", #20
                    d.GAS,d.NON.GAS,"")) %>% #21
  mutate(Strat=ordered(Strat, 1:64, c("Male","Female","Interaction", #1
                                     ">65 years","≤65 years","Interaction", #2
                                     "<35 years","35-65 years","Interaction", #3
                                     "Underweight","Overweight","Obesity","Interaction", #4
                                     "Never","Previous-Current","Interaction", #5
                                     "Diabetes","Without Diabetes","Interaction", #6
                                     "Hypertension","Without Hypertension","Interaction", #7
                                     "Use","Non-User","Interaction", #8
                                     "Use","Non-User","Interaction", #9
                                     "Use","Non-User","Interaction", #10
                                     "≥100 mg/dL","<100 mg/dL","Interaction", #11
                                     "<150 mg/dL","≥150 mg/dL","Interaction", #12
                                     "≤70 mg/dL",">70 mg/dL","Interaction", #13
                                     "<130 mg/dL","≥130 mg/dL","Interaction", #14
                                     "≥p75","<p75","Interaction", #15
                                     "≥3.0 g/L","<3.0 g/L","Interaction", #16
                                     "≥p75","<p75","Interaction", #17
                                     "≥p75","<p75","Interaction", #18
                                     "<p25","≥p25","Interaction", #19
                                     "≥p75 ","<p75 ","Interaction", #21
                                     "≥p75 ","<p75 ","Interaction"))) %>%  #22
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))


## Build forestplot ##
fnfp1 <- fpDrawNormalCI; fnfp2 <- fpDrawDiamondCI; gp1 <- gpar(fill="#6699CC"); gp2 <- gpar(fill="#994455"); 

#Label of Figure
d.labs <- c("Sex",NA,NA,NA, #1
            "Age",NA,NA,NA, #2
            "Age at CVD",NA,NA,NA, #3
            "BMI Categories",NA,NA,NA,NA, #4
            "Smoking Status",NA,NA,NA, #5
            "Diabetes",NA,NA,NA, #6
            "Arterial Hypertension",NA,NA,NA, #7
            "Statin",NA,NA,NA, #8
            "Anticoagulants Use",NA,NA,NA, #9
            "Aspirin Use",NA,NA,NA, #10
            "Glucose",NA,NA,NA, #11
            "Triglycerides",NA,NA,NA, #12
            "LDL-C",NA,NA,NA, #13
            "Non-HDL-C",NA,NA,NA, #14
            "APO-B",NA,NA,NA, #15
            "C-Reactive Protein",NA,NA,NA, #16
            "HOMA-IR",NA,NA,NA, #17
            "ADIPO-IR",NA,NA,NA, #18
            "Adiponectin",NA,NA,NA, #19
            "Total Adiposse Tissue",NA,NA,NA,#21
            "Subcutaneous Adiposse Tissue",NA,NA,NA)#22


p.inter.val<- round(rbind(NA,NA,NA,summ1c$coefficients[3,5],
                          NA,NA,NA,summ2c$coefficients[3,5],
                          NA,NA,NA,summ3c$coefficients[3,5],
                          NA,NA,NA,NA,summ4d$coefficients[3,5],
                          NA,NA,NA,summ5c$coefficients[3,5],
                          NA,NA,NA,summ6c$coefficients[3,5],
                          NA,NA,NA,summ7c$coefficients[3,5],
                          NA,NA,NA,summ8c$coefficients[3,5],
                          NA,NA,NA,summ9c$coefficients[3,5],
                          NA,NA,NA,summ10c$coefficients[3,5],
                          NA,NA,NA,summ11c$coefficients[3,5],
                          NA,NA,NA,summ12c$coefficients[3,5],
                          NA,NA,NA,summ13c$coefficients[3,5],
                          NA,NA,NA,summ14c$coefficients[3,5],
                          NA,NA,NA,summ15c$coefficients[3,5],
                          NA,NA,NA,summ16c$coefficients[3,5],
                          NA,NA,NA,summ17c$coefficients[3,5],
                          NA,NA,NA,summ18c$coefficients[3,5],
                          NA,NA,NA,summ19c$coefficients[3,5],
                          NA,NA,NA,summ20c$coefficients[3,5],
                          NA,NA,NA,summ21c$coefficients[3,5]),3)

rbind(NA, F.TAB_WWI[1:3,], #1
      NA, F.TAB_WWI[4:6,], #2
      NA, F.TAB_WWI[7:9,], #3
      NA, F.TAB_WWI[10:13,], #4
      NA, F.TAB_WWI[14:16,], #5
      NA, F.TAB_WWI[17:19,], #6
      NA, F.TAB_WWI[20:22,], #8 
      NA, F.TAB_WWI[23:25,], #9
      NA, F.TAB_WWI[26:28,], #10
      NA, F.TAB_WWI[29:31,], #11
      NA, F.TAB_WWI[32:34,], #12
      NA, F.TAB_WWI[35:37,], #13
      NA, F.TAB_WWI[38:40,], #14
      NA, F.TAB_WWI[41:43,], #15
      NA, F.TAB_WWI[44:46,], #16
      NA, F.TAB_WWI[47:49,], #17
      NA, F.TAB_WWI[50:52,], #18
      NA, F.TAB_WWI[53:55,], #19
      NA, F.TAB_WWI[56:58,], #20
      NA, F.TAB_WWI[59:61,], #21
      NA, F.TAB_WWI[62:64,]) %>% #21
  cbind(p.inter.val)%>%
  mutate(Cause=d.labs)%>%
  forestplot::forestplot(
    xticks = c(-0.015, 0.028),
    col = fpColors(lines = "black", zero = "black"),
    labeltext = c("Cause", "Strat", "Deaths", "HR","p.inter.val"), 
    graph.pos = 4, xlog = T)%>% 
  fp_add_header(Cause="Stratification", 
                Strat = c("Categories"), 
                Deaths = "MACE", 
                HR = c("HR (95% CI)"),
                p.inter.val = c("p for interaction")) %>%
  fp_set_style(txt_gp = fpTxtGp(label = gpar(fontfamily = "Arial"), 
                                ticks = gpar(cex = 1), 
                                xlab  = gpar(cex=1.20, fontface="bold.italic"),
                                summary = gpar(fontfamily = "Arial", fontface="bold.italic")), align="c", summary = "#7A316F") %>% 
  fp_set_zebra_style("gray90","white") -> forest1

png("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/VAT - MACE/Figure4.png",
    width=35, height=36, units = "cm", res = 300); forest1
dev.off(); print(plot(1))


#####Analysis: Incidence of Fatal and Non-Fatal MACE Events by VAT Tertiles and KM Curvers (Supplementary Figure 1A-C)#####

#NON-FATAL-MACE
#Overall Non-Fatal MACE
rates1.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group) %>% dplyr::filter(seguimiento_anos>=0) %>% filter(MACE_CAT %in% c(0,1)) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group) %>% dplyr::filter(seguimiento_anos>=0) %>% filter(MACE_CAT %in% c(0,1)) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
sup.tab.0_NF_MACE <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

#Lower-Tertile
rates1.1.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,1)) %>% dplyr::filter(gav_terciles_neav_2==1,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.1.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,1)) %>% dplyr::filter(gav_terciles_neav_2==1, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
rates1A <- rbind(rates1.1.M, rates1.1.F) %>% un.ag.st.CI

#Middle-Tertile
rates1.2.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,1)) %>% dplyr::filter(gav_terciles_neav_2==2,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.2.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,1)) %>% dplyr::filter(gav_terciles_neav_2==2, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
rates1B <- rbind(rates1.2.M, rates1.2.F) %>% un.ag.st.CI

#Upper-Tertile
rates1.3.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,1)) %>% dplyr::filter(gav_terciles_neav_2==3,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.3.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,1)) %>% dplyr::filter(gav_terciles_neav_2==3, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
rates1C <- rbind(rates1.3.M, rates1.3.F) %>% un.ag.st.CI

#Merging Datasets
Figure1A.df.1.1<-round(rbind(rates1A,rates1B,rates1C),1)
Figure1A.df.1.1$group<-c("Lower-Tertile\n(<140 cm²)",
                         "Middle-Tertile\n(140-194 cm²)", 
                         "Upper-Tertile\n(>194 cm²)")
#FATAL-MACE
#Overall Fatal MACE
rates1.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group) %>% dplyr::filter(seguimiento_anos>=0) %>% filter(MACE_CAT %in% c(0,2)) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group) %>% dplyr::filter(seguimiento_anos>=0) %>% filter(MACE_CAT %in% c(0,2)) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
sup.tab.0_F_MACE <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

#Lower-Tertile
rates1.1.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,2)) %>% dplyr::filter(gav_terciles_neav_2==1,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.1.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,2)) %>% dplyr::filter(gav_terciles_neav_2==1, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
rates1A <- rbind(rates1.1.M, rates1.1.F) %>% un.ag.st.CI

#Middle-Tertile
rates1.2.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,2)) %>% dplyr::filter(gav_terciles_neav_2==2,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.2.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,2)) %>% dplyr::filter(gav_terciles_neav_2==2, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
rates1B <- rbind(rates1.2.M, rates1.2.F) %>% un.ag.st.CI

#Upper-Tertile
rates1.3.M <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,2)) %>% dplyr::filter(gav_terciles_neav_2==3,seguimiento_anos>=0) %>% filter(sexo==1) %>% group_by(MACE_CAT) %>% un.ag.st
rates1.3.F <- base %>% select(seguimiento_anos,MACE_CAT,sexo,age_group,gav_terciles_neav_2) %>% filter(MACE_CAT %in% c(0,2)) %>% dplyr::filter(gav_terciles_neav_2==3, seguimiento_anos>=0) %>% filter(sexo==0) %>% group_by(MACE_CAT) %>% un.ag.st
rates1C <- rbind(rates1.3.M, rates1.3.F) %>% un.ag.st.CI

#Merging Datasets
Figure1A.df.1.2<-round(rbind(rates1A,rates1B,rates1C),1)
Figure1A.df.1.2$group<-c("Lower-Tertile\n(<140 cm²)",
                                  "Middle-Tertile\n(140-194 cm²)", 
                                  "Upper-Tertile\n(>194 cm²)")

#Figures
#Non-Fatal-MACE
Sup.Figure1A<-Figure1A.df.1.1 %>%
  ggplot( aes(x=group, y=lambda,fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2,
                position=position_dodge(.9)) +
  theme_classic()+
  xlab("")+
  ylab("Events per 1,000 Persons-Years\n(Age and Sex Standarized)")+
  scale_fill_manual(values = c("gray80","gray50","gray10"))+
  geom_text(
    aes(label = paste0(lambda,"","\n","(",Lower,"-",Upper,")"), y = lambda + 10),
    position = position_dodge(0.9),
    vjust = .5)+
  scale_y_continuous(limits = c(-.1,50))+
  labs(title = "Incidence Rate of Non-Fatal MACE",
       fill="VAT-Tertiles")+
  annotate("text",x = 1.5,y = 45,label = paste0("Overall Rate: ",sup.tab.0_NF_MACE$lambda,"",""," (",sup.tab.0_NF_MACE$Lower,"-",sup.tab.0_NF_MACE$Upper,")"))

#Fatal-MACE
Sup.Figure1C<-Figure1A.df.1.2 %>%
  ggplot( aes(x=group, y=lambda,fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2,
                position=position_dodge(.9)) +
  theme_classic()+
  xlab("")+
  ylab("Events per 1,000 Persons-Years\n(Age and Sex Standarized)")+
  scale_fill_manual(values = c("gray80","gray50","gray10"))+
  geom_text(
    aes(label = paste0(lambda,"","\n","(",Lower,"-",Upper,")"), y = lambda + 10),
    position = position_dodge(0.9),
    vjust = .5)+
  scale_y_continuous(limits = c(-.1,50))+
  labs(title = "Incidence Rate of Fatal MACE",
       fill="VAT-Tertiles")+
  annotate("text",x = 1.5,y = 45,label = paste0("Overall Rate: ",sup.tab.0_F_MACE$lambda,"",""," (",sup.tab.0_F_MACE$Lower,"-",sup.tab.0_F_MACE$Upper,")"))

#####Analysis: Kaplan-Meier Survival Curves For Fatal and Non-Fatal MACE (Supplementary Figure 1B-D)######

#Non Fatal MACE KM
mod1_km<-survfit(Surv(seguimiento_anos, mace_no_fatal) ~factor(gav_terciles_neav), data = base)
SF.KM_fig1<-ggsurvplot(mod1_km, data = base, size = 1,palette =c("gray80","gray50","gray10"),conf.int = T,
                    risk.table = T,
                    ggtheme = theme_classic(),
                    xlab="Time of Follow-Up, (Years)",
                    ylab="Survival Probability, (%)",
                    legend.labs = c("Lower-Tertile\n(<140 cm\u00b2)",
                                    "Middle-Tertile\n(140-194 cm\u00b2)",
                                    "Upper-Tertile\n(>194 cm\u00b2)"),
                    xlim = c(0,5),
                    ylim= c(0,1.0),
                    break.y.by= c(0.1),
                    break.x.by= c(1),
                    pval = TRUE, 
                    pval.method = TRUE,
                    log.rank.weights = "1", 
                    pval.method.size = 3,
                    pval.coord = c(1, 0.95),
                    pval.method.coord = c(1, 0.965))

Sup.Figure1B <-SF.KM_fig1 + theme_survminer(base_size = 8,
                                     base_family = "Arial",
                                     font.x = c(8, "plain" ), 
                                     font.y = c(8, "plain"),
                                     font.main = c(8, "plain"),
                                     font.caption = c(8, "plain"), 
                                     font.legend = c(8, "plain"),
                                     font.tickslab = c(8, "plain"))

Sup.Figure1B_TOP<-Sup.Figure1B$plot+
  ggbreak::scale_y_break(breaks = c(0.15,0.85),scales = 5,space = 0.45)+
  scale_x_continuous(limits = c(0,5))+
  theme_survminer(base_size = 9,base_family = "Arial")+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position = "top",
        axis.text.y.right =element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())+
  labs(title="Survival Probability for Non-Fatal MACE")

Sup.Figure1B_DOWN<-Sup.Figure1B$table

Sup.Figure1B<-Sup.Figure1B_TOP / Sup.Figure1B_DOWN + 
  plot_layout(heights =  c(2, 0.5))

#FATAL MACE
mod1_km<-survfit(Surv(seguimiento_anos, mace_fatal) ~factor(gav_terciles_neav), data = base)
SF.KM_fig2<-ggsurvplot(mod1_km, data = base, size = 1,palette =c("gray80","gray50","gray10"),conf.int = T,
                       risk.table = T,
                       ggtheme = theme_classic(),
                       xlab="Time of Follow-Up, (Years)",
                       ylab="Survival Probability, (%)",
                       legend.labs = c("Lower-Tertile\n(<140 cm\u00b2)",
                                       "Middle-Tertile\n(140-194 cm\u00b2)",
                                       "Upper-Tertile\n(>194 cm\u00b2)"),
                       xlim = c(0,5),
                       ylim= c(0,1.0),
                       break.y.by= c(0.1),
                       break.x.by= c(1),
                       pval = TRUE, 
                       pval.method = TRUE,
                       log.rank.weights = "1", 
                       pval.method.size = 3,
                       pval.coord = c(1, 0.95),
                       pval.method.coord = c(1, 0.965))

Sup.Figure1D <-SF.KM_fig2 + theme_survminer(base_size = 8,
                                            base_family = "Arial",
                                            font.x = c(8, "plain" ), 
                                            font.y = c(8, "plain"),
                                            font.main = c(8, "plain"),
                                            font.caption = c(8, "plain"), 
                                            font.legend = c(8, "plain"),
                                            font.tickslab = c(8, "plain"))

Sup.Figure1D_TOP<-Sup.Figure1D$plot+
  ggbreak::scale_y_break(breaks = c(0.15,0.85),scales = 5,space = 0.45)+
  scale_x_continuous(limits = c(0,5))+
  theme_survminer(base_size = 9,base_family = "Arial")+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position = "top",
        axis.text.y.right =element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())+
  labs(title="Survival Probability for Fatal MACE")

Sup.Figure1D_DOWN<-Sup.Figure1D$table

Sup.Figure1D<-Sup.Figure1D_TOP / Sup.Figure1D_DOWN + 
  plot_layout(heights =  c(2, 0.5))


Sup.Figure1<-ggarrange(Sup.Figure1A,Sup.Figure1B,
          Sup.Figure1C,Sup.Figure1D,
          ncol = 2,nrow = 2,labels = c(LETTERS[1:4]))

ggsave(Sup.Figure1,
       filename = "Supplementary_Figure_1.png", 
       bg = "white",
       width = 40, 
       height = 32,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)


#####Analysis: Interaction with VAT and Cardiometabolic Risk Factors for Fatal and Non-Fatal MACE (Supplementary Figure 4)#####

#Non_FATAL
#Sex 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(sexo==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(sexo==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*sexo, data=base_gea_lexit_NF_MACE)
summ1a <- summary(c1); summ1b <- summary(c2);summ1c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(sexo) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.MALE <- deaths$Deaths[deaths$sexo==1]
d.FEMALE <- deaths$Deaths[deaths$sexo==0]

#Age 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(edad_65==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(edad_65==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*edad_65, data=base_gea_lexit_NF_MACE)
summ2a <- summary(c1); summ2b <- summary(c2);summ2c <- summary(c3)

## Deaths by age definition ##
deaths <- base %>% group_by(edad_65) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.AGE.65 <- deaths$Deaths[deaths$edad_65==1]
d.AGE.non65 <- deaths$Deaths[deaths$edad_65==0]

#Age at CVD 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(AGE_CVD_CAT==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(AGE_CVD_CAT==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*AGE_CVD_CAT, data=base_gea_lexit_NF_MACE)
summ3a <- summary(c1); summ3b <- summary(c2);summ3c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(AGE_CVD_CAT) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.CVD.EDAD.35<- deaths$Deaths[deaths$AGE_CVD_CAT==1]
d.CVD.EDAD.35.65 <- deaths$Deaths[deaths$AGE_CVD_CAT==0]

#BMI Categories
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(BMI_CAT==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(BMI_CAT==2))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(BMI_CAT==3))
c4<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*BMI_CAT, data=base_gea_lexit_NF_MACE)
summ4a <- summary(c1); summ4b <- summary(c2);summ4c <- summary(c3);summ4d <- summary(c4)

## Deaths by sex definition ##
deaths <- base %>% group_by(BMI_CAT) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NORMO<- deaths$Deaths[deaths$BMI_CAT==1]
d.SOBRE <- deaths$Deaths[deaths$BMI_CAT==2]
d.OBES <- deaths$Deaths[deaths$BMI_CAT==3]

#Smoking Status
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(tabaquismo_cat=="Never"))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(tabaquismo_cat %in% c("Previous","Current")))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*(tabaquismo_cat=="Never"), data=base_gea_lexit_NF_MACE)
summ5a <- summary(c1); summ5b <- summary(c2);summ5c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(tabaquismo_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NEVER <- deaths$Deaths[deaths$tabaquismo_cat=="Never"]
d.PRE.CURR <- sum(deaths$Deaths[deaths$tabaquismo_cat=="Previous" | deaths$tabaquismo_cat=="Current"])

#Diabetes
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(diabetes==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(diabetes==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*diabetes, data=base_gea_lexit_NF_MACE)
summ6a <- summary(c1); summ6b <- summary(c2);summ6c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(diabetes) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.DIABETES <- deaths$Deaths[deaths$diabetes==1]
d.NON.DIABETES <- deaths$Deaths[deaths$diabetes==0]

#Hypertension Status
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(hta==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(hta==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*hta, data=base_gea_lexit_NF_MACE)
summ7a <- summary(c1); summ7b <- summary(c2);summ7c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(hta) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.HTA <- deaths$Deaths[deaths$hta==1]
d.NON.HTA <- deaths$Deaths[deaths$hta==0]

#Statins Use
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(Estatinas_REC==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(Estatinas_REC==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*Estatinas_REC, data=base_gea_lexit_NF_MACE)
summ8a <- summary(c1); summ8b <- summary(c2);summ8c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(Estatinas_REC) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ESTATINA <- deaths$Deaths[deaths$Estatinas_REC==1]
d.NON.ESTATINA <- deaths$Deaths[deaths$Estatinas_REC==0]

#Anticoagulation
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(anticoagulantes_orales==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(anticoagulantes_orales==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*anticoagulantes_orales, data=base_gea_lexit_NF_MACE)
summ9a <- summary(c1); summ9b <- summary(c2);summ9c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(anticoagulantes_orales) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ANTICOAG <- deaths$Deaths[deaths$anticoagulantes_orales==1]
d.NON.ANTICOAG <- deaths$Deaths[deaths$anticoagulantes_orales==0]

#Aspirin
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(aspirina==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(aspirina==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*aspirina, data=base_gea_lexit_NF_MACE)
summ10a <- summary(c1); summ10b <- summary(c2);summ10c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(aspirina) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ASA <- deaths$Deaths[deaths$aspirina==1]
d.NON.ASA <- deaths$Deaths[deaths$aspirina==0]

#Glucose
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(gluc_100==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(gluc_100==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*gluc_100, data=base_gea_lexit_NF_MACE)
summ11a <- summary(c1); summ11b <- summary(c2);summ11c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(gluc_100) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GLUC <- deaths$Deaths[deaths$gluc_100==1]
d.NON.GLUC <- deaths$Deaths[deaths$gluc_100==0]

#Tryglycerides
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(tg_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(tg_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*tg_cat, data=base_gea_lexit_NF_MACE)
summ12a <- summary(c1); summ12b <- summary(c2);summ12c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(tg_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.TG <- deaths$Deaths[deaths$tg_cat==1]
d.NON.TG <- deaths$Deaths[deaths$tg_cat==0]


#LDL
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(cldl_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(cldl_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*cldl_cat, data=base_gea_lexit_NF_MACE)
summ13a <- summary(c1); summ13b <- summary(c2);summ13c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(cldl_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.LDL <- deaths$Deaths[deaths$cldl_cat==1]
d.NON.LDL <- deaths$Deaths[deaths$cldl_cat==0]


#NON-HDL
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(nhdl_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(nhdl_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*nhdl_cat, data=base_gea_lexit_NF_MACE)
summ14a <- summary(c1); summ14b <- summary(c2);summ14c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(nhdl_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NHDL <- deaths$Deaths[deaths$nhdl_cat==1]
d.NON.NHDL <- deaths$Deaths[deaths$nhdl_cat==0]

#APO-B p75 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(apob_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(apob_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*apob_p75_neav, data=base_gea_lexit_NF_MACE)
summ15a <- summary(c1); summ15b <- summary(c2); summ15c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(apob_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.APOB <- deaths$Deaths[deaths$apob_p75_neav==1]
d.NON.APOB <- deaths$Deaths[deaths$apob_p75_neav==0]


#CRP
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(pcr_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(pcr_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*pcr_cat, data=base_gea_lexit_NF_MACE)
summ16a <- summary(c1); summ16b <- summary(c2);summ16c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(pcr_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.PCR <- deaths$Deaths[deaths$pcr_cat==1]
d.NON.PCR <- deaths$Deaths[deaths$pcr_cat==0]


#HOMA-2IR
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(homa_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(homa_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*homa_p75_neav, data=base_gea_lexit_NF_MACE)
summ17a <- summary(c1); summ17b <- summary(c2);summ17c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(homa_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.HOMA <- deaths$Deaths[deaths$homa_p75_neav==1]
d.NON.HOMA <- deaths$Deaths[deaths$homa_p75_neav==0]

#RITA
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(rita_p25_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(rita_p25_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*rita_p25_neav, data=base_gea_lexit_NF_MACE)
summ18a <- summary(c1); summ18b <- summary(c2); summ18c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(rita_p25_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.RITA <- deaths$Deaths[deaths$rita_p25_neav==1]
d.NON.RITA <- deaths$Deaths[deaths$rita_p25_neav==0]


#Adiponectin
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(adipo_p25_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(adipo_p25_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*adipo_p25_neav, data=base_gea_lexit_NF_MACE)
summ19a <- summary(c1); summ19b <- summary(c2); summ19c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(adipo_p25_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ADIPONEC <- deaths$Deaths[deaths$adipo_p25_neav==1]
d.NON.ADIPONEC <- deaths$Deaths[deaths$adipo_p25_neav==0]


#GAT p>75
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(GAT_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(GAT_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*GAT_p75_neav, data=base_gea_lexit_NF_MACE)
summ20a <- summary(c1); summ20b <- summary(c2); summ20c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(GAT_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GAT <- deaths$Deaths[deaths$GAT_p75_neav==1]
d.NON.GAT <- deaths$Deaths[deaths$GAT_p75_neav==0]

#GAS p>75
#Stratification

c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(GAS_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_NF_MACE%>%dplyr::filter(GAS_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*GAS_p75_neav, data=base_gea_lexit_NF_MACE)
summ21a <- summary(c1); summ21b <- summary(c2); summ21c <- summary(c3)


## Deaths by sex definition ##
deaths <- base %>% group_by(GAS_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GAS <- deaths$Deaths[deaths$GAS_p75_neav==1]
d.NON.GAS <- deaths$Deaths[deaths$GAS_p75_neav==0]

#Dataset Merging
F.TAB_WWI <- rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4]),c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4]),c(summ1c$coefficients[3,2], summ1c$conf.int[3,3:4]),
  c(summ2a$coefficients[1,2], summ2a$conf.int[1,3:4]),c(summ2b$coefficients[1,2], summ2b$conf.int[1,3:4]),c(summ2c$coefficients[3,2], c(1,1)),
  c(summ3a$coefficients[1,2], summ3a$conf.int[1,3:4]),c(summ3b$coefficients[1,2], summ3b$conf.int[1,3:4]),c(summ3c$coefficients[3,2], summ3c$conf.int[3,3:4]),
  c(summ4a$coefficients[1,2], summ4a$conf.int[1,3:4]),c(summ4b$coefficients[1,2], summ4b$conf.int[1,3:4]),c(summ4c$coefficients[1,2], summ4c$conf.int[1,3:4]),c(summ4d$coefficients[3,2], summ4d$conf.int[3,3:4]),
  c(summ5a$coefficients[1,2], summ5a$conf.int[1,3:4]),c(summ5b$coefficients[1,2], summ5b$conf.int[1,3:4]),c(summ5c$coefficients[3,2], summ5c$conf.int[3,3:4]),
  c(summ6a$coefficients[1,2], summ6a$conf.int[1,3:4]),c(summ6b$coefficients[1,2], summ6b$conf.int[1,3:4]),c(summ6c$coefficients[3,2], summ6c$conf.int[3,3:4]),
  c(summ7a$coefficients[1,2], summ7a$conf.int[1,3:4]),c(summ7b$coefficients[1,2], summ7b$conf.int[1,3:4]),c(summ7c$coefficients[3,2], summ7c$conf.int[3,3:4]),
  c(summ8a$coefficients[1,2], summ8a$conf.int[1,3:4]),c(summ8b$coefficients[1,2], summ8b$conf.int[1,3:4]),c(summ8c$coefficients[3,2], summ8c$conf.int[3,3:4]),
  c(summ9a$coefficients[1,2], summ9a$conf.int[1,3:4]),c(summ9b$coefficients[1,2], summ9b$conf.int[1,3:4]),c(summ9c$coefficients[3,2], summ9c$conf.int[3,3:4]),
  c(summ10a$coefficients[1,2], summ10a$conf.int[1,3:4]),c(summ10b$coefficients[1,2], summ10b$conf.int[1,3:4]),c(summ10c$coefficients[3,2], summ10c$conf.int[3,3:4]),
  c(summ11a$coefficients[1,2], summ11a$conf.int[1,3:4]),c(summ11b$coefficients[1,2], summ11b$conf.int[1,3:4]),c(summ11c$coefficients[3,2], summ11c$conf.int[3,3:4]),
  c(summ12a$coefficients[1,2], summ12a$conf.int[1,3:4]),c(summ12b$coefficients[1,2], summ12b$conf.int[1,3:4]),c(summ12c$coefficients[3,2], summ12c$conf.int[3,3:4]),
  c(summ13a$coefficients[1,2], summ13a$conf.int[1,3:4]),c(summ13b$coefficients[1,2], summ13b$conf.int[1,3:4]),c(summ13c$coefficients[3,2], summ13c$conf.int[3,3:4]),
  c(summ14a$coefficients[1,2], summ14a$conf.int[1,3:4]),c(summ14b$coefficients[1,2], summ14b$conf.int[1,3:4]),c(summ14c$coefficients[3,2], summ14c$conf.int[3,3:4]),
  c(summ15a$coefficients[1,2], summ15a$conf.int[1,3:4]),c(summ15b$coefficients[1,2], summ15b$conf.int[1,3:4]),c(summ15c$coefficients[3,2], summ15c$conf.int[3,3:4]),
  c(summ16a$coefficients[1,2], summ16a$conf.int[1,3:4]),c(summ16b$coefficients[1,2], summ16b$conf.int[1,3:4]),c(summ16c$coefficients[3,2], summ16c$conf.int[3,3:4]),
  c(summ17a$coefficients[1,2], summ17a$conf.int[1,3:4]),c(summ17b$coefficients[1,2], summ17b$conf.int[1,3:4]),c(summ17c$coefficients[3,2], summ17c$conf.int[3,3:4]),
  c(summ18a$coefficients[1,2], summ18a$conf.int[1,3:4]),c(summ18b$coefficients[1,2], summ18b$conf.int[1,3:4]),c(summ18c$coefficients[3,2], summ18c$conf.int[3,3:4]),
  c(summ19a$coefficients[1,2], summ19a$conf.int[1,3:4]),c(summ19b$coefficients[1,2], summ19b$conf.int[1,3:4]),c(summ19c$coefficients[3,2], summ19c$conf.int[3,3:4]),
  c(summ20a$coefficients[1,2], summ20a$conf.int[1,3:4]),c(summ20b$coefficients[1,2], summ20b$conf.int[1,3:4]),c(summ20c$coefficients[3,2], summ20c$conf.int[3,3:4]),
  c(summ21a$coefficients[1,2], summ21a$conf.int[1,3:4]),c(summ21b$coefficients[1,2], summ21b$conf.int[1,3:4]),c(summ21c$coefficients[3,2], summ21c$conf.int[3,3:4])) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper")) %>% 
  mutate("Cause"=c(rep(1,3),
                   rep(2,3),
                   rep(3,3),
                   rep(4,4),
                   rep(5,3),
                   rep(6,3),
                   rep(7,3),
                   rep(8,3),
                   rep(9,3),
                   rep(10,3),
                   rep(11,3),
                   rep(12,3),
                   rep(13,3),
                   rep(14,3),
                   rep(15,3),
                   rep(16,3),
                   rep(17,3),
                   rep(18,3),
                   rep(19,3),
                   rep(20,3),
                   rep(21,3)), 
         "Strat"=rep(1:64,1)) %>% 
  mutate("Deaths"=c(d.MALE,d.FEMALE,"", #1
                    d.AGE.65,d.AGE.non65,"", #2
                    d.CVD.EDAD.35,d.CVD.EDAD.35.65,"", #3
                    d.NORMO,d.SOBRE,d.OBES,"", #4
                    d.NEVER,d.PRE.CURR,"", #5
                    d.DIABETES,d.NON.DIABETES,"", #6
                    d.HTA,d.NON.HTA,"", #7
                    d.ESTATINA,d.NON.ESTATINA,"", #8
                    d.ANTICOAG,d.NON.ANTICOAG,"", #9
                    d.ASA,d.NON.ASA,"", #10
                    d.LDL,d.NON.LDL,"", #11
                    d.GLUC,d.NON.GLUC,"", #12
                    d.TG,d.NON.TG,"", #13
                    d.NHDL,d.NON.NHDL,"", #14
                    d.APOB,d.NON.APOB,"", #15
                    d.PCR,d.NON.PCR,"", #16
                    d.HOMA,d.NON.HOMA,"", #17
                    d.RITA,d.NON.RITA,"", #18
                    d.ADIPONEC,d.NON.ADIPONEC,"", #19
                    d.GAT,d.NON.GAT,"", #20
                    d.GAS,d.NON.GAS,"")) %>% #21
  mutate(Strat=ordered(Strat, 1:64, c("Male","Female","Interaction", #1
                                      ">65 years","≤65 years","Interaction", #2
                                      "<35 years","35-65 years","Interaction", #3
                                      "Underweight","Overweight","Obesity","Interaction", #4
                                      "Never","Previous-Current","Interaction", #5
                                      "Diabetes","Without Diabetes","Interaction", #6
                                      "Hypertension","Without Hypertension","Interaction", #7
                                      "Use","Non-User","Interaction", #8
                                      "Use","Non-User","Interaction", #9
                                      "Use","Non-User","Interaction", #10
                                      "≥100 mg/dl","<100 mg/dl","Interaction", #11
                                      "<150 mg/dl","≥150 mg/dl","Interaction", #12
                                      "≤70 mg/dl",">70 mg/dl","Interaction", #13
                                      "<130 mg/dl","≥130 mg/dl","Interaction", #14
                                      "≥p75","<p75","Interaction", #15
                                      "≥3.0 pg/ml","<3.0 pg/ml","Interaction", #16
                                      "≥p75","<p75","Interaction", #17
                                      "≥p75","<p75","Interaction", #18
                                      "<p25","≥p25","Interaction", #19
                                      "≥p75 ","<p75 ","Interaction", #21
                                      "≥p75 ","<p75 ","Interaction"))) %>%  #22
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))

F.TAB_WWI$mean[is.na(F.TAB_WWI$mean)]<-1
F.TAB_WWI$lower[is.na(F.TAB_WWI$lower)]<-1
F.TAB_WWI$upper[is.na(F.TAB_WWI$upper)]<-1

F.TAB_WWI$mean[is.infinite(F.TAB_WWI$mean)]<-1
F.TAB_WWI$lower[is.infinite(F.TAB_WWI$lower)]<-1
F.TAB_WWI$upper[is.infinite(F.TAB_WWI$upper)]<-1

F.TAB_WWI$mean[F.TAB_WWI$mean==0]<-1
F.TAB_WWI$lower[F.TAB_WWI$lower==0]<-1
F.TAB_WWI$upper[F.TAB_WWI$upper==0]<-1

F.TAB_WWI$upper<-replace(F.TAB_WWI$upper, F.TAB_WWI$upper>10, 1)
F.TAB_WWI$lower<-replace(F.TAB_WWI$lower, F.TAB_WWI$lower>10, 1)
F.TAB_WWI$mean<-replace(F.TAB_WWI$mean, F.TAB_WWI$mean>10, 1)

F.TAB_WWI$upper<-replace(F.TAB_WWI$upper, F.TAB_WWI$upper<0.5, 1)
F.TAB_WWI$lower<-replace(F.TAB_WWI$lower, F.TAB_WWI$lower<0.5, 1)
F.TAB_WWI$mean<-replace(F.TAB_WWI$mean, F.TAB_WWI$mean<0.5, 1)

## Build forestplot ##
fnfp1 <- fpDrawNormalCI; fnfp2 <- fpDrawDiamondCI; gp1 <- gpar(fill="#6699CC"); gp2 <- gpar(fill="#994455"); 

#Label of Figure
d.labs <- c("Sex",NA,NA,NA, #1
            "Age",NA,NA,NA, #2
            "Age at CVD",NA,NA,NA, #3
            "BMI Categories",NA,NA,NA,NA, #4
            "Smoking Status",NA,NA,NA, #5
            "Diabetes",NA,NA,NA, #6
            "Arterial Hypertension",NA,NA,NA, #7
            "Statin",NA,NA,NA, #8
            "Anticoagulants",NA,NA,NA, #9
            "Aspirin",NA,NA,NA, #10
            "Glucose",NA,NA,NA, #11
            "Triglycerides Goals",NA,NA,NA, #12
            "LDL Goals",NA,NA,NA, #13
            "Non-HDL-C Goals",NA,NA,NA, #14
            "APO-B Levels",NA,NA,NA, #15
            "C-Reactive Protein",NA,NA,NA, #16
            "HOMA-IR",NA,NA,NA, #17
            "ADIPO-IR",NA,NA,NA, #18
            "Adiponectin",NA,NA,NA, #19
            "Total Adiposse Tissue",NA,NA,NA,#21
            "Subcutaneous Adiposse Tissue",NA,NA,NA)#22


p.inter.val<- round(rbind(NA,NA,NA,summ1c$coefficients[3,5],
                          NA,NA,NA,summ2c$coefficients[3,5],
                          NA,NA,NA,summ3c$coefficients[3,5],
                          NA,NA,NA,NA,summ4d$coefficients[3,5],
                          NA,NA,NA,summ5c$coefficients[3,5],
                          NA,NA,NA,summ6c$coefficients[3,5],
                          NA,NA,NA,summ7c$coefficients[3,5],
                          NA,NA,NA,summ8c$coefficients[3,5],
                          NA,NA,NA,summ9c$coefficients[3,5],
                          NA,NA,NA,summ10c$coefficients[3,5],
                          NA,NA,NA,summ11c$coefficients[3,5],
                          NA,NA,NA,summ12c$coefficients[3,5],
                          NA,NA,NA,summ13c$coefficients[3,5],
                          NA,NA,NA,summ14c$coefficients[3,5],
                          NA,NA,NA,summ15c$coefficients[3,5],
                          NA,NA,NA,summ16c$coefficients[3,5],
                          NA,NA,NA,summ17c$coefficients[3,5],
                          NA,NA,NA,summ18c$coefficients[3,5],
                          NA,NA,NA,summ19c$coefficients[3,5],
                          NA,NA,NA,summ20c$coefficients[3,5],
                          NA,NA,NA,summ21c$coefficients[3,5]),3)

rbind(NA, F.TAB_WWI[1:3,], #1
      NA, F.TAB_WWI[4:6,], #2
      NA, F.TAB_WWI[7:9,], #3
      NA, F.TAB_WWI[10:13,], #4
      NA, F.TAB_WWI[14:16,], #5
      NA, F.TAB_WWI[17:19,], #6
      NA, F.TAB_WWI[20:22,], #8 
      NA, F.TAB_WWI[23:25,], #9
      NA, F.TAB_WWI[26:28,], #10
      NA, F.TAB_WWI[29:31,], #11
      NA, F.TAB_WWI[32:34,], #12
      NA, F.TAB_WWI[35:37,], #13
      NA, F.TAB_WWI[38:40,], #14
      NA, F.TAB_WWI[41:43,], #15
      NA, F.TAB_WWI[44:46,], #16
      NA, F.TAB_WWI[47:49,], #17
      NA, F.TAB_WWI[50:52,], #18
      NA, F.TAB_WWI[53:55,], #19
      NA, F.TAB_WWI[56:58,], #20
      NA, F.TAB_WWI[59:61,], #21
      NA, F.TAB_WWI[62:64,]) %>% #21
  cbind(p.inter.val)%>%
  mutate(Cause=d.labs)%>%
  forestplot::forestplot(
    xticks = c(-0.015, 0.028),
    col = fpColors(lines = "black", zero = "black"),
    labeltext = c("Cause", "Strat", "Deaths", "HR","p.inter.val"), 
    graph.pos = 4, xlog = T)%>% 
  fp_add_header(Cause="Stratification", 
                Strat = c("Categories"), 
                Deaths = "Non-Fatal MACE Events", 
                HR = c("HR (95% CI)"),
                p.inter.val = c("p for interaction"))%>%
  fp_set_style(txt_gp = fpTxtGp(label = gpar(fontfamily = "Arial"), 
                                ticks = gpar(cex = 1), 
                                xlab  = gpar(cex=1.20, fontface="bold.italic"),
                                summary = gpar(fontfamily = "Arial", fontface="bold.italic")), align="c", summary = "#7A316F") %>% 
  fp_set_zebra_style("gray90","white") -> sup.forest1

png("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/VAT - MACE/Supplementary_Figure_2A.png",
    width=35, height=32.5, units = "cm", res = 300);sup.forest1
dev.off(); print(plot(1))

#Fatal MACE
#Sex 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(sexo==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(sexo==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*sexo, data=base_gea_lexit_F_MACE)
summ1a <- summary(c1); summ1b <- summary(c2);summ1c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(sexo) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.MALE <- deaths$Deaths[deaths$sexo==1]
d.FEMALE <- deaths$Deaths[deaths$sexo==0]

#Age 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(edad_65==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(edad_65==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*edad_65, data=base_gea_lexit_F_MACE)
summ2a <- summary(c1); summ2b <- summary(c2);summ2c <- summary(c3)

## Deaths by age definition ##
deaths <- base %>% group_by(edad_65) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.AGE.65 <- deaths$Deaths[deaths$edad_65==1]
d.AGE.non65 <- deaths$Deaths[deaths$edad_65==0]

#Age at CVD 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(AGE_CVD_CAT==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(AGE_CVD_CAT==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*AGE_CVD_CAT, data=base_gea_lexit_F_MACE)
summ3a <- summary(c1); summ3b <- summary(c2);summ3c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(AGE_CVD_CAT) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.CVD.EDAD.35<- deaths$Deaths[deaths$AGE_CVD_CAT==1]
d.CVD.EDAD.35.65 <- deaths$Deaths[deaths$AGE_CVD_CAT==0]

#BMI Categories
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(BMI_CAT==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(BMI_CAT==2))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(BMI_CAT==3))
c4<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*BMI_CAT, data=base_gea_lexit_F_MACE)
summ4a <- summary(c1); summ4b <- summary(c2);summ4c <- summary(c3);summ4d <- summary(c4)

## Deaths by sex definition ##
deaths <- base %>% group_by(BMI_CAT) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NORMO<- deaths$Deaths[deaths$BMI_CAT==1]
d.SOBRE <- deaths$Deaths[deaths$BMI_CAT==2]
d.OBES <- deaths$Deaths[deaths$BMI_CAT==3]

#Smoking Status
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(tabaquismo_cat=="Never"))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(tabaquismo_cat %in% c("Previous","Current")))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*(tabaquismo_cat=="Never"), data=base_gea_lexit_F_MACE)
summ5a <- summary(c1); summ5b <- summary(c2);summ5c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(tabaquismo_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NEVER <- deaths$Deaths[deaths$tabaquismo_cat=="Never"]
d.PRE.CURR <- sum(deaths$Deaths[deaths$tabaquismo_cat=="Previous" | deaths$tabaquismo_cat=="Current"])

#Diabetes
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(diabetes==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(diabetes==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*diabetes, data=base_gea_lexit_F_MACE)
summ6a <- summary(c1); summ6b <- summary(c2);summ6c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(diabetes) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.DIABETES <- deaths$Deaths[deaths$diabetes==1]
d.NON.DIABETES <- deaths$Deaths[deaths$diabetes==0]

#Hypertension Status
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(hta==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(hta==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*hta, data=base_gea_lexit_F_MACE)
summ7a <- summary(c1); summ7b <- summary(c2);summ7c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(hta) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.HTA <- deaths$Deaths[deaths$hta==1]
d.NON.HTA <- deaths$Deaths[deaths$hta==0]

#Statins Use
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(Estatinas_REC==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(Estatinas_REC==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*Estatinas_REC, data=base_gea_lexit_F_MACE)
summ8a <- summary(c1); summ8b <- summary(c2);summ8c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(Estatinas_REC) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ESTATINA <- deaths$Deaths[deaths$Estatinas_REC==1]
d.NON.ESTATINA <- deaths$Deaths[deaths$Estatinas_REC==0]

#Anticoagulation
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(anticoagulantes_orales==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(anticoagulantes_orales==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*anticoagulantes_orales, data=base_gea_lexit_F_MACE)
summ9a <- summary(c1); summ9b <- summary(c2);summ9c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(anticoagulantes_orales) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ANTICOAG <- deaths$Deaths[deaths$anticoagulantes_orales==1]
d.NON.ANTICOAG <- deaths$Deaths[deaths$anticoagulantes_orales==0]

#Aspirin
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(aspirina==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(aspirina==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*aspirina, data=base_gea_lexit_F_MACE)
summ10a <- summary(c1); summ10b <- summary(c2);summ10c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(aspirina) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ASA <- deaths$Deaths[deaths$aspirina==1]
d.NON.ASA <- deaths$Deaths[deaths$aspirina==0]

#Glucose
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(gluc_100==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(gluc_100==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*gluc_100, data=base_gea_lexit_F_MACE)
summ11a <- summary(c1); summ11b <- summary(c2);summ11c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(gluc_100) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GLUC <- deaths$Deaths[deaths$gluc_100==1]
d.NON.GLUC <- deaths$Deaths[deaths$gluc_100==0]

#Tryglycerides
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(tg_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(tg_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*tg_cat, data=base_gea_lexit_F_MACE)
summ12a <- summary(c1); summ12b <- summary(c2);summ12c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(tg_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.TG <- deaths$Deaths[deaths$tg_cat==1]
d.NON.TG <- deaths$Deaths[deaths$tg_cat==0]


#LDL
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(cldl_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(cldl_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*cldl_cat, data=base_gea_lexit_F_MACE)
summ13a <- summary(c1); summ13b <- summary(c2);summ13c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(cldl_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.LDL <- deaths$Deaths[deaths$cldl_cat==1]
d.NON.LDL <- deaths$Deaths[deaths$cldl_cat==0]


#NON-HDL
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(nhdl_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(nhdl_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*nhdl_cat, data=base_gea_lexit_F_MACE)
summ14a <- summary(c1); summ14b <- summary(c2);summ14c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(nhdl_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.NHDL <- deaths$Deaths[deaths$nhdl_cat==1]
d.NON.NHDL <- deaths$Deaths[deaths$nhdl_cat==0]

#APO-B p75 
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(apob_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(apob_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*apob_p75_neav, data=base_gea_lexit_F_MACE)
summ15a <- summary(c1); summ15b <- summary(c2); summ15c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(apob_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.APOB <- deaths$Deaths[deaths$apob_p75_neav==1]
d.NON.APOB <- deaths$Deaths[deaths$apob_p75_neav==0]


#CRP
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(pcr_cat==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(pcr_cat==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*pcr_cat, data=base_gea_lexit_F_MACE)
summ16a <- summary(c1); summ16b <- summary(c2);summ16c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(pcr_cat) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.PCR <- deaths$Deaths[deaths$pcr_cat==1]
d.NON.PCR <- deaths$Deaths[deaths$pcr_cat==0]


#HOMA-2IR
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(homa_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(homa_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*homa_p75_neav, data=base_gea_lexit_F_MACE)
summ17a <- summary(c1); summ17b <- summary(c2);summ17c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(homa_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.HOMA <- deaths$Deaths[deaths$homa_p75_neav==1]
d.NON.HOMA <- deaths$Deaths[deaths$homa_p75_neav==0]

#RITA
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(rita_p25_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(rita_p25_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*rita_p25_neav, data=base_gea_lexit_F_MACE)
summ18a <- summary(c1); summ18b <- summary(c2); summ18c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(rita_p25_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.RITA <- deaths$Deaths[deaths$rita_p25_neav==1]
d.NON.RITA <- deaths$Deaths[deaths$rita_p25_neav==0]


#Adiponectin
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(adipo_p25_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(adipo_p25_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*adipo_p25_neav, data=base_gea_lexit_F_MACE)
summ19a <- summary(c1); summ19b <- summary(c2); summ19c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(adipo_p25_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.ADIPONEC <- deaths$Deaths[deaths$adipo_p25_neav==1]
d.NON.ADIPONEC <- deaths$Deaths[deaths$adipo_p25_neav==0]


#GAT p>75
#Stratification
c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(GAT_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(GAT_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*GAT_p75_neav, data=base_gea_lexit_F_MACE)
summ20a <- summary(c1); summ20b <- summary(c2); summ20c <- summary(c3)

## Deaths by sex definition ##
deaths <- base %>% group_by(GAT_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GAT <- deaths$Deaths[deaths$GAT_p75_neav==1]
d.NON.GAT <- deaths$Deaths[deaths$GAT_p75_neav==0]

#GAS p>75
#Stratification

c1<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(GAS_p75_neav==1))
c2<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav, data=base_gea_lexit_F_MACE%>%dplyr::filter(GAS_p75_neav==0))
c3<-coxph(Surv(time_at_entry, time_at_exit, exit_status) ~ gav*GAS_p75_neav, data=base_gea_lexit_F_MACE)
summ21a <- summary(c1); summ21b <- summary(c2); summ21c <- summary(c3)


## Deaths by sex definition ##
deaths <- base %>% group_by(GAS_p75_neav) %>% summarise(
  "MACE"=sum(morbimortalidad_mace_2)) %>%
  pivot_longer(!1, names_to = "Cause", values_to = "Deaths")
d.GAS <- deaths$Deaths[deaths$GAS_p75_neav==1]
d.NON.GAS <- deaths$Deaths[deaths$GAS_p75_neav==0]

#Dataset Merging
F.TAB_WWI <- rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4]),c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4]), c(summ1c$coefficients[3,2], summ1c$conf.int[3,3:4]),
  c(summ2a$coefficients[1,2], summ2a$conf.int[1,3:4]),c(summ2b$coefficients[1,2], summ2b$conf.int[1,3:4]),c(summ2c$coefficients[3,2], summ2c$conf.int[3,3:4]),
  c(summ3a$coefficients[1,2], summ3a$conf.int[1,3:4]),c(summ3b$coefficients[1,2], summ3b$conf.int[1,3:4]),c(summ3c$coefficients[3,2], summ3c$conf.int[3,3:4]),
  c(summ4a$coefficients[1,2], summ4a$conf.int[1,3:4]),c(summ4b$coefficients[1,2], summ4b$conf.int[1,3:4]),c(summ4c$coefficients[1,2], summ4c$conf.int[1,3:4]),c(summ4d$coefficients[3,2], summ4d$conf.int[3,3:4]),
  c(summ5a$coefficients[1,2], summ5a$conf.int[1,3:4]),c(summ5b$coefficients[1,2], summ5b$conf.int[1,3:4]),c(summ5c$coefficients[3,2], summ5c$conf.int[3,3:4]),
  c(summ6a$coefficients[1,2], summ6a$conf.int[1,3:4]),c(summ6b$coefficients[1,2], summ6b$conf.int[1,3:4]),c(summ6c$coefficients[3,2], summ6c$conf.int[3,3:4]),
  c(summ7a$coefficients[1,2], summ7a$conf.int[1,3:4]),c(summ7b$coefficients[1,2], summ7b$conf.int[1,3:4]),c(summ7c$coefficients[3,2], summ7c$conf.int[3,3:4]),
  c(summ8a$coefficients[1,2], summ8a$conf.int[1,3:4]),c(summ8b$coefficients[1,2], summ8b$conf.int[1,3:4]),c(summ8c$coefficients[3,2], summ8c$conf.int[3,3:4]),
  c(summ9a$coefficients[1,2], summ9a$conf.int[1,3:4]),c(summ9b$coefficients[1,2], summ9b$conf.int[1,3:4]),c(summ9c$coefficients[3,2], summ9c$conf.int[3,3:4]),
  c(summ10a$coefficients[1,2], summ10a$conf.int[1,3:4]),c(summ10b$coefficients[1,2], summ10b$conf.int[1,3:4]),c(summ10c$coefficients[3,2], summ10c$conf.int[3,3:4]),
  c(summ11a$coefficients[1,2], summ11a$conf.int[1,3:4]),c(summ11b$coefficients[1,2], summ11b$conf.int[1,3:4]),c(summ11c$coefficients[3,2], summ11c$conf.int[3,3:4]),
  c(summ12a$coefficients[1,2], summ12a$conf.int[1,3:4]),c(summ12b$coefficients[1,2], summ12b$conf.int[1,3:4]),c(summ12c$coefficients[3,2], summ12c$conf.int[3,3:4]),
  c(summ13a$coefficients[1,2], summ13a$conf.int[1,3:4]),c(summ13b$coefficients[1,2], summ13b$conf.int[1,3:4]),c(summ13c$coefficients[3,2], summ13c$conf.int[3,3:4]),
  c(summ14a$coefficients[1,2], summ14a$conf.int[1,3:4]),c(summ14b$coefficients[1,2], summ14b$conf.int[1,3:4]),c(summ14c$coefficients[3,2], summ14c$conf.int[3,3:4]),
  c(summ15a$coefficients[1,2], summ15a$conf.int[1,3:4]),c(summ15b$coefficients[1,2], summ15b$conf.int[1,3:4]),c(summ15c$coefficients[3,2], summ15c$conf.int[3,3:4]),
  c(summ16a$coefficients[1,2], summ16a$conf.int[1,3:4]),c(summ16b$coefficients[1,2], summ16b$conf.int[1,3:4]),c(summ16c$coefficients[3,2], summ16c$conf.int[3,3:4]),
  c(summ17a$coefficients[1,2], summ17a$conf.int[1,3:4]),c(summ17b$coefficients[1,2], summ17b$conf.int[1,3:4]),c(summ17c$coefficients[3,2], summ17c$conf.int[3,3:4]),
  c(summ18a$coefficients[1,2], summ18a$conf.int[1,3:4]),c(summ18b$coefficients[1,2], summ18b$conf.int[1,3:4]),c(summ18c$coefficients[3,2], summ18c$conf.int[3,3:4]),
  c(summ19a$coefficients[1,2], summ19a$conf.int[1,3:4]),c(summ19b$coefficients[1,2], summ19b$conf.int[1,3:4]),c(summ19c$coefficients[3,2], summ19c$conf.int[3,3:4]),
  c(summ20a$coefficients[1,2], summ20a$conf.int[1,3:4]),c(summ20b$coefficients[1,2], summ20b$conf.int[1,3:4]),c(summ20c$coefficients[3,2], summ20c$conf.int[3,3:4]),
  c(summ21a$coefficients[1,2], summ21a$conf.int[1,3:4]),c(summ21b$coefficients[1,2], summ21b$conf.int[1,3:4]),c(summ21c$coefficients[3,2], summ21c$conf.int[3,3:4])) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper")) %>% 
  mutate("Cause"=c(rep(1,3),
                   rep(2,3),
                   rep(3,3),
                   rep(4,4),
                   rep(5,3),
                   rep(6,3),
                   rep(7,3),
                   rep(8,3),
                   rep(9,3),
                   rep(10,3),
                   rep(11,3),
                   rep(12,3),
                   rep(13,3),
                   rep(14,3),
                   rep(15,3),
                   rep(16,3),
                   rep(17,3),
                   rep(18,3),
                   rep(19,3),
                   rep(20,3),
                   rep(21,3)), 
         "Strat"=rep(1:64,1)) %>% 
  mutate("Deaths"=c(d.MALE,d.FEMALE,"", #1
                    d.AGE.65,d.AGE.non65,"", #2
                    d.CVD.EDAD.35,d.CVD.EDAD.35.65,"", #3
                    d.NORMO,d.SOBRE,d.OBES,"", #4
                    d.NEVER,d.PRE.CURR,"", #5
                    d.DIABETES,d.NON.DIABETES,"", #6
                    d.HTA,d.NON.HTA,"", #7
                    d.ESTATINA,d.NON.ESTATINA,"", #8
                    d.ANTICOAG,d.NON.ANTICOAG,"", #9
                    d.ASA,d.NON.ASA,"", #10
                    d.LDL,d.NON.LDL,"", #11
                    d.GLUC,d.NON.GLUC,"", #12
                    d.TG,d.NON.TG,"", #13
                    d.NHDL,d.NON.NHDL,"", #14
                    d.APOB,d.NON.APOB,"", #15
                    d.PCR,d.NON.PCR,"", #16
                    d.HOMA,d.NON.HOMA,"", #17
                    d.RITA,d.NON.RITA,"", #18
                    d.ADIPONEC,d.NON.ADIPONEC,"", #19
                    d.GAT,d.NON.GAT,"", #20
                    d.GAS,d.NON.GAS,"")) %>% #21
  mutate(Strat=ordered(Strat, 1:64, c("Male","Female","Interaction", #1
                                      ">65 years","≤65 years","Interaction", #2
                                      "<35 years","35-65 years","Interaction", #3
                                      "Underweight","Overweight","Obesity","Interaction", #4
                                      "Never","Previous-Current","Interaction", #5
                                      "Diabetes","Without Diabetes","Interaction", #6
                                      "Hypertension","Without Hypertension","Interaction", #7
                                      "Use","Non-User","Interaction", #8
                                      "Use","Non-User","Interaction", #9
                                      "Use","Non-User","Interaction", #10
                                      "≥100 mg/dl","<100 mg/dl","Interaction", #11
                                      "<150 mg/dl","≥150 mg/dl","Interaction", #12
                                      "≤70 mg/dl",">70 mg/dl","Interaction", #13
                                      "<130 mg/dl","≥130 mg/dl","Interaction", #14
                                      "≥p75","<p75","Interaction", #15
                                      "≥3.0 pg/ml","<3.0 pg/ml","Interaction", #16
                                      "≥p75","<p75","Interaction", #17
                                      "≥p75","<p75","Interaction", #18
                                      "<p25","≥p25","Interaction", #19
                                      "≥p75 ","<p75 ","Interaction", #21
                                      "≥p75 ","<p75 ","Interaction"))) %>%  #22
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))

F.TAB_WWI$mean[is.na(F.TAB_WWI$mean)]<-1
F.TAB_WWI$lower[is.na(F.TAB_WWI$lower)]<-1
F.TAB_WWI$upper[is.na(F.TAB_WWI$upper)]<-1

F.TAB_WWI$mean[is.infinite(F.TAB_WWI$mean)]<-1
F.TAB_WWI$lower[is.infinite(F.TAB_WWI$lower)]<-1
F.TAB_WWI$upper[is.infinite(F.TAB_WWI$upper)]<-1

F.TAB_WWI$mean[F.TAB_WWI$mean==0]<-1
F.TAB_WWI$lower[F.TAB_WWI$lower==0]<-1
F.TAB_WWI$upper[F.TAB_WWI$upper==0]<-1

F.TAB_WWI$upper<-replace(F.TAB_WWI$upper, F.TAB_WWI$upper>10, 1)
F.TAB_WWI$lower<-replace(F.TAB_WWI$lower, F.TAB_WWI$lower>10, 1)
F.TAB_WWI$mean<-replace(F.TAB_WWI$mean, F.TAB_WWI$mean>10, 1)

F.TAB_WWI$upper<-replace(F.TAB_WWI$upper, F.TAB_WWI$upper<0.5, 1)
F.TAB_WWI$lower<-replace(F.TAB_WWI$lower, F.TAB_WWI$lower<0.5, 1)
F.TAB_WWI$mean<-replace(F.TAB_WWI$mean, F.TAB_WWI$mean<0.5, 1)


## Build forestplot ##
fnfp1 <- fpDrawNormalCI; fnfp2 <- fpDrawDiamondCI; gp1 <- gpar(fill="#6699CC"); gp2 <- gpar(fill="#994455"); 

#Label of Figure
d.labs <- c("Sex",NA,NA,NA, #1
            "Age",NA,NA,NA, #2
            "Age at CVD",NA,NA,NA, #3
            "BMI Categories",NA,NA,NA,NA, #4
            "Smoking Status",NA,NA,NA, #5
            "Diabetes",NA,NA,NA, #6
            "Arterial Hypertension",NA,NA,NA, #7
            "Statin",NA,NA,NA, #8
            "Anticoagulants",NA,NA,NA, #9
            "Aspirin",NA,NA,NA, #10
            "Glucose",NA,NA,NA, #11
            "Triglycerides Goals",NA,NA,NA, #12
            "LDL Goals",NA,NA,NA, #13
            "Non-HDL-C Goals",NA,NA,NA, #14
            "APO-B Levels",NA,NA,NA, #15
            "C-Reactive Protein",NA,NA,NA, #16
            "HOMA-IR",NA,NA,NA, #17
            "ADIPO-IR",NA,NA,NA, #18
            "Adiponectin",NA,NA,NA, #19
            "Total Adiposse Tissue",NA,NA,NA,#21
            "Subcutaneous Adiposse Tissue",NA,NA,NA)#22


p.inter.val<- round(rbind(NA,NA,NA,summ1c$coefficients[3,5],
                          NA,NA,NA,summ2c$coefficients[3,5],
                          NA,NA,NA,summ3c$coefficients[3,5],
                          NA,NA,NA,NA,summ4d$coefficients[3,5],
                          NA,NA,NA,summ5c$coefficients[3,5],
                          NA,NA,NA,summ6c$coefficients[3,5],
                          NA,NA,NA,summ7c$coefficients[3,5],
                          NA,NA,NA,summ8c$coefficients[3,5],
                          NA,NA,NA,summ9c$coefficients[3,5],
                          NA,NA,NA,summ10c$coefficients[3,5],
                          NA,NA,NA,summ11c$coefficients[3,5],
                          NA,NA,NA,summ12c$coefficients[3,5],
                          NA,NA,NA,summ13c$coefficients[3,5],
                          NA,NA,NA,summ14c$coefficients[3,5],
                          NA,NA,NA,summ15c$coefficients[3,5],
                          NA,NA,NA,summ16c$coefficients[3,5],
                          NA,NA,NA,summ17c$coefficients[3,5],
                          NA,NA,NA,summ18c$coefficients[3,5],
                          NA,NA,NA,summ19c$coefficients[3,5],
                          NA,NA,NA,summ20c$coefficients[3,5],
                          NA,NA,NA,summ21c$coefficients[3,5]),3)

rbind(NA, F.TAB_WWI[1:3,], #1
      NA, F.TAB_WWI[4:6,], #2
      NA, F.TAB_WWI[7:9,], #3
      NA, F.TAB_WWI[10:13,], #4
      NA, F.TAB_WWI[14:16,], #5
      NA, F.TAB_WWI[17:19,], #6
      NA, F.TAB_WWI[20:22,], #8 
      NA, F.TAB_WWI[23:25,], #9
      NA, F.TAB_WWI[26:28,], #10
      NA, F.TAB_WWI[29:31,], #11
      NA, F.TAB_WWI[32:34,], #12
      NA, F.TAB_WWI[35:37,], #13
      NA, F.TAB_WWI[38:40,], #14
      NA, F.TAB_WWI[41:43,], #15
      NA, F.TAB_WWI[44:46,], #16
      NA, F.TAB_WWI[47:49,], #17
      NA, F.TAB_WWI[50:52,], #18
      NA, F.TAB_WWI[53:55,], #19
      NA, F.TAB_WWI[56:58,], #20
      NA, F.TAB_WWI[59:61,], #21
      NA, F.TAB_WWI[62:64,]) %>% #21
  cbind(p.inter.val)%>%
  mutate(Cause=d.labs)%>%
  forestplot::forestplot(
    xticks = c(-0.015, 0.028),
    col = fpColors(lines = "black", zero = "black"),
    labeltext = c("Cause", "Strat", "Deaths", "HR","p.inter.val"), 
    graph.pos = 4, xlog = T)%>% 
  fp_add_header(Cause="Stratification", 
                Strat = c("Categories"), 
                Deaths = "Fatal MACE Events", 
                HR = c("HR (95% CI)"),
                p.inter.val = c("p for interaction")) %>%
  fp_set_style(txt_gp = fpTxtGp(label = gpar(fontfamily = "Arial"), 
                                ticks = gpar(cex = 1), 
                                xlab  = gpar(cex=1.20, fontface="bold.italic"),
                                summary = gpar(fontfamily = "Arial", fontface="bold.italic")), align="c", summary = "#7A316F") %>% 
  fp_set_zebra_style("gray90","white") -> sup.forest2


png("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INCICh/VAT - MACE/Supplementary_Figure_2B.png",
    width=35, height=32.5, units = "cm", res = 300);sup.forest2
dev.off(); print(plot(1))


#####Analysis: Incidence Rate for MACE and its components (Text)####

table(base$tipo_mace,base$mace_no_fatal,useNA = "always")

chisq.test(table(base$tipo_mace==6,base$gav_terciles_neav))


base %>% 
  dplyr::select(gav_terciles_neav,mace_no_fatal,
                tipo_mace)%>%
  tbl_summary(by = gav_terciles_neav,
              missing = "ifany")%>%
  bold_labels()%>%
  add_overall()%>%
  add_p()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()

#####Analysis: MACE Event Breakdown (Supplementary Table 2)#####


#MUERTES REC
base$muerte_string_REC<-NULL
base$muerte_string_REC[base$muerte_string=="IAM"]<-1
base$muerte_string_REC[base$muerte_string=="Arritmia"]<-2
base$muerte_string_REC[base$muerte_string=="FA"]<-2
base$muerte_string_REC[base$muerte_string=="HF"]<-3
base$muerte_string_REC[base$muerte_string=="HF"]<-3
base$muerte_string_REC[base$muerte_string=="EVC"]<-4
base$muerte_string_REC[base$muerte_string=="Stroke"]<-4
base$muerte_string_REC[base$muerte_string=="Fallecio"]<-5
base$muerte_string_REC[base$muerte_string=="fallecio"]<-5
base$muerte_string_REC[base$muerte_string=="fallecio"]<-5
base$muerte_string_REC[base$muerte_string=="Enfisema"]<-5
base$muerte_string_REC[base$muerte_string=="Neumonia"]<-5
base$muerte_string_REC[base$muerte_string=="Sepsis"]<-5
base$muerte_string_REC[base$muerte_string=="Shock"]<-5
base$muerte_string_REC[base$muerte_string=="Trauma"]<-5
base$muerte_string_REC<-na.tools::na.replace(base$muerte_string_REC,0)
base$muerte_string_REC[base$muerte_string_REC==5 & base$mace_fatal==0]<-0
base$muerte_string_REC[base$muerte_string_REC==0 & base$mace_fatal==1]<-1
base$muerte_string_REC[base$muerte_string_REC==0]<-NA

base %>% 
  dplyr::select(gav_terciles_neav,morbimortalidad_mace,
                mace_no_fatal,
                mace_fatal,
                muerte_string_REC)%>%
  tbl_summary(by = gav_terciles_neav,
              missing = "ifany")%>%
  bold_labels()%>%
  add_overall()%>%
  add_p()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Supplementary_Table_2.docx")





#####Analysis: Fully Adjusted Model Association with VAT and MACE (Supplementary Table 3)#####

##OVERALL MACE###
#VAT Grs
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav, data=base_gea_lexit);summary(m0)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit);summary(m1)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit);summary(m2)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.1<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 cm\u00b2"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (Scale)
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav), data=base_gea_lexit);summary(m0)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit);summary(m1)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit);summary(m2)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.2<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 Standard Deviation\n(66.5 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (p75)
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2), data=base_gea_lexit)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.3<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Middle-Tertile\n(140-194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")


Sup.tab.df.4<-rbind(
  c(summ1a$coefficients[2,2], summ1a$conf.int[2,3:4],round(summ1a$coefficients[2,5],4)), 
  c(summ1b$coefficients[2,2], summ1b$conf.int[2,3:4],round(summ1b$coefficients[2,5],4)), 
  c(summ1c$coefficients[2,2], summ1c$conf.int[2,3:4],round(summ1c$coefficients[2,5],4)),
  c(summ1d$coefficients[2,2], summ1d$conf.int[2,3:4],round(summ1d$coefficients[2,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Upper-Tertile\n(>194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

Sup.tab.1.OVERALL.MACE<-rbind(Sup.tab.df.1,Sup.tab.df.2,Sup.tab.df.3,Sup.tab.df.4)


##Non_Fatal MACE
#VAT Grs
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav, data=base_gea_lexit_NF_MACE);summary(m0)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_NF_MACE);summary(m1)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_NF_MACE);summary(m2)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_NF_MACE);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.1<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 cm\u00b2"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (Scale)
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav), data=base_gea_lexit_NF_MACE);summary(m0)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_NF_MACE);summary(m1)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_NF_MACE);summary(m2)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_NF_MACE);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.2<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 Standard Deviation\n(66.5 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (p75)
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2), data=base_gea_lexit_NF_MACE)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_NF_MACE)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_NF_MACE)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_NF_MACE)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.3<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Middle-Tertile\n(140-194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")


Sup.tab.df.4<-rbind(
  c(summ1a$coefficients[2,2], summ1a$conf.int[2,3:4],round(summ1a$coefficients[2,5],4)), 
  c(summ1b$coefficients[2,2], summ1b$conf.int[2,3:4],round(summ1b$coefficients[2,5],4)), 
  c(summ1c$coefficients[2,2], summ1c$conf.int[2,3:4],round(summ1c$coefficients[2,5],4)),
  c(summ1d$coefficients[2,2], summ1d$conf.int[2,3:4],round(summ1d$coefficients[2,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Upper-Tertile\n(>194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

Sup.tab.1.Non.Fatal.MACE<-rbind(Sup.tab.df.1,Sup.tab.df.2,Sup.tab.df.3,Sup.tab.df.4)

#Fatal MACE
#VAT Grs
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav, data=base_gea_lexit_F_MACE);summary(m0)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_F_MACE);summary(m1)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_F_MACE);summary(m2)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_F_MACE);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.1<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 cm\u00b2"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (Scale)
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav), data=base_gea_lexit_F_MACE);summary(m0)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_F_MACE);summary(m1)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_F_MACE);summary(m2)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(gav)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_F_MACE);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.2<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 Standard Deviation\n(66.5 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (p75)
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2), data=base_gea_lexit_F_MACE)
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_F_MACE)
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_F_MACE)
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+strata(age.cat)+strata(sexo)+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_F_MACE)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.3<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,5],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,5],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,5],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Middle-Tertile\n(140-194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")


Sup.tab.df.4<-rbind(
  c(summ1a$coefficients[2,2], summ1a$conf.int[2,3:4],round(summ1a$coefficients[2,5],4)), 
  c(summ1b$coefficients[2,2], summ1b$conf.int[2,3:4],round(summ1b$coefficients[2,5],4)), 
  c(summ1c$coefficients[2,2], summ1c$conf.int[2,3:4],round(summ1c$coefficients[2,5],4)),
  c(summ1d$coefficients[2,2], summ1d$conf.int[2,3:4],round(summ1d$coefficients[2,5],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Upper-Tertile\n(>194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

Sup.tab.1.Fatal.MACE<-rbind(Sup.tab.df.1,Sup.tab.df.2,Sup.tab.df.3,Sup.tab.df.4)
Sup.tab.1<-rbind(Sup.tab.1.OVERALL.MACE,Sup.tab.1.Non.Fatal.MACE,Sup.tab.1.Fatal.MACE)
Table1_Flex<-flextable::align(flextable::flextable(Sup.tab.1,cwidth=4),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table1_Flex,path="Supplementary_Table_3.docx")


#####Analysis: VIF Models for Cox Models (Supplementary Figure 4)#####

#Overall MACE
#Continous Models
#Unadjusted
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo), data=base_gea_lexit);summary(m1)
survminer::ggcoxzph(cox.zph(m0))[1]
summary(m0)
BIC(m0)

#Model 1
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit);summary(m1)
survminer::ggcoxzph(cox.zph(m1))[1]
summary(m1)
BIC(m1)
Sup.Fig.3.1.1<-as.data.frame(car::vif(m1))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 2")


#Model 2
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit);summary(m2)
summary(m2)
BIC(m2)
survminer::ggcoxzph(cox.zph(m2))[1]
Sup.Fig.3.1.2<-as.data.frame(car::vif(m2))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High Apoliprotein-B"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 3")

#Model 3
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit);summary(m3)
summary(m3)
BIC(m3)
survminer::ggcoxzph(cox.zph(m3))[1]
Sup.Fig.3.1.3<-as.data.frame(car::vif(m3))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High-Apoliprotein-B",
                  "High-HOMA-IR",
                  "High-ADIPO-IR",
                  "Low-Adiponectin",
                  "High-SAT"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 4")

#Categorical
#Unadjusted
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo), data=base_gea_lexit);summary(m1)
survminer::ggcoxzph(cox.zph(m0))[1]
summary(m0)
BIC(m0)

#Model 1
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit);summary(m1)
survminer::ggcoxzph(cox.zph(m1))[1]
summary(m1)
BIC(m1)
Sup.Fig.3.1.4<-as.data.frame(car::vif(m1))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 2")

#Model 2
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit);summary(m2)
summary(m2)
BIC(m2)
survminer::ggcoxzph(cox.zph(m2))[1]
Sup.Fig.3.1.5<-as.data.frame(car::vif(m2))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High Apoliprotein-B"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 3")

#Model 3
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit);summary(m3)
summary(m3)
BIC(m3)
survminer::ggcoxzph(cox.zph(m3))[1]
Sup.Fig.3.1.6<-as.data.frame(car::vif(m3))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High-Apoliprotein-B",
                  "High-HOMA-IR",
                  "High-ADIPO-IR",
                  "Low-Adiponectin",
                  "High-SAT"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 4")

#Non-Fatal MACE
#Overall MACE
#Continous Models
#Unadjusted
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+strata(age.cat)+strata(sexo), data=base_gea_lexit_NF_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m0))[1]
summary(m0)
BIC(m0)

#Model 1
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_NF_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m1))[1]
summary(m1)
BIC(m1)
Sup.Fig.3.2.1<-as.data.frame(car::vif(m1))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 2")


#Model 2
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_NF_MACE);summary(m2)
summary(m2)
BIC(m2)
survminer::ggcoxzph(cox.zph(m2))[1]
Sup.Fig.3.2.2<-as.data.frame(car::vif(m2))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High Apoliprotein-B"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 3")

#Model 3
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_NF_MACE);summary(m3)
summary(m3)
BIC(m3)
survminer::ggcoxzph(cox.zph(m3))[1]
Sup.Fig.3.2.3<-as.data.frame(car::vif(m3))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High-Apoliprotein-B",
                  "High-HOMA-IR",
                  "High-ADIPO-IR",
                  "Low-Adiponectin",
                  "High-SAT"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 4")

#Categorical
#Unadjusted
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+strata(age.cat)+strata(sexo), data=base_gea_lexit_NF_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m0))[1]
summary(m0)
BIC(m0)

#Model 1
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_NF_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m1))[1]
summary(m1)
BIC(m1)
Sup.Fig.3.2.4<-as.data.frame(car::vif(m1))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 2")

#Model 2
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_NF_MACE);summary(m2)
summary(m2)
BIC(m2)
survminer::ggcoxzph(cox.zph(m2))[1]
Sup.Fig.3.2.5<-as.data.frame(car::vif(m2))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High Apoliprotein-B"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 3")

#Model 3
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_NF_MACE);summary(m3)
summary(m3)
BIC(m3)
survminer::ggcoxzph(cox.zph(m3))[1]
Sup.Fig.3.2.6<-as.data.frame(car::vif(m3))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High-Apoliprotein-B",
                  "High-HOMA-IR",
                  "High-ADIPO-IR",
                  "Low-Adiponectin",
                  "High-SAT"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 4")

#Fatal MACE
#Overall MACE
#Continous Models
#Unadjusted
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav, data=base_gea_lexit_F_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m0))[1]
summary(m0)
BIC(m0)

#Model 1
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_F_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m1))[1]
summary(m1)
BIC(m1)
Sup.Fig.3.3.1<-as.data.frame(car::vif(m1))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 2")


#Model 2
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_F_MACE);summary(m2)
summary(m2)
BIC(m2)
survminer::ggcoxzph(cox.zph(m2))[1]
Sup.Fig.3.3.2<-as.data.frame(car::vif(m2))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High Apoliprotein-B"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 3")

#Model 3
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~gav+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_F_MACE);summary(m3)
summary(m3)
BIC(m3)
survminer::ggcoxzph(cox.zph(m3))[1]
Sup.Fig.3.3.3<-as.data.frame(car::vif(m3))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High-Apoliprotein-B",
                  "High-HOMA-IR",
                  "High-ADIPO-IR",
                  "Low-Adiponectin",
                  "High-SAT"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Continuous Model\nModel 4")

#Categorical
#Unadjusted
m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2), data=base_gea_lexit_F_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m0))[1]
summary(m0)
BIC(m0)

#Model 1
m1<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+ BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina, data=base_gea_lexit_F_MACE);summary(m1)
survminer::ggcoxzph(cox.zph(m1))[1]
summary(m1)
BIC(m1)
Sup.Fig.3.3.4<-as.data.frame(car::vif(m1))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 2")

#Model 2
m2<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base_gea_lexit_F_MACE);summary(m2)
summary(m2)
BIC(m2)
survminer::ggcoxzph(cox.zph(m2))[1]
Sup.Fig.3.3.5<-as.data.frame(car::vif(m2))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High Apoliprotein-B"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 3")

#Model 3
m3<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~factor(gav_terciles_neav_2)+edad+sexo+seguimiento_anos+edad_eac+
            BMI_CAT+tabaquismo_cat+dm+hta+Estatinas_REC+anticoagulantes_orales+aspirina+
            gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            rita_p25_neav+adipo_p25_neav+gasp75, data=base_gea_lexit_F_MACE);summary(m3)
summary(m3)
BIC(m3)
survminer::ggcoxzph(cox.zph(m3))[1]
Sup.Fig.3.3.6<-as.data.frame(car::vif(m3))%>%
  janitor::clean_names()%>%
  cbind("names"=c("VAT",
                  "Age",
                  "Sex",
                  "Follow-up Time",
                  "Age at First pCAD",
                  "BMI Categories",
                  "Smoking Status",
                  "Diabetes",
                  "Hypertension",
                  "Statin Use",
                  "Anticoagulations Use",
                  "Aspirin Use",
                  "High Glucose",
                  "High Tryglycerides",
                  "High-LDL-C",
                  "Low-HDL-C",
                  "Non-HDL-C",
                  "High-Apoliprotein-B",
                  "High-ADIPO-IR",
                  "Low-Adiponectin",
                  "High-SAT"))%>%
  ggplot(aes(names,gvif_1_2_df))+
  geom_bar(stat='identity')+
  coord_flip()+
  scale_y_continuous(limits = c(0,5))+ 
  geom_hline(yintercept = 4,col="red",size=1.5)+
  xlab("Variables")+
  ylab("Variance Inflation Factor")+
  theme_classic()+
  ggtitle("Categorical Model\nModel 4")

Sup.Fig.3A.1<-ggarrange(Sup.Fig.3.1.1,Sup.Fig.3.1.2,Sup.Fig.3.1.3,ncol = 1,nrow = 3)
Sup.Fig.3A.2<-ggarrange(Sup.Fig.3.1.4,Sup.Fig.3.1.5,Sup.Fig.3.1.6,ncol = 1,nrow = 3)

Sup.Fig.3B.1<-ggarrange(Sup.Fig.3.2.1,Sup.Fig.3.2.2,Sup.Fig.3.2.3,ncol = 1,nrow = 3)
Sup.Fig.3B.2<-ggarrange(Sup.Fig.3.2.4,Sup.Fig.3.2.5,Sup.Fig.3.2.6,ncol = 1,nrow = 3)

Sup.Fig.3C.1<-ggarrange(Sup.Fig.3.3.1,Sup.Fig.3.3.2,Sup.Fig.3.3.3,ncol = 1,nrow = 3)
Sup.Fig.3C.2<-ggarrange(Sup.Fig.3.3.4,Sup.Fig.3.3.5,Sup.Fig.3.3.6,ncol = 1,nrow = 3)

Sup.Fig.3<-ggarrange(Sup.Fig.3A.1,Sup.Fig.3A.2,
          Sup.Fig.3B.1,Sup.Fig.3B.2,
          Sup.Fig.3C.1,Sup.Fig.3C.2,
          ncol = 6,nrow = 1)

ggsave(Sup.Fig.3,
       filename = "Supplementary_Figure_3.png", 
       bg = "white",
       width = 70, 
       height = 25,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)

#####Analysis: Propensity Matched Scoring######

set.seed(123)
match<-matchit(morbimortalidad_mace_2~edad+sexo+edad_eac+dm+hta+Estatinas_REC, data=base,method = "full", estimand = "ATT")
summary(match)
base2<-match.data(match,distance = "ps")
nrow(base2)
#plot(match)


#####Analysis: Descriptive characteristics of the study sample (Propensity Matched Base) (Supplementary Table 5)#####

base2 %>% 
  dplyr::select(gav_terciles_neav,
                edad,sexo,seguimiento_anos,edad_eac,
                imc,BMI_CAT,
                tabaquismo_cat,dm,hta,
                Estatinas_REC,anticoagulantes_orales,aspirina,
                glucosa,gluc_100,
                tg,tg_cat,
                ct,
                cldl,cldl_cat,
                chdl,chdl_cat,
                cno_hdl,nhdl_cat,
                apo_b,apob_p75_neav,
                pcr,pcr_cat,
                homa,homa_p75_neav,
                rita,rita_p25_neav,
                adiponectina,adipo_p25_neav,
                gas,gasp75,
                gat,gatp75,
                gav,
                morbimortalidad_mace,mace_no_fatal,mace_fatal)%>%
  tbl_summary(by = gav_terciles_neav,
              missing = "ifany")%>%
  bold_labels()%>%
  add_overall()%>%
  add_p()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Supplementary_Table_5.docx")


#####Analysis: Cox-Models to Evaluate the Association with VAT and MACE (Propensity Matched Base) (Supplementary Table 6)#####

##OVERALL MACE###
#VAT Grs
m0<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~gav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m0)
m1<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~gav+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m1)
m2<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~gav+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m2)
m3<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~gav+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.1<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 cm\u00b2"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (Scale)
m0<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~scale(gav), data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m0)
m1<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~scale(gav)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m1)
m2<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~scale(gav)+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m2)
m3<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~scale(gav)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.2<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 Standard Deviation\n(66.5 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (p75)
m0<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~factor(gav_terciles_neav_2), data=base2,robust = TRUE, weights = weights, cluster = subclass)
m1<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~factor(gav_terciles_neav_2)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass)
m2<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~factor(gav_terciles_neav_2)+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass)
m3<-coxph(Surv(seguimiento_anos,morbimortalidad_mace_2)~factor(gav_terciles_neav_2)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.3<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Middle-Tertile\n(140-194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")


Sup.tab.df.4<-rbind(
  c(summ1a$coefficients[2,2], summ1a$conf.int[2,3:4],round(summ1a$coefficients[2,6],4)), 
  c(summ1b$coefficients[2,2], summ1b$conf.int[2,3:4],round(summ1b$coefficients[2,6],4)), 
  c(summ1c$coefficients[2,2], summ1c$conf.int[2,3:4],round(summ1c$coefficients[2,6],4)),
  c(summ1d$coefficients[2,2], summ1d$conf.int[2,3:4],round(summ1d$coefficients[2,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Upper-Tertile\n(>194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

Sup.tab.1.OVERALL.MACE<-rbind(Sup.tab.df.1,Sup.tab.df.2,Sup.tab.df.3,Sup.tab.df.4)


##Non_Fatal MACE
#VAT Grs
m0<-coxph(Surv(seguimiento_anos,mace_no_fatal)~gav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m0)
m1<-coxph(Surv(seguimiento_anos,mace_no_fatal)~gav+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m1)
m2<-coxph(Surv(seguimiento_anos,mace_no_fatal)~gav+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m2)
m3<-coxph(Surv(seguimiento_anos,mace_no_fatal)~gav+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.1<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 cm\u00b2"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (Scale)
m0<-coxph(Surv(seguimiento_anos,mace_no_fatal)~scale(gav), data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m0)
m1<-coxph(Surv(seguimiento_anos,mace_no_fatal)~scale(gav)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m1)
m2<-coxph(Surv(seguimiento_anos,mace_no_fatal)~scale(gav)+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m2)
m3<-coxph(Surv(seguimiento_anos,mace_no_fatal)~scale(gav)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.2<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 Standard Deviation\n(66.5 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (p75)
m0<-coxph(Surv(seguimiento_anos,mace_no_fatal)~factor(gav_terciles_neav_2), data=base2,robust = TRUE, weights = weights, cluster = subclass)
m1<-coxph(Surv(seguimiento_anos,mace_no_fatal)~factor(gav_terciles_neav_2)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass)
m2<-coxph(Surv(seguimiento_anos,mace_no_fatal)~factor(gav_terciles_neav_2)+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass)
m3<-coxph(Surv(seguimiento_anos,mace_no_fatal)~factor(gav_terciles_neav_2)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.3<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Middle-Tertile\n(140-194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")


Sup.tab.df.4<-rbind(
  c(summ1a$coefficients[2,2], summ1a$conf.int[2,3:4],round(summ1a$coefficients[2,6],4)), 
  c(summ1b$coefficients[2,2], summ1b$conf.int[2,3:4],round(summ1b$coefficients[2,6],4)), 
  c(summ1c$coefficients[2,2], summ1c$conf.int[2,3:4],round(summ1c$coefficients[2,6],4)),
  c(summ1d$coefficients[2,2], summ1d$conf.int[2,3:4],round(summ1d$coefficients[2,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Upper-Tertile\n(>194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

Sup.tab.1.Non.Fatal.MACE<-rbind(Sup.tab.df.1,Sup.tab.df.2,Sup.tab.df.3,Sup.tab.df.4)

#Fatal MACE
#VAT Grs
m0<-coxph(Surv(seguimiento_anos,mace_fatal)~gav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m0)
m1<-coxph(Surv(seguimiento_anos,mace_fatal)~gav+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m1)
m2<-coxph(Surv(seguimiento_anos,mace_fatal)~gav+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m2)
m3<-coxph(Surv(seguimiento_anos,mace_fatal)~gav+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.1<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 cm\u00b2"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (Scale)
m0<-coxph(Surv(seguimiento_anos,mace_fatal)~scale(gav), data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m0)
m1<-coxph(Surv(seguimiento_anos,mace_fatal)~scale(gav)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m1)
m2<-coxph(Surv(seguimiento_anos,mace_fatal)~scale(gav)+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m2)
m3<-coxph(Surv(seguimiento_anos,mace_fatal)~scale(gav)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass);summary(m3)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.2<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("+1 Standard Deviation\n(66.5 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

#VAT (p75)
m0<-coxph(Surv(seguimiento_anos,mace_fatal)~factor(gav_terciles_neav_2), data=base2,robust = TRUE, weights = weights, cluster = subclass)
m1<-coxph(Surv(seguimiento_anos,mace_fatal)~factor(gav_terciles_neav_2)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav, data=base2,robust = TRUE, weights = weights, cluster = subclass)
m2<-coxph(Surv(seguimiento_anos,mace_fatal)~factor(gav_terciles_neav_2)+homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass)
m3<-coxph(Surv(seguimiento_anos,mace_fatal)~factor(gav_terciles_neav_2)+gluc_100+tg_cat+cldl_cat+chdl_cat+nhdl_cat+apob_p75_neav+
            homa_p75_neav+rita_p25_neav+adipo_p25_neav+gasp75, data=base2,robust = TRUE, weights = weights, cluster = subclass)


summ1a <- summary(m0); summ1b <- summary(m1); summ1c <- summary(m2); summ1d <- summary(m3)

Sup.tab.df.3<-rbind(
  c(summ1a$coefficients[1,2], summ1a$conf.int[1,3:4],round(summ1a$coefficients[1,6],4)), 
  c(summ1b$coefficients[1,2], summ1b$conf.int[1,3:4],round(summ1b$coefficients[1,6],4)), 
  c(summ1c$coefficients[1,2], summ1c$conf.int[1,3:4],round(summ1c$coefficients[1,6],4)),
  c(summ1d$coefficients[1,2], summ1d$conf.int[1,3:4],round(summ1d$coefficients[1,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Middle-Tertile\n(140-194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")


Sup.tab.df.4<-rbind(
  c(summ1a$coefficients[2,2], summ1a$conf.int[2,3:4],round(summ1a$coefficients[2,6],4)), 
  c(summ1b$coefficients[2,2], summ1b$conf.int[2,3:4],round(summ1b$coefficients[2,6],4)), 
  c(summ1c$coefficients[2,2], summ1c$conf.int[2,3:4],round(summ1c$coefficients[2,6],4)),
  c(summ1d$coefficients[2,2], summ1d$conf.int[2,3:4],round(summ1d$coefficients[2,6],4))) %>%
  as.data.frame() %>% `colnames<-`(c("mean", "lower", "upper","p-value")) %>% 
  mutate("Cause"=c("Upper-Tertile\n(>194 cm\u00b2)"),
         "Model"=c(paste0("Model ",rep(1:4,1)))) %>% 
  arrange(Cause) %>% 
  mutate("HR"=paste0(sprintf("%#.3f", mean), " (",sprintf("%#.3f", lower), "-", sprintf("%#.3f", upper), ")"))%>%
  dplyr::select(Model,Cause,HR,"p-value")

Sup.tab.1.Fatal.MACE<-rbind(Sup.tab.df.1,Sup.tab.df.2,Sup.tab.df.3,Sup.tab.df.4)
Sup.tab.1<-rbind(Sup.tab.1.OVERALL.MACE,Sup.tab.1.Non.Fatal.MACE,Sup.tab.1.Fatal.MACE)
Table1_Flex<-flextable::align(flextable::flextable(Sup.tab.1,cwidth=4),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table1_Flex,path="Supplementary_Table_6.docx")




##Cox Proportional Hazard Regression Models
#Univariate Association
base_gea_lexit$ldl_no_hdl<-base_gea_lexit$cldl/base_gea_lexit$cno_hdl
base_gea_lexit$apo_b
m1<-lm(apo_b~ldl_no_hdl,data=base_gea_lexit)
base_gea_lexit$res_apo_b<-m1$residuals
scale(m1$residuals)

m0<-coxph(Surv(time_at_entry, time_at_exit, exit_status)~scale(res_apo_b)+ldl_no_hdl+strata(age.cat)+strata(sexo), data=base_gea_lexit);summary(m0)
cox.zph(m0); summary(m0); BIC(m0)


