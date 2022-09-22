
/*******************************************************************************
/*******************************************************************************
/*******************************************************************************
prepare variables on omega-3 for Jorick

MLVS FFQ (trans_avg omega3_avg omega6_avg ala_avg epa_avg dpa_avg dha_avg omega3_noala_avg)
  trans_ffq1   a_trn07_fs_ffq1           Num     8  Antilog Total Trans, gm, Sacks, 2007
  omega3_ffq1  a_pfn307_fs_ffq1          Num     8  Antilog Total Omega 3, gms, Sacks, 2007
  omega6_ffq1  a_pfn607_fs_ffq1          Num     8  Antilog Total Omega 6, gm, Sacks, 2007
  ala_ffq1     a_pfa183n3c07_fs_ffq1     Num     8  Antilog Alpha Linolenic Acid, gms, Sacks, 2007
  epa_ffq1     a_pf205n3c07_fs_ffq1      Num     8  Antilog Eicosapentaenoic EPA fatty acid, gms, Sacks P10, 2007
  dpa_ffq1     a_pf225n3c07_fs_ffq1      Num     8  Antilog Docosapentaenoic fatty acid, DPA, gms, Sacks P13, 2007
  dha_ffq1     a_pf226n3c07_fs_ffq1      Num     8  Antilog Docosahexaenoic fatty acid, DHA, gms, Sacks P14, 2007
  omega3_noala_ffq1 total omega-3 excluding ALA

MLVS DDR
note: no total omega-6, no ALA
  a_f205_fs_dr_w1avg    Num    8 energy-adjusted one week average of pufa 20:5 (eicosapentaenoic acid [epa]) - g, from both food and supp.
  a_f226_fs_dr_w1a      Num    8 energy-adjusted one week average of pufa 22:6 (docosahexaenoic acid [dha]) - g, from both food and supp.
  a_p22_5_fs_dr_w1avg   Num    8 energy-adjusted one week average of pufa 22:5 (docosapentaenoic acid [dpa]) - g, from both food and supp.
  a_trn07_fo_dr_w1avg   Num    8 energy-adjusted one week average of total trans-fatty acids (trans) - g, from food
  a_omega3_fs_dr_w1avg8 Num    8 energy-adjusted one week average of omega-3 fatty acids - g, from both food and supp.

cohort FFQ (ala10v epa10v dha10v dpa10v trans10v omega610v omega310v omega3_noala10v)
86-98
  ALA      EPA      DHA      DPA      total trans   total w3  total w6
  f183s86a f205s86a f226s86a          trnss86a                n6s86a
  f183s90a f205s90a f226s90a          trnss90a                n6s90a
  f183s94a f205s94a f226s94a          trnss94a                n6s94a
  a183098a pf20598a pf22698a pf22598a trn0098a      pfn3098a  n6s0098a

02-10
  pfa183n3c0202a pf205n3c0202a pf226n3c0202a pf225n3c0202a trn0202a pfn30202a n60202a
  pfa183n3c0706a pf205n3c0706a pf226n3c0706a pf225n3c0706a trn0706a pfn30706a n60706a
  pfa183n3c1110a pf205n3c1110a pf226n3c1110a pf225n3c1110a trn1110a pfn31110a n61110a

********************************************************************************
********************************************************************************
*******************************************************************************/

/*FOR EVERY PROJECT*/
libname hpfsfmt '/proj/hpsass/hpsas00/formats';
*libname nhsfmt '/proj/nhsass/nhsas00/formats';
*filename nhstools '/proj/nhsass/nhsas00/nhstools/sasautos/';
filename hpstools '/proj/hpsass/hpsas00/hpstools/sasautos';
filename channing '/usr/local/channing/sasautos/';
options mautosource sasautos=(channing nhstools hpstools);
options fmtsearch=(nhsfmt hpfsfmt) nofmterr nocenter nonumber nodate formdlim=' ';
options mautosource sasautos=(channing /*nhstools*/ hpstools);
options nocenter linesize=130 pagesize=150;


libname DDR '/proj/polvss/polvs00/MLVS/data_for_analysis/diet_7ddr/sas_data';
libname NTS '/proj/polvss/polvs00/MLVS/data_for_analysis/diet_7ddr/sas_data';
libname FFQ '/proj/polvss/polvs00/MLVS/data_for_analysis/diet_ffq/sas_data';
libname STOOLQQ '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_questionnaire/sas_data';
libname SPECIES '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/taxonomy';
libname DATE '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_sample_details/sas_data';
libname ECD '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_analysis/sas_data/june2017/ecdna';
libname ECR '/udd/nhwma/mlvs/reference/201710_galeb/non-pathway/created_with_jlp_script';
libname crp '/udd/nhwma/mlvs/pipeline/crp';
libname blood '/proj/polvss/polvs00/MLVS/data_for_analysis/blood_analysis/sas_data/';

%include '/proj/polvss/polvs00/MLVS/data_for_analysis/diet_7ddr/formats/NDS_format.sas';
%include '/proj/polvss/polvs00/MLVS/data_for_analysis/fecal_questionnaire/formats/fecal_qu_format.sas';

/*******************************************************************************
/*******************************************************************************
/*******************************************************************************
*                              READ IN VARIABLES                               *
********************************************************************************
*******************************************************************************/
%hp_der (keep=id agemlvs dbmy09 age86 rtmnyr86 rtmnyr87 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96
rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 smoke86 smoke88 smoke90 smoke92 smoke94
smoke96 smoke98 smoke00 smoke02 smoke04 smoke06 wt86 wt88 wt90 wt92 wt94 wt96 wt98 wt00
wt02 wt04 wt06 bmi86 bmi88 bmi90 bmi92 bmi94 bmi96 bmi98 bmi00 bmi02 bmi04 bmi06
hbp86 hbp88 hbp90 db86 db88 db90f chol86 chol88 chol90 mi86 mi88 mi90 cabg86
cabg88 cabg90 ang86 ang88 ang90 str86 str88 str90 can86 can88 can90 height);
agemlvs=(1345 /*31 December 2012*/ - dbmy09)/12;
run;

%hp_der_2 (keep=id rtmnyr08 rtmnyr10 smoke08 smoke10 wt08 wt10 bmi08 bmi10);
run;

%hp86 (keep=id /*no mvt*/ asp86 motrn86 tag86 alcf86 seuro86 scand86 ocauc86 afric86 asian86 oanc86);
run;

%hp88 (keep=id mvt88 asp88 motrn88 tag88 alc88);
run;

%hp90 (keep=id /*no mvt*/ asp90 motrn90 tag90); /*no alc as # of days, only #years with >3/day*/
run;

/*after 1990, start keeping diagnoses again since they weren't in hp_der*/

%hp92 (keep=id mvt92 hbp92 db92 chol92 mi92 ang92 cabg92 str92 canc92 asp92
  motrn92 tag92); /*alc asked in unusual ways again: largest # drinks/more than 50 in life*/
run;

%hp94 (keep=id /*no mvt*/ hbp94 db94 chol94 mi94 ang94 cabg94 str94 canc94
  asp94 motrn94 tag94 /*no alc*/);
run;

%hp96 (keep=id mvt96 hbp96 db96 chol96 mi96 ang96 cabg96 str96 canc96 nasp96
  motrn96 tag96);
run;

%hp98 (keep=id /*mvt*/ hbp98 db98 chol98 mi98 ang98 cabg98 str98 canc98 aspw98
  motrn98 nsaid98 /*new*/ tag98);
run;


%hp00 (keep=id mvt00 hbp00 db00 chol00 mi00 ang00 cabg00 strk00 canc00 asp00
  motrn00 nsaid00 tag00);
run;

%hp02 (keep=id mvt02 hbp02 db02 chol02 mi02 ang02 cabg02 strk02 canc02 asp02
  mtrn02 /*no nsaid*/ tag02);
run;

/*start of PPI querying*/
%hp04 (keep=id flag04 mvt04 hbp04 db04 chol04 mi04 ang04 cabg04 strk04 canc04 asp04
  mtrn04 tag04 pril04);
run;

%hp06 (keep=id flag06 mvt06 hbp06 db06 chol06 mi06 ang06 cabg06 strk06 canc06 asp06
  mtrn06 tag06 pril06);
run;

%hp08 (keep=id flag08 mvt08 hbp08 db08 chol08 mi08 ang08 cabg08 strk08 canc08 asp08
  mtrn08 tag08 pril08);
run;

%hp10 (keep=id flag10 mvt10 hbp10 db10 chol10 cabg10 mi10 ang10 strk10 canc10 asp10
  mtrn10 tag10 pril10);

%hp12 (keep=id smk12 wt12 ncig12);
run;


/*proc print;var id id1; where in_key=1;endsas;*/


data hpfsvars;
  merge hp_der hp_der_2 hp86 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04 hp06 hp08 hp10 hp12;
  by id;
  if first.id;
run;


proc sort; by id;
run;

/*physical activity*/
%hmet8616(keep=id act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16);


proc means data=hmet8616 n nmiss min max;
run;


data hp_mets;
   set hmet8616;
array acts{*} act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16;

do i=1 to DIM(acts);
   if acts{i}=998 or acts{i}=999 then acts{i}=.;
end;

do k=DIM(acts) to 2 by -1;
    if acts{k} eq . and acts{k-1} ne . then acts{k}=acts{k-1};
end;
run;

data mets;
    set hp_mets;
run;

proc sort data=mets out=mets;
  by id;
run;

/*
proc means data=mets n nmiss mean median min max;
run;
*/


/*nutrients/alcohol*/
%h86_nts(keep= id calor86n alco86n);
  run;
%h90_nts(keep= id calor90n alco90n);
  run;
%h94_nts(keep= id calor94n alco94n);
  run;
%h98_nts(keep= id calor98n alco98n);
  run;
%h02_nts(keep= id calor02n alco02n);
  run;
%h06_nts(keep= id calor06n alco06n);
  run;
%h10_nts(keep= id calor10n alco10n);
  run;

%h86_ant(keep= id aofib86a frtaf86a ceraf86a vegaf86a f183s86a f205s86a f226s86a          trnss86a n6s86a);
  run;
%h90_ant(keep= id aofib90a frtaf90a ceraf90a vegaf90a f183s90a f205s90a f226s90a          trnss90a n6s90a);
  run;
%h94_ant(keep= id aofib94a frtaf94a ceraf94a vegaf94a f183s94a f205s94a f226s94a          trnss94a n6s94a);
  run;
%h98_ant(keep= id aofib98a frtaf98a ceraf98a vegaf98a a183098a pf20598a pf22698a pf22598a trn0098a pfn3098a n6s0098a);
  run;
%h02_ant(keep= id aofib02a frtaf02a ceraf02a vegaf02a pfa183n3c0202a pf205n3c0202a pf226n3c0202a pf225n3c0202a trn0202a pfn30202a n60202a);
  run;
%h06_ant(keep= id aofib06a frtaf06a ceraf06a vegaf06a pfa183n3c0706a pf205n3c0706a pf226n3c0706a pf225n3c0706a trn0706a pfn30706a n60706a);
  run;
%h10_ant(keep= id aofib10a frtaf10a ceraf10a vegaf10a pfa183n3c1110a pf205n3c1110a pf226n3c1110a pf225n3c1110a trn1110a pfn31110a n61110a);
  run;



data nutrients;
  merge h86_nts h90_nts h94_nts h98_nts h02_nts h06_nts h10_nts h86_ant h90_ant h94_ant h98_ant h02_ant h06_ant h10_ant;
  by id;
  if first.id;
run;


proc sort data=nutrients out=nutrients;
  by id;
run;


/**********************************
/**********************************
*             MLVS 7DDR           *
***********************************
**********************************/

/*
DDR nutrients
Only chose the ones that were previously targeted for the "nutrient" study
This recoding step is because the nutrient variables were renamed at some point so instead of changing
all the downstream code, just change it here.
*/



proc contents data = NTS.dr_nts_wk1_mean_e_adjusted;
run;

data ntsw1;
  set NTS.dr_nts_wk1_mean_e_adjusted;
  a_vpro_fo_dr_w1avg        =         a_vprot_fo_dr_w1avg;
  a_star_fo_dr_w1avg        =         a_st_fo_dr_w1avg;
  a_vk_fo_dr_w1avg          =         a_vitk_fo_dr_w1avg;
  a_sfa_fo_dr_w1avg         =         a_satfat_fo_dr_w1avg;
  a_sfa_fs_dr_w1avg         =         a_satfat_fs_dr_w1avg;
  a_mfa_fo_dr_w1avg         =         a_monfat_fo_dr_w1avg;
  a_mfa_fs_dr_w1avg         =         a_monfat_fs_dr_w1avg;
  a_pfa_fo_dr_w1avg         =         a_poly_fo_dr_w1avg;
  a_pfa_fs_dr_w1avg         =         a_poly_fs_dr_w1avg;
  a_vd_fo_dr_w1avg          =         a_vitd_fo_dr_w1avg;
  a_vd_fs_dr_w1avg          =         a_vitd_fs_dr_w1avg;

run;


data ntsw1;
  set ntsw1
  (keep =

  /*PROTEIN*/
  a_pro_fo_dr_w1avg
  a_pro_fs_dr_w1avg
  a_vpro_fo_dr_w1avg
  a_aprot_fo_dr_w1avg

  /*FAT*/
  a_fat_fo_dr_w1avg
  a_fat_fs_dr_w1avg

  /*CARBOHYDRATES*/
  a_carbo_fo_dr_w1avg
  a_carbo_fs_dr_w1avg

  /*SUCROSE*/
  a_sucr_fo_dr_w1avg
  a_sucr_fs_dr_w1avg

  /*FRUCTOSE*/
  a_fruct_fo_dr_w1avg
  a_fruct_fs_dr_w1avg

  /*LACTOSE*/
  a_lact_fo_dr_w1avg

  /*STARCH*/
  a_star_fo_dr_w1avg

  /*GLUCOSE*/
  a_glu_fo_dr_w1avg
  a_glu_fs_dr_w1avg

  /*FIBER*/
  a_aofib_fo_dr_w1avg
  a_aofib_fs_dr_w1avg

  /*pectin fiber*/
  a_pect_fo_dr_w1avg

  /*insoluble fiber*/
  a_ifib_fo_dr_w1avg /*from foods*/
  a_ifib_fs_dr_w1avg /*from foods and supplements*/

  /*soluble fiber*/
  a_wsdf_fo_dr_w1avg /*from foods*/

  /*CALCIUM FROM FOOD*/
  a_calc_fo_dr_w1avg

  /*CALCIUM FROM BOTH*/
  a_calc_fs_dr_w1avg

  /*IRON TOTAL (NOT JUST HEME)*/
  a_iron_fo_dr_w1avg
  a_iron_fs_dr_w1avg

  /*FOLATE*/
  a_dfe_fo_dr_w1avg
  a_dfe_fs_dr_w1avg
  a_fdfol_fo_dr_w1avg
  a_fol98_fo_dr_w1avg
  a_folic_fo_dr_w1avg
  a_folic_fs_dr_w1avg

  /*VITAMIN B12*/
  a_b12_fo_dr_w1avg
  a_b12_fs_dr_w1avg

  /*CHOLINE*/
  a_choline_fo_dr_w1avg
  a_choline_fs_dr_w1avg

  /*VITAMIN K*/
  a_vk_fo_dr_w1avg

  /*SATURATED FATTY ACIDS*/
  a_sfa_fo_dr_w1avg
  a_sfa_fs_dr_w1avg

  /*MONOUNSATURATED FATTY ACIDS*/
  a_mfa_fo_dr_w1avg
  a_mfa_fs_dr_w1avg

  /*POLY UNSATURATED FATTY ACIDS*/
  a_pfa_fo_dr_w1avg
  a_pfa_fs_dr_w1avg

  /*ALCOHOL*/
  a_alco_fo_dr_w1avg

  /*VITAMIN D*/
  a_vd_fo_dr_w1avg
  a_vd_fs_dr_w1avg

  /*omega-3*/
  a_omega3_fs_dr_w1avg
  a_f205_fs_dr_w1avg /*EPA*/
  a_p22_5_fs_dr_w1avg  /*DPA*/
  a_f226_fs_dr_w1avg /*DHA*/

  /*TRANS FAT*/
  a_trn07_fo_dr_w1avg


  id)
  ;
  proc sort;
    by id;
run;

/* proc contents data = ntsw1;
run; */

data ntsw2;
  set NTS.dr_nts_wk2_mean_e_adjusted;
  a_vpro_fo_dr_w2avg        =         a_vprot_fo_dr_w2avg;
  a_star_fo_dr_w2avg        =         a_st_fo_dr_w2avg;
  a_vk_fo_dr_w2avg          =         a_vitk_fo_dr_w2avg;
  a_sfa_fo_dr_w2avg         =         a_satfat_fo_dr_w2avg;
  a_sfa_fs_dr_w2avg         =         a_satfat_fs_dr_w2avg;
  a_mfa_fo_dr_w2avg         =         a_monfat_fo_dr_w2avg;
  a_mfa_fs_dr_w2avg         =         a_monfat_fs_dr_w2avg;
  a_pfa_fo_dr_w2avg         =         a_poly_fo_dr_w2avg;
  a_pfa_fs_dr_w2avg         =         a_poly_fs_dr_w2avg;
  a_vd_fo_dr_w2avg          =         a_vitd_fo_dr_w2avg;
  a_vd_fs_dr_w2avg          =         a_vitd_fs_dr_w2avg;
run;

data ntsw2;
  set ntsw2

  (keep =
  /*PROTEIN*/
  a_pro_fo_dr_w2avg
  a_pro_fs_dr_w2avg
  a_vpro_fo_dr_w2avg
  a_aprot_fo_dr_w2avg


  /*FAT*/
  a_fat_fo_dr_w2avg
  a_fat_fs_dr_w2avg

  /*CARBOHYDRATES*/
  a_carbo_fo_dr_w2avg
  a_carbo_fs_dr_w2avg

  /*SUCROSE*/
  a_sucr_fo_dr_w2avg
  a_sucr_fs_dr_w2avg

  /*FRUCTOSE*/
  a_fruct_fo_dr_w2avg
  a_fruct_fs_dr_w2avg

  /*LACTOSE*/
  a_lact_fo_dr_w2avg

  /*STARCH*/
  a_star_fo_dr_w2avg

  /*GLUCOSE*/
  a_glu_fo_dr_w2avg
  a_glu_fs_dr_w2avg

  /*FIBER*/
  a_aofib_fo_dr_w2avg
  a_aofib_fs_dr_w2avg

  /*pectin*/
  a_pect_fo_dr_w2avg

  /*insoluble fiber*/
  a_ifib_fo_dr_w2avg /*from foods*/
  a_ifib_fs_dr_w2avg /*from foods and supplements*/

  /*soluble fiber*/
  a_wsdf_fo_dr_w2avg /*from foods*/

  /*CALCIUM FROM FOOD*/
  a_calc_fo_dr_w2avg

  /*CALCIUM FROM BOTH*/
  a_calc_fs_dr_w2avg

  /*IRON TOTAL (NOT JUST HEME)*/
  a_iron_fo_dr_w2avg
  a_iron_fs_dr_w2avg

  /*FOLATE*/
  a_dfe_fo_dr_w2avg
  a_dfe_fs_dr_w2avg
  a_fdfol_fo_dr_w2avg
  a_fol98_fo_dr_w2avg
  a_folic_fo_dr_w2avg
  a_folic_fs_dr_w2avg

  /*VITAMIN B12*/
  a_b12_fo_dr_w2avg
  a_b12_fs_dr_w2avg

  /*CHOLINE*/
  a_choline_fo_dr_w2avg
  a_choline_fs_dr_w2avg

  /*VITAMIN K*/
  a_vk_fo_dr_w2avg

  /*SATURATED FATTY ACIDS*/
  a_sfa_fo_dr_w2avg
  a_sfa_fs_dr_w2avg

  /*MONOUNSATURATED FATTY ACIDS*/
  a_mfa_fo_dr_w2avg
  a_mfa_fs_dr_w2avg

  /*POLY UNSATURATED FATTY ACIDS*/
  a_pfa_fo_dr_w2avg
  a_pfa_fs_dr_w2avg

  /*ALCOHOL*/
  a_alco_fo_dr_w2avg

  /*VITAMIN D*/
  a_vd_fo_dr_w2avg
  a_vd_fs_dr_w2avg

  /*omega-3*/
  a_omega3_fs_dr_w2avg

  a_f205_fs_dr_w2avg /*EPA*/
  a_p22_5_fs_dr_w2avg /*DPA*/
  a_f226_fs_dr_w2avg /*DHA*/

  /*TRANS FAT*/
  a_trn07_fo_dr_w2avg

  id)
  ;

  proc sort;
    by id;
run;

/* proc contents data = ntsw2;
run; */

data ntsw12;
  merge ntsw1 ntsw2;
  by id;
  if first.id;
run;



proc sort data=ntsw12 out=ntsw12;
  by id;
run;

/*eliminate people in the MLVS who did not give stool*/
data ntsw12;
  set ntsw12;
  if id ne .;
run;


/*calorie counts from the unadjusted files*/
data energy_w1;
  set NTS.dr_nts_wk1_mean
  (keep = id calor_fs_dr_w1avg);
  proc sort;
    by id;
run;

data energy_w2;
  set NTS.dr_nts_wk2_mean
  (keep = id calor_fs_dr_w2avg);
  proc sort;
    by id;
run;

data energyw12;
  merge energy_w1 energy_w2;
  by id;
  if first.id;
run;


proc sort data=energyw12 out=energyw12;
  by id;
run;


/*eliminate people in the MLVS who did not give stool*/
data energyw12;
  set energyw12;
  if id ne .;
run;


/**********************************
/**********************************
*             MLVS FFQs           *
***********************************
**********************************/

/*REPLACE PATHS ONCE MLVS FFQS CLEANED*/

/*week 1 FFQ-nutrients*/

data ffq1;
	set FFQ.mlvsffq1_noexcl_ant;
	proc sort; by id;
run;

proc contents data=ffq1;
run;


/*renaming step because var names changed at some point*/
data ffq1;
  set ffq1;
  aofib             =         a_aofib_fs_ffq1;
  aprot             =         a_aprot_fs_ffq1;
  vprot             =         a_vprot_fs_ffq1;
  prot              =         a_prot_fs_ffq1;
  b1_sulf           =         a_b1_sulf_fs_ffq1;
  b1_sulf_wo        =         a_b1_sulf_fo_ffq1;
  biot_sulf         =         a_biot_sulf_fs_ffq1;
  cys_sulf          =         a_cys_sulf_fs_ffq1;
  meth_sulf         =         a_meth_sulf_fs_ffq1;
  prot_sulf         =         a_prot_sulf_fs_ffq1;
  sulfur            =         a_sulfur_fs_ffq1;
  ala               =         a_pfa183n3c07_fs_ffq1;
  epa               =         a_pf205n3c07_fs_ffq1;
  dpa               =         a_pf225n3c07_fs_ffq1;
  dha               =         a_pf226n3c07_fs_ffq1;
  omega3            =         a_pfn307_fs_ffq1; /*Total Omega 3*/
  omega6            =         a_pfn607_fs_ffq1; /*Total omega 6*/
  trans             =         a_trn07_fs_ffq1; /*total trans*/
run;


data ffq1;
  set ffq1
  (keep=
  id aofib aprot vprot prot b1_sulf b1_sulf_wo biot_sulf cys_sulf
  meth_sulf prot_sulf sulfur
  ala epa dpa dha omega3 omega6 trans)
  end=_end_;

  array oldsulf1{*}

  aofib aprot vprot prot b1_sulf b1_sulf_wo biot_sulf cys_sulf
  meth_sulf prot_sulf sulfur
  ala epa dpa dha omega3 omega6 trans;

  array newsulf1{*}

  aofib_ffq1 aprot_ffq1 vprot_ffq1 prot_ffq1 b1_sulf_ffq1 b1_sulf_wo_ffq1
  biot_sulf_ffq1 cys_sulf_ffq1 meth_sulf_ffq1 prot_sulf_ffq1 sulfur_ffq1
  ala_ffq1 epa_ffq1 dpa_ffq1 dha_ffq1 omega3_ffq1 omega6_ffq1 trans_ffq1;

  do i=1 to dim(oldsulf1);
    newsulf1{i}=oldsulf1{i};
    end;

  nonprot_sulf_ffq1=(sulfur_ffq1)-(prot_sulf_ffq1);
  omega3_noala_ffq1=omega3_ffq1-ala_ffq1;

run;



data ffq1;
  set ffq1 (keep=id aofib_ffq1 aprot_ffq1 vprot_ffq1 prot_ffq1 prot_sulf_ffq1 nonprot_sulf_ffq1 sulfur_ffq1 ala_ffq1 epa_ffq1 dpa_ffq1 dha_ffq1 omega3_ffq1 omega6_ffq1 trans_ffq1 omega3_noala_ffq1)
  end=_end_;

run;


data nts_ffq1;
  set FFQ.mlvsffq1_noexcl_ant;
  proc sort; by id;
run;

proc contents data = nts_ffq1;
run;

data nts_ffq1;
  set nts_ffq1;
  aprot             =           a_aprot_fs_ffq1;
  dprot             =           a_dprot_fs_ffq1;
  prot              =           a_prot_fs_ffq1;
  vprot             =           a_vprot_fs_ffq1;
  afat              =           a_afat_fs_ffq1;
  dfat              =           a_dfat_fs_ffq1;
  tfat              =           a_tfat_fs_ffq1;
  vfat              =           a_vfat_fs_ffq1;
  carbo             =           a_carbo_fs_ffq1;
  sucr              =           a_sucr_fs_ffq1;
  fruct             =           a_fruct_fs_ffq1;
  lact              =           a_lact_fs_ffq1;
  st                =           a_st_fs_ffq1;
  glu               =           a_glu_fs_ffq1;
  aofib             =           a_aofib_fs_ffq1;
  calc_wo           =           a_calc_fo_ffq1;
  calc              =           a_calc_fs_ffq1;
  iron              =           a_iron_fo_ffq1;
  iron_wo           =           a_iron_fs_ffq1;
  heme              =           a_heme_fs_ffq1;
  fdfol             =           a_fdfol_fs_ffq1;
  fol98             =           a_fol98_fs_ffq1;
  fol98_wo          =           a_fol98_fo_ffq1;
  folic             =           a_folic_fs_ffq1;
  dfe               =           a_dfe_fs_ffq1;
  b12               =           a_b12_fs_ffq1;
  b12_wo            =           a_b12_fo_ffq1;
  choline           =           a_choline_fs_ffq1;
  choline_wo        =           a_choline_fo_ffq1;
  dvitk             =           a_dvitk_fs_ffq1;
  vitk              =           a_vitk_fs_ffq1;
  vitk2_m4          =           a_vitk2_m4_fs_ffq1;
  vitk_wo           =           a_vitk_fo_ffq1;
  satfat            =           a_satfat_fs_ffq1;
  sft07             =           a_sft07_fs_ffq1;
  mft07             =           a_mft07_fs_ffq1;
  monfat            =           a_monfat_fs_ffq1;
  ply07             =           a_ply07_fs_ffq1;
  poly              =           a_poly_fs_ffq1;
  alco              =           a_alco_fs_ffq1;
  vitd              =           a_vitd_fs_ffq1;
  vitd_wo           =           a_vitd_fo_ffq1;
  dvitd             =           a_dvitd_fs_ffq1;
  glcsin            =           a_glcsin_fs_ffq1;
run;

data nts_ffq1;
  set nts_ffq1 (keep =
    id

    /*PROTEIN*/
    aprot
    dprot
    prot
    vprot

    /*FAT*/
    afat
    dfat
    vfat

    /*CARBOHYDRATES*/
    carbo

    /*SUCROSE*/
    sucr

    /*FRUCTOSE*/
    fruct

    /*LACTOSE*/
    lact

    /*STARCH*/
    st

    /*GLUCOSE*/
    glu

    /*FIBER*/
    aofib

    /*CALCIUM FROM FOOD*/
    calc_wo

    /*CALCIUM FROM BOTH*/
    calc

    /*IRON TOTAL (NOT JUST HEME)*/
    iron
    iron_wo

    /*IRON HEME*/
    heme

    /*FOLATE*/
    fdfol
    fol98
    fol98_wo
    folic
    dfe

    /*VITAMIN B12*/
    b12
    b12_wo

    /*CHOLINE*/
    choline
    choline_wo

    /*VITAMIN K*/
    dvitk
    vitk
    vitk2_m4
    vitk_wo

    /*SATURATED FATTY ACIDS*/
    satfat
    sft07

    /*MONOUNSATURATED FATTY ACIDS*/
    mft07
    monfat

    /*POLY UNSATURATED FATTY ACIDS*/
    ply07
    poly


    /*ALCOHOL*/
    alco

    /*VITAMIN D*/
    vitd
    vitd_wo
    dvitd

    /*GLUCOSINOLATE*/
    glcsin

    );

    array oldnts1{*}

    /*PROTEIN*/
    aprot
    dprot
    prot
    vprot

    /*FAT*/
    afat
    dfat
    tfat
    vfat

    /*CARBOHYDRATES*/
    carbo

    /*SUCROSE*/
    sucr

    /*FRUCTOSE*/
    fruct

    /*LACTOSE*/
    lact

    /*STARCH*/
    st

    /*GLUCOSE*/
    glu

    /*FIBER*/
    aofib

    /*CALCIUM FROM FOOD*/
    calc_wo

    /*CALCIUM FROM BOTH*/
    calc

    /*IRON TOTAL (NOT JUST HEME)*/
    iron
    iron_wo

    /*IRON HEME*/
    heme

    /*FOLATE*/
    fdfol
    fol98
    fol98_wo
    folic
    dfe

    /*VITAMIN B12*/
    b12
    b12_wo

    /*CHOLINE*/
    choline
    choline_wo

    /*VITAMIN K*/
    dvitk
    vitk
    vitk2_m4
    vitk_wo

    /*SATURATED FATTY ACIDS*/
    satfat
    sft07

    /*MONOUNSATURATED FATTY ACIDS*/
    mft07
    monfat

    /*POLY UNSATURATED FATTY ACIDS*/
    ply07
    poly

    /*ALCOHOL*/
    alco

    /*VITAMIN D*/
    vitd
    vitd_wo
    dvitd

    /*GLUCOSINOLATE*/
    glcsin

    ;

    array newnts1{*}

    /*PROTEIN*/
    aprot_ffq1
    dprot_ffq1
    prot_ffq1
    vprot_ffq1

    /*FAT*/
    afat_ffq1
    dfat_ffq1
    tfat_ffq1
    vfat_ffq1

    /*CARBOHYDRATES*/
    carbo_ffq1

    /*SUCROSE*/
    sucr_ffq1

    /*FRUCTOSE*/
    fruct_ffq1

    /*LACTOSE*/
    lact_ffq1

    /*STARCH*/
    st_ffq1

    /*GLUCOSE*/
    glu_ffq1

    /*FIBER*/
    aofib_ffq1

    /*CALCIUM FROM FOOD*/
    calc_wo_ffq1

    /*CALCIUM FROM BOTH*/
    calc_ffq1

    /*IRON TOTAL (NOT JUST HEME)*/
    iron_ffq1
    iron_wo_ffq1

    /*IRON HEME*/
    heme_ffq1

    /*FOLATE*/
    fdfol_ffq1
    fol98_ffq1
    fol98_wo_ffq1
    folic_ffq1
    dfe_ffq1

    /*VITAMIN B12*/
    b12_ffq1
    b12_wo_ffq1

    /*CHOLINE*/
    choline_ffq1
    choline_wo_ffq1

    /*VITAMIN K*/
    dvitk_ffq1
    vitk_ffq1
    vitk2_m4_ffq1
    vitk_wo_ffq1

    /*SATURATED FATTY ACIDS*/
    satfat_ffq1
    sft07_ffq1

    /*MONOUNSATURATED FATTY ACIDS*/
    mft07_ffq1
    monfat_ffq1

    /*POLY UNSATURATED FATTY ACIDS*/
    ply07_ffq1
    poly_ffq1

    /*ALCOHOL*/
    alco_ffq1

    /*VITAMIN D*/
    vitd_ffq1
    vitd_wo_ffq1
    dvitd_ffq1

    /*GLUCOSINOLATE*/
    glcsin_ffq1
    ;

    do i=1 to dim(oldnts1);
      newnts1{i}=oldnts1{i};
    end;

run;

data energy_ffq1;
  set FFQ.mlvsffq1_noexcl_nts;
  calor_ffq1= calor_fs_ffq1;
run;

data energy_ffq1;
  set energy_ffq1 (keep =
    id
    calor_ffq1
    );
run;

proc sort data=energy_ffq1;
  by id;
run;


/*week 2 FFQ*/
data ffq2;
	set FFQ.mlvsffq2_noexcl_ant;
	proc sort; by id;
run;

/*renaming step because var names changed at some point*/

data ffq2;
  set ffq2;
  aofib             =         a_aofib_fs_ffq2;
  aprot             =         a_aprot_fs_ffq2;
  vprot             =         a_vprot_fs_ffq2;
  prot              =         a_prot_fs_ffq2;
  b1_sulf           =         a_b1_sulf_fs_ffq2;
  b1_sulf_wo        =         a_b1_sulf_fo_ffq2;
  biot_sulf         =         a_biot_sulf_fs_ffq2;
  cys_sulf          =         a_cys_sulf_fs_ffq2;
  meth_sulf         =         a_meth_sulf_fs_ffq2;
  prot_sulf         =         a_prot_sulf_fs_ffq2;
  sulfur            =         a_sulfur_fs_ffq2;
  ala               =         a_pfa183n3c07_fs_ffq2;
  epa               =         a_pf205n3c07_fs_ffq2;
  dpa               =         a_pf225n3c07_fs_ffq2;
  dha               =         a_pf226n3c07_fs_ffq2;
  omega3            =         a_pfn307_fs_ffq2; /*Total Omega 3*/
  omega6            =         a_pfn607_fs_ffq2; /*Total omega 6*/
  trans             =         a_trn07_fs_ffq2; /*total trans*/
run;

/*sulfur*/

data ffq2;
  set ffq2 (keep=id aofib aprot vprot prot b1_sulf b1_sulf_wo biot_sulf cys_sulf
  meth_sulf prot_sulf sulfur
  ala epa dpa dha omega3 omega6 trans)
  end=_end_;

  array oldsulf2{*}

  aofib aprot vprot prot b1_sulf b1_sulf_wo biot_sulf cys_sulf
  meth_sulf prot_sulf sulfur
  ala epa dpa dha omega3 omega6 trans;

  array newsulf2{*}

  aofib_ffq2 aprot_ffq2 vprot_ffq2 prot_ffq2 b1_sulf_ffq2 b1_sulf_wo_ffq2
  biot_sulf_ffq2 cys_sulf_ffq2 meth_sulf_ffq2 prot_sulf_ffq2 sulfur_ffq2
  ala_ffq2 epa_ffq2 dpa_ffq2 dha_ffq2 omega3_ffq2 omega6_ffq2 trans_ffq2;

  do i=1 to dim(oldsulf2);
    newsulf2{i}=oldsulf2{i};
    end;

  nonprot_sulf_ffq2=(sulfur_ffq2)-(prot_sulf_ffq2);

  omega3_noala_ffq2=omega3_ffq2-ala_ffq2;

run;


data ffq2;
    set ffq2 (keep=id aofib_ffq2 aprot_ffq2 vprot_ffq2 prot_ffq2 prot_sulf_ffq2 nonprot_sulf_ffq2 sulfur_ffq2 ala_ffq2 epa_ffq2 dpa_ffq2 dha_ffq2 omega3_ffq2 omega6_ffq2 trans_ffq2 omega3_noala_ffq2)
    end=_end_;
run;

data nts_ffq2;
  set FFQ.mlvsffq2_noexcl_ant;
  aprot             =           a_aprot_fs_ffq2;
  dprot             =           a_dprot_fs_ffq2;
  prot              =           a_prot_fs_ffq2;
  vprot             =           a_vprot_fs_ffq2;
  afat              =           a_afat_fs_ffq2;
  dfat              =           a_dfat_fs_ffq2;
  tfat              =           a_tfat_fs_ffq2;
  vfat              =           a_vfat_fs_ffq2;
  carbo             =           a_carbo_fs_ffq2;
  sucr              =           a_sucr_fs_ffq2;
  fruct             =           a_fruct_fs_ffq2;
  lact              =           a_lact_fs_ffq2;
  st                =           a_st_fs_ffq2;
  glu               =           a_glu_fs_ffq2;
  aofib             =           a_aofib_fs_ffq2;
  calc_wo           =           a_calc_fo_ffq2;
  calc              =           a_calc_fs_ffq2;
  iron              =           a_iron_fo_ffq2;
  iron_wo           =           a_iron_fs_ffq2;
  heme              =           a_heme_fs_ffq2;
  fdfol             =           a_fdfol_fs_ffq2;
  fol98             =           a_fol98_fs_ffq2;
  fol98_wo          =           a_fol98_fo_ffq2;
  folic             =           a_folic_fs_ffq2;
  dfe               =           a_dfe_fs_ffq2;
  b12               =           a_b12_fs_ffq2;
  b12_wo            =           a_b12_fo_ffq2;
  choline           =           a_choline_fs_ffq2;
  choline_wo        =           a_choline_fo_ffq2;
  dvitk             =           a_dvitk_fs_ffq2;
  vitk              =           a_vitk_fs_ffq2;
  vitk2_m4          =           a_vitk2_m4_fs_ffq2;
  vitk_wo           =           a_vitk_fo_ffq2;
  satfat            =           a_satfat_fs_ffq2;
  sft07             =           a_sft07_fs_ffq2;
  mft07             =           a_mft07_fs_ffq2;
  monfat            =           a_monfat_fs_ffq2;
  ply07             =           a_ply07_fs_ffq2;
  poly              =           a_poly_fs_ffq2;
  alco              =           a_alco_fs_ffq2;
  vitd              =           a_vitd_fs_ffq2;
  vitd_wo           =           a_vitd_fo_ffq2;
  dvitd             =           a_dvitd_fs_ffq2;
  glcsin            =           a_glcsin_fs_ffq2;
  proc sort; by id;
run;

data nts_ffq2;
  set nts_ffq2 (keep =
    id
    /*PROTEIN*/
    aprot
    dprot
    prot
    vprot

    /*FAT*/
    afat
    dfat
    tfat
    vfat

    /*CARBOHYDRATES*/
    carbo

    /*SUCROSE*/
    sucr

    /*FRUCTOSE*/
    fruct

    /*LACTOSE*/
    lact

    /*STARCH*/
    st

    /*GLUCOSE*/
    glu

    /*FIBER*/
    aofib

    /*CALCIUM FROM FOOD*/
    calc_wo

    /*CALCIUM FROM BOTH*/
    calc

    /*IRON TOTAL (NOT JUST HEME)*/
    iron
    iron_wo

    /*IRON HEME*/
    heme

    /*FOLATE*/
    fdfol
    fol98
    fol98_wo
    folic
    dfe

    /*VITAMIN B12*/
    b12
    b12_wo

    /*CHOLINE*/
    choline
    choline_wo

    /*VITAMIN K*/
    dvitk
    vitk
    vitk2_m4
    vitk_wo

    /*SATURATED FATTY ACIDS*/
    satfat
    sft07

    /*MONOUNSATURATED FATTY ACIDS*/
    mft07
    monfat

    /*POLY UNSATURATED FATTY ACIDS*/
    ply07
    poly

    /*ALCOHOL*/
    alco

    /*VITAMIN D*/
    vitd
    vitd_wo
    dvitd

    /*GLUCOSINOLATE*/
    glcsin

    );

    array oldnts2{*}

    /*PROTEIN*/
    aprot
    dprot
    prot
    vprot

    /*FAT*/
    afat
    dfat
    tfat
    vfat

    /*CARBOHYDRATES*/
    carbo

    /*SUCROSE*/
    sucr

    /*FRUCTOSE*/
    fruct

    /*LACTOSE*/
    lact

    /*STARCH*/
    st

    /*GLUCOSE*/
    glu

    /*FIBER*/
    aofib

    /*CALCIUM FROM FOOD*/
    calc_wo

    /*CALCIUM FROM BOTH*/
    calc

    /*IRON TOTAL (NOT JUST HEME)*/
    iron
    iron_wo

    /*IRON HEME*/
    heme

    /*FOLATE*/
    fdfol
    fol98
    fol98_wo
    folic
    dfe

    /*VITAMIN B12*/
    b12
    b12_wo

    /*CHOLINE*/
    choline
    choline_wo

    /*VITAMIN K*/
    dvitk
    vitk
    vitk2_m4
    vitk_wo

    /*SATURATED FATTY ACIDS*/
    satfat
    sft07

    /*MONOUNSATURATED FATTY ACIDS*/
    mft07
    monfat

    /*POLY UNSATURATED FATTY ACIDS*/
    ply07
    poly

    /*ALCOHOL*/
    alco

    /*VITAMIN D*/
    vitd
    vitd_wo
    dvitd

    /*GLUCOSINOLATE*/
    glcsin
    ;

    array newnts2{*}

    /*PROTEIN*/
    aprot_ffq2
    dprot_ffq2
    prot_ffq2
    vprot_ffq2

    /*FAT*/
    afat_ffq2
    dfat_ffq2
    tfat_ffq2
    vfat_ffq2

    /*CARBOHYDRATES*/
    carbo_ffq2

    /*SUCROSE*/
    sucr_ffq2

    /*FRUCTOSE*/
    fruct_ffq2

    /*LACTOSE*/
    lact_ffq2

    /*STARCH*/
    st_ffq2

    /*GLUCOSE*/
    glu_ffq2

    /*FIBER*/
    aofib_ffq2

    /*CALCIUM FROM FOOD*/
    calc_wo_ffq2

    /*CALCIUM FROM BOTH*/
    calc_ffq2

    /*IRON TOTAL (NOT JUST HEME)*/
    iron_ffq2
    iron_wo_ffq2

    /*IRON HEME*/
    heme_ffq2

    /*FOLATE*/
    fdfol_ffq2
    fol98_ffq2
    fol98_wo_ffq2
    folic_ffq2
    dfe_ffq2

    /*VITAMIN B12*/
    b12_ffq2
    b12_wo_ffq2

    /*CHOLINE*/
    choline_ffq2
    choline_wo_ffq2

    /*VITAMIN K*/
    dvitk_ffq2
    vitk_ffq2
    vitk2_m4_ffq2
    vitk_wo_ffq2

    /*SATURATED FATTY ACIDS*/
    satfat_ffq2
    sft07_ffq2

    /*MONOUNSATURATED FATTY ACIDS*/
    mft07_ffq2
    monfat_ffq2

    /*POLY UNSATURATED FATTY ACIDS*/
    ply07_ffq2
    poly_ffq2

    /*ALCOHOL*/
    alco_ffq2

    /*VITAMIN D*/
    vitd_ffq2
    vitd_wo_ffq2
    dvitd_ffq2

    /*GLUCOSINOLATE*/
    glcsin_ffq2
    ;

    do i=1 to dim(oldnts2);
      newnts2{i}=oldnts2{i};
    end;

run;

data energy_ffq2;
  set FFQ.mlvsffq2_noexcl_nts;
  calor_ffq2 = calor_fs_ffq2;
run;

data energy_ffq2;
  set energy_ffq2 (keep =
    id
    calor_ffq2
    );
run;

proc sort data=energy_ffq2;
  by id;
run;

data energyffq12;
    merge energy_ffq1 energy_ffq2
    end=_end_;
  by id;
  if first.id;
run;


data energyffq12;
  set energyffq12;
  if id ne .;
run;

proc sort data=energyffq12 out=energyffq12;
  by id;
run;


data ffq12;
    merge ffq1 ffq2
    end=_end_;
  by id;
  if first.id;
run;


data ffq12;
  set ffq12;
  if id ne .;
run;

proc sort data=ffq12 out=ffq12;
  by id;
run;



data ntsffq12;
    merge nts_ffq1 nts_ffq2
    end=_end_;
  by id;
  if first.id;
run;



data ntsffq12;
  set ntsffq12;
  if id ne .;
run;

proc sort data=ntsffq12 out=ntsffq12;
  by id;
run;




/*ffq 1 raw foods*/
data raw1;
  set FFQ.mlvsffq1_raw;
  proc sort; by id;
run;

proc contents data=raw1;
run;


data raw1;
  set raw1
  end=_end_;

  array old1{*} /*need to rename raw foods so they are unique*/

  a_j_ffq1 apple_ffq1 apricot_ffq1 avocado_ffq1 bacon_ffq1 ban_ffq1 beans_ffq1 beer_ffq1 fr_fish_kids_ffq1 blue_ffq1 beef_ffq1
  sand_bf_ham_ffq1 br_rice_ffq1 broc_ffq1 brusl_ffq1 spread_bu_ffq1 cabb_ffq1 cakehr_ffq1 cant_ffq1 caul_ffq1 coke_ffq1 carrot_c_ffq1 candy_nuts_ffq1
  candy_ffq1 celery_ffq1 cold_cer_ffq1 cer_ffq1 chow_ffq1 chix_sk_ffq1 chix_no_ffq1 ckd_cer_ffq1 coff_ffq1 cof_wht_ffq1 coox_brn_rf_ffq1 coox_brn_home_ffq1 coox_brn_ffq1
  corn_ffq1 cot_ch_ffq1 crax_ffq1 cream_ffq1 cr_ch_ffq1 spin_ckd_ffq1 chix_dog_ffq1 chix_no_sand_ffq1 tuna_ffq1 decaf_ffq1 choc_dark_ffq1 tea_decaf_ffq1
  dk_fish_ffq1 hotdog_ffq1 donut_ffq1 o_v_ffq1 eggs_ffq1 eggs_omega_ffq1 zuke_ffq1 eng_muff_ffq1 ff_pot_ffq1 grfrt_ffq1 peppers_ffq1 h2o_ffq1	hamb_ffq1
  xtrlean_hamburg_ffq1 ice_cr_ffq1 ice_let_ffq1 jam_ffq1 kale_ffq1 catsup_ffq1 beer_lite_ffq1 dietsoda_caf_ffq1 liq_ffq1	liver_ffq1 chix_liver_ffq1 mayo_d_ffq1 margarine_ffq1
  mayo_ffq1 choc_ffq1 mix_veg_ffq1 milk2_ffq1 muff_ffq1 oatmeal_bran_ffq1 oat_bran_ffq1 oth_fish_ffq1 o_j_ca_d_ffq1 o_j_ffq1 dietsoda_nocaf_ffq1 onions_ffq1 onions1_ffq1
  oth_nuts_ffq1 proc_mts_ffq1 orang_ffq1 yel_sqs_ffq1 oth_ch_ffq1 oth_f_j_ffq1 soda_nocaf_ffq1 pasta_ffq1 p_bu_ffq1 pancak_ffq1 pot_chip_ffq1 peaches_ffq1
  peas_ffq1 pie_comm_ffq1 pizza_ffq1 yog_plain_ffq1 pork_ffq1 nuts_ffq1 popc_ff_ffq1 popc_ffq1 pot_ffq1 pretzel_ffq1 prun_ffq1 prun_j_ffq1 punch_ffq1
  raisgrp_ffq1 carrot_r_ffq1 rom_let_ffq1 spin_raw_ffq1 r_wine_ffq1 rye_br_ffq1 salsa_ffq1 st_beans_ffq1 bologna_ffq1 yogurt_frozen_ffq1 shrimp_ckd_ffq1 skim_kids_ffq1 soymilk_fort_ffq1 s_roll_h_ffq1
  straw_ffq1 yog_ffq1 tea_ffq1 tofu_ffq1 tom_j_ffq1 tom_ffq1 tortillas_ffq1 tom_s_ffq1 dk_br_ffq1 wh_br_ffq1 milk_ffq1 walnuts_ffq1 wh_rice_ffq1 w_wine_ffq1 swt_pot_ffq1;

  /*
  List from other studies before formatting to MLVS style as above --keep for reference
  q2appl07d q2apric07d q2avo07d q2bacon07d q2ban07d q2bean07d q2beer07d q2bfsh07d q2blueb07d q2bmain07d
  q2bmix07d q2brice07d q2brocc07d q2bruss07d q2sbu07d q2cabb07d q2cakh07d q2cant07d q2cauli07d q2cola07d q2ccar07d q2cdyw07d
  q2cdywo07d q2cel07d q2cer07d q2cerbr07d q2chowd07d q2chwi07d q2chwo07d q2ckcer07d q2coff07d q2cofwh07d q2coknf07d q2cokh07d q2cokr07d
  q2corn07d q2cotch07d q2crack07d q2cream07d q2crmch07d q2cspin07d q2ctdog07d q2chksa07d q2ctuna07d  q2decaf07d q2dchoc07d q2dtea07d
  q2dkfsh07d q2dog07d q2donut07d  q2dress07d q2egg07d q2eggom07d q2eggpl07d q2engl07d q2fries07d q2grfr07d q2grpep07d q2h2o07d q2hamb07d
  q2hambl07d q2icecr07d q2ilett07d q2jam07d q2kale07d q2ketch07d q2lbeer07d q2lcbar07d q2liq07d q2livb07d q2livc07d q2lmayo07d q2marg07d
  q2mayo07d q2mchoc07d q2mixv07d q2m1or207d q2muff07d q2oat07d q2oatbr07d q2ofish07d q2ojca07d q2oj07d q2lcnoc07d q2oniog07d q2oniov07d
  q2onut07d ****q2bkoliv07d********** q2procm07d q2oran07d q2osqua07d q2otch07d q2othj07d q2otsug07d q2pasta07d q2pbut07d q2pcake07d q2pchip07d q2peach07d
  q2peas07d q2pieh07d q2pizza07d  q2plyog07d q2pmain07d q2pnut07d q2ffpop07d  q2popc07d q2pot07d q2pretz07d q2prune07d q2prunj07d q2punch07d
  q2rais07d q2rcar07d q2rlett07d q2rspin07d q2rwine07d  q2ryebr07d q2salsa07d q2sbean07 pmsan07d q2sherb07d q2shrim07d q2skim07d q2soy07d q2srolr07d
  q2straw07d q2flyog07d q2tea07d q2tofu07d q2toj0y q2tom07d q2tort07d q2tosau0y q2dkbr07d q2whbr07d q2whole07d q2wnut07d q2wrice07d q2wwine07d q2yam07d;
  */

  array new1{*}

  ajffq1	applffq1	apricffq1	avoffq1	baconffq1	banffq1	beanffq1	beerffq1	bfshffq1	bluebffq1	bmainffq1
  bmixffq1	briceffq1	broccffq1	brussffq1	butffq1	cabbffq1	cakeffq1	cantffq1	cauliffq1	cbvcfsffq1	ccarffq1	cdywffq1
  cdywoffq1	celffq1	cerffq1	cerbrffq1	chowdffq1	chwiffq1	chwoffq1	ckcerffq1	coffffq1	cofwhffq1	cokffffq1	cokhffq1	cokrffq1
  cornffq1	cotchffq1	crackotffq1	creamffq1	crmchffq1	cspinffq1	ctdogffq1	ctsanffq1	ctunaffq1	dcafffq1	dchocffq1	dcteaffq1
  dkfshffq1	dogffq1	donutffq1	dressffq1	eggffq1	eggomffq1	eggplffq1	englffq1	friesffq1	grfrjffq1	grpepffq1	h2offq1	hambffq1
  hamblffq1	icecrffq1	ilettffq1	jamffq1	kaleffq1	ketchffq1	lbeerffq1	lcbcfffq1	liqffq1	livbffq1	livcffq1	lmayoffq1	margffq1
  mayoffq1	mchocffq1	mixvffq1	mlk12ffq1	muffffq1	oatffq1	oatbrffq1	ofishffq1	ojcaffq1	ojregffq1	olcncffq1	oniogffq1	oniovffq1
  onutffq1	opromffq1	oranffq1	osquaffq1	otchffq1	othjffq1	otsugffq1	pastaffq1	pbutffq1	pcakeffq1	pchipffq1	peachffq1
  peasffq1	piehrffq1	pizzaffq1	plyogffq1	pmainffq1	pnutffq1	popclffq1	popcrffq1	potaffq1	pretzffq1	pruneffq1	prunjffq1	punchffq1
  raisffq1	rcarffq1	rlettffq1	rspinffq1	rwineffq1	ryeffq1	salsaffq1	sbeanffq1	sbolffq1	sherbffq1	shrimffq1	skimffq1	soymffq1	srollffq1
  strawffq1	sweyogffq1	teaffq1	tofuffq1	tojffq1	tomffq1	tortffq1	tosauffq1	wgrbrffq1	whbrffq1	wholeffq1	wnutffq1	wriceffq1	wwineffq1	yamffq1;

  do i=1 to dim(old1);
  	new1{i}=old1{i};
  end;


/*convert to servings per day*/
 do i=1 to DIM(new1);
   if new1{i}=0 then new1{i}=0; /*never, <1/month*/
    else if new1{i}=1 then new1{i}=0.07; /*1-3/month*/
    else if new1{i}=2 then new1{i}=0.14; /*1/week*/
    else if new1{i}=3  then new1{i}=0.43; /*2-4/week*/
    else if new1{i}=4  then new1{i}=0.79; /*5-6/week*/
    else if new1{i}=5  then new1{i}=1; /*1/day*/
    else if new1{i}=6  then new1{i}=2.5; /*2-3/day*/
    else if new1{i}=7  then new1{i}=4.5; /*4-5/day*/
    else if new1{i}=8  then new1{i}=6; /*6+/day*/
    else if new1{i}=9  then new1{i}=0; /*PT*/
    else new1{i}=0; /*missing*/
  end;

/*coding for livers - beef, calf, pork, chicken*/
  array livffq1{2} livcffq1 livbffq1;
     do i=1 to 2;
      if livffq1{i} in (.,1,6) then livffq1{i}=0; /*missing, never, passthrough*/
        else if livffq1{i} in (2) then livffq1{i}=0.01; /*>1/month*/
        else if livffq1{i} in (3) then livffq1{i}=0.03; /*1/month*/
        else if livffq1{i} in (4) then livffq1{i}=0.08; /*2-3/month*/
        else if livffq1{i} in (5) then livffq1{i}=0.14; /*1+/week*/
     end;

/*distinguish between whole grain and refined grain breakfast cereals*/
  if cerbrffq1 in
          (29,43,84,96,111,112,114,115,119,127,128,129,131,132,133,134,146,147,
          148,151,152,154,155,157,164,165,166,167,176,184,185,186,188,0,4,11,12,
          13,15,16,17,19,20,21,22,23,25,30,31,33,35,42,45,46,59,60,61,62,63,64,
          65,66,68,74,75,76,79,80,81,82,83,85,86,88,95,97,98,101,109,113,123,
          124,149,150,153,158,170,173,177,37,38,47,73,130,182)
     then rcerffq1 = cerffq1;
     else wcerffq1 = cerffq1;

  if rcerffq1=. then rcerffq1 = 0;
  if wcerffq1=. then wcerffq1 = 0;

  promeatffq1 = sum(dogffq1, baconffq1, ctdogffq1, opromffq1, sbolffq1);
  redmeatffq1 = sum(hambffq1, bmixffq1, bmainffq1, pmainffq1, hamblffq1);
  orgmeatffq1 = sum(livcffq1, livbffq1);
  fishffq1    = sum(ctunaffq1, dkfshffq1, ofishffq1, shrimffq1, bfshffq1);
  poultffq1   = sum(chwoffq1, chwiffq1, ctsanffq1);
  eggsffq1    = sum(eggffq1, eggomffq1);
  butterffq1  = butffq1;
  margffq1    = margffq1;
  lowdaiffq1  = sum(skimffq1, sherbffq1, cotchffq1, sweyogffq1, plyogffq1, mlk12ffq1);
  highdaiffq1 = sum(wholeffq1, creamffq1, crmchffq1, icecrffq1, otchffq1);
  wineffq1    = sum(wwineffq1, rwineffq1);
  liquorffq1  = liqffq1;
  beerffq1    = sum(beerffq1, lbeerffq1);
  teaffq1     = sum(teaffq1, dcteaffq1);
  coffeeffq1  = sum(dcafffq1, coffffq1);
  fruitffq1   = sum(raisffq1, oranffq1, pruneffq1, banffq1, cantffq1, applffq1, strawffq1, bluebffq1, peachffq1, avoffq1, apricffq1);
  frujuffq1   = sum(othjffq1, ajffq1, grfrjffq1, ojcaffq1, prunjffq1, ojregffq1);
  cruvegffq1  = sum(broccffq1, cabbffq1, kaleffq1, cauliffq1, brussffq1);
  yelvegffq1  = sum(ccarffq1, rcarffq1, osquaffq1, yamffq1);
  tomatoffq1  = sum(tomffq1, tojffq1, tosauffq1);
  leafvegffq1 = sum(rspinffq1, cspinffq1, ilettffq1, rlettffq1);
  legumeffq1  = sum(tofuffq1, sbeanffq1, peasffq1, beanffq1, soymffq1);
  othvegffq1  = sum(cornffq1, mixvffq1, eggplffq1, celffq1, grpepffq1, oniovffq1, oniogffq1);
  potatoffq1  = potaffq1;
  frenchffq1  = friesffq1;
  wholegffq1  = sum(briceffq1, oatffq1, ckcerffq1, oatbrffq1, ryeffq1, wgrbrffq1);
  refingffq1  = sum(whbrffq1, wriceffq1, englffq1, muffffq1, pastaffq1, pcakeffq1, tortffq1);
  pizzaffq1   = pizzaffq1;
  sugdrkffq1  = sum(punchffq1, otsugffq1, cbvcfsffq1);
  lowdrkffq1  = sum(lcbcfffq1, olcncffq1);
  snackffq1   = sum(pchipffq1, pretzffq1, popclffq1, popcrffq1, crackotffq1);
  nutsffq1    = sum(pbutffq1, pnutffq1, onutffq1, wnutffq1);
  mayoffq1    = sum(mayoffq1, lmayoffq1);
  dressffq1   = dressffq1;
  crmsoupffq1 = chowdffq1;
  sweetsffq1  = sum(cdywffq1, cdywoffq1, cokhffq1, cokrffq1, donutffq1, cakeffq1, mchocffq1, dchocffq1, jamffq1, srollffq1, cokffffq1, piehrffq1);
  condimffq1  = sum(cofwhffq1, ketchffq1, salsaffq1);
run;

/* proc export data=raw1 REPLACE
	outfile="/udd/nhlng/sulfurmicrobiome/pipeline/exposure/raw1b.txt"; */

/* removed compared to others' analysis because I couldn't find equivalent MLVS groups: garlic, olive oil, and "other soup" */

proc factor data=raw1 rotate=varimax mineigen=1.5 fuzz=0.15 nfactor=3 scree reorder score
            out=out1; var promeatffq1 redmeatffq1 orgmeatffq1 fishffq1
            poultffq1 eggsffq1 butterffq1 margffq1 lowdaiffq1 highdaiffq1 h2offq1 liquorffq1
            wineffq1 beerffq1 teaffq1 coffeeffq1 fruitffq1 frujuffq1 cruvegffq1
            yelvegffq1 tomatoffq1 leafvegffq1 legumeffq1 othvegffq1 potatoffq1
            frenchffq1 wholegffq1 refingffq1 pizzaffq1 snackffq1 nutsffq1
            sugdrkffq1 lowdrkffq1 mayoffq1 crmsoupffq1
            sweetsffq1 condimffq1 dressffq1;

data out1;
  set out1;

proc rank groups=5 out=factorffq1;
  var factor1 factor2 factor3;
  ranks f1q_ffq1 f2q_ffq1 f3q_ffq1;
run;

proc sort data=factorffq1;
  by id;
run;

data factorffq1;
  set factorffq1;

if factor1=. then f1_ffq1=.;
  else f1_ffq1=factor1;
if factor2=. then f2_ffq1=.;
  else f2_ffq1=factor2;
if factor3=. then f3_ffq1=.;
  else f3_ffq1=factor3;
run;


/*ffq 2 2 raw foods*/
data raw2;
  set FFQ.mlvsffq2_raw;
  proc sort; by id;
run;

/* proc contents data=raw2;
run; */


data raw2;
  set raw2
  end=_end;
  array old2{*} /*need to rename raw food1s so they are unique*/

  a_j_ffq2 apple_ffq2 apricot_ffq2 avocado_ffq2 bacon_ffq2 ban_ffq2 beans_ffq2 beer_ffq2 fr_fish_kids_ffq2 blue_ffq2 beef_ffq2
  sand_bf_ham_ffq2 br_rice_ffq2 broc_ffq2 brusl_ffq2 spread_bu_ffq2 cabb_ffq2 cakehr_ffq2 cant_ffq2 caul_ffq2 coke_ffq2 carrot_c_ffq2 candy_nuts_ffq2
  candy_ffq2 celery_ffq2 cold_cer_ffq2 cer_ffq2 chow_ffq2 chix_sk_ffq2 chix_no_ffq2 ckd_cer_ffq2 coff_ffq2 cof_wht_ffq2 coox_brn_rf_ffq2 coox_brn_home_ffq2 coox_brn_ffq2
  corn_ffq2 cot_ch_ffq2 crax_ffq2 cream_ffq2 cr_ch_ffq2 spin_ckd_ffq2 chix_dog_ffq2 chix_no_sand_ffq2 tuna_ffq2 decaf_ffq2 choc_dark_ffq2 tea_decaf_ffq2
  dk_fish_ffq2 hotdog_ffq2 donut_ffq2 o_v_ffq2 eggs_ffq2 eggs_omega_ffq2 zuke_ffq2 eng_muff_ffq2 ff_pot_ffq2 grfrt_ffq2 peppers_ffq2 h2o_ffq2	hamb_ffq2
  xtrlean_hamburg_ffq2 ice_cr_ffq2 ice_let_ffq2 jam_ffq2 kale_ffq2 catsup_ffq2 beer_lite_ffq2 dietsoda_caf_ffq2 liq_ffq2	liver_ffq2 chix_liver_ffq2 mayo_d_ffq2 margarine_ffq2
  mayo_ffq2 choc_ffq2 mix_veg_ffq2 milk2_ffq2 muff_ffq2 oatmeal_bran_ffq2 oat_bran_ffq2 oth_fish_ffq2 o_j_ca_d_ffq2 o_j_ffq2 dietsoda_nocaf_ffq2 onions_ffq2 onions1_ffq2
  oth_nuts_ffq2 proc_mts_ffq2 orang_ffq2 yel_sqs_ffq2 oth_ch_ffq2 oth_f_j_ffq2 soda_nocaf_ffq2 pasta_ffq2 p_bu_ffq2 pancak_ffq2 pot_chip_ffq2 peaches_ffq2
  peas_ffq2 pie_comm_ffq2 pizza_ffq2 yog_plain_ffq2 pork_ffq2 nuts_ffq2 popc_ff_ffq2 popc_ffq2 pot_ffq2 pretzel_ffq2 prun_ffq2 prun_j_ffq2 punch_ffq2
  raisgrp_ffq2 carrot_r_ffq2 rom_let_ffq2 spin_raw_ffq2 r_wine_ffq2 rye_br_ffq2 salsa_ffq2 st_beans_ffq2 bologna_ffq2 yogurt_frozen_ffq2 shrimp_ckd_ffq2 skim_kids_ffq2 soymilk_fort_ffq2 s_roll_h_ffq2
  straw_ffq2 yog_ffq2 tea_ffq2 tofu_ffq2 tom_j_ffq2 tom_ffq2 tortillas_ffq2 tom_s_ffq2 dk_br_ffq2 wh_br_ffq2 milk_ffq2 walnuts_ffq2 wh_rice_ffq2 w_wine_ffq2 swt_pot_ffq2;

  /*
  List from other studies before formatting to MLVS style as above --keep for reference
  q2appl07d q2apric07d q2avo07d q2bacon07d q2ban07d q2bean07d q2beer07d q2bfsh07d q2blueb07d q2bmain07d
  q2bmix07d q2brice07d q2brocc07d q2bruss07d q2sbu07d q2cabb07d q2cakh07d q2cant07d q2cauli07d q2cola07d q2ccar07d q2cdyw07d
  q2cdywo07d q2cel07d q2cer07d q2cerbr07d q2chowd07d q2chwi07d q2chwo07d q2ckcer07d q2coff07d q2cofwh07d q2coknf07d q2cokh07d q2cokr07d
  q2corn07d q2cotch07d q2crack07d q2cream07d q2crmch07d q2cspin07d q2ctdog07d q2chksa07d q2ctuna07d  q2decaf07d q2dchoc07d q2dtea07d
  q2dkfsh07d q2dog07d q2donut07d  q2dress07d q2egg07d q2eggom07d q2eggpl07d q2engl07d q2fries07d q2grfr07d q2grpep07d q2h2o07d q2hamb07d
  q2hambl07d q2icecr07d q2ilett07d q2jam07d q2kale07d q2ketch07d q2lbeer07d q2lcbar07d q2liq07d q2livb07d q2livc07d q2lmayo07d q2marg07d
  q2mayo07d q2mchoc07d q2mixv07d q2m1or207d q2muff07d q2oat07d q2oatbr07d q2ofish07d q2ojca07d q2oj07d q2lcnoc07d q2oniog07d q2oniov07d
  q2onut07d q2bkoliv07d q2procm07d q2oran07d q2osqua07d q2otch07d q2othj07d q2otsug07d q2pasta07d q2pbut07d q2pcake07d q2pchip07d q2peach07d
  q2peas07d q2pieh07d q2pizza07d  q2plyog07d q2pmain07d q2pnut07d q2ffpop07d  q2popc07d q2pot07d q2pretz07d q2prune07d q2prunj07d q2punch07d
  q2rais07d q2rcar07d q2rlett07d q2rspin07d q2rwine07d  q2ryebr07d q2salsa07d q2sbean07 pmsan07d q2sherb07d q2shrim07d q2skim07d q2soy07d q2srolr07d
  q2straw07d q2flyog07d q2tea07d q2tofu07d q2toj0y q2tom07d q2tort07d q2tosau0y q2dkbr07d q2whbr07d q2whole07d q2wnut07d q2wrice07d q2wwine07d q2yam07d;
  */

  array new2{*}

  ajffq2	applffq2	apricffq2	avoffq2	baconffq2	banffq2	beanffq2	beerffq2	bfshffq2	bluebffq2	bmainffq2
  bmixffq2	briceffq2	broccffq2	brussffq2	butffq2	cabbffq2	cakeffq2	cantffq2	cauliffq2	cbvcfsffq2	ccarffq2	cdywffq2
  cdywoffq2	celffq2	cerffq2	cerbrffq2	chowdffq2	chwiffq2	chwoffq2	ckcerffq2	coffffq2	cofwhffq2	cokffffq2	cokhffq2	cokrffq2
  cornffq2	cotchffq2	crackotffq2	creamffq2	crmchffq2	cspinffq2	ctdogffq2	ctsanffq2	ctunaffq2	dcafffq2	dchocffq2	dcteaffq2
  dkfshffq2	dogffq2	donutffq2	dressffq2	eggffq2	eggomffq2	eggplffq2	englffq2	friesffq2	grfrjffq2	grpepffq2	h2offq2	hambffq2
  hamblffq2	icecrffq2	ilettffq2	jamffq2	kaleffq2	ketchffq2	lbeerffq2	lcbcfffq2	liqffq2	livbffq2	livcffq2	lmayoffq2	margffq2
  mayoffq2	mchocffq2	mixvffq2	mlk12ffq2	muffffq2	oatffq2	oatbrffq2	ofishffq2	ojcaffq2	ojregffq2	olcncffq2	oniogffq2	oniovffq2
  onutffq2		opromffq2	oranffq2	osquaffq2	otchffq2	othjffq2	otsugffq2	pastaffq2	pbutffq2	pcakeffq2	pchipffq2	peachffq2
  peasffq2	piehrffq2	pizzaffq2	plyogffq2	pmainffq2	pnutffq2	popclffq2	popcrffq2	potaffq2	pretzffq2	pruneffq2	prunjffq2	punchffq2
  raisffq2	rcarffq2	rlettffq2	rspinffq2	rwineffq2	ryeffq2	salsaffq2	sbeanffq2	sbolffq2	sherbffq2	shrimffq2	skimffq2	soymffq2	srollffq2
  strawffq2	sweyogffq2	teaffq2	tofuffq2	tojffq2	tomffq2	tortffq2	tosauffq2	wgrbrffq2	whbrffq2	wholeffq2	wnutffq2	wriceffq2	wwineffq2	yamffq2;

  do i=1 to dim(old2);
    new2{i}=old2{i};
    end;

  /*convert to servings per day*/
  do i=1 to DIM(new2);
   if new2{i}=0 then new2{i}=0; /*never, <1/month*/
    else if new2{i}=1 then new2{i}=0.07; /*1-3/month*/
    else if new2{i}=2 then new2{i}=0.14; /*1/week*/
    else if new2{i}=3  then new2{i}=0.43; /*2-4/week*/
    else if new2{i}=4  then new2{i}=0.79; /*5-6/week*/
    else if new2{i}=5  then new2{i}=1; /*1/day*/
    else if new2{i}=6  then new2{i}=2.5; /*2-3/day*/
    else if new2{i}=7  then new2{i}=4.5; /*4-5/day*/
    else if new2{i}=8  then new2{i}=6; /*6+/day*/
    else if new2{i}=9  then new2{i}=0; /*PT*/
    else new2{i}=0; /*missing*/
  end;

  /*coding for livers - beef, calf, pork, chicken*/
  array livffq2{2} livcffq2 livbffq2;
     do i=1 to 2;
      if livffq2{i} in (.,1,6) then livffq2{i}=0; /*missing, never, passthrough*/
        else if livffq2{i} in (2) then livffq2{i}=0.01; /*>1/month*/
        else if livffq2{i} in (3) then livffq2{i}=0.03; /*1/month*/
        else if livffq2{i} in (4) then livffq2{i}=0.08; /*2-3/month*/
        else if livffq2{i} in (5) then livffq2{i}=0.14; /*1+/week*/
     end;

  /*distinguish between whole grain and refined grain breakfast cereals*/
  if cerbrffq2 in
          (29,43,84,96,111,112,114,115,119,127,128,129,131,132,133,134,146,147,
          148,151,152,154,155,157,164,165,166,167,176,184,185,186,188,0,4,11,12,
          13,15,16,17,19,20,21,22,23,25,30,31,33,35,42,45,46,59,60,61,62,63,64,
          65,66,68,74,75,76,79,80,81,82,83,85,86,88,95,97,98,101,109,113,123,
          124,149,150,153,158,170,173,177,37,38,47,73,130,182)
     then rcerffq2 = cerffq2;
     else wcerffq2 = cerffq2;

  if rcerffq2=. then rcerffq2 = 0;
  if wcerffq2=. then wcerffq2 = 0;


  promeatffq2 = sum(dogffq2, baconffq2, ctdogffq2, opromffq2, sbolffq2);
  redmeatffq2 = sum(hambffq2, bmixffq2, bmainffq2, pmainffq2, hamblffq2);
  orgmeatffq2 = sum(livcffq2, livbffq2);
  fishffq2    = sum(ctunaffq2, dkfshffq2, ofishffq2, shrimffq2, bfshffq2);
  poultffq2   = sum(chwoffq2, chwiffq2, ctsanffq2);
  eggsffq2    = sum(eggffq2, eggomffq2);
  butterffq2  = butffq2;
  margffq2    = margffq2;
  lowdaiffq2  = sum(skimffq2, sherbffq2, cotchffq2, sweyogffq2, plyogffq2, mlk12ffq2);
  highdaiffq2 = sum(wholeffq2, creamffq2, crmchffq2, icecrffq2, otchffq2);
  wineffq2    = sum(wwineffq2, rwineffq2);
  liquorffq2  = liqffq2;
  beerffq2    = sum(beerffq2, lbeerffq2);
  teaffq2     = sum(teaffq2, dcteaffq2);
  coffeeffq2  = sum(dcafffq2, coffffq2);
  fruitffq2   = sum(raisffq2, oranffq2, pruneffq2, banffq2, cantffq2, applffq2, strawffq2, bluebffq2, peachffq2, avoffq2, apricffq2);
  frujuffq2   = sum(othjffq2, ajffq2, grfrjffq2, ojcaffq2, prunjffq2, ojregffq2);
  cruvegffq2  = sum(broccffq2, cabbffq2, kaleffq2, cauliffq2, brussffq2);
  yelvegffq2  = sum(ccarffq2, rcarffq2, osquaffq2, yamffq2);
  tomatoffq2  = sum(tomffq2, tojffq2, tosauffq2);
  leafvegffq2 = sum(rspinffq2, cspinffq2, ilettffq2, rlettffq2);
  legumeffq2  = sum(tofuffq2, sbeanffq2, peasffq2, beanffq2, soymffq2);
  othvegffq2  = sum(cornffq2, mixvffq2, eggplffq2, celffq2, grpepffq2, oniovffq2, oniogffq2);
  potatoffq2  = potaffq2;
  frenchffq2  = friesffq2;
  wholegffq2  = sum(briceffq2, oatffq2, ckcerffq2, oatbrffq2, ryeffq2, wgrbrffq2);
  refingffq2  = sum(whbrffq2, wriceffq2, englffq2, muffffq2, pastaffq2, pcakeffq2, tortffq2);
  pizzaffq2   = pizzaffq2;
  sugdrkffq2  = sum(punchffq2, otsugffq2, cbvcfsffq2);
  lowdrkffq2  = sum(lcbcfffq2, olcncffq2);
  snackffq2   = sum(pchipffq2, pretzffq2, popclffq2, popcrffq2, crackotffq2);
  nutsffq2    = sum(pbutffq2, pnutffq2, onutffq2, wnutffq2);
  mayoffq2    = sum(mayoffq2, lmayoffq2);
  dressffq2   = dressffq2;
  crmsoupffq2 = chowdffq2;
  sweetsffq2  = sum(cdywffq2, cdywoffq2, cokhffq2, cokrffq2, donutffq2, cakeffq2, mchocffq2, dchocffq2, jamffq2, srollffq2, cokffffq2, piehrffq2);
  condimffq2  = sum(cofwhffq2, ketchffq2, salsaffq2);

run;



proc factor data=raw2 rotate=varimax mineigen=1.5 fuzz=0.15 nfactor=3 scree reorder
            out=out2; var promeatffq2 redmeatffq2 orgmeatffq2 fishffq2
            poultffq2 eggsffq2 butterffq2 margffq2 lowdaiffq2 highdaiffq2 h2offq2 liquorffq2
            wineffq2 beerffq2 teaffq2 coffeeffq2 fruitffq2 frujuffq2 cruvegffq2
            yelvegffq2 tomatoffq2 leafvegffq2 legumeffq2 othvegffq2 potatoffq2
            frenchffq2 wholegffq2 refingffq2 pizzaffq2 snackffq2 nutsffq2
            sugdrkffq2 lowdrkffq2 mayoffq2 crmsoupffq2
            sweetsffq2 condimffq2 dressffq2 ;

data out2;
  set out2;

proc rank groups=5 out=factorffq2;
  var factor1 factor2 factor3;
  ranks f1q_ffq2 f2q_ffq2 f3q_ffq2;
run;

proc sort data=factorffq2;
  by id;
run;

data factorffq2;
  set factorffq2;

if factor1=. then f1_ffq2=.;
  else f1_ffq2=factor1;
if factor2=. then f2_ffq2=.;
  else f2_ffq2=factor2;
if factor3=. then f3_ffq2=.;
  else f3_ffq2=factor3;
run;


data raw12;
  merge raw1 raw2;
  by id;
run;


data raw12;
  set raw12;
  if id ne .;
run;

proc sort data=raw12 out=raw12;
  by id;
run;

data out12;
  merge out1 out2;
  by id;
run;



data out12;
  set out12;
  if id ne .;
run;

proc sort data=out12 out=out12;
  by id;
run;

data factorffq12;
  merge factorffq1 factorffq2;
  by id;
run;


data factorffq12;
  set factorffq12;
  if id ne .;
run;

proc sort data=factorffq12 out=factorffq12;
  by id;
run;


/**********************************
/**********************************
*          MLVS FECAL QQs         *
***********************************
**********************************/

data stoolqq1;
	set STOOLQQ.fecal_qu1;
	proc sort; by id;
run;

data stoolqq2;
	set STOOLQQ.fecal_qu2;
	proc sort; by id;
run;

data stoolqq12;
  merge stoolqq1 stoolqq2;
  by id;
  if first.id;
run;


proc sort data=stoolqq12 out=stoolqq12;
  by id;
run;


/**********************************
/**********************************
*          MLVS BLOOD -- CRP
"below LLOQ (<0.01 mg/L)"
"above upper detect limit (>20 mg/dL).

All those that were below LOD were flagged with -99 because lab marked these sameple as <0.10.
Those that were above LOD the value marked of 20 was kept.
***********************************
**********************************/

data blood1;
set blood.mlvs_plasma1;
if crp_plasma1=-99 then crplow1=1;
else crplow1=0;
if crp_plasma1=20 then crphigh1=1;
else crphigh1=0;
if crplow1=1 then crp_plasma1=.;
if crphigh1=1 then crp_plasma1=.;
run;
proc sort; by id;
run;

data blood2;
set blood.mlvs_plasma2;
if crp_plasma2=-99 then crplow2=1;
else crplow2=0;
if crp_plasma2=20 then crphigh2=1;
else crphigh2=0;
if crplow2=1 then crp_plasma2=.;
if crphigh2=1 then crp_plasma2=.;
proc sort; by id;
run;

data blood12;
merge blood1 blood2;
by id;
run;


proc sort data=blood12 out=blood12;
by id;
run;

proc sort data=hpfsvars nodupkey; by id;
proc sort data=mets nodupkey; by id;
proc sort data=nutrients nodupkey; by id;
proc sort data=ntsw12 nodupkey; by id;
proc sort data=ntsffq12  nodupkey; by id;
proc sort data=energyw12 nodupkey; by id;
proc sort data=ffq12 nodupkey; by id;
proc sort data=stoolqq12 nodupkey; by id;
proc sort data=blood12 nodupkey; by id;
proc sort data=raw12 nodupkey; by id;
proc sort data=factorffq12 nodupkey; by id;
proc sort data=energyffq12 nodupkey; by id;



data mlvs_exposure;
    merge hpfsvars (in=inder)
    mets
    nutrients
    ntsw12
    ntsffq12
    energyw12
    ffq12
    stoolqq12
    blood12
    raw12
    /*factorw12*/
    factorffq12
    energyffq12
    /*stooldates*/
  end=_end_;
  by id;
  exrec=1;
  if first.id and inder then exrec=0;

/*******************************************************************************
********************************************************************************
********************************************************************************

  note the change in included IDs here--exporting alias IDs to conform to data
  practices for data external to capecod (to avoid using "identifying" cohort id
  in microbiome pipelines on Odyseey)

  remove colectomy: 706007

********************************************************************************
********************************************************************************
********************************************************************************
********************************************************************************
*/

   if id in (
    715382 132386 502046 601436 180684 710438 713121 123028 169147 169345
    102256 711184 121224 187624 103088 123024 720066 172920 707531 193642
    187253 702133 144188 201963 703563 119015 126618 152051 184788 152745
    164624 143834 301604 191681 700427 714588 105722 711703 715980 175247
    300095 146735 404604 202899 /*706007*/ 132814 507513 703402 407665 132081
    718156 306356 702790 125037 710188 160922 303967 708498 709987 136027
    507589 152868 504477 191141 178944 721551 711501 132653 307486 115583
    718509 184002 157602 188525 191061 154751 719199 186591 185402 706964
    128705 714001	713674 406577	506358 186852	180711 147583	704675 124977
    300380 401383 138640 151803	191958 106513	721829 702077	170822 705560
    191864 700221	171806 605111	111397 154567	501268 309938	138800 170764
    113594 600566	719721 157657	186998 181015	188877 169446	189098 604307
    707455 134614	135950 501318	119467 139085	124734 158115	310074 705793
    189423 503282	703470 194507	152531 185012 153169 505514	722520 186423
    307836 186679	700654 177106	134863 708042	105720 721801	721352 174812
    709517 132333	600174 135092	408039 186911 193858 715520	194659 702107
    700634 136266	709059 718848	192790 709780	702658 169386	718461 190849
    164524 185274 106580 700499	704045 192035	118842 167006	167377 603938
    153186 601743	704791 713527	191297 718054	704685 145396	504446 506472
    201172 712059 700243 718325	172663 185178	186404 113632	147046 135063
    500102 508298	115918 102191	189421 131692	169192 308560 135400 600662
    708007 167548	193391 304397	712536 505643 715370 718068	605099 139232
    703264 702746	716130 194575	704195 704710	132636 700742	504433 300523
    706038 202747 149534 301165	176569 150347	127956 711104	101944 114965
    169313 400160	124752 149708	702416 156855	709381 185855	309876 700375
    700770 180754 194558 189208	201956 710751	713146 704904	134957 186616
    717849 400954	702031 301517	143977 717211	710651 132675 186884 715886
    406855 195395	722625 186389	134964 408015	722236 195016	136240 160746
    706639 712861	307648 507367	188818 700528 150178 177243	307293 704369
    403521 714515	601775 718329	307838 405972	180686 714418	721863 710661
    193481 115879	402496 124363	138509 712661	190351 174450    );


 * if inder;
  if first.id;
run;


/*******************************************************************************
/*******************************************************************************
*                                 MAIN EXPOSURES                               *
********************************************************************************
*******************************************************************************/

/*********************************************************
* Cumulative intake of long-term fiber and calorie intake*
***********************************/
data mlvs_exposure;
  set mlvs_exposure
  end=_end_;



array nut {7, 10}
aofib86a frtaf86a ceraf86a vegaf86a calor86n f183s86a       f205s86a      f226s86a      trnss86a  n6s86a
aofib90a frtaf90a ceraf90a vegaf90a calor90n f183s90a       f205s90a      f226s90a      trnss90a  n6s90a
aofib94a frtaf94a ceraf94a vegaf94a calor94n f183s94a       f205s94a      f226s94a      trnss94a  n6s94a
aofib98a frtaf98a ceraf98a vegaf98a calor98n a183098a       pf20598a      pf22698a      trn0098a  n6s0098a
aofib02a frtaf02a ceraf02a vegaf02a calor02n pfa183n3c0202a pf205n3c0202a pf226n3c0202a trn0202a  n60202a
aofib06a frtaf06a ceraf06a vegaf06a calor06n pfa183n3c0706a pf205n3c0706a pf226n3c0706a trn0706a  n60706a
aofib10a frtaf10a ceraf10a vegaf10a calor10n pfa183n3c1110a pf205n3c1110a pf226n3c1110a trn1110a  n61110a
;
array nutv {7, 10}
aofib86v frtaf86v ceraf86v vegaf86v calor86v ala86v         epa86v        dha86v        trans86v  omega686v
aofib90v frtaf90v ceraf90v vegaf90v calor90v ala90v         epa90v        dha90v        trans90v  omega690v
aofib94v frtaf94v ceraf94v vegaf94v calor94v ala94v         epa94v        dha94v        trans94v  omega694v
aofib98v frtaf98v ceraf98v vegaf98v calor98v ala98v         epa98v        dha98v        trans98v  omega698v
aofib02v frtaf02v ceraf02v vegaf02v calor02v ala02v         epa02v        dha02v        trans02v  omega602v
aofib06v frtaf06v ceraf06v vegaf06v calor06v ala06v         epa06v        dha06v        trans06v  omega606v
aofib10v frtaf10v ceraf10v vegaf10v calor10v ala10v         epa10v        dha10v        trans10v  omega610v
;
do i=2 to 7;
  do j=1 to 10;
  if nut{i,j}<0 then nut{i,j}=nut{i-1,j};
  end;
end;


do j=1 to 10;
    nutv{1,j}=nut{1,j};
    do i=2 to 7;
        sumvar=0;
        n=0;
        do k=1 to i;
            if (nut{k,j} ne .) then do;
                n=n+1;
                sumvar=sumvar+nut{k,j};
            end;
        end;
    if n=0 then nutv{i,j}=nut{1,j};
    else nutv{i,j}=sumvar/n;
    end;
end;

/*DPA was not available until 1998*/
/*omega3 excluding ALA*/
omega3_noala98a=pfn3098a-a183098a;
omega3_noala02a=pfn30202a-pfa183n3c0202a;
omega3_noala06a=pfn30706a-pfa183n3c0706a;
omega3_noala10a=pfn31110a-pfa183n3c1110a;

array dpa {4, 3}
pf22598a       pfn3098a    omega3_noala98a
pf225n3c0202a  pfn30202a   omega3_noala02a
pf225n3c0706a  pfn30706a   omega3_noala06a
pf225n3c1110a  pfn31110a   omega3_noala10a

;
array dpav {4, 3}
dpa98v         omega398v   omega3_noala98v
dpa02v         omega302v   omega3_noala02v
dpa06v         omega306v   omega3_noala06v
dpa10v         omega310v   omega3_noala10v
;
do i=2 to 4;
  do j=1 to 3;
  if dpa{i,j}<0 then dpa{i,j}=dpa{i-1,j};
  end;
end;


do j=1 to 3;
    dpav{1,j}=dpa{1,j};
    do i=2 to 4;
        sumvar=0;
        n=0;
        do k=1 to i;
            if (dpa{k,j} ne .) then do;
                n=n+1;
                sumvar=sumvar+dpa{k,j};
            end;
        end;
    if n=0 then dpav{i,j}=dpa{1,j};
    else dpav{i,j}=sumvar/n;
    end;
end;

run;




/**********************************
*     Average of two MLVS FFQs        *
***********************************/


data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

  aofib_avg       = mean(aofib_ffq1, aofib_ffq2);
  aprot_avg       = mean(aprot_ffq1, aprot_ffq2);
  vprot_avg       = mean(vprot_ffq1, vprot_ffq2);
  prot_avg        = mean(prot_ffq1, prot_ffq2);
  prot_sulf_avg   = mean(prot_sulf_ffq1, prot_sulf_ffq2);
  sulfur_avg      = mean(sulfur_ffq1, sulfur_ffq2);
  nonprot_sulf_avg = mean(nonprot_sulf_ffq1, nonprot_sulf_ffq2);

  /*omega-3*/
  ala_avg         = mean(ala_ffq1, ala_ffq2);
  epa_avg         = mean(epa_ffq1, epa_ffq2);
  dpa_avg         = mean(dpa_ffq1, dpa_ffq2);
  dha_avg         = mean(dha_ffq1, dha_ffq2);
  trans_avg       = mean(trans_ffq1, trans_ffq2);
  omega3_avg      = mean(omega3_ffq1, omega3_ffq2);
  omega6_avg      = mean(omega6_ffq1, omega6_ffq2);
  omega3_noala_avg = mean(omega3_noala_ffq1, omega3_noala_ffq2);

run;


proc rank data=mlvs_exposure out=mlvs_exposure groups=5;
  var agemlvs aprot_avg vprot_avg prot_avg prot_sulf_avg sulfur_avg nonprot_sulf_avg;
  ranks agemlvsq aprot_avgq vprot_avgq prot_avgq prot_sulf_avgq sulfur_avgq nonprot_sulf_avgq;
run;




/**********************************
*           NUTRIENTS             *
***********************************/
data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

  calor_avg=mean (calor_ffq1, calor_ffq2);

  aprot_avg=mean (aprot_ffq1, aprot_ffq2);
  dprot_avg=mean (dprot_ffq1, dprot_ffq2);
  prot_avg=mean (prot_ffq1, prot_ffq2);
  vprot_avg=mean (vprot_ffq1, vprot_ffq2);
  afat_avg=mean (afat_ffq1, afat_ffq2);
  dfat_avg=mean (dfat_ffq1, dfat_ffq2);
  tfat_avg=mean (tfat_ffq1, tfat_ffq2);
  vfat_avg=mean (vfat_ffq1, vfat_ffq2);
  carbo_avg=mean (carbo_ffq1, carbo_ffq2);

  /*SUCROSE*/
  sucr_avg=mean (sucr_ffq1, sucr_ffq2);

  /*FRUCTOSE*/
  fruct_avg=mean (fruct_ffq1, fruct_ffq2);

  /*LACTOSE*/
  lact_avg=mean (lact_ffq1, lact_ffq2);

  /*STARCH*/
  st_avg=mean (st_ffq1, st_ffq2);

  /*GLUCOSE*/
  glu_avg=mean (glu_ffq1, glu_ffq2);

  /*FIBER*/
  aofib_avg=mean (aofib_ffq1, aofib_ffq2);

  /*CALCIUM FROM FOOD*/
  calc_wo_avg=mean (calc_wo_ffq1, calc_wo_ffq2);

  /*CALCIUM FROM BOTH*/
  calc_avg=mean (calc_ffq1, calc_ffq2);

  /*IRON TOTAL (NOT JUST HEME)*/
  iron_avg=mean (iron_ffq1, iron_ffq2);
  iron_wo_avg=mean (iron_wo_ffq1, iron_wo_ffq2);

  /*IRON HEME*/
  heme_avg=mean (heme_ffq1, heme_ffq2);

  /*FOLATE*/
  fdfol_avg=mean (fdfol_ffq1, fdfol_ffq2);
  fol98_avg=mean (fol98_ffq1, fol98_ffq2);
  fol98_wo_avg=mean (fol98_wo_ffq1, fol98_wo_ffq2);
  folic_avg=mean (folic_ffq1, folic_ffq2);
  dfe_avg=mean (dfe_ffq1, dfe_ffq2);

  /*VITAMIN B12*/
  b12_avg=mean (b12_ffq1, b12_ffq2);
  b12_wo_avg=mean (b12_wo_ffq1, b12_wo_ffq2);

  /*CHOLINE*/
  choline_avg=mean (choline_ffq1, choline_ffq2);
  choline_wo_avg=mean (choline_wo_ffq1, choline_wo_ffq2);

  /*VITAMIN K*/
  dvitk_avg=mean (dvitk_ffq1, dvitk_ffq2);
  vitk_avg=mean (vitk_ffq1, vitk_ffq2);
  vitk2_m4_avg=mean (vitk2_m4_ffq1, vitk2_m4_ffq2);
  vitk_wo_avg=mean (vitk_wo_ffq1, vitk_wo_ffq2);

  /*SATURATED FATTY ACIDS*/
  satfat_avg=mean (satfat_ffq1, satfat_ffq2);
  sft07_avg=mean (sft07_ffq1, sft07_ffq2);

  /*MONOUNSATURATED FATTY ACIDS*/
  mft07_avg=mean (mft07_ffq1, mft07_ffq2);
  monfat_avg=mean (monfat_ffq1, monfat_ffq2);

  /*POLY UNSATURATED FATTY ACIDS*/
  ply07_avg=mean (ply07_ffq1, ply07_ffq2);
  poly_avg=mean (poly_ffq1, poly_ffq2);

  /*ALCOHOL*/
  alco_avg=mean (alco_ffq1, alco_ffq2);

  /*VITAMIN D*/
  vitd_avg=mean (vitd_ffq1, vitd_ffq2);
  vitd_wo_avg=mean (vitd_wo_ffq1, vitd_wo_ffq2);
  dvitd_avg=mean (dvitd_ffq1, dvitd_ffq2);

  /*GLUCOSINOLATE*/
  glcsin_avg=mean (glcsin_ffq1, glcsin_ffq2);

run;



/**********************************
***********************************
*        DIETARY PATTERNS         *
***********************************
***********************************/

/****************
*  SULFUR RISK  *
*****************/

/*calculate average intake per food group*/

data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

  promeatavg = mean(promeatffq1, promeatffq2);
  redmeatavg = mean(redmeatffq1, redmeatffq2);
  orgmeatavg = mean(orgmeatffq1, orgmeatffq2);
  fishavg    = mean(fishffq1, fishffq2);
  poultavg   = mean(poultffq1, poultffq2);
  eggsavg    = mean(eggsffq1, eggsffq2);
  butteravg  = mean(butterffq1, butterffq2);
  margavg    = mean(margffq1, margffq2);
  lowdaiavg  = mean(lowdaiffq1, lowdaiffq2);
  highdaiavg = mean(highdaiffq1, highdaiffq2);
  wineavg    = mean(wineffq1, wineffq2);
  liquoravg  = mean(liquorffq1, liquorffq2);
  beeravg    = mean(beerffq1, beerffq2);
  teaavg     = mean(teaffq1, teaffq2);
  coffeeavg  = mean(coffeeffq1, coffeeffq2);
  fruitavg   = mean(fruitffq1, fruitffq2);
  frujuavg   = mean(frujuffq1, frujuffq2);
  cruvegavg  = mean(cruvegffq1, cruvegffq2);
  yelvegavg  = mean(yelvegffq1, yelvegffq2);
  tomatoavg  = mean(tomatoffq1, tomatoffq2);
  leafvegavg = mean(leafvegffq1, leafvegffq2);
  legumeavg  = mean(legumeffq1, legumeffq2);
  othvegavg  = mean(othvegffq1, othvegffq2);
  potatoavg  = mean(potatoffq1, potatoffq2);
  frenchavg  = mean(frenchffq1, frenchffq2);
  wholegavg  = mean(wholegffq1, wholegffq2);
  refingavg  = mean(refingffq1, refingffq2);
  pizzaavg   = mean(pizzaffq1, pizzaffq2);
  sugdrkavg  = mean(sugdrkffq1, sugdrkffq2);
  lowdrkavg  = mean(lowdrkffq1, lowdrkffq2);
  snackavg   = mean(snackffq1, snackffq2);
  nutsavg    = mean(nutsffq1, nutsffq2);
  mayoavg    = mean(mayoffq1, mayoffq2);
  dressavg   = mean(dressffq1, dressffq2);
  crmsoupavg = mean(crmsoupffq1, crmsoupffq2);
  sweetsavg  = mean(sweetsffq1, sweetsffq2);
  condimavg  = mean(condimffq1, condimffq2);

  label promeatavg = 'processed meats';
  label redmeatavg = 'red meats';
  label orgmeatavg = 'organ meats';
  label fishavg    = 'fish';
  label poultavg   = 'poultry';
  label eggsavg    = 'eggs';
  label butteravg  = 'butter';
  label margavg    = 'margarine';
  label lowdaiavg  = 'low fat dairy';
  label highdaiavg = 'high fat dairy';
  label wineavg    = 'wine';
  label liquoravg  = 'liquor';
  label beeravg    = 'beer';
  label teaavg     = 'tea';
  label coffeeavg  = 'coffee';
  label fruitavg   = 'fruits';
  label frujuavg   = 'fruit juice';
  label cruvegavg  = 'cruciferous veggies';
  label yelvegavg  = 'yellow veggies';
  label tomatoavg  = 'tomatoes';
  label leafvegavg = 'leafy veggies';
  label legumeavg  = 'legumes';
  label othvegavg  = 'other veggies';
  label potatoavg  = 'potatoes';
  label frenchavg  = 'french fries';
  label wholegavg  = 'whole grains';
  label refingavg  = 'refined grains';
  label pizzaavg   = 'pizza';
  label sugdrkavg  = 'sugar beverages';
  label lowdrkavg  = 'diet beverages';
  label snackavg   = 'snacks';
  label nutsavg    = 'nuts';
  label mayoavg    = 'mayo';
  label dressavg   = 'dressing';
  label crmsoupavg = 'cream soup';
  label sweetsavg  = 'sweets';
  label condimavg  = 'condiments';

run;


/****************
*WESTERN/PRUDENT*
*****************/

data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

west_avg  = mean(f2_ffq1, f2_ffq2);
prud_avg  = mean(f1_ffq1, f1_ffq2);


proc rank groups=5 out=factoravg;
  var west_avg prud_avg /*west_w1 west_w2 prud_w1 prud_w2*/;
  ranks westq_avg prudq_avg /*westq_w1 westq_w2 prudq_w1 prudq_w2*/;
run;

data mlvs_exposure;
  merge mlvs_exposure factoravg;
  by id;
run;


/**********************************
/**********************************
*           COVARIATES            *
***********************************
**********************************/

/*age calculated near derived vars above (only at single point)*/

/*cumulative average BMI*/

data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

*BMI in 2012;
if height>0 and wt12>0 then bmi12=(wt12*0.45359237)/((height*25.4/1000)*(height*25.4/1000));

  array bmia{*} bmi86 bmi88 bmi90 bmi92 bmi94 bmi96 bmi98 bmi00 bmi02 bmi04 bmi06 bmi08 bmi10 bmi12;
  array bmicum{*} bmic86 bmic88 bmic90 bmic92 bmic94 bmic96 bmic98 bmic00 bmic02 bmic04 bmic06 bmic08 bmic10 bmic12;

  do i=1 to DIM(bmia);
  if bmia{i}=0 then bmia{i}=.;
  end;

/*
   do k=DIM(bmia) to 2 by -1;
    if bmia{k} eq . and bmia{k-1} ne . then bmia{k}=bmia{k-1};
  end;
*/

  bmic86=bmi86;
  bmic88=mean(bmi86, bmi88);
  bmic90=mean(bmi86, bmi88, bmi90);
  bmic92=mean(bmi86, bmi88, bmi90, bmi92);
  bmic94=mean(bmi86, bmi88, bmi90, bmi92, bmi94);
  bmic96=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96);
  bmic98=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98);
  bmic00=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98, bmi00);
  bmic02=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98, bmi00, bmi02);
  bmic04=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98, bmi00, bmi02, bmi04);
  bmic06=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98, bmi00, bmi02, bmi04, bmi06);
  bmic08=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98, bmi00, bmi02, bmi04, bmi06, bmi08);
  bmic10=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98, bmi00, bmi02, bmi04, bmi06, bmi08, bmi10);
  bmic12=mean(bmi86, bmi88, bmi90, bmi92, bmi94, bmi96, bmi98, bmi00, bmi02, bmi04, bmi06, bmi08, bmi10, bmi12);


  if bmi12<18.5 then bmi12cat=1;
    else if 18.5  ge bmi12 lt 25 then bmi12cat=2;
    else if 25 ge bmi12 lt 30 then bmi12cat=3;
    else if 30 ge bmi12 then bmi12cat=4;
  if bmi12=. then bmi12cat=.; /*no one is missing*/

  bmicatmlvs=bmi12cat;


  /*smoking*/
  /*Wenjie Ma updated on 03/2018: smoking variables in 2012*/
  if smk12=2 then smoke12=4;                      /* if currently smoke yes */
  else if smk12=1 then do;
          if smoke10=1 then smoke12=1;            /* if not now & never */
          else if smoke10=2 then smoke12=2;       /* if not now & past yes */
          else if smoke10=3 then smoke12=3;       /* if not now & past uncertain*/
          else if smoke10=4 then smoke12=2;       /* if not now & yes in 08 */
          else if smoke10=5 then do;
                  if smoke08=2 | smoke08=4 |
                     smoke06=2 | smoke06=4 |
                     smoke04=2 | smoke04=4 |
                     smoke02=2 | smoke02=4 |
                     smoke00=2 | smoke00=4 |
                     smoke98=2 | smoke98=4 |
                     smoke96=2 | smoke96=4 |
                     smoke94=2 | smoke94=4 |
                     smoke92=2 | smoke92=4 |
                     smoke90=2 | smoke90=4 |
                     smoke88=2 | smoke88=4 |
                     smoke86=2 | smoke86=4 then smoke12=2;
                  else smoke12=3;
          end;
  end;
  else if smk12=3 then smoke12=5;                 /* if currently missing status*/
  else smoke12=5;


  if smoke12=4 then smkcurrent=1; else smkcurrent=0;

  /***ethnicity data****/
  if seuro86=1 then ethnic =1;
  else if scand86=1 then ethnic =1;
  else if ocauc86=1 then ethnic =1;
  else if afric86=1 then ethnic =2;
  else if asian86=1 then ethnic =3;
  else if oanc86=1  then ethnic =4;
  else ethnic=.;
  %indic3(vbl=ethnic, prefix=ethnic, reflev=1, min=2, max=4, missing=., usemiss=1,
  label1='causcasian',
  label2='african',
  label3='asian',
  label4='other');
  if ethnic>3 then eth3g=3;
  else             eth3g=ethnic;
  %indic3(vbl=eth3g, prefix=eth3g, reflev=1, min=2, max=3, missing=., usemiss=1,
  label1='causcasian',
  label2='african',
  label3='asian & others');
  if eth3g=1 then white=1; else white=0; if ethnic=. then white=.;
  %indic3(vbl=white, prefix=white, reflev=0, min=1, max=1, missing=., usemiss=1);

run;

proc rank data=mlvs_exposure out=mlvs_exposure groups=5;
  var bmicatmlvs;
  ranks bmicatmlvsq;
run;


data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

/*ppi*/

  if flag04 eq 1 then do;
    if pril04=1 then ppi04=1;
      else ppi04=0;
    if tag04=1 then h2ra04=1;
      else h2ra04=0;
  end;

  if flag06 eq 1 then do;
    if pril06=1 then ppi06=1;
      else ppi06=0;
    if tag06=1 then h2ra06=1;
      else h2ra06=0;
  end;

  if flag08 eq 1 then do;
    if pril08=1 then ppi08=1;
      else ppi08=0;
    if tag08=1 then h2ra08=1;
      else h2ra08=0;
  end;

  if flag10 eq 1 then do;
    if pril10=1 then ppi10=1;
      else ppi10=0;
    if tag10=1 then h2ra10=1;
      else h2ra10=0;
  end;

  array ppi_cur{4} ppi04 ppi06 ppi08 ppi10;

  do k=4 to 2 by -1;
    if ppi_cur{k}=. and ppi_cur{k-1} ne . then ppi_cur{k}=ppi_cur{k-1};
  end;

	ppi=ppi_cur{k};
	%indic3(vbl=ppi, prefix=ppi, min=1, max=1, reflev=0, missing=., usemiss=0, label1='current PPI use');

run;

/*******************************************************************************
/*******************************************************************************
/*******************************************************************************

                        STOOL QUESTIONNAIRE RESPONSES
       D/W DR. HUTTENHOWER RE: INCLUSION OF WHICH METADATUM TO INCLUDE

******************************TO BE INCLUDED************************************

16 acid_med_     Num    8 ACID_MED_2MO_.             BEST32.  11) In the past 2 months, have you used any acid reducing
2mo_qu1                                                    medications chronically (>1 week), including
                   Omeprazole,Protonix,Esomeprazole,Lansoprazole,Dexlansoprazole,
                   Ranitidine,Famotidine,Nizatidine ,Cimetidine?
56 acid_med_     Num    8 ACID_MED_2MO_.             BEST32.  11) In the past 2 months, have you used any acid reducing
2mo_qu2                                                    medications chronically (>1 week), including
                   Omeprazole,Protonix,Esomeprazole,Lansoprazole,Dexlansoprazole,
                   Ranitidine,Famotidine,Nizatidine ,Cimetidine?
*****
15 alc_qu1       Num    8 ALC_.                      BEST32.  10) How often do you consume alcoholic beverages?
55 alc_qu2       Num    8 ALC_.                      BEST32.  10) How often do you consume alcoholic beverages?
*****
4 ant_12mo_qu1  Num    8 ANT_12MO_.                 BEST32.  1) In the past 12 months, please check if you have received any
of the following medications in pill or intravenous (IV) form (DO
NOT INCLUDE INHALERS OR TOPICAL CREAMS): (choice=antibiotics)
44 ant_12mo_qu2  Num    8 ANT_12MO_.                 BEST32.  1) In the past 12 months, please check if you have received any
of the following medications in pill or intravenous (IV) form (DO
NOT INCLUDE INHALERS OR TOPICAL CREAMS): (choice=antibiotics)
*****
7 colsc_2mo_qu1 Num    8 COLSC_2MO_.                BEST32.  2) In the past 2 months, have you undergone a colonoscopy
                   or other procedure requiring bowel preparation?
47 colsc_2mo_qu2 Num    8 COLSC_2MO_.                BEST32.  2) In the past 2 months, have you undergone a colonoscopy
                   or other procedure requiring bowel preparation?
*****
14 dietpref_qu1  Num    8 DIETPREF_.                 BEST32.  9) What are your dietary preferences with respect to meat?
54 dietpref_qu2  Num    8 DIETPREF_.                 BEST32.  9) What are your dietary preferences with respect to meat?
*****
12 probio_       Num    8 PROBIO_2MO_.               BEST32.  7) In the past 2 months, have you consumed any probiotics
2mo_qu1                                                    (other than yogurt; see below) at least once per week?
52 probio_       Num    8 PROBIO_2MO_.               BEST32.  7) In the past 2 months, have you consumed any probiotics
2mo_qu2
*****
40 stool_type_   Num    8                                     3) Based on the attached Bristol chart, what was
batch1_qu1                                                 the appearance of the stool at the site you collected
                   the sample to put into the tube (SET 1)?
41 stool_type_   Num    8                                     3) Based on the attached Bristol chart, what was
batch2_qu1                                                 the appearance of the stool at the site you collected
                   the sample to put into the tube (SET 2)?
80 stool_type_   Num    8                                     3) Based on the attached Bristol chart, what was
batch3_qu2                                                 the appearance of the stool at the site you collected
                   the sample to put into the tube (SET 1)?
81 stool_type_   Num    8                                     3) Based on the attached Bristol chart, what was
batch4_qu2                                                 the appearance of the stool at the site you collected
                   the sample to put into the tube (SET 2)?
*****
30 tthhealt_qu1  Num    8 TTHHEALT_.                 BEST32.  18) Overall, how would you rate the health of your teeth and gums?
70 tthhealt_qu2  Num    8 TTHHEALT_.                 BEST32.  18) Overall, how would you rate the health of your teeth and gums?
*****
13 yog_2mo_qu1   Num    8 YOG_2MO_.                  BEST32.  8) In the past 2 months, how often have you consumed
                   yogurt or other foods containing active bacterial
                   cultures (kefir, sauerkraut, etc.)?
53 yog_2mo_qu2   Num    8 YOG_2MO_.                  BEST32.  8) In the past 2 months, how often have you consumed
                   yogurt or other foods containing active bacterial
                   cultures (kefir, sauerkraut, etc.)?


*********************OMITTED DUE TO NOT ENOUGH VARIATION************************

11 adiarr_       Num    8 ADIARR_2MO_.               BEST32.  6) In the past 2 months, have you had any acute diarrhea illnesses?
2mo_qu1
51 adiarr_       Num    8 ADIARR_2MO_.               BEST32.  6) In the past 2 months, have you had any acute diarrhea illnesses?
2mo_qu2
*****
22 alcmwash_qu1  Num    8 ALCMWASH_.                 BEST32.  Alcohol-based mouthwash
62 alcmwash_qu2  Num    8 ALCMWASH_.                 BEST32.  Alcohol-based mouthwash
*****
17 bile_med_     Num    8 BILE_MED_2MO_.             BEST32.  12) In the past 2 months, have you used any medications
2mo_qu1                                                    modifying bile acid production chronically (>1 week),
                   includingCholestyramine, cholestipol, colesevelam,
                   chenodeoxycholic acids, or ursodeoxycholic acid?
57 bile_med_     Num    8 BILE_MED_2MO_.             BEST32.  12) In the past 2 months, have you used any medications
2mo_qu2                                                    modifying bile acid production chronically (>1 week),
                   includingCholestyramine, cholestipol, colesevelam,
                   chenodeoxycholic acids, or ursodeoxycholic acid?
*****
39 brief_        Num    8 BRIEF_LIFESTYLE_QUES_V_1_. BEST32.  Complete?
lifestyle_
ques_v_1_qu1
79 brief_        Num    8 BRIEF_LIFESTYLE_QUES_V_1_. BEST32.  Complete?
lifestyle_
ques_v_1_qu2
******
9 cdiarr_       Num    8 CDIARR_2MO_.               BEST32.  4) In the past 2 months, have you had chronic diarrhea?
2mo_qu1
49 cdiarr_       Num    8 CDIARR_2MO_.               BEST32.  4) In the past 2 months, have you had chronic diarrhea?
2mo_qu2
*****
5 chm_12mo_qu1  Num    8 CHM_12MO_.                 BEST32.  1) In the past 12 months, please check if you have received any
                   of the following medications in pill or intravenous (IV) form (DO
                   NOT INCLUDE INHALERS OR TOPICAL CREAMS): (choice=chemotherapy)
45 chm_12mo_qu2  Num    8 CHM_12MO_.                 BEST32.  1) In the past 12 months, please check if you have received any
                   of the following medications in pill or intravenous (IV) form (DO
                   NOT INCLUDE INHALERS OR TOPICAL CREAMS): (choice=chemotherapy)
*****
2 complete_qu1  Num    8 COMPLETE_.                 BEST32.  Complete
42 complete_qu2  Num    8 COMPLETE_.                 BEST32.  Complete
*****
8 ctscan_       Num    8 CTSCAN_2MO_.               BEST32.  3) In the past 2 months, have you ingested
2mo_qu1                                                    a contrast agent for a CT scan or x-ray?
48 ctscan_       Num    8 CTSCAN_2MO_.               BEST32.  3) In the past 2 months, have you ingested
2mo_qu2                                                    a contrast agent for a CT scan or x-ray?
*****
18 dentist_qu1   Num    8 DENTIST_.                  BEST32.  13) The last time you saw a dentist (generalist, prosthodontist,
                   periodontist, etc.) was within the past:
58 dentist_qu2   Num    8 DENTIST_.                  BEST32.  13) The last time you saw a dentist (generalist, prosthodontist,
                   periodontist, etc.) was within the past:
*****
28 denture_      Num    8 DENTURE_CLEAN_.            BEST32.  16) If you wear dentures, how often do you clean them?
clean_qu1
68 denture_      Num    8 DENTURE_CLEAN_.            BEST32.  16) If you wear dentures, how often do you clean them?
clean_qu2
*****
21 etbrush_qu1   Num    8 ETBRUSH_.                  BEST32.  Electric toothbrush
61 etbrush_qu2   Num    8 ETBRUSH_.                  BEST32.  Electric toothbrush
******
24 floss_qu1     Num    8 FLOSS_.                    BEST32.  Floss
64 floss_qu2     Num    8 FLOSS_.                    BEST32.  Floss
******
35 freq_floss_   Num    8 BEST12.                    BEST32.  23) Aside from brushing your teeth with a toothbrush,
7days_qu1                                                  in the last 7 days, how many times did you use dental
                   floss or any other device to clean between your teeth?
75 freq_floss_   Num    8 BEST12.                    BEST32.  23) Aside from brushing your teeth with a toothbrush,
7days_qu2                                                  in the last 7 days, how many times did you use dental
                   floss or any other device to clean between your teeth?
******
36 freq_wash_    Num    8 BEST12.                    BEST32.  24) Aside from brushing your teeth with a toothbrush,
7days_qu1                                                  in the last 7 days, how many times did you use
                   mouthwash or any other dental rinse product?
76 freq_wash_    Num    8 BEST12.                    BEST32.  24) Aside from brushing your teeth with a toothbrush,
7days_qu2                                                  in the last 7 days, how many times did you use
                   mouthwash or any other dental rinse product?
*****
29 gumdis_qu1    Num    8 GUMDIS_.                   BEST32.  Gum disease (or periodontal disease) is a common problem
                   with the mouth. People with gum disease often have swollen,
                   receding, sore, or infected gums, and/or loose teeth.
                   17) Do you think you might have gum disease?
69 gumdis_qu2    Num    8 GUMDIS_.                   BEST32.  Gum disease (or periodontal disease) is a common problem
                   with the mouth. People with gum disease often have swollen,
                   receding, sore, or infected gums, and/or loose teeth.
                   17) Do you think you might have gum disease?
*****
10 hsp_2mo_qu1   Num    8 HSP_2MO_.                  BEST32.  5) In the past 2 months, have you been hospitalized for any reason?
50 hsp_2mo_qu2   Num    8 HSP_2MO_.                  BEST32.  5) In the past 2 months, have you been hospitalized for any reason?
*****
19 hygenist_qu1  Num    8 HYGENIST_.                 BEST32.  14) The last time you had a professional dental cleaning
                   (dentist, hygienist, etc.) was within the past:
59 hygenist_qu2  Num    8 HYGENIST_.                 BEST32.  14) The last time you had a professional dental cleaning
                   (dentist, hygienist, etc.) was within the past:
*****
1 id            Num    8 BEST12.                    BEST32.
*****
6 ims_12mon_qu1 Num    8 IMS_12MON_.                BEST32.  1) In the past 12 months, please check if you have received
                   any of the following medications in pill or intravenous
                   (IV) form (DO NOT INCLUDE INHALERS OR TOPICAL CREAMS):
                   (choice=immunosuppressants (e.g. oral corticosteroids))
46 ims_12mon_qu2 Num    8 IMS_12MON_.                BEST32.  1) In the past 12 months, please check if you have received
                   any of the following medications in pill or intravenous
                   (IV) form (DO NOT INCLUDE INHALERS OR TOPICAL CREAMS):
                   (choice=immunosuppressants (e.g. oral corticosteroids))
*****
33 lostbone_qu1  Num    8 LOSTBONE_.                 BEST32.  21) Have you ever been told by a dental professional
                   that you lost bone around your teeth?
73 lostbone_qu2  Num    8 LOSTBONE_.                 BEST32.  21) Have you ever been told by a dental professional
                   that you lost bone around your teeth?
*****
20 mtbrush_qu1   Num    8 MTBRUSH_.                  BEST32.  Manual toothbrush
60 mtbrush_qu2   Num    8 MTBRUSH_.                  BEST32.  Manual toothbrush
******
23 nonalcmwash_  Num    8 NONALCMWASH_.              BEST32.  Non-alcoholic mouthwash
qu1
63 nonalcmwash_  Num    8 NONALCMWASH_.              BEST32.  Non-alcoholic mouthwash
qu2
*****
37 perd_qu1      Num    8 PERD_.                     BEST32.  26) Have you been previously diagnosed with
                   periodontal disease with bone loss?
77 perd_qu2      Num    8 PERD_.                     BEST32.  26) Have you been previously diagnosed with
                   periodontal disease with bone loss?
*****                                                  (other than yogurt; see below) at least once per week?
3 sampleid_qu1  Num    8 SAMPLEID_.                 BEST32.  Sample ID
43 sampleid_qu2  Num    8 SAMPLEID_.                 BEST32.  Sample ID
*****
38 teeth_qu1     Num    8 TEETH_.                    BEST32.  25) How many natural teeth do you have?
78 teeth_qu2     Num    8 TEETH_.                    BEST32.  25) How many natural teeth do you have?
*****
26 tongueclean_  Num    8 TONGUECLEAN_.              BEST32.  Tongue cleaner
qu1
66 tongueclean_  Num    8 TONGUECLEAN_.              BEST32.  Tongue cleaner
qu2
*****
32 tooth_loose_  Num    8 TOOTH_LOOSE_.              BEST32.  20) Have you ever had any teeth become
qu1                                                        loose on their own without an injury?
72 tooth_loose_  Num    8 TOOTH_LOOSE_.              BEST32.  20) Have you ever had any teeth become
qu2                                                        loose on their own without an injury?
******
27 toothwhite_   Num    8 TOOTHWHITE_.               BEST32.  Tooth-whiteners
qu1
67 toothwhite_   Num    8 TOOTHWHITE_.               BEST32.  Tooth-whiteners
qu2
******
31 trt_gum_qu1   Num    8 TRT_GUM_.                  BEST32.  19) Have you ever had treatment for gum disease, such as
                   scaling and root planing, sometimes called deep cleaning?
71 trt_gum_qu2   Num    8 TRT_GUM_.                  BEST32.  19) Have you ever had treatment for gum disease, such as
                   scaling and root planing, sometimes called deep cleaning?
*****
34 tth_notright_ Num    8 TTH_NOTRIGHT_.             BEST32.  22) During the past 3 months, have you noticed
qu1                                                        a tooth that doesn't look quite right?
74 tth_notright_ Num    8 TTH_NOTRIGHT_.             BEST32.  22) During the past 3 months, have you noticed
qu2                                                        a tooth that doesn't look quite right?
*****
25 waterpick_qu1 Num    8 WATERPICK_.                BEST32.  Water-based pick/jet
65 waterpick_qu2 Num    8 WATERPICK_.                BEST32.  Water-based pick/jet


********************************************************************************
*******************************************************************************/

data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

  /*self-reported acid suppression*/

  if acid_med_2mo_qu1 or acid_med_2mo_qu2 eq 1 then acid_avg=1;
    else acid_avg=0;

  if acid_med_2mo_qu1 eq 1 then acid_w1=1;
    else acid_w1=0;

  if acid_med_2mo_qu2 eq 1 then acid_w2=1;
    else acid_w2=0;

  /* alcohol 1	Never 2	Rarely 3	1-6 times a week 4	Daily 5	More than daily */

  alc_avg = max (alc_qu1, alc_qu2);

  alc_w1 = alc_qu1;
  alc_w2 = alc_qu2;

  /*antibiotics in the past 12 months, any=1*/

  if ant_12mo_qu1 or ant_12mo_qu2 eq 1 then abx_avg=1;
    else abx_avg=0;

  abx_w1=ant_12mo_qu1;

  abx_w2=ant_12mo_qu2;

  /*colo prep*/
  if colsc_2mo_qu1 or colsc_2mo_qu2 eq 1 then prep_avg=1;
    else prep_avg=0;

  prep_w1=colsc_2mo_qu1;
  prep_w2=colsc_2mo_qu2;

  /*
  self-diet
  1	Standard diet
  2	Standard diet with poultry and/or fish (no red meat)
  3	Vegetarian (no meat)
  4	Vegan (no meat, dairy, or animal products)
  */

  if dietpref_qu1 or dietpref_qu2 eq 1 then selfdiet_avg=0;
    else selfdiet_avg=1; /*lump all others as non-standard diet*/

  if dietpref_qu1 then selfdiet_w1=0;
    else selfdiet_w1=1;

  if dietpref_qu2 then selfdiet_w2=0;
    else selfdiet_w2=1;

  /*probiotic*/
  if probio_2mo_qu1 or probio_2mo_qu2 eq 1 then probx_avg=1;
    else probx_avg=0;

  if probio_2mo_qu1 or probio_2mo_qu2 eq 1 then probx_avg=1;
    else probx_avg=0;

  probx_w1=probio_2mo_qu1;
  probx_w2=probio_2mo_qu2;

  /*bristol stool type: grade HIGHEST stool, as in loosest, out of a given week*/
  bristol5=stool_type_batch1_qu1;
  bristol6=stool_type_batch2_qu1;
  bristol7=stool_type_batch3_qu2;
  bristol8=stool_type_batch4_qu2;

  if bristol5 in (3,4,5) then bristolcat5=0; /*normal stool*/
    else if bristol5 in (1,2) then bristolcat5=1; /*hard stool*/
    else if bristol5 in (6,7) then bristolcat5=2; /*loose stool*/
    else bristolcat5=.; /*missing*/

  if bristol6 in (3,4,5) then bristolcat6=0; /*normal stool*/
    else if bristol6 in (1,2) then bristolcat6=1; /*hard stool*/
    else if bristol6 in (6,7) then bristolcat6=2; /*loose stool*/
    else bristolcat6=.; /*missing*/

  if bristol7 in (3,4,5) then bristolcat7=0; /*normal stool*/
    else if bristol7 in (1,2) then bristolcat7=1; /*hard stool*/
    else if bristol7 in (6,7) then bristolcat7=2; /*loose stool*/
    else bristolcat7=.; /*missing*/

  if bristol8 in (3,4,5) then bristolcat8=0; /*normal stool*/
    else if bristol8 in (1,2) then bristolcat8=1; /*hard stool*/
    else if bristol8 in (6,7) then bristolcat8=2; /*loose stool*/
    else bristolcat8=.; /*missing*/

  bristol_w1=mean(bristol5, bristol6);
  if 3<=bristol_w1<=5 then bristolcat_w1=0; /*normal stool*/
    else if 1<=bristol_w1<3 then bristolcat_w1=1; /*hard stool*/
    else if 5<bristol_w1<=7 then bristolcat_w1=2; /*loose stool*/
    else bristolcat_w1=.; /*missing*/

  bristol_w2=mean(bristol7, bristol8);
  if 3<=bristol_w2<=5 then bristolcat_w2=0; /*normal stool*/
    else if 1<=bristol_w2<3 then bristolcat_w2=1; /*hard stool*/
    else if 5<bristol_w2<=7 then bristolcat_w2=2; /*loose stool*/
    else bristolcat_w2=.; /*missing*/

  bristol_avg=mean(bristol5, bristol6, bristol7, bristol8);
  if 3<=bristol_avg<=5 then bristolcat_avg=0; /*normal stool*/
    else if 1<=bristol_avg<3 then bristolcat_avg=1; /*hard stool*/
    else if 5<bristol_avg<=7 then bristolcat_avg=2; /*loose stool*/
    else bristolcat_avg=.; /*missing*/


  /*oral health 1Excellent 2	Very good 3	Good 4	Fair 5	Poor 6	Don't know*/
  if tthhealt_qu1 or tthhealt_qu2 eq 1 then oral_avg=0; /*ref=excellent oral health*/
    else oral_avg=1; /*all others*/

  if tthhealt_qu1 eq 1 then oral_w1=0; /*ref=excellent oral health*/
    else oral_w1=1; /*all others*/

  if tthhealt_qu2 eq 1 then oral_w2=0; /*ref=excellent oral health*/
    else oral_w2=1; /*all others*/

  /*yogurt 1	Never 2	Rarely 3	1-6 times a week 4	Daily 5	More than daily*/
  if yog_2mo_qu1 or yog_2mo_qu2 in (3,4,5) then yogurt_avg=1; /* >weekly yogurt */
    else yogurt_avg=0;

  if yog_2mo_qu1 in (3,4,5) then yogurt_w1=1; /* >weekly yogurt */
    else yogurt_w1=0;

  if yog_2mo_qu2 in (3,4,5) then yogurt_w2=1; /* >weekly yogurt */
    else yogurt_w2=0;


run;


/**********************************
*         CRP and lipids        *
***********************************/

data mlvs_exposure;
  set mlvs_exposure
  end=_end_;

/*if crp_plasma1=. then delete;*/
if crp_plasma1 NE . and crp_plasma2 NE . then crpchg=crp_plasma2-crp_plasma1;
crp_avg=mean(crp_plasma1, crp_plasma2);


if crp_plasma1 NE . and CRP_plasma1 NE 0 then logcrp_plasma1=log(crp_plasma1);
if crp_plasma2 NE . and CRP_plasma2 NE 0 then logcrp_plasma2=log(crp_plasma2);

if 0 < crp_plasma1 < 1 then crp_plasma1cat=1;
else if 1 <= crp_plasma1 < 3 then crp_plasma1cat=2;
else if 3 <= crp_plasma1 then crp_plasma1cat=3;
if crp_plasma1 eq . then crp_plasma1cat=.;

if 0 < crp_plasma2 < 1 then crp_plasma2cat=1;
else if 1 <= crp_plasma2 < 3 then crp_plasma2cat=2;
else if 3<= crp_plasma2  then crp_plasma2cat=3;
if crp_plasma2 eq . then crp_plasma2cat=.;

if hdlc_plasma1 NE . and hdlc_plasma1 NE 0 then tchdl_plasma1=tc_plasma1/hdlc_plasma1;
if hdlc_plasma2 NE . and hdlc_plasma2 NE 0 then tchdl_plasma2=tc_plasma2/hdlc_plasma2;

if tchdl_plasma1 NE . and tchdl_plasma1 NE 0 then logtchdl_plasma1=log(tchdl_plasma1);
if tchdl_plasma2 NE . and tchdl_plasma2 NE 0 then logtchdl_plasma2=log(tchdl_plasma2);

if tc_plasma1 NE . and tc_plasma1 NE 0 then logtc_plasma1=log(tc_plasma1);
if tc_plasma2 NE . and tc_plasma2 NE 0 then logtc_plasma2=log(tc_plasma2);

if hdlc_plasma1 NE . and hdlc_plasma1 NE 0 then loghdl_plasma1=log(hdlc_plasma1);
if hdlc_plasma2 NE . and hdlc_plasma2 NE 0 then loghdl_plasma2=log(hdlc_plasma2);

if tg_plasma1 NE . and tg_plasma1 NE 0 then logtg_plasma1=log(tg_plasma1);
if tg_plasma2 NE . and tg_plasma2 NE 0 then logtg_plasma2=log(tg_plasma2);

run;

proc means data=mlvs_exposure n nmiss mean median min max;
var crp_plasma1 crp_plasma2 crp_avg crpchg logcrp_plasma1 logcrp_plasma2
tc_plasma1 tc_plasma2 hdlc_plasma1 hdlc_plasma2 tchdl_plasma1 tchdl_plasma2 tg_plasma1 tg_plasma2 logtchdl_plasma1 logtchdl_plasma2
logtc_plasma1 logtc_plasma2 logtg_plasma1 logtg_plasma2 loghdl_plasma1 loghdl_plasma2;
run;


proc freq data=mlvs_exposure;
tables crp_plasma1cat crp_plasma2cat;
run;



/*******************************************************************************
/*******************************************************************************
              KEPT VARIABLES/EXPORTING TXTs FOR MICROBIOME PIPELINES

********************************************************************************
*******************************************************************************/

  /*export all*/

  /*
  legacy variables i don't need. need to limit export because
  sas will truncate when header + data exceeds the maximum length of a SAS statement (32767)
  */

  data mlvs_exposure;
    set mlvs_exposure (drop = act00 act02 act04 act06 act86 act88 act90 act92 act94 act96 act98 age86 alco02n
    alco06n alco86n alco90n alco94n alco98n ang00 ang02 ang04 ang06 ang08 ang10
    ang86 ang88 ang90 ang92 ang94 ang96 ang98 asp00 asp02 asp04 asp06 asp08 asp86
    asp88 asp90 asp92 asp94 aspw98 bmi00 bmi02 bmi04 bmi06 bmi08 bmi86 bmi88 bmi90
    bmi92 bmi94 bmi96 bmi98 cabg00 cabg02 cabg04 cabg06 cabg08 cabg10 cabg86 cabg88
    cabg90 cabg92 cabg94 cabg96 cabg98 can86 can88 can90 canc00 canc02 canc04 canc06
    canc08 canc10 canc92 canc94 canc96 canc98 db00 db02 db04 db06 db08 db10 db86
    db88 db90f db92 db94 db96 db98 dbmy09 hbp00 hbp02 hbp04 hbp06 hbp08 hbp10 hbp86
    hbp88 hbp90 hbp92 hbp94 hbp96 hbp98 mi00 mi02 mi04 mi06 mi08 mi10 mi86 mi88 mi90
    mi92 mi94 mi96 mi98 motrn00 motrn86 motrn88 motrn90 motrn92 motrn94 motrn96
    motrn98 mvt00 mvt02 mvt04 mvt06 mvt08 mvt10 mvt88 mvt92 mvt96 pril04 pril06
    pril08 pril10 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr86
    rtmnyr87 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 smoke00 smoke02
    smoke04 smoke06 smoke08 smoke86 smoke88 smoke90 smoke92 smoke94 smoke96 smoke98
    str86 str88 str90 str92 str94 str96 str98 strk00 strk02 strk04 strk06 strk08
    strk10 tag00 tag02 tag04 tag06 tag08 tag10 tag86 tag88 tag90 tag92 tag94 tag96
    tag98 wt00 wt02 wt04 wt06 wt08 wt86 wt88 wt90 wt92 wt94 wt96 wt98
    chol86	chol88	chol90 alcf86	alc88	chol92	chol94	chol96	nasp96	chol98	nsaid98
    chol00	nsaid00	chol02	mtrn02	chol04	mtrn04	flag04 chol06	mtrn06	flag06	chol08	mtrn08	flag08
    mulfrq_ffq1 brndpt_ffq1 centrum_silver_ffq1 theragran_m_ffq1 centrum_ffq1 oneaday_essent_ffq1 cvs_mv_fe_ffq1 cvs_mvm_ffq1
    krkland_dly_multv_ffq1 kirk_pak_ffq1 nm_ess_bal_ffq1 double_x_ffq1 onesource_multv_ffq1 mm_2000_hb_mvm_ffq1
    shaklee_vitalea_ffq1 trader_joe_ffq1 mvwithiron_ffq1 mvwithmin_ffq1 other_pak_ffq1 seniorvit_ffq1 gnc_ultra_wom_ffq1
    mel_vita_ffq1 oneaday_within_ffq1 mvwomen_ffq1 cvs_mv_ffq1 flintstones_compl_ffq1 occuvite_ffq1 icaps_plus_ffq1
    oneaday_wgt_ffq1 stresstab_ffq1 stuart_natal_ffq1 mulfrq_ffq2 brndpt_ffq2 centrum_silver_ffq2 theragran_m_ffq2 centrum_ffq2
    oneaday_essent_ffq2 cvs_mv_fe_ffq2 cvs_mvm_ffq2 krkland_dly_multv_ffq2 kirk_pak_ffq2 nm_ess_bal_ffq2 double_x_ffq2
    onesource_multv_ffq2 mm_2000_hb_mvm_ffq2 shaklee_vitalea_ffq2 trader_joe_ffq2 mvwithiron_ffq2 mvwithmin_ffq2 other_pak_ffq2
    seniorvit_ffq2 gnc_ultra_wom_ffq2 mel_vita_ffq2 oneaday_within_ffq2 mvwomen_ffq2 cvs_mv_ffq2 flintstones_compl_ffq2
    occuvite_ffq2 icaps_plus_ffq2 oneaday_wgt_ffq2 stresstab_ffq2 stuart_natal_ffq2
    seuro86 scand86 ocauc86 afric86 asian86 oanc86
    );
  run;


  proc export data=mlvs_exposure REPLACE
	outfile="mlvs_exposure_for_jorick.txt";

proc means data=mlvs_exposure n nmiss mean min max;
  var aofib10a aofib10v aofib_avg aofib_ffq1 aofib_ffq2
      a_aofib_fs_dr_w1avg a_aofib_fs_dr_w2avg a_aofib_fo_dr_w1avg a_aofib_fo_dr_w2avg
      a_pect_fo_dr_w1avg a_pect_fo_dr_w2avg
      calor10n calor10v calor_avg calor_ffq1 calor_ffq2 calor_fs_dr_w1avg calor_fs_dr_w2avg
      frtaf10a frtaf10v ceraf10a ceraf10v vegaf10a vegaf10v
      ala_ffq1 ala_ffq2 ala_avg epa_ffq1 epa_ffq2 epa_avg dha_ffq1 dha_ffq2 dha_avg dpa_ffq1 dpa_ffq2 dpa_avg
      ala86v ala90v ala94v ala98v ala02v ala06v ala10v
      epa86v epa90v epa94v epa98v epa02v epa06v epa10v
      dha86v dha90v dha94v dha98v dha02v dha06v dha10v
                                  dpa02v dpa06v dpa10v
      a_omega3_fs_dr_w1avg a_omega3_fs_dr_w2avg

      ;
run;

endsas;

proc means data=mlvs_exposure n nmiss mean min max;
  var bristol_w1 bristol_w2 bristol_avg;
run;

proc freq data=mlvs_exposure;
  tables bristolcat5 bristolcat6 bristolcat7 bristolcat8 bristolcat_w1 bristolcat_w2 ;
run;

proc corr data=mlvs_exposure spearman;
  var aofib10v ceraf10v frtaf10v vegaf10v a_aofib_fs_dr_w1avg a_aofib_fs_dr_w2avg a_pect_fo_dr_w1avg a_pect_fo_dr_w2avg;
run;

proc corr data=mlvs_exposure spearman;
  var a_aofib_fs_dr_w1avg a_aofib_fs_dr_w2avg a_pect_fo_dr_w1avg a_pect_fo_dr_w2avg a_ifib_fs_dr_w1avg a_ifib_fs_dr_w2avg a_wsdf_fo_dr_w1avg a_wsdf_fo_dr_w2avg;
run;


proc means data=mlvs_exposure n nmiss mean min max;
  var bmi12 act12;
run;

endsas;


/*******************************************************************************
/*******************************************************************************
*                           Table 1                     *
********************************************************************************
*******************************************************************************/


data table1;
  set mlvs_exposure;
  keep id agemlvs white bmi12 alc_w1 alc_w2 alc_avg aofib10a aofib10v ceraf10v frtaf10v vegaf10v aofib_avg aofib_ffq1 aofib_ffq2 a_aofib_fs_dr_w1avg a_aofib_fs_dr_w2avg
                  calor10n calor10v calor_avg calor_ffq1 calor_ffq2 calor_fs_dr_w1avg calor_fs_dr_w2avg
                  crp_plasma1 crp_plasma2 tc_plasma1 tc_plasma2 hdlc_plasma1 hdlc_plasma2 tchdl_plasma1 tchdl_plasma2 tg_plasma1 tg_plasma2
                  bmic12 act12 smkcurrent abx_w1 abx_w2 abx_avg probx_w1 probx_w2 probx_avg bristolcat_w1 bristolcat_w2 bristolcat_avg prep_w1 prep_w2 prep_avg;
bristolcat_w1=bristolcat_w1+1;
bristolcat_w2=bristolcat_w2+1;
bristolcat_avg=bristolcat_avg+1;
run;

proc format;
value alc
      1='Never drinking'
      2='Rarely'
      3='1-6 times/wk'
      4='Daily'
      5='More than daily';
value bristolcat
      1='Normal stool'
      2='Hard stool'
      3='Loose stool';
run;

data table1;
set table1;
format alc_w1 alc.;
format alc_w2 alc.;
format alc_avg alc.;
format bristolcat_w1 bristolcat.;
format bristolcat_w2 bristolcat.;
format bristolcat_avg bristolcat.;
run;

proc rank data=table1 groups=4 out=table1;
 var a_aofib_fs_dr_w1avg a_aofib_fs_dr_w2avg aofib10v aofib_avg;
 ranks a_aofib_fs_dr_w1avgq a_aofib_fs_dr_w2avgq aofib10vq aofib_avgq;
run;

proc freq data=table1;
  tables alc_w1 alc_w2 alc_avg;
run;


%table1(data=table1,
        exposure=  a_aofib_fs_dr_w1avgq,
        varlist  = agemlvs bmi12 act12 alc_w1 alc_w2 calor_fs_dr_w1avg a_aofib_fs_dr_w1avg crp_plasma1 crp_plasma2 tc_plasma1 tc_plasma2 hdlc_plasma1 hdlc_plasma2 tchdl_plasma1 tchdl_plasma2 tg_plasma1 tg_plasma2 a_aofib_fs_dr_w1avgq
                  aofib_ffq1 aofib10v smkcurrent abx_w1 abx_w2 abx_avg  probx_w1 probx_w2 bristolcat_w1 bristolcat_w2 ceraf10v frtaf10v vegaf10v,
        cat      = smkcurrent abx_w1 abx_w2 abx_avg probx_w1 probx_w2,
        poly =  alc_w1 alc_w2 bristolcat_w1 bristolcat_w2,
        ageadj = F,
        rtftitle =  Baseline characteristics of participants according to quartiles of total fiber intake in MLVS,
        landscape=  F,
        file     =  result.table1.mlvs.aofibq,
        dec      =  1,
        uselbl   =  F);
run;


%table1(data=table1,
        noexp =  t,
        varlist  = agemlvs white bmi12 act12 alc_w1 alc_w2 alc_avg calor_fs_dr_w1avg a_aofib_fs_dr_w1avg crp_plasma1 crp_plasma2 tc_plasma1 tc_plasma2 hdlc_plasma1 hdlc_plasma2 tchdl_plasma1 tchdl_plasma2 tg_plasma1 tg_plasma2 a_aofib_fs_dr_w1avgq
                  aofib_ffq1 aofib10v smkcurrent abx_w1 abx_w2 abx_avg probx_w1 probx_w2 probx_avg bristolcat_w1 bristolcat_w2 bristolcat_avg prep_w1 prep_w2 prep_avg ceraf10v frtaf10v vegaf10v,
        cat      = white smkcurrent abx_w1 abx_w2 abx_avg probx_w1 probx_w2 probx_avg prep_w1 prep_w2 prep_avg,
        poly =  alc_w1 alc_w2 alc_avg bristolcat_w1 bristolcat_w2 bristolcat_avg,
        ageadj = F,
        rtftitle =  Baseline characteristics of participants in MLVS,
        landscape=  F,
        file     =  result.table1.mlvs,
        dec      =  1,
        uselbl   =  F);
run;


%macro glm(outcome);
proc glm data = table1;
*class a_aofib_fs_dr_w1avgq;
model &outcome = a_aofib_fs_dr_w1avgq;
run;
%mend;

%glm(agemlvs);
%glm(bmi12);
%glm(act12);
%glm(calor_fs_dr_w1avg);
%glm(a_aofib_fs_dr_w1avg);
%glm(crp_plasma1);
%glm(tc_plasma1);
%glm(hdlc_plasma1);
%glm(tg_plasma1);
%glm(ceraf10v);
%glm(frtaf10v);
%glm(vegaf10v);


%macro cat(outcome);
data one;
    set table1;
    if &outcome eq . then delete;
  run;

proc freq data=one;
  tables a_aofib_fs_dr_w1avgq*&outcome/cmh;
run;
%mend;

%cat(alc_w1);
%cat(abx_w1);
%cat(probx_w1);
%cat(bristolcat_w1);








/*******************************************************************************
/*******************************************************************************
*                           PROC DATASETS/FINAL OUTPUT                         *
********************************************************************************
*******************************************************************************/

proc datasets nolist;
	delete hp_der hp_der_2 hp86 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04
	hp06 hp08 hp10 hpfsvars hp_mets mets h86_nts h90_nts h94_nts h98_nts h02_nts
  h06_nts h10_nts nutrients ffq1 ffq2 ffq12 raw2 raw2 raw22 stoolqq1 stoolqq2
  stoolqq12 blood12
  out1 out2
  factorffq1 factorffq2 factorffq12 energyffq12
  ntsw12
  energyw12
run;
