/*** 1. Load Nutrient Intake Data  ***/

libname myxpt xport "/home/u63546026/sasuser.v94/STA_3013/DR1TOT_H.xpt";
proc copy in=myxpt out=work;
run;

libname myxpt xport "/home/u63546026/sasuser.v94/STA_3013/DR2TOT_H.xpt";
proc copy in=myxpt out=work;
run;

/*** 2. Load Cognitive Function Data ***/

libname myxpt xport "/home/u63546026/sasuser.v94/STA_3013/CFQ_H.xpt";
proc copy in=myxpt out=work;
run;


/*** 3. Factor Analysis on Nutrients ***/

/* a. Sort both days by ID */
proc sort data=work.dr1tot_h; by SEQN; run;
proc sort data=work.dr2tot_h; by SEQN; run;

/* b. Merge and Average using Arrays */
data work.nutrients_averaged;
    merge work.dr1tot_h(rename=(
            DR1TKCAL=d1_1  DR1TPROT=d1_2  DR1TCARB=d1_3  DR1TSUGR=d1_4  
            DR1TFIBE=d1_5  DR1TTFAT=d1_6  DR1TSFAT=d1_7  DR1TMFAT=d1_8  
            DR1TPFAT=d1_9  DR1TCHOL=d1_10 DR1TATOC=d1_11 DR1TVARA=d1_12 
            DR1TVB6=d1_13  DR1TFOLA=d1_14 DR1TCHL=d1_15  DR1TVB12=d1_16 
            DR1TVC=d1_17   DR1TVD=d1_18   DR1TVK=d1_19   DR1TCALC=d1_20 
            DR1TMAGN=d1_21 DR1TIRON=d1_22 DR1TZINC=d1_23 DR1TSODI=d1_24 
            DR1TPOTA=d1_25 DR1TCAFF=d1_26))
            
          work.dr2tot_h(rename=(
            DR2TKCAL=d2_1  DR2TPROT=d2_2  DR2TCARB=d2_3  DR2TSUGR=d2_4  
            DR2TFIBE=d2_5  DR2TTFAT=d2_6  DR2TSFAT=d2_7  DR2TMFAT=d2_8  
            DR2TPFAT=d2_9  DR2TCHOL=d2_10 DR2TATOC=d2_11 DR2TVARA=d2_12 
            DR2TVB6=d2_13  DR2TFOLA=d2_14 DR2TCHL=d2_15  DR2TVB12=d2_16 
            DR2TVC=d2_17   DR2TVD=d2_18   DR2TVK=d2_19   DR2TCALC=d2_20 
            DR2TMAGN=d2_21 DR2TIRON=d2_22 DR2TZINC=d2_23 DR2TSODI=d2_24 
            DR2TPOTA=d2_25 DR2TCAFF=d2_26));
    by SEQN;

/* c. Define arrays for Day 1, Day 2, and the new Averaged variables */
    array day1[26] d1_1-d1_26;
    array day2[26] d2_1-d2_26;
    array avg[26]  KCAL PROT CARB SUGR FIBE TFAT 
                   SFAT MFAT PFAT CHOL ATOC VARA 
                   VB6  FOLA CHL  VB12 VC   VD   
                   VK   CALC MAGN IRON ZINC SODI 
                   POTA CAFF;

    do i = 1 to 26;
        avg[i] = mean(day1[i], day2[i]);
    end;

    drop d1_: d2_: i;
    
label
    PROT = "Protein (gm)"
    CARB = "Carbohydrate (gm)"
    SUGR = "Total sugars (gm)"
    FIBE = "Dietary fiber (gm)"
    TFAT = "Total fat (gm)"
    SFAT = "Total saturated fatty acids (gm)"
    MFAT = "Total monounsaturated fatty acids (gm)"
    PFAT = "Total polyunsaturated fatty acids (gm)"
    CHOL = "Cholesterol (mg)"
    ATOC = "Vitamin E (mg)"
    VARA = "Vitamin A (mcg RAE)"
    VB6  = "Vitamin B6 (mg)"
    FOLA = "Total folate (mcg)"
    CHL  = "Total choline (mg)"
    VB12 = "Vitamin B12 (mcg)"
    VC   = "Vitamin C (mg)"
    VD   = "Vitamin D (mcg)"
    VK   = "Vitamin K (mcg)"
    CALC = "Calcium (mg)"
    MAGN = "Magnesium (mg)"
    IRON = "Iron (mg)"
    ZINC = "Zinc (mg)"
    SODI = "Sodium (mg)"
    POTA = "Potassium (mg)"
    CAFF = "Caffeine (mg)";
    
run;

/* d. Run Proc Factor on averaged variables */

options ls=75 ps=65 nodate nonumber;

/* Note on Data Cleaning: 
   Approximately 13% of observations (1,278 out of 9,813) were excluded 
   due to missing values in the 26 nutrient variables. The remaining 
   sample size (n=8,535) remains highly robust for factor extraction.
*/

proc factor data=work.nutrients_averaged 
        method=prin 
        res 
        nfact=6 
        scree 
        rotate=varimax
        score
        out=factor_scores
        prefix=Factor;
        
    var PROT CARB SUGR FIBE TFAT 
        SFAT MFAT PFAT CHOL ATOC VARA 
        VB6  FOLA CHL  VB12 VC   VD   
        VK   CALC MAGN IRON ZINC SODI 
        POTA CAFF;
run;



/*** 4. Standardize Cognitive Scores and Create Composite ***/

/* a. Standardize scores */

proc standard data=work.cfq_h mean=0 std=1 out=cfq_std;
    var CFDCST1 CFDCST2 CFDCST3 CFDCSR CFDAST CFDDS;
run;

/* b. Convert into T scores */

data cfq_std;
    set cfq_std;
    
    Z_Cognitive_Score = mean(of CFDCST1 CFDCST2 CFDCST3 CFDCSR CFDAST CFDDS);
    
    Cognitive_T_Score = (Z_Cognitive_Score * 10) + 50;

    label Cognitive_T_Score = "Global Cognitive Function (T-Score: Mean=50, SD=10)";
run;


/*** 5. Merge Cognitive Scores with Factor Scores ***/

proc sort data=factor_scores; by SEQN; run;
proc sort data=cfq_std; by SEQN; run;
proc sort data=work.nutrients_averaged; by SEQN; run;

data merged;
    merge cfq_std (in=a) factor_scores (in=b)
          work.nutrients_averaged(keep=SEQN KCAL);
    by SEQN;
    if nmiss(Cognitive_T_Score, KCAL, Factor1, Factor2, Factor3, Factor4, Factor5, Factor6) = 0;
run;



/*** 6. Label factors for descriptive regression output ***/


data final_analysis;
    set merged;
    
    label 
        Factor1 = "Animal-based / protein & fat pattern"
        Factor2 = "Whole Grains & Legumes"
        Factor3 = "Dairy & Fortified foods"
        Factor4 = "Fruits & Vegetables"
        Factor5 = "Sugar & Refined Carbs"
        Factor6 = "Caffeine & Stimulants";
run;


/*** 7. Run regression using labeled factors ***/
proc reg data=final_analysis;
    model Cognitive_T_Score = Factor1 Factor2 Factor3 Factor4 Factor5 Factor6 KCAL / stb;
    title "The Impact of Dietary Patterns on Global Cognitive Function (T-Scores)";
run;
quit;