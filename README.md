# sparse-nowcasting
codes for project with Christian on nowcasting with sparse factor models
- -------------------------------------------------------------------------------

repo structure

* /data
* /estim
* /eval

local file management structure

C:/Users/ * / Disseration/sparse nowcasting
  * /data 
    * /in: 
      - FRED-MD vintages in folder 
      - vintages of US GDP (vintagesGDP.xls) 
      - Philly Fed BOS data (bos_historyx.xls
      - raw German data vintages (datassets_orig.mat)that are transformed into datasetsGER.mat 
      - vintages of German GDP (Vintages_GDP_GER.xls)
    * /out: 
      - datasetsUS.mat
      - datasetsGER.mat
      - truegdpGER.mat (first, second, final realizations of German GDP from 2006Q1 until 2018Q4)
      - truegdpUS.mat (first, second, final realizations of German GDP from 2000Q1 until 2018Q4)
    * /IRIS_Tbx: [IRIS toolbox](<https://iris.igpmn.org/>)
  * /estim
    * ineff_facs.Rda: data frame containing the inefficiency factors of the different models (calculated using the script in the repo /estim/calc_inefffac.R
    * plot_ineff_facs.pdf: figure plotting the inefficiency factors (output from the repo script /estim/plot_ineff_fac.R
  * /eval
    * /GER
      - /benchmark: B-AR(1) forecasts of GDP corresponding to the vintages in ../*/datasetsGER.mat
      - /level rec/.. result files
    * /US
      - /benchmark: B-AR(1) forecasts of GDP corresponding to the vintages in ../*/datasetsUS.mat 
      - /level rec/.. result files
    * /tables (linked): formatted evaluation tables that are linked to tables folders in /GER and /US
    * /latex_tables: tables in /GER and /US but in .tex format
  * /documentation (this folder contains files that arent used locally but stored for documentation purposes)
    * /PH_GER: mat files with predictive densities for German GDP
    * /PH_US: mat files with predictive densities for US GDP
    * /results eval mat files: evaluation mat files for different model specs (surveys in levels vs diffs; first, second release; Np = {1, 3}
