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
    * /PH_GER: mat-files of model outputs  
  * /eval
    * /GER
      - /benchmark: B-AR(1) forecasts of GDP corresponding to the vintages in ../*/datasetsGER.mat 
    * /US
      - /benchmark: B-AR(1) forecasts of GDP corresponding to the vintages in ../*/datasetsUS.mat 
