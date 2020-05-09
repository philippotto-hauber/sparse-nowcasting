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
    * /GER (copied from HPC folder)
       - /level rec/.. result files
    * /US (copied from HPC folder)
       - /level rec/.. result files
    * /tables (linked): formatted evaluation tables that are linked to tables folders in /GER and /US
    * /latex_tables: tables in /GER and /US but in .tex format
    
HPC file management
/sfs/fs2/work-sh1/swwiw663/eval-hpc
 * /benchmark GER: benchmark B-AR models that were formerly in the respective country folder
 * /benchmark US: benchmark B-AR models that were formerly in the respective country folder
 * trueGDPUS.mat, truegdpGER.mat: GDP growth realizations (copied from data folder)
 * /GER 
       - /level rec/.. result files
 * /US 
       - /level rec/.. result files
      
 

