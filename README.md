Replication Code for: "Optimal longevity of a dynasty"
======================================================
Authors: Satoshi Nakano and Kazuhiko Nishimura
Journal: Homo Oeconomicus

Description:
This repository contains the R script necessary to replicate the theoretical 
figures (Figure 2 through Figure 6) presented in the manuscript "Optimal 
longevity of a dynasty". The code numerically solves the finite-horizon 
dynamic programming problem for the critical-level utilitarian (CLU) 
framework under both the AK setting and the Zero Discounting (ZD) setting.

File Structure:
- Replication_Code.R : The master R script that generates all figures.
- README.txt         : This documentation file.

System Requirements & Dependencies:
The code was written and tested in R. To run the script, the following 
packages must be installed:
  - ggplot2
  - dplyr
  - DescTools  (for Gini index calculations)
  - gglorenz   (for Lorenz curve plotting)

Instructions:
1. Open "Replication_Code.R" in your R environment (e.g., RStudio).
2. Ensure that the required packages are installed. You can uncomment 
   the 'install.packages()' lines at the top of the script if needed.
3. Run the script. The code is divided into clear sections corresponding 
   to the figures in the paper.
4. The script will automatically generate and save the figures as .eps 
   files in your working directory.

Contact:
For any questions regarding the code, please contact Kazuhiko Nishimura 
at nishimura@lets.chukyo-u.ac.jp.
