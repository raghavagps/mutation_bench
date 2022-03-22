# **Prediction of high-risk cancer patients using mutation profiles**
Benchmarking of mutation calling techniques by developing classification and regresion prediction models to predict the high-risk cancer patients.
## Introduction
In this method, a user can predict the high-risk cancer patients using their mutation profiles in the form of Variant Calling Format (VCF) and Mutation Annotation Format (MAF) derived using four widely used mutation calling techniques, such as, MuTect2, MuSE, VarScan2, and SomaticSniper. A comparison can be made between the formats or between the techniques for same formats. In case of classification, the whole dataset is divided into 80:20 ratio, where 80% is used for training purpose, on which five fold cross-validation is appiled while training the model, whereas 20% dataset is used as testing dataset, which is kept to test the trained model. On the other hand, in case of regression, five-fold cross validation is applied on whole dataset. This method provides the following seven files as output:

1. **Classification results:** This file contains the performance measures for the seven different classifiers such as Decision tree (DT), Support Vector Classifier (SVC), Random Forest (RF), XGBoost (XGB), Gaussian Naive Bayes (GNB), Logistic Regression (LR), and k-nearest neighbors (KN). The perfomance of each classifiers is measured in terms of sensitivity (Sens), specificity (Spec), accuracy (Acc), Area Under the Receiver Operating Characteristic (AUC) F1-score (F1), Kappa, and Matthews Correlation Coefficient (MCC), for training (tr) and testing (te) dataset.

2. **Regression results:** This file contains the performance measures for the seven different regressors such as Random Forest (RFR), Ridge (RID), Lasso (LAS), Decision Tree (DTR), Elastic Net (ENT), Linear Regression (LR), and Support Vector Regression (SVR). The performance for each regressor is calculated in terms of mean absolute error (MAE), root mean-square error (RMSE), R2, Hazard Ratio (HR), and p-value.

3. **Top10 Correlation results:** This file contains the correlation results between number of mutations/gene/sample and overall survival time. Gene name, Correlation coefficents and p-value is reported for top10 genes based on their correlation coefficients. Classification and regression models are developed using mutations/gene/sample values from these top-10 genes.

4. **Correlation results:** This file contains the correlation results between number of mutations/gene/sample and overall survival time for all the genes sorted in order of coefficents.

5. **Mutations per sample per gene file:** This file reports the number of mutations/gene/sample along with overall survival time (OS.time) and overall status (OS).

6. **Best Classification Model:** This is model file which user can use to make the survival group prediction such as High-/low-risk group for unknown samples based on top-10 genes. Model with the highest AUROC will be saved.

7. **Best Regression Model:** This is model file which user can use to make the survival time prediction for unknown samples based on top-10 genes. Model with the highest HR value will be saved.

## Standalone
The Standalone version of this mthod is written in python3 and following libraries are necessary for the successful run:
- scikit-learn
- Pandas
- Numpy
- rpy2
- tqdm

## Important Note
In order to run the provided example, please download the following files:
1. Database file required by annovar to map the coordinates with gene names. [**Click Me**](https://github.com/sumeetpatiyal/muthrp/raw/main/humandb.zip)
2. Example VCF files. [**Click Me**](https://github.com/sumeetpatiyal/muthrp/raw/main/test.zip)
3. Download these files and unzip them.

## Minimum USAGE
To know about the available option for the stanadlone, type the following command:
```
python3 muthrp.py -h
```
To run the example, type the following command:
```
python3 -W ignore muthrp.py -i test/ -t MUTECT2 -f VCF -s gdc_sample_sheet.tsv -c clinical_data.tsv
```
This will provide the five output files as afore-mentioned. It will use other parameters by default. It will save the output in .csv format with files string "mutation_based_results.csv" as the suffix.

## Full Usage
```
usage: muthrp.py [-h]
		-i INPUT
		-t {VARSCAN2,MUTECT2,MUSE,SOMATICSNIPER}
		-f {VCF,MAF}
		-s SAMPLE 
		[-o OUTPUT]
		[-d DATABASE]
		[-c CLINICAL]
```
```
Please provide following arguments for the sucessful run

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, -I INPUT, --input INPUT
                        Input Directory: Please provide the path of directory containing either VCF or MAF files.
  -t {VARSCAN2,MUTECT2,MUSE,SOMATICSNIPER}, -T {VARSCAN2,MUTECT2,MUSE,SOMATICSNIPER}, --tech {VARSCAN2,MUTECT2,MUSE,SOMATICSNIPER}
                        Techniques: Please provide the techniques to be consider from the available options. By default, its MUTECT2
  -f {VCF,MAF}, -F {VCF,MAF}, --format {VCF,MAF}
                        File Type Formats: VCF: Variat Calling Format; MAF: Mutation Annotation Format.
  -s SAMPLE, -S SAMPLE, --sample SAMPLE
                        Sample Data File: Please provide the sample file containing the file IDs of patients to map with the TCGA IDs.
  -o OUTPUT, -O OUTPUT, --output OUTPUT
                        Output: This would be the prefix that will be added in the output filename.
  -d DATABASE, -D DATABASE, --database DATABASE
                        Database: Please provide the path the database required by annovar to map the coordinates with gene names.
  -c CLINICAL, -C CLINICAL, --clinical CLINICAL
                        Clinical Data File: Please provide the file containing the clinical information of patients with OS and OS.time to calculate HR.
```

**Input File:** This argument takes the path of directory containing VCF or MAF files.

**Output File:** This is the string which would be incorporated in the output filenames stored in the .csv format. 

**Technique**: User is allowed to choose one of four techniques provided in the help.

**Sample file**: This file comprises the information of the sample and their case IDs, which is used to trace the clinical data.

**Clinical Data**: This file comprises of clinical information of patients.

**Database**: This folder contains the files to get the gene names by mapping VCF file coordinates on it.

**Format**: User is allowed to choose between two formats of mutations, such as VCF or MAF.

Package Files
=============
It contain following files, brief descript of these files given below

README.md                       		: This file provide information about this package

muthrp.py                       		: Main python program

test                            		: This folder contains the VCF files for test run.

humandb                         		: This folder contain the files for gene mapping required by annovar script.

clinical_data.tsv               		: This file contains the clinical data i.e. OS and OS.time for each patient.

gdc_sample_sheet.tsv            		: This file contains the information of samples.

convert2annovar.pl              		: perl script to convert"genotype calling" format into ANNOVAR format.

annotate_variation.pl                		: perl script for annotate variations.

Mutations_gene_sample_MUTECT2_VCF_test_run.csv	: This is the example output which reports mutations/gene/sample.

Correlation_MUTECT2_VCF_test_run.csv 		: This is the example output reporting correlation between mutations/gene and OS_time.

Top10_Correlated_genes_MUTECT2_VCF_test_run.csv : This is the example output exhibits top-10 genes used to train classification and regression models.

Classification_MUTECT2_VCF_test_run.csv		: This is example output for classification results.

Regression_MUTECT2_VCF_test_run.csv		: This is example output for performanc of regression models.

LR_Mutect2_VCF_Classification.pkl		: This is the pickle file for best classification model trained on top10 genes.

ENT_Mutect2_VCF_Regression.pkl			: This is the pickle file for best regression model trained on top10 genes.
