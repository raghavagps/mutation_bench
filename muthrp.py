import glob
import os
import re
import sys
import argparse
import warnings
import zipfile
import pickle
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import pearsonr
from math import *
from sklearn.metrics import roc_auc_score
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import f1_score
from sklearn.tree import DecisionTreeClassifier as DT
from sklearn.ensemble import RandomForestClassifier as RF
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier as KN
from sklearn.model_selection import KFold
from sklearn.naive_bayes import GaussianNB as GNB
from sklearn.linear_model import LogisticRegression as LR
from xgboost import XGBClassifier as XGB
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor as RFR
from sklearn.tree import DecisionTreeRegressor as DTR
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from argparse import RawTextHelpFormatter
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
pandas2ri.activate()
os.environ['KMP_DUPLICATE_LIB_OK']='True'
warnings.filterwarnings("ignore", category=FutureWarning)
parser = argparse.ArgumentParser(description='Please provide following arguments for the sucessful run',formatter_class=RawTextHelpFormatter)
parser.add_argument("-i","-I","--input", type=str, required=True, help="Input Directory: Please provide the path of directory containing either VCF or MAF files.")
parser.add_argument("-t","-T", "--tech",type=str.upper, required=True, choices = ['VARSCAN2','MUTECT2','MUSE','SOMATICSNIPER'],help="Techniques: Please provide the techniques to be consider from the available options. By default, its MUTECT2")
parser.add_argument("-f","-F", "--format",type=str.upper, required=True, choices = ['VCF','MAF'],help="File Type Formats: VCF: Variat Calling Format; MAF: Mutation Annotation Format.")
parser.add_argument("-s","-S", "--sample",type=str, required=True ,help="Sample Data File: Please provide the sample file containing the file IDs of patients to map with the TCGA IDs.")
parser.add_argument("-o","-O", "--output",type=str, help="Output: This would be the prefix that will be added in the output filename.")
parser.add_argument("-d","-D", "--database",type=str ,help="Database: Please provide the path the database required by annovar to map the coordinates with gene names.")
parser.add_argument("-c","-C", "--clinical",type=str ,help="Clinical Data File: Please provide the file containing the clinical information of patients with OS and OS.time to calculate HR.")
args = parser.parse_args()
def annovar_run(filename1,humdb):
    print("Running Annovar. Please wait.....")
    for j in tqdm(filename1):
        if not glob.glob(j+'/*.new.vcf4'):
            k = j.split('/')[-1]
            os.system('gunzip '+j+'/*.vcf.gz')
            os.system('perl convert2annovar.pl '+j+'/*.vcf -format vcf4old > '+j+'/'+k+'.new.vcf4')
            os.system('perl annotate_variation.pl -geneanno '+j+'/'+k+'.new.vcf4 -buildver hg19 '+humdb)
def Clean_names(gene_name):
    if re.search('\(.*', gene_name):
        pos = re.search('\(.*', gene_name).start()
        return gene_name[:pos]
    else:
        return gene_name
def final_genes(filename1):
    ss = []
    uu = []
    print("Getting whole set of genes. Please wait....")
    for i in tqdm(filename1):
        k = i.split('/')[-1]
        dd = pd.read_csv(i+'/'+k+'.new.vcf4.variant_function',sep="\t",header=None)
        cc = []
        for j in range(0,len(dd)):
            cc.extend(dd[1][j].split(","))
        for k in range(0,len(cc)):
            ss.extend(cc[k].split(";"))
    df2 = pd.DataFrame(ss)
    df2[0] = df2[0].apply(Clean_names)
    df3 = pd.DataFrame(df2[0].unique())
    for i in df3[0]:
        if i.startswith('NM_') == False:
            uu.append(i)
    df4 = pd.DataFrame(uu)
    df4 = df4.loc[df4[0]!='None']
    return df4
def gene2matrix_vcf(filename1):
    print("Generating matrix comprising number of mutations/gene/samples from VCF file. Please wait.....")
    df5 = final_genes(filename1)
    nn = []
    pp = []
    for i in tqdm(filename1):
        kkk = i.split('/')[-1]
        dd = pd.read_csv(i+'/'+kkk+'.new.vcf4.variant_function',sep="\t",header=None)
        cc = []
        ss = []
        ee = []
        for j in range(0,len(dd)):
            cc.extend(dd[1][j].split(","))
        for k in range(0,len(cc)):
            ss.extend(cc[k].split(";"))
        df2 = pd.DataFrame(ss)
        df2[0] = df2[0].apply(Clean_names)
        for l in range(0,len(df5)):
            ee.append(len(df2[0].loc[df2[0]==df5[0][l]]))
        nn.append(ee)
        pp.append(kkk)
    df6 = pd.concat([pd.DataFrame(pp),pd.DataFrame(nn)],axis=1)
    mm = df5[0].tolist()
    mm.insert(0,'File ID')
    df6.columns = mm
    return df6
def gene2matrix_maf(filename1,tech):
    df_3 = pd.DataFrame()
    for j in filename1:
        if not glob.glob(j+'/*.maf'):
            k = j.split('/')[-1]
            os.system('gunzip '+j+'/*.gz')
        if 'varscan' in glob.glob(j+'/*.maf')[0]:
            mct = 'VARSCAN2'
        elif 'mutect' in glob.glob(j+'/*.maf')[0]:
            mct = 'MUTECT2'
        elif 'muse' in glob.glob(j+'/*.maf')[0]:
            mct = 'MUSE'
        elif 'somaticsniper' in glob.glob(j+'/*.maf')[0]:
            mct = 'SOMATICSNIPER'
        with open(glob.glob(j+'/*.maf')[0]) as f:
            fob = f.readlines()
        jc = []
        for i in fob:
            if i.startswith('#') == False:
                jc.append(i)
        cc = []
        for i in jc:
            cc.append(i.split('\t'))
        xx = pd.DataFrame(cc)
        xx.columns = xx.loc[0]
        xx = xx.loc[1:]
        xx.reset_index(drop=True, inplace=True)
        xy = xx[['Hugo_Symbol','Tumor_Sample_Barcode']]
        xy['MCT'] = mct
        df_3 = df_3.append(xy)
    df_3.reset_index(drop=True,inplace=True)
    df2 = df_3.loc[df_3.MCT==tech].reset_index(drop=True)
    df_1 = pd.DataFrame(df2['Hugo_Symbol'].apply(Clean_names))
    df1 = pd.DataFrame(df_1['Hugo_Symbol'].unique())
    genes = list(df1[0])
    cc = []
    ee = []
    for i in tqdm(df2['Tumor_Sample_Barcode'].unique()):
        cc.append([])
        dd = []
        ee.append(i[0:12])
        for j in genes:
            cc[len(cc)-1].append(len(df2.loc[(df2['Hugo_Symbol']==j) & (df2['Tumor_Sample_Barcode']==i)]))
    df3 = pd.concat([pd.DataFrame(ee),pd.DataFrame(cc)],axis=1)
    genes.insert(0,'bcr_patient_barcode')
    df3.columns=genes
    return df3
def mct_anno(filename1):
    print("Dividing files based on mutation calling techniques. Please wait.....")
    tt = []
    ii = []
    for pp in tqdm(filename1):
        k = pp.split('/')[-1]
        if any(fname.endswith('.vep.vcf') for fname in os.listdir(pp)):
            with open(pp+'/'+k+'.vep.vcf') as f:
                lines = f.readlines()
            for i in lines:
                if 'somatic_mutation_calling_workflow' in i:
                    mct = i.split(',')
                    for j in mct:
                        if 'Name=' in j:
                            tt.append(j.split('=')[1].upper())
                            ii.append(k)
        else:
            tt.append('NG')
            ii.append(k)
    df = pd.concat([pd.DataFrame(ii),pd.DataFrame(tt)],axis=1)
    df.columns = ['File ID','MCT']
    df = df.loc[df.MCT!='NG'].reset_index(drop=True)
    return df
def correlation(file1):
    gene = []
    coe = []
    pv = []
    print("Calculating correlation. Please wait.....")
    for i in tqdm(file1.columns[2:]):
        gene.append(i)
        coe.append(pearsonr(file1[i],file1['OS.time'])[0])
        pv.append(pearsonr(file1[i],file1['OS.time'])[1])
    df_cor = pd.concat([pd.DataFrame(gene),pd.DataFrame(coe),pd.DataFrame(pv)],axis=1)
    df_cor.columns = ['Gene','Coefficient','pvalue']
    top10 = df_cor.loc[df_cor.pvalue<0.05].sort_values('Coefficient')[:10].reset_index(drop=True)
    complete = df_cor.loc[df_cor.pvalue<0.05].sort_values('Coefficient').reset_index(drop=True)
    return top10,complete
def performance(y_actual, y_hat, thr):
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    spec = 0
    mcc = 0
    sens = 0
    acc = 0
    i = 0
    #print (len(y_hat),len(y_actual))
    while i < len(y_hat):
        if y_actual[i]== 1 and y_hat[i] >= thr:
            TP += 1
        elif y_actual[i]== 0 and y_hat[i] >= thr:
            FP += 1
        elif y_actual[i] == 0 and y_hat[i] < thr:
            TN += 1
        elif y_actual[i]==1 and y_hat[i] < thr:
            FN += 1
        i += 1
    binder=TP+FN
    nonb=TN+FP
    total=TP+TN+FP+FN
    Pred = list(map(lambda x: 1 if x >= thr else 0, y_hat))
    if binder!=0:
        sens=(TP/binder)*100
    else:
        sens == 0
    if nonb!=0:
        spec=TN/(nonb)*100
    else:
        spec == 0
    acc=((TP+TN)/total)*100
    FPR=100 - spec
    f1 = 2*TP/((2*TP)+FP+FN)
    F1 = f1_score(y_actual, Pred,zero_division=0)
    auc1=roc_auc_score(y_actual, y_hat)
    kappa = cohen_kappa_score(Pred,y_actual)
    if ((TP+FN)*(TP+FP)*(TN+FP)*(TN+FN)) != 0:
        mcc=(TP*TN-FP*FN)/((TN+FN)*(TP+FN)*(TN+FP)*(TP+FP))**0.5
    else:
        mcc=0
    return(sens,spec,acc,auc1,F1,kappa,mcc)
def classif(file1,file2,m_learn):
    load1 = file1
    load2 = load1['label']
    load3 = file2
    X_test_1 = load3.drop(['label'], axis=1)
    y_test_1 = load3['label']
    X = load1.drop(['label'], axis=1) 
    y = load2
    merge = pd.concat([X,y], axis = 1).reset_index(drop = True)
    if m_learn == 'DT' :
        clf = DT(random_state=42,class_weight='balanced')
        clf.fit(X,y.values.ravel())
    if m_learn == 'RF' :
        clf = RF(random_state=42,class_weight='balanced')
        clf.fit(X,y.values.ravel())
    if m_learn == 'ET' :
        clf = ET(random_state=42,class_weight='balanced')
        clf.fit(X,y.values.ravel())
    if m_learn == 'SVC' :
        clf = SVC(random_state=42,probability=True)
        clf.fit(X,y.values.ravel())
    if m_learn == 'KN' :
        clf = KN()
        clf.fit(X,y.values.ravel())
    if m_learn == 'XGB' :
        clf = XGB(random_state=42,class_weight='balanced')
        clf.fit(X,y.values.ravel())
    if m_learn == 'GNB' :
        clf = GNB()
        clf.fit(X, y.values.ravel())
    if m_learn == 'LR' :
        clf = LR(random_state=42,class_weight='balanced')
        clf.fit(X, y.values.ravel())
    kfold = KFold(5,shuffle=True,random_state=42)
    x2 = pd.DataFrame()
    x3 = pd.DataFrame()
    x4 = pd.DataFrame()
    x_pred = []
    x_true = []
    for train_index, test_index in kfold.split(merge):
        x1 = pd.DataFrame(test_index)
        x2 = x2.append(x1, ignore_index = True)
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        clf.fit(X_train, y_train.values.ravel())
        y_true1, y_pred1 = y_test, clf.predict(X_test)
        y_p_score=clf.predict_proba(X_test)
        y_t=y_true1
        y_p=y_pred1
        y_p_s1=y_p_score[:,-1]
        x_pred = pd.DataFrame(y_p_s1)
        x_true = pd.DataFrame(y_t)
        x3 = x3.append(x_pred, ignore_index = True)
        x4 = x4.append(x_true, ignore_index = True)
    x2["pred"] = x3
    x2["true"] = x4
    predictions = pd.DataFrame()
    y_true11, y_pred11 = y_test_1, clf.predict(X_test_1)
    y_p_score1=clf.predict_proba(X_test_1)
    y_t1=y_true11
    y_p1=y_pred11
    y_p_s11=y_p_score1[:,-1]
    predictions['ID'] = X_test_1.index
    predictions['pred'] = list(y_p_s11)
    predictions['true'] = list(y_t1)
    return x2,predictions
def regres(file1,m_learn):
    ro.r('''
              f <- function(train) {
                       set.seed(42)
                       library(survival)
                       library(ranger)
                       surv_obj_HLA <- Surv(time = train$true,event = train$OS)
                       cox = coxph(surv_obj_HLA ~ train$pred < median(train$pred), data = train)
                       cbind(coef(summary(cox))[2],coef(summary(cox))[5])
               }
               ''')
    r_f = ro.globalenv['f']
    load1 = file1
    load2 = load1['OS.time']*0.0329
    X = load1.drop(['OS.time',"OS"], axis=1)
    y = load2
    merge = pd.concat([X,y], axis = 1).reset_index(drop = True)
    if m_learn == 'DTR' :
        clf = DTR(random_state=42)
        clf.fit(X,y.values.ravel())
    if m_learn == 'RFR' :
        clf = RFR(random_state=42)
        clf.fit(X,y.values.ravel())
    if m_learn == 'SVR' :
        clf = SVR()
        clf.fit(X,y.values.ravel())
    if m_learn == 'LAS' :
        clf = Lasso(random_state=42)
        clf.fit(X,y.values.ravel())
    if m_learn == 'RID' :
        clf = Ridge(random_state=42)
        clf.fit(X,y.values.ravel())
    if m_learn == 'LR' :
        clf = LinearRegression()
        clf.fit(X,y.values.ravel())
    if m_learn == 'ENT' :
        clf = ElasticNet(random_state=42)
        clf.fit(X,y.values.ravel())
    kfold = KFold(5,shuffle=True,random_state=42)
    x2 = pd.DataFrame()
    x3 = pd.DataFrame()
    x4 = pd.DataFrame()
    x_pred = []
    x_true = []
    MAE = []
    RMSE = []
    R2 = []
    for train_index, test_index in kfold.split(merge):
        x1 = pd.DataFrame(test_index)
        x2 = x2.append(x1, ignore_index = True)
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        clf.fit(X_train, y_train.values.ravel())
        y_true1, y_pred1 = y_test, clf.predict(X_test)
        y_t=y_true1
        y_p=y_pred1
        x_pred = pd.DataFrame(y_p)
        x_true = pd.DataFrame(y_t)
        x3 = x3.append(x_pred, ignore_index = True)
        x4 = x4.append(x_true, ignore_index = True)
        rmse = np.sqrt(mean_squared_error(y_t,y_p))
        mae = mean_absolute_error(y_t,y_p)
        r2 = r2_score(y_t,y_p)
        MAE.append(mae)
        RMSE.append(rmse)
        R2.append(r2)
    x2["pred"] = x3
    x2["true"] = x4
    x2['OS'] = [load1['OS'][x2[0][p]] for p in range(0,len(x2[0]))]
    with localconverter (ro.default_converter+pandas2ri.converter):
        r_from_pd_df_A = ro.conversion.py2rpy(x2)
    T1_A=(r_f(r_from_pd_df_A))
    cc_N = []
    gg_N = []
    hh_N = []
    cc_N.append(sum(MAE)/len(MAE))
    gg_N.append(sum(RMSE)/len(RMSE))
    hh_N.append(sum(R2)/len(R2))
    df = pd.DataFrame()
    df["MAE"] = cc_N
    df["RMSE"] = gg_N
    df["R2"] = hh_N
    df['HR'] = [T1_A[0][0]]
    df['p-value'] = [T1_A[0][1]]
    df[["MAE","RMSE","R2","HR"]] = df[["MAE","RMSE","R2","HR"]].round(3)
    return df
def full_run(clin,samp,loca,hudb,tech,ftyp):
    filedir = glob.glob(loca+'/*')
    filename = []
    for i in filedir:
        if os.path.isdir(i) == True:
            filename.append(i)
    df_sample = pd.read_csv(samp, sep='\t')
    df_clincal = pd.read_csv(clin, sep="\t")
    df_clincal.dropna(inplace=True)
    if ftyp == 'VCF':
        annovar_run(filename,hudb)
        df_mct = mct_anno(filename)
        final_matrix = gene2matrix_vcf(filename)
        df1 = pd.merge(df_mct,final_matrix,on='File ID')
        df_mut = df1.loc[df1.MCT==tech].reset_index(drop=True)
        df_mut_sample = pd.merge(df_sample[['File ID','Case ID']],df_mut,on='File ID')
        df_mut_sample['bcr_patient_barcode'] = [df_mut_sample['Case ID'][i].split(',')[0] for i in range(len(df_mut_sample))]
        df_mut_sample_clincal = pd.merge(df_clincal,df_mut_sample,on='bcr_patient_barcode')
        df_mut_sample_clincal.iloc[:,:7].to_csv(filename[0].split('/')[0]+'/Reference_VCF_file.csv', index=None)
        df_mut_sample_clincal.drop(columns= ['bcr_patient_barcode','type','File ID','Case ID','MCT'],inplace=True)
    else:
        final_matrix = gene2matrix_maf(filename,tech)
        df_mut_sample_clincal = pd.merge(df_clinical,final_matrix,on='bcr_patient_barcode')
        df_mut_sample_clincal.drop(columns= ['bcr_patient_barcode','type','MCT'],inplace=True)
    df_c,df_complete = correlation(df_mut_sample_clincal)
    colnames = df_c.Gene.tolist()
    colnames.append('OS')
    colnames.append('OS.time')
    df_class = pd.DataFrame()
    df_reg = pd.DataFrame()
    zz = df_mut_sample_clincal.iloc[:,:]
    df_class = df_class.append(zz)
    df_class = zz[colnames]
    df_class['label'] = [1 if df_class['OS.time'][i]<df_class['OS.time'].median() else 0 for i in range(len(df_class))]
    df_class.drop(columns=['OS','OS.time'],inplace=True)
    df_reg = df_mut_sample_clincal[colnames]
    X_train_C, X_test_C, y_train_C, y_test_C = train_test_split(df_class.iloc[:,:-1], df_class.iloc[:,-1], test_size=0.2, random_state=42)
    train_C, test_C = pd.concat([X_train_C,y_train_C],axis=1).reset_index(drop=True), pd.concat([X_test_C,y_test_C],axis=1).reset_index(drop=True)
    ss = []
    nn = []
    df4 = pd.DataFrame()
    print("Generating classification models. Please wait")
    for q in tqdm(['DT','RF','LR','XGB','KN','GNB','SVC']):
        yy = []
        aa = []
        ll = []
        df_1, df_2 = classif(train_C,test_C,q)
        nn.append(q)
        for qq in range(100,-1,-1):
            num = qq/100
            yy.append(performance(df_1['true'],df_1['pred'],num))
            ll.append(performance(df_2['true'],df_2['pred'],num))
            aa.append(num)
        df3 = pd.concat([pd.DataFrame(yy),pd.DataFrame(ll)],axis=1).reset_index(drop=True)
        df3["Threshold"] = aa
        df3 = df3.round(3)
        df3.columns = ["Sens_tr","Spec_tr","Accuracy_tr","AUC_tr","F1_tr","Kappa_tr","MCC_tr","Sens_te","Spec_te","Accuracy_te","AUC_te","F1_te","Kappa_te","MCC_te","Threshold"]
        df3 = df3[["Threshold","Sens_tr","Spec_tr","Accuracy_tr","AUC_tr","F1_tr","Kappa_tr","MCC_tr","Sens_te","Spec_te","Accuracy_te","AUC_te","F1_te","Kappa_te","MCC_te"]]
        df3["min"] = abs(df3['Sens_tr']-df3["Spec_tr"])
        df4 = df4.append(df3.loc[df3["min"] == min(df3['min'])].reset_index(drop=True).loc[0])
    df4["Classifier"] = nn
    df4 = df4[["Classifier","Sens_tr","Spec_tr","Accuracy_tr","AUC_tr","F1_tr","Kappa_tr","MCC_tr","Sens_te","Spec_te","Accuracy_te","AUC_te","F1_te","Kappa_te","MCC_te"]].reset_index(drop=True)
    ss_p = []
    df_4 = pd.DataFrame()
    print("Generating regression models. Please wait")
    for q_r in tqdm(['DTR','RFR','LR','LAS','RID','ENT','SVR']):
        ss_p.append(q_r)
        df_1 = regres(df_reg,q_r)
        df_4 = df_4.append(df_1)
    df_4['Regressor'] = ss_p
    df_4 = df_4[["Regressor","MAE","RMSE","R2","HR","p-value"]].reset_index(drop=True)
    return df4,df_4,df_c,df_mut_sample_clincal,df_complete
def model_save_C(file1,file2,file3,tech,ftype):
    df1 = file1
    df2 = df1.sort_values('AUC_tr',ascending=False).reset_index(drop=True)
    df3 = file2
    xx = df3.Gene.tolist()
    xx.append('label')
    df4 = file3
    df4['label'] = [1 if df4['OS.time'][i]<df4['OS.time'].median() else 0 for i in range(len(df4))]
    df5 = df4[xx]
    if df2['Classifier'][0] == 'DT' :
        clf = DT(random_state=42,class_weight='balanced')
    if df2['Classifier'][0] == 'RF' :
        clf = RF(random_state=42,class_weight='balanced')
    if df2['Classifier'][0] == 'ET' :
        clf = ET(random_state=42,class_weight='balanced')
    if df2['Classifier'][0] == 'SVC' :
        clf = SVC(random_state=42,probability=True)
    if df2['Classifier'][0] == 'KN' :
        clf = KN()
    if df2['Classifier'][0] == 'XGB' :
        clf = XGB(random_state=42,class_weight='balanced')
    if df2['Classifier'][0] == 'GNB' :
        clf = GNB()
    if df2['Classifier'][0] == 'LR' :
        clf = LR(random_state=42,class_weight='balanced')
    clf.fit(df5.iloc[:,:-1],df5.iloc[:,-1])
    pickle.dump(clf,open(df2['Classifier'][0]+'_'+tech+'_'+ftype+'_Classification.pkl','wb'))
def model_save_R(file1,file2,file3,tech,ftype):
    df1 = file1
    df2 = df1.sort_values('HR',ascending=False).reset_index(drop=True)
    df3 = file2
    xx = df3.Gene.tolist()
    xx.append('OS.time')
    df4 = file3
    df5 = df4[xx]
    if df2['Regressor'][0] == 'DTR' :
        clf = DTR(random_state=42)
    if df2['Regressor'][0] == 'RFR' :
        clf = RFR(random_state=42)
    if df2['Regressor'][0] == 'SVR' :
        clf = SVR()
    if df2['Regressor'][0] == 'LAS' :
        clf = Lasso(random_state=42)
    if df2['Regressor'][0] == 'RID' :
        clf = Ridge(random_state=42)
    if df2['Regressor'][0] == 'LR' :
        clf = LinearRegression()
    if df2['Regressor'][0] == 'ENT' :
        clf = ElasticNet(random_state=42)
    clf.fit(df5.iloc[:,:-1],df5.iloc[:,-1])
    pickle.dump(clf,open(df2['Regressor'][0]+'_'+tech+'_'+ftype+'_Regression.pkl','wb'))
################Arguments##############
in_file = args.input
technique = args.tech
formatf = args.format
samplef = args.sample
if args.output == None:
    o_file = 'mutation_based_results'
else:
    o_file = args.output
if args.database == None:
    db_file = 'humandb/'
else:
    db_file = args.database
if args.clinical == None:
    cl_file = 'clincal_data.tsv'
else:
    cl_file = args.clinical
#######################################
if formatf == 'VCF':
    print("====================================================================================================================================")
    print('Summary of Parameters for VCF based files:''\n')
    print('Input directory:',in_file,'; Output filename contains:',o_file,)
    print('Technique Chosen:',technique,'File Format:',formatf,'; Sample file:', samplef,)
    print('Database Path:',db_file,'; Clinical data file:', cl_file,)
    print("=====================================================================================================================================")
elif formatf == 'MAF':
    print("====================================================================================================================================")
    print('Summary of Parameters for MAF based files:''\n')
    print('Input directory:',in_file,'; Output filename contains:',o_file,)
    print('Technique Chosen:',technique,'File Format:',formatf,)
    print('; Sample file:', samplef,'; Clinical data file:', cl_file,)
    print("=====================================================================================================================================")
############################################Complete Run####################################################
if os.path.isdir('humandb') == False:
    with zipfile.ZipFile('./humandb.zip', 'r') as zip_ref:
        zip_ref.extractall('.')
else:
    pass
df1, df2, df3, df4, df5 = full_run(cl_file,samplef,in_file,db_file,technique,formatf)
model_save_C(df1,df3,df4,technique,formatf)
model_save_R(df2,df3,df4,technique,formatf)
df1.to_csv("Classification_"+technique+"_"+formatf+"_"+o_file+'.csv',index=None)
df2.to_csv("Regression_"+technique+"_"+formatf+"_"+o_file+'.csv',index=None)
df3.to_csv("Top10_Correlated_genes_"+technique+"_"+formatf+"_"+o_file+'.csv',index=None)
df5.to_csv("Correlation_"+technique+"_"+formatf+"_"+o_file+'.csv',index=None)
df4.to_csv("Mutations_gene_sample_"+technique+"_"+formatf+"_"+o_file+'.csv',index=None)
