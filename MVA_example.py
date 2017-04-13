import numpy as np
import root_numpy as rnp
import sys
from copy import deepcopy
import os
import scipy.stats as stats
import pickle
import sklearn.metrics as skmetric
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score

from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, TGraph, TCanvas, TMultiGraph, TLegend
from array import array

trainFraction = 0.6
batch_size = 128

MVAScoreHists={}
MVAPRED = {}

def NomrmalizeInput(inputarray):
    for i in range(inputarray.shape[1]):
        inputarray[:,i] = inputarray[:,i]-min(inputarray[:,i])
        inputarray[:,i] = inputarray[:,i] / np.std(inputarray[:,i])
    

def buildBDT(dataset, channel, period, btag, MVALeaveOut, doTest = False, splitIndex = None):

    max_depth = 1
    n_estimators = 300
    if channel == "2cen":
        n_estimators = 200
    if channel == "4cen":
        n_estimators = 300

    if doTest:
        max_depth = 1
        n_estimators = 20

    print ("building bdt")
    X_train = dataset['X_train']
    y_train = dataset['y_train']
    weights_train = dataset['weights_train']

    if splitIndex != None:

        X_train = dataset['X_train_splits'][splitIndex]
        y_train = dataset['y_train_splits'][splitIndex]
        weights_train = dataset['weights_train_splits'][splitIndex]

        print "splitting"
        print "X_train", type(X_train), type(dataset['X_train_splits'][splitIndex])
        print "y_train", type(y_train), type(dataset['y_train_splits']), type(dataset['y_train_splits'][splitIndex])
        print "weights_train", type(weights_train), type(dataset['weights_train_splits']), type(dataset['weights_train_splits'][splitIndex])

    print "no split"
    print "X_train", type(X_train)
    print "y_train", type(y_train)
    print "weights_train", type(weights_train)
            
    bdt_discrete = AdaBoostClassifier( DecisionTreeClassifier(max_depth=max_depth), n_estimators=n_estimators, algorithm="SAMME.R", learning_rate=0.2)
    bdt_discrete.fit(X_train, y_train, weights_train)
    
    filename = "BDT"+channel+period+btag+'.sav'
    if MVALeaveOut != "No":
        filename = "BDT"+channel+period+btag+"_No"+MVALeaveOut+'.sav'

    if splitIndex == None:
        pickle.dump(bdt_discrete, open(filename, 'wb'))

    print ("finished building bdt")
    return bdt_discrete


def PredictModel(model, modelname, dataset, sysdataset, multi_class = False):
    print ("predicting ", modelname)
        
    X_test = dataset['X_test']
    y_test = dataset['y_test']
    X_train = dataset['X_train']
    y_train = dataset['y_train']
    ggF    = dataset['ggF']
    Zbb_EWK    = dataset['Zbb_EWK']
    Zbb_QCD    = dataset['Zbb_QCD']
    data_0tag    = dataset['data_0tag']

    ggF_full = dataset['ggF_full']
    VBF_full = dataset['VBF_full']
    data_0tag_full = dataset['data_0tag_full']
    data_2tag_full = dataset['data_2tag_full']
    Zbb_EWK_full = dataset['Zbb_EWK_full']
    Zbb_QCD_full = dataset['Zbb_QCD_full']
    
    pred_train = None
    pred_test = None
    pred_ggF = None
    pred_ggF_full = None
    pred_Zbb_EWK = None
    pred_Zbb_EWK_full = None
    pred_Zbb_QCD = None
    pred_Zbb_QCD_full = None
    pred_VBF_full = None
    pred_0tag_full = None
    pred_2tag_full = None

    histup = 1
    histdn = 0

    if modelname == "BDT":
        
        pred_train = model.decision_function( X_train )
        pred_test = model.decision_function( X_test )
        print "finished predicting X"
        pred_ggF = model.decision_function( ggF )
        pred_Zbb_EWK = model.decision_function( Zbb_EWK )
        pred_Zbb_QCD = model.decision_function( Zbb_QCD )
        pred_data_0tag = model.decision_function( data_0tag )
        print "finished predicting 0tag"
        pred_VBF_full = model.decision_function( VBF_full )
        pred_ggF_full = model.decision_function( ggF_full )
        pred_Zbb_EWK_full = model.decision_function( Zbb_EWK_full )
        pred_Zbb_QCD_full = model.decision_function( Zbb_QCD_full )
        print "finished predicting MC"
        pred_0tag_full = model.decision_function( data_0tag_full )
        pred_2tag_full = model.decision_function( data_2tag_full )
        print "finished predicting data"

        histup = 0.1
        histdn = -0.2

        for key in sysdataset:
            #tmp_pred = model.predict_log_proba(sysdataset[key])
            #MVAPRED[modelname+"_"+key+"_full"] = tmp_pred[:,0]/tmp_pred[:,1]
            #print key, sysdataset[key].shape

            tmp_pred = model.decision_function(sysdataset[key])
            MVAPRED[modelname+"_"+key+"_full"] = tmp_pred
            
        print ('BDT importance', model.feature_importances_)
        
    if not multi_class:
        pred_sig_train =  pred_train[ y_train==True ]
        pred_bkg_train =  pred_train[ y_train==False ]
        pred_sig_test =  pred_test[ y_test==True ]
        pred_bkg_test =  pred_test[ y_test==False ]
    else:
        pred_sig_train =  pred_train[ y_train[:,0]==0 ]
        pred_bkg_train =  pred_train[ y_train[:,0]==1 ]
        pred_sig_test =  pred_test[ y_test[:,0]==0 ]
        pred_bkg_test =  pred_test[ y_test[:,0]==1 ]
        
    MVAPRED[modelname+"_sig_train"] = pred_sig_train
    MVAPRED[modelname+"_bkg_train"] = pred_bkg_train
    MVAPRED[modelname+"_sig_test"] = pred_sig_test
    MVAPRED[modelname+"_bkg_test"] = pred_bkg_test

    MVAPRED[modelname+"_ggF_full"] = pred_ggF_full
    MVAPRED[modelname+"_VBF_full"] = pred_VBF_full
    MVAPRED[modelname+"_Zbb_QCD_full"] = pred_Zbb_QCD_full
    MVAPRED[modelname+"_Zbb_EWK_full"] = pred_Zbb_EWK_full

    MVAPRED[modelname+"_VBF_ggF_full"] = np.hstack((pred_VBF_full,pred_ggF_full))
    MVAPRED[modelname+"_data_0tag_full"] = pred_0tag_full
    MVAPRED[modelname+"_data_2tag_full"] = pred_2tag_full
    MVAPRED[modelname+"_ggF"] = pred_ggF
    MVAPRED[modelname+"_Zbb_QCD"] = pred_Zbb_QCD
    MVAPRED[modelname+"_Zbb_EWK"] = pred_Zbb_EWK

    print "data training ", pred_bkg_train.shape[0]
    print "VBF training ",  pred_sig_train.shape[0]
    
    hist_array_sig_train = np.histogram(pred_sig_train, 30, (histdn, histup))[0]
    hist_array_bkg_train = np.histogram(pred_bkg_train, 30, (histdn, histup))[0]
    hist_array_sig_test = np.histogram(pred_sig_test, 30, (histdn, histup))[0]
    hist_array_bkg_test = np.histogram(pred_bkg_test, 30, (histdn, histup))[0]
    hist_array_ggF      = np.histogram(pred_ggF, 30, (histdn, histup))[0]
    hist_array_Zbb_QCD  = np.histogram(pred_Zbb_QCD, 30, (histdn, histup))[0]
    hist_array_Zbb_EWK  = np.histogram(pred_Zbb_EWK, 30, (histdn, histup))[0]
    hist_array_data_0tag= np.histogram(pred_data_0tag, 30, (histdn, histup))[0]
    
    Hist_sig_train = TH1D(modelname+"_Hist_sig_train", modelname+"_Hist_sig_train", 30, histdn, histup)
    Hist_sig_test  = TH1D(modelname+"_Hist_sig_test", modelname+"_Hist_sig_test", 30, histdn, histup)
    Hist_bkg_train = TH1D(modelname+"_Hist_bkg_train", modelname+"_Hist_bkg_train", 30, histdn, histup)
    Hist_bkg_test  = TH1D(modelname+"_Hist_bkg_test", modelname+"_Hist_bkg_test", 30, histdn, histup)
    Hist_ggF       = TH1D(modelname+"_Hist_ggF", modelname+"_Hist_ggF", 30, histdn, histup)
    Hist_Zbb_QCD   = TH1D(modelname+"_Hist_Zbb_QCD", modelname+"_Hist_Zbb_QCD", 30, histdn, histup)
    Hist_Zbb_EWK   = TH1D(modelname+"_Hist_Zbb_EWK", modelname+"_Hist_Zbb_EWK", 30, histdn, histup)
    Hist_data_0tag = TH1D(modelname+"_Hist_data_0tag", modelname+"_Hist_data_0tag", 30, histdn, histup)
    
    print ("hist", Hist_sig_train)
    
    rnp.array2hist(hist_array_sig_test, Hist_sig_test)
    rnp.array2hist(hist_array_sig_train, Hist_sig_train)
    rnp.array2hist(hist_array_bkg_test, Hist_bkg_test)
    rnp.array2hist(hist_array_bkg_train, Hist_bkg_train)
    rnp.array2hist(hist_array_ggF, Hist_ggF)
    rnp.array2hist(hist_array_Zbb_EWK, Hist_Zbb_EWK)
    rnp.array2hist(hist_array_Zbb_QCD, Hist_Zbb_QCD)
    rnp.array2hist(hist_array_data_0tag, Hist_data_0tag)
    
    MVAScoreHists[modelname+"_sig_train"] =  Hist_sig_train
    MVAScoreHists[modelname+"_sig_test"]  =  Hist_sig_test
    MVAScoreHists[modelname+"_bkg_train"] =  Hist_bkg_train
    MVAScoreHists[modelname+"_bkg_test"]  =  Hist_bkg_test
    
    Hist_vbf = deepcopy(Hist_sig_train)
    Hist_vbf.Add(Hist_sig_test)
    MVAScoreHists[modelname+"_VBF"]       =  Hist_vbf

    Hist_bkg = deepcopy(Hist_bkg_train)
    Hist_bkg.Add(Hist_bkg_test)
    MVAScoreHists[modelname+"_data_2tag"]       =  Hist_bkg

    MVAScoreHists[modelname+"_ggF"]       =  Hist_ggF
    MVAScoreHists[modelname+"_Zbb_EWK"]   =  Hist_Zbb_EWK
    MVAScoreHists[modelname+"_Zbb_QCD"]   =  Hist_Zbb_QCD
    MVAScoreHists[modelname+"_data_0tag"] =  Hist_data_0tag
