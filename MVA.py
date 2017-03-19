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

#from keras.preprocessing import sequence
#from keras.optimizers import SGD, RMSprop, Adagrad
#from keras.utils import np_utils
#from keras.models import Sequential
#from keras.layers.core import Dense, Dropout, Activation, Merge, Flatten, TimeDistributedDense, Masking, Lambda
#from keras.layers.embeddings import Embedding
#from keras.layers.recurrent import LSTM, GRU
#from keras.models import model_from_json

from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, TGraph, TCanvas, TMultiGraph, TLegend
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
from array import array

trainFraction = 0.6
batch_size = 128

MVAScoreHists={}
MVAPRED = {}

def NomrmalizeInput(inputarray):
    for i in range(inputarray.shape[1]):
        inputarray[:,i] = inputarray[:,i]-min(inputarray[:,i])
        inputarray[:,i] = inputarray[:,i] / np.std(inputarray[:,i])
    

def makeData(normalize_input=True, multi_class=False, data_2tag_data=None, data_0tag_data=None, VBF_data=None, ggF_data=None, Zbb_QCD_data=None, Zbb_EWK_data=None, sys_data = None ):   #Zbb_QCD_data_0tag=None, Zbb_EWK_data_0tag=None, 

    data_2tag_data_weight = data_2tag_data[0,:]/float( sum( data_2tag_data[0,:]) ) ## normalize data weight to 1/N Data events
    print "data weight shape", data_2tag_data_weight.shape
    data_2tag_data = np.vstack( (data_2tag_data, np.zeros(data_2tag_data.shape[1])) )
    data_2tag_data_full = data_2tag_data[1:-1, :].transpose()

    data_0tag_data = np.vstack( (data_0tag_data, np.ones(data_0tag_data.shape[1])*3 ) )
    data_0tag_data_full = data_0tag_data[1:-1, :].transpose()

    ggF_data = np.vstack( (ggF_data, np.ones(ggF_data.shape[1])*2) )
    ggF_data_full = ggF_data[1:-1, :].transpose()

    VBF_data_weight = VBF_data[0,:] /float(sum( VBF_data[0,:]) ) ## normalize VBF weight to 1/ sum of event weights
    VBF_data = np.vstack( (VBF_data, np.ones(VBF_data.shape[1])*1) )
    VBF_data_full = VBF_data[1:-1, :].transpose()

    Zbb_EWK_data = np.vstack( (Zbb_EWK_data, np.ones(Zbb_EWK_data.shape[1])*4) )
    Zbb_EWK_data_full = Zbb_EWK_data[1:-1, :].transpose()

    Zbb_QCD_data = np.vstack( (Zbb_QCD_data, np.ones(Zbb_QCD_data.shape[1])*5) )
    Zbb_QCD_data_full = Zbb_QCD_data[1:-1, :].transpose()

    all_data = np.hstack( (data_2tag_data, VBF_data) )
    all_data_weight = np.hstack( (data_2tag_data_weight, VBF_data_weight) )

    print all_data_weight

    print ("ggF shape", ggF_data.shape)
    print ("VBF shape", VBF_data.shape)
    print ("Zbb_EWK shape", Zbb_EWK_data.shape)
    print ("Zbb_QCD shape", Zbb_QCD_data.shape)
    print ("data 2tag shape", data_2tag_data.shape)
    print ("data 0tag shape", data_0tag_data.shape)
    print ("all_data shape", all_data.shape)
    print ("all_data_weight shape", all_data_weight.shape)

    print all_data_weight

    ## randomize the array
    random_index = np.random.permutation(range(all_data.shape[1]))
    all_data_weight = all_data_weight[random_index]
    all_data = all_data[ :, random_index]

    ##event weight cut
    all_data_weight = all_data_weight[all_data[0,:]!=0 ]
    all_data = all_data[ :, all_data[0,:]!=0 ]
    ggF_data = ggF_data[ :, ggF_data[0,:]!=0 ]
    Zbb_EWK_data = Zbb_EWK_data[ :, Zbb_EWK_data[0,:]!=0 ]
    Zbb_QCD_data = Zbb_QCD_data[ :, Zbb_QCD_data[0,:]!=0 ]
    data_0tag_data = data_0tag_data[ :, data_0tag_data[0,:]!=0]

    y = all_data[-1,:]==1
    X = all_data[1:-1, :]
    X = X.transpose()
    evtweight = all_data[0, :]

    ##normalize inputs
    if normalize_input:
        for i in range(X.shape[1]):
            X[:,i] = X[:,i] - min(X[:,i])
            X[:,i] = X[:,i]/ np.std(X[:,i])

    weights_train, weights_test = np.split( all_data_weight, [ int(trainFraction*X.shape[0]) ] )
    X_train, X_test = np.split( X, [ int(trainFraction*X.shape[0]) ] )
    y_train, y_test = np.split( y, [ int(trainFraction*X.shape[0]) ] )
    #weights_train, weights_test = np.split( evtweight, [ int(trainFraction*X.shape[0]) ] )

    ggF_data = ggF_data[1:-1].transpose()
    Zbb_EWK_data = Zbb_EWK_data[1:-1].transpose()
    Zbb_QCD_data = Zbb_QCD_data[1:-1].transpose()
    data_0tag_data = data_0tag_data[1:-1].transpose()

    ##normalize inputs
    if normalize_input:
        NomrmalizeInput(ggF_data)
        NomrmalizeInput(Zbb_QCD_data)
        NomrmalizeInput(Zbb_EWK_data)
        NomrmalizeInput(data_0tag_data)
        NomrmalizeInput(data_0tag_data_full)
        NomrmalizeInput(data_2tag_data_full)
        NomrmalizeInput(ggF_data_full)
        NomrmalizeInput(VBF_data_full)
        NomrmalizeInput(Zbb_QCD_data_full)
        NomrmalizeInput(Zbb_EWK_data_full)

    dataset = {
        "X": X,
        "X_train": X_train,
        "X_test": X_test,
        "y": y,
        "y_train": y_train,
        "y_test": y_test,
        "weights_train": weights_train,
        "weights_test": weights_test,
        "ggF": ggF_data,
        "Zbb_QCD": Zbb_QCD_data,
        "Zbb_EWK": Zbb_EWK_data,
        "data_0tag": data_0tag_data,
        "ggF_full": ggF_data_full,
        "VBF_full": VBF_data_full,
        "Zbb_QCD_full": Zbb_QCD_data_full,
        "Zbb_EWK_full": Zbb_EWK_data_full,
        "data_0tag_full": data_0tag_data_full,
        "data_2tag_full": data_2tag_data_full
	}

    
    sys_dataset = {}

    ## treating sys
    for key in sys_data:
        thisarray = sys_data[key][1: ,:].transpose()
        #thisarray = thisarray[:, thisarray[0,:]!=0]

        if normalize_input:
            NomrmalizeInput(thisarray)
        sys_dataset[key] = thisarray

    return dataset, sys_dataset


def buildModel(dataset, multi_class = False):

    print ("Building Model ...")
    nb_epoch = 5

    X_train = dataset['X_train']
    y_train = dataset['y_train']

    model = Sequential()
    model.add(Dense(128, input_dim=(X_train.shape[1])))
    model.add(Activation('relu'))
    model.add(Dropout(0.4))

    model.add(Dense(64))
    model.add(Activation('relu'))
    model.add(Dropout(0.4))

    if multi_class:
        model.add(Dense(3))
        model.add(Activation('softmax'))
    else:
        model.add(Dense(1))
        model.add(Activation('sigmoid'))

    if multi_class:
        model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=["accuracy"])
    else:
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=["accuracy"])

    print ("Finish Compilation")
    print ("Train...")

    history = model.fit( X_train , y_train, batch_size=batch_size, nb_epoch=nb_epoch, validation_split=0.25, shuffle = True)

    print ("Finished fitting")
    return (model, history)


def buildBDT(dataset, channel, period, btag, doTest = False):

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
    
    bdt_discrete = AdaBoostClassifier( DecisionTreeClassifier(max_depth=max_depth), n_estimators=n_estimators, algorithm="SAMME.R", learning_rate=0.2)
    bdt_discrete.fit(X_train, y_train, weights_train)
    
    filename = "BDT"+channel+period+btag+'.sav'

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
        histup = -1
        histdn = 1
        
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
