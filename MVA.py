import numpy as np
import root_numpy as rnp
import sys
from copy import deepcopy
import os
import scipy.stats as stats
import sklearn.metrics as skmetric
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier

from keras.preprocessing import sequence
from keras.optimizers import SGD, RMSprop, Adagrad
from keras.utils import np_utils
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Merge, Flatten, TimeDistributedDense, Masking, Lambda
from keras.layers.embeddings import Embedding
from keras.layers.recurrent import LSTM, GRU
from keras.models import model_from_json

from ROOT import gDirectory, TLorentzVector, TVector3, TH1D, TH2D, TFile, TGraph, TCanvas, TMultiGraph, TLegend
from HistoLab import HistoTool, CopyHist, StackHists, GetPlot, ScaleHist, Project2D, CopyHists, project1D, NormHists, AddHistPairs, ResetBins, GetHistMedian, GetHistIQR, GetRatios, GetAveRatios, GetGausMean, CalSysVariation, CalSysDiff
from array import array


trainFraction = 0.5
batch_size = 128


MVAScoreHists={}
MVAPRED = {}

def NomrmalizeInput(inputarray):
    for i in range(inputarray.shape[1]):
        inputarray[:,i] = inputarray[:,i]-min(inputarray[:,i])
        inputarray[:,i] = inputarray[:,i] / np.std(inputarray[:,i])
    

def makeData(normalize_input=True, multi_class=False, data_2tag_data=None, data_0tag_data=None, VBF_data=None, ggF_data=None, Zbb_QCD_data=None, Zbb_EWK_data=None, Zbb_QCD_data_0tag=None, Zbb_EWK_data_0tag=None, sys_data = None ):


    data_2tag_data = np.vstack( (data_2tag_data, np.zeros(data_2tag_data.shape[1])) )
    data_2tag_data_full = data_2tag_data[1:-1, :].transpose()

    data_0tag_data = np.vstack( (data_0tag_data, np.ones(data_0tag_data.shape[1])*3 ) )
    data_0tag_data_full = data_0tag_data[1:-1, :].transpose()

    ggF_data = np.vstack( (ggF_data, np.ones(ggF_data.shape[1])*2) )
    ggF_data_full = ggF_data[1:-1, :].transpose()

    VBF_data = np.vstack( (VBF_data, np.ones(VBF_data.shape[1])*1) )
    VBF_data_full = VBF_data[1:-1, :].transpose()

    Zbb_EWK_data = np.vstack( (Zbb_EWK_data, np.ones(Zbb_EWK_data.shape[1])*4) )
    Zbb_EWK_data_full = Zbb_EWK_data[1:-1, :].transpose()

    Zbb_QCD_data = np.vstack( (Zbb_QCD_data, np.ones(Zbb_QCD_data.shape[1])*5) )
    Zbb_QCD_data_full = Zbb_QCD_data[1:-1, :].transpose()

    Zbb_EWK_data_0tag = np.vstack( (Zbb_EWK_data_0tag, np.ones(Zbb_EWK_data_0tag.shape[1])*4) )
    Zbb_EWK_data_0tag_full = Zbb_EWK_data_0tag[1:-1, :].transpose()

    Zbb_QCD_data_0tag = np.vstack( (Zbb_QCD_data_0tag, np.ones(Zbb_QCD_data_0tag.shape[1])*5) )
    Zbb_QCD_data_0tag_full = Zbb_QCD_data_0tag[1:-1, :].transpose()

    all_data = np.hstack( (data_0tag_data, VBF_data) )
    all_data = np.hstack( (all_data, ggF_data) )
    all_data = np.hstack( (all_data, data_2tag_data) )

    print ("ggF shape", ggF_data.shape)
    print ("VBF shape", VBF_data.shape)
    print ("Zbb_EWK shape", Zbb_EWK_data.shape)
    print ("Zbb_QCD shape", Zbb_QCD_data.shape)
    print ("data 2tag shape", data_2tag_data.shape)
    print ("data 0tag shape", data_0tag_data.shape)
    print ("all_data shape", all_data.shape)

    ## randomize the array
    random_index = np.random.permutation(range(all_data.shape[1]))
    all_data = all_data[ :, random_index]

    ##train on 2tag 
    all_data = all_data[ :, all_data[-1,:]!=3 ]
    all_data = all_data[ :, all_data[-1,:]!=4 ]
    all_data = all_data[ :, all_data[-1,:]!=5 ]
    ##train on VBF only
    if not multi_class:
        all_data = all_data[ :, all_data[-1,:]!=2 ]

    ##event weight cut
    all_data = all_data[ :, all_data[0,:]!=0 ]
    ggF_data = ggF_data[ :, ggF_data[0,:]!=0 ]
    Zbb_EWK_data = Zbb_EWK_data[ :, Zbb_EWK_data[0,:]!=0 ]
    Zbb_QCD_data = Zbb_QCD_data[ :, Zbb_QCD_data[0,:]!=0 ]
    Zbb_EWK_data_0tag = Zbb_EWK_data_0tag[ :, Zbb_EWK_data_0tag[0,:]!=0 ]
    Zbb_QCD_data_0tag = Zbb_QCD_data_0tag[ :, Zbb_QCD_data_0tag[0,:]!=0 ]
    data_0tag_data = data_0tag_data[ :, data_0tag_data[0,:]!=0]

    y = all_data[-1,:]==1
    X = all_data[1:-1, :]
    X = X.transpose()
    evtweight = all_data[0, :]

    if multi_class:
        y = np.ndarray(shape =(all_data.shape[1],3), dtype=float)
        y[:, 0] = (all_data[-1, :]==0).transpose()
        y[:, 1] = (all_data[-1, :]==1).transpose()
        y[:, 2] = (all_data[-1, :]==2).transpose()

    ##normalize inputs
    if normalize_input:
        for i in range(X.shape[1]):
            X[:,i] = X[:,i] - min(X[:,i])
            X[:,i] = X[:,i]/ np.std(X[:,i])

    X_train, X_test = np.split( X, [ int(trainFraction*X.shape[0]) ] )
    y_train, y_test = np.split( y, [ int(trainFraction*X.shape[0]) ] )
    weights_train, weights_test = np.split( evtweight, [ int(trainFraction*X.shape[0]) ] )

    ggF_data = ggF_data[1:-1].transpose()
    Zbb_EWK_data = Zbb_EWK_data[1:-1].transpose()
    Zbb_QCD_data = Zbb_QCD_data[1:-1].transpose()
    Zbb_EWK_data_0tag = Zbb_EWK_data_0tag[1:-1].transpose()
    Zbb_QCD_data_0tag = Zbb_QCD_data_0tag[1:-1].transpose()
    data_0tag_data = data_0tag_data[1:-1].transpose()

    ##normalize inputs
    if normalize_input:
        NomrmalizeInput(ggF_data)
        NomrmalizeInput(Zbb_QCD_data)
        NomrmalizeInput(Zbb_EWK_data)
        NomrmalizeInput(Zbb_QCD_data_0tag)
        NomrmalizeInput(Zbb_EWK_data_0tag)
        NomrmalizeInput(data_0tag_data)
        NomrmalizeInput(data_0tag_data_full)
        NomrmalizeInput(data_2tag_data_full)
        NomrmalizeInput(ggF_data_full)
        NomrmalizeInput(VBF_data_full)
        NomrmalizeInput(Zbb_QCD_data_full)
        NomrmalizeInput(Zbb_EWK_data_full)
        NomrmalizeInput(Zbb_QCD_data_0tag_full)
        NomrmalizeInput(Zbb_EWK_data_0tag_full)


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
        "Zbb_QCD_0tag": Zbb_QCD_data_0tag,
        "Zbb_EWK_0tag": Zbb_EWK_data_0tag,
        "data_0tag": data_0tag_data,
        "ggF_full": ggF_data_full,
        "VBF_full": VBF_data_full,
        "Zbb_QCD_full": Zbb_QCD_data_full,
        "Zbb_EWK_full": Zbb_EWK_data_full,
        "Zbb_QCD_0tag_full": Zbb_QCD_data_0tag_full,
        "Zbb_EWK_0tag_full": Zbb_EWK_data_0tag_full,
        "data_0tag_full": data_0tag_data_full,
        "data_2tag_full": data_2tag_data_full
	}

    
    sys_dataset = {}

    ## treating sys
    for key in sys_data:
        print ("sys ", key)
        thisarray = sys_data[key][1:,:].transpose()
        thisarray = thisarray[:, thisarray[0,:]!=0]
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


def buildBDT(dataset, doTest = False):

    max_depth = 1
    n_estimators = 100
    if doTest:
        max_depth = 1
        n_estimators = 20

    print ("building bdt")
    X_train = dataset['X_train']
    y_train = dataset['y_train']
    
    bdt_discrete = AdaBoostClassifier( DecisionTreeClassifier(max_depth=max_depth), n_estimators=n_estimators,  learning_rate=1,  algorithm="SAMME.R")
    bdt_discrete.fit(X_train, y_train)

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
    Zbb_EWK_0tag    = dataset['Zbb_EWK_0tag']
    Zbb_QCD_0tag    = dataset['Zbb_QCD_0tag']
    data_0tag    = dataset['data_0tag']

    ggF_full = dataset['ggF_full']
    VBF_full = dataset['VBF_full']
    data_0tag_full = dataset['data_0tag_full']
    data_2tag_full = dataset['data_2tag_full']
    Zbb_EWK_full = dataset['Zbb_EWK_full']
    Zbb_QCD_full = dataset['Zbb_QCD_full']
    Zbb_EWK_0tag_full = dataset['Zbb_EWK_0tag_full']
    Zbb_QCD_0tag_full = dataset['Zbb_QCD_0tag_full']
    
    pred_train = None
    pred_test = None
    pred_ggF = None
    pred_ggF_full = None
    pred_Zbb_EWK = None
    pred_Zbb_EWK_0tag = None
    pred_Zbb_EWK_full = None
    pred_Zbb_EWK_0tag_full = None
    pred_Zbb_QCD = None
    pred_Zbb_QCD_0tag = None
    pred_Zbb_QCD_full = None
    pred_Zbb_QCD_0tag_full = None
    pred_VBF_full = None
    pred_0tag_full = None
    pred_2tag_full = None

    histup = 1
    histdn = 0

    if "NN" in modelname:

        pred_train = model.predict( X_train )
        pred_test = model.predict( X_test )
        pred_ggF = model.predict( ggF )
        pred_Zbb_EWK = model.predict( Zbb_EWK )
        pred_Zbb_QCD = model.predict( Zbb_QCD )
        pred_Zbb_EWK_0tag = model.predict( Zbb_EWK_0tag )
        pred_Zbb_QCD_0tag = model.predict( Zbb_QCD_0tag )

        pred_data_0tag = model.predict( data_0tag )
        
        pred_VBF_full = model.predict( VBF_full )
        pred_ggF_full = model.predict( ggF_full )
        pred_Zbb_EWK_0tag_full = model.predict( Zbb_EWK_0tag_full )
        pred_Zbb_QCD_0tag_full = model.predict( Zbb_QCD_0tag_full )
        pred_Zbb_EWK_full = model.predict( Zbb_EWK_full )
        pred_Zbb_QCD_full = model.predict( Zbb_QCD_full )
        pred_0tag_full = model.predict( data_0tag_full )
        pred_2tag_full = model.predict( data_2tag_full )

        if True:
            pred_train =  pred_train[:,0]
            pred_test =   pred_test[:,0]
            pred_ggF =    pred_ggF[:,0]
            pred_Zbb_EWK =    pred_Zbb_EWK[:,0]
            pred_Zbb_QCD =    pred_Zbb_QCD[:,0]
            pred_Zbb_EWK_0tag =    pred_Zbb_EWK_0tag[:,0]
            pred_Zbb_QCD_0tag =    pred_Zbb_QCD_0tag[:,0]
            pred_data_0tag = pred_data_0tag[:,0]
            
            pred_VBF_full = pred_VBF_full[:,0]
            pred_ggF_full = pred_ggF_full[:,0]
            pred_Zbb_EWK_full =  pred_Zbb_EWK_full[:,0]
            pred_Zbb_QCD_full =  pred_Zbb_QCD_full[:,0]
            pred_Zbb_EWK_0tag_full =  pred_Zbb_EWK_0tag_full[:,0]
            pred_Zbb_QCD_0tag_full =  pred_Zbb_QCD_0tag_full[:,0]
            pred_0tag_full = pred_0tag_full[:,0]
            pred_2tag_full = pred_2tag_full[:,0]
            
            
#            for key in sysdataset:
#                tmp_pred = model.predict(sysdataset[key])
#                MVAPRED[modelname+"_"+key+"_full"] = tmp_pred[:,0]

                    
    elif modelname == "BDT":
        histup = 1.1
        histdn = 0.85
        
        pred_train = model.predict_log_proba( X_train )
        pred_train = pred_train[:,0]/(pred_train[:,1])
        
        pred_test = model.predict_log_proba( X_test )
        pred_test = pred_test[:,0]/(pred_test[:,1])
        
        pred_ggF = model.predict_log_proba( ggF )
        pred_ggF = pred_ggF[:,0]/(pred_ggF[:,1])

        print "Zbb contribution", len(Zbb_EWK)

        pred_Zbb_EWK = model.predict_log_proba( Zbb_EWK )
        pred_Zbb_EWK = pred_Zbb_EWK[:,0]/(pred_Zbb_EWK[:,1])

        pred_Zbb_QCD = model.predict_log_proba( Zbb_QCD )
        pred_Zbb_QCD = pred_Zbb_QCD[:,0]/(pred_Zbb_QCD[:,1])

        pred_Zbb_EWK_0tag = model.predict_log_proba( Zbb_EWK_0tag )
        pred_Zbb_EWK_0tag = pred_Zbb_EWK_0tag[:,0]/(pred_Zbb_EWK_0tag[:,1])

        pred_Zbb_QCD_0tag = model.predict_log_proba( Zbb_QCD_0tag )
        pred_Zbb_QCD_0tag = pred_Zbb_QCD_0tag[:,0]/(pred_Zbb_QCD_0tag[:,1])
        
        pred_data_0tag = model.predict_log_proba( data_0tag )
        pred_data_0tag = pred_data_0tag[:,0]/(pred_data_0tag[:,1])

        pred_VBF_full = model.predict_log_proba( VBF_full )
        pred_VBF_full = pred_VBF_full[:,0]/pred_VBF_full[:,1]

        pred_ggF_full = model.predict_log_proba( ggF_full )
        pred_ggF_full = pred_ggF_full[:,0]/pred_ggF_full[:,1]

        pred_Zbb_EWK_full = model.predict_log_proba( Zbb_EWK_full )
        pred_Zbb_EWK_full = pred_Zbb_EWK_full[:,0]/(pred_Zbb_EWK_full[:,1])

        pred_Zbb_QCD_full = model.predict_log_proba( Zbb_QCD_full )
        pred_Zbb_QCD_full = pred_Zbb_QCD_full[:,0]/(pred_Zbb_QCD_full[:,1])

        pred_Zbb_EWK_0tag_full = model.predict_log_proba( Zbb_EWK_0tag_full )
        pred_Zbb_EWK_0tag_full = pred_Zbb_EWK_0tag_full[:,0]/(pred_Zbb_EWK_0tag_full[:,1])

        pred_Zbb_QCD_0tag_full = model.predict_log_proba( Zbb_QCD_0tag_full )
        pred_Zbb_QCD_0tag_full = pred_Zbb_QCD_0tag_full[:,0]/(pred_Zbb_QCD_0tag_full[:,1])

        pred_0tag_full = model.predict_log_proba( data_0tag_full )
        pred_0tag_full = pred_0tag_full[:,0]/pred_0tag_full[:,1]

        pred_2tag_full = model.predict_log_proba( data_2tag_full )
        pred_2tag_full = pred_2tag_full[:,0]/pred_2tag_full[:,1]

        for key in sysdataset:
            tmp_pred = model.predict_log_proba(sysdataset[key])
            MVAPRED[modelname+"_"+key+"_full"] = tmp_pred[:,0]/tmp_pred[:,1]
            
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

    MVAPRED[modelname+"_Zbb_QCD_0tag_full"] = pred_Zbb_QCD_0tag_full
    MVAPRED[modelname+"_Zbb_EWK_0tag_full"] = pred_Zbb_EWK_0tag_full

    MVAPRED[modelname+"_VBF_ggF_full"] = np.hstack((pred_VBF_full,pred_ggF_full))
    MVAPRED[modelname+"_data_0tag_full"] = pred_0tag_full
    MVAPRED[modelname+"_data_2tag_full"] = pred_2tag_full
    MVAPRED[modelname+"_ggF"] = pred_ggF
    MVAPRED[modelname+"_Zbb_QCD"] = pred_Zbb_QCD
    MVAPRED[modelname+"_Zbb_EWK"] = pred_Zbb_EWK
    MVAPRED[modelname+"_Zbb_QCD_0tag"] = pred_Zbb_QCD_0tag
    MVAPRED[modelname+"_Zbb_EWK_0tag"] = pred_Zbb_EWK_0tag

    
    hist_array_sig_train = np.histogram(pred_sig_train, 100, (histdn, histup))[0]
    hist_array_bkg_train = np.histogram(pred_bkg_train, 100, (histdn, histup))[0]
    hist_array_sig_test = np.histogram(pred_sig_test, 100, (histdn, histup))[0]
    hist_array_bkg_test = np.histogram(pred_bkg_test, 100, (histdn, histup))[0]
    hist_array_ggF      = np.histogram(pred_ggF, 100, (histdn, histup))[0]
    hist_array_Zbb_QCD  = np.histogram(pred_Zbb_QCD, 100, (histdn, histup))[0]
    hist_array_Zbb_EWK  = np.histogram(pred_Zbb_EWK, 100, (histdn, histup))[0]
    hist_array_data_0tag= np.histogram(pred_data_0tag, 100, (histdn, histup))[0]
    
    Hist_sig_train = TH1D(modelname+"_Hist_sig_train", modelname+"_Hist_sig_train", 100, histdn, histup)
    Hist_sig_test  = TH1D(modelname+"_Hist_sig_test", modelname+"_Hist_sig_test", 100, histdn, histup)
    Hist_bkg_train = TH1D(modelname+"_Hist_bkg_train", modelname+"_Hist_bkg_train", 100, histdn, histup)
    Hist_bkg_test  = TH1D(modelname+"_Hist_bkg_test", modelname+"_Hist_bkg_test", 100, histdn, histup)
    Hist_ggF       = TH1D(modelname+"_Hist_ggF", modelname+"_Hist_ggF", 100, histdn, histup)
    Hist_Zbb_QCD   = TH1D(modelname+"_Hist_Zbb_QCD", modelname+"_Hist_Zbb_QCD", 100, histdn, histup)
    Hist_Zbb_EWK   = TH1D(modelname+"_Hist_Zbb_EWK", modelname+"_Hist_Zbb_EWK", 100, histdn, histup)
    Hist_data_0tag = TH1D(modelname+"_Hist_data_0tag", modelname+"_Hist_data_0tag", 100, histdn, histup)
    
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
    MVAScoreHists[modelname+"_ggF"]       =  Hist_ggF
    MVAScoreHists[modelname+"_Zbb_EWK"]   =  Hist_Zbb_EWK
    MVAScoreHists[modelname+"_Zbb_QCD"]   =  Hist_Zbb_QCD
    MVAScoreHists[modelname+"_data_0tag"] =  Hist_data_0tag
