%% June 2015
% Jorge Luis Guevara Diaz
% jorge.jorjasso@gmail.com
% strategy boostraping I (with stratified sampling and replacement) (alto bias, solo existe un conjunto de test)
% -----------------
% 0 generar indices para training&validation and test.
% 1 crear 2*N conjuntos boostrap a partir de los conjuntos training&validation and test.
% 2 Usar training y validation for model selection
% 3 test in an independent test set
% obs: error la media de los 100 errores, acceso a la varianza
% pensar en usar la matrix de confunsion

% strategy boostraping II (with stratified sampling and replacement)(alta varianza, existen muchas versiones del conjunto de test)
% -----------------
% 0 separar cada classe en un conjunto.
% 1, Respetando los porcentajes de cada classe generar N conjuntos de
% training&valitation and test
% 2 Usar training y validation for model selection
% 3 test in an independent test set
% obs: error la media de los 100 errores, acceso a la varianza
% pensar en usar la matrix de confunsion
addpath ./SVM-KM/
addpath ./kernels/

% reading dataset
fid = fopen('long-4.dat');
C = textscan(fid, '[%n %n] [%n %n] [%n %n] [%n %n] {%n %s %*[^\n]','delimiter', ',','headerlines',9);
X1=[C{1,1},C{1,3},C{1,5},C{1,7}];
X2=[C{1,2},C{1,4},C{1,6},C{1,8}];
y1=[C{1,9}];

%preprocessing: eliminate samples belonging to two classes at the same time
[X1N,X2N,y1N]=preprocessing(X1,X2,y1,C{1,10});

%preprocessing: scaling data

%%scale between [-1,1]
X=[X1N;X2N];
X=bsxfun(@rdivide, 2*bsxfun(@minus, X,min(X)), (max(X)-min(X)))-1;
[m,~]=size(X1N);
[p,~]=size(X2N);
X1N=X(1:m,:);
X2N=X(m+1:m+p,:);

% counting statistics of delete data, and % of samples by class
% format: [nroDataOriginal, nroDataAfterProcessing, %ClassOne ]
statisticsPerClass=[length(y1) length(y1N) sum(y1N==1)*100/length(y1N)];

%random sampling, without replacement using  stratified sampling
nBoostrap=1000;
kernelList=[1:25];
statistics = cell(nBoostrap,25,3);% nroBoostrap x nroKernels x {CVstat, trainingStat, testStat}
for nB=1:nBoostrap
    CVP=cvpartition(y1N,'holdout',0.25); % 25% for test set
    trIdx = CVP.training; % training indices
    teIdx = CVP.test; % test indices
    
    % data set for training and validation TrVd = training&validation
    %---------------------------------------------------------------
    X1TrVd=X1N(trIdx,:);  X2TrVd=X2N(trIdx,:); y1TrVd=y1N(trIdx,:);
    
    % Divide ...TrVd sets in training and validation sets
    kFold=7;
    CVP=cvpartition(y1TrVd,'k',kFold);% split trIdx in training and validation indices
    % the reason behind 7 in cvpartition is because with a number higher
    % than 7 some folds do not contain points from all the groups
    Training=cell(kFold,3); %  kfold x (interval left matrix, interval right matrix, labels)
    Validation=cell(kFold,3);
    for i = 1:CVP.NumTestSets
        trIdx_MS = CVP.training(i); % training index in the model selection step
        vdIdx_MS = CVP.test(i); % validation index in the model selection step
        Training{i,1}=X1TrVd(trIdx_MS,:);
        Training{i,2}=X2TrVd(trIdx_MS,:);
        Training{i,3}=y1TrVd(trIdx_MS);
        
        Validation{i,1}=X1TrVd(vdIdx_MS,:);
        Validation{i,2}=X2TrVd(vdIdx_MS,:);
        Validation{i,3}=y1TrVd(vdIdx_MS);
    end
    %test data set
    %---------------------------------------------------------------
    Test=cell(1,3);
    Test{1,1}=X1N(teIdx,:); Test{1,2}=X2N(teIdx,:); Test{1,3}=y1N(teIdx,:);
    
    wholeTraining=cell(1,3);
    wholeTraining{1,1}=X1N(trIdx,:); wholeTraining{1,2}=X2N(trIdx,:); wholeTraining{1,3}=y1N(trIdx,:);
    % SVM Stuff
    %---------------------------------------------------------------
    
    kernelList=[1 3:25];
    %kernelList=[22];
    % FOR each kernel do
    for kernelOption = kernelList
        % Model selection using training and validation sets
        %[bestC,bestGamma, accuracy, varianza, nrosvs,grid_accuracy,grid_errorA,grid_errorN,grid_var,grid_nrosvs,grid_AUC]=gridSearch(Training, Validation, kFold,kernelOption);
        [bestC,bestGamma, accuracy, std, nrosvs, Err_RateA,Err_RateN,AUC]=gridSearch(Training, Validation, kFold,kernelOption);
        
        %cross validation statistics:
        %----------------------------
        statistics{nB,kernelOption,1}= [bestC,bestGamma, accuracy, std, nrosvs, 100-Err_RateA,100-Err_RateN,AUC];

        
        %Training in the whole training data set (training + validation)
        %with the best parameters                
        [b,alph1,pos, nSV, Err_Rate, Err_RateA,Err_RateN,~,AUC,stat]=getStatisticsSVM(wholeTraining, wholeTraining,kernelOption, bestGamma,bestC);
        %training set statistics
        %----------------------------
        
        statistics{nB,kernelOption,2}=[nSV, 100-Err_Rate, 100-Err_RateA,100-Err_RateN,AUC,stat];
        
        nSV=length(pos);
        
        %Testing        
        if (kernelOption == 1 ||kernelOption == 2 || kernelOption == 3)                        
            X=(wholeTraining{1}(pos,:)+wholeTraining{2}(pos,:))/2;
            Z=(Test{1}+Test{2})/2;
        else
            X=cell(1,3);% only considering the support vectors set
            X{1,1}=wholeTraining{1}(pos,:); X{1,2}=wholeTraining{2}(pos,:); X{1,3}=wholeTraining{3}(pos,:);
            Z=Test;
        end
        y=wholeTraining{3};
        yTest= Test{1,3};
        V=getKernel(kernelOption,X,Z,bestGamma);% prediction over test dataset
        ypred = V'*(y(pos).*alph1) + b;
        
        % Test statistics 
        %----------------------------
        [Err_Rate, Err_RateA,Err_RateN,ConfMat,AUC,stat]=computeError(yTest,ypred);        
        statistics{nB,kernelOption,3}=[100-Err_Rate, 100-Err_RateA,100-Err_RateN,AUC,stat];
        %---------------
        
        [kernelOption, accuracy, nrosvs, nSV, Err_Rate, Err_RateA,Err_RateN,AUC,stat];
        
    end
    display('boostraping set nro');
    nB
    
end

save expLong4.mat
