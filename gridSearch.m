%function [bestC,bestGamma, accuracy, varianza, nrosvs,grid_accuracy,grid_errorA,grid_errorN,grid_var,grid_nrosvs,grid_AUC]=gridSearch(Training, Validation, kFold,kernelOption)
function [bestC,bestGamma, accuracy, standarDev, nrosvs, Err_RateA,Err_RateN,AUC]=gridSearch(Training, Validation, kFold,kernelOption)


accuracy_history=zeros(1,kFold);
accuracy_history_errorA=zeros(1,kFold);
accuracy_history_errorN=zeros(1,kFold);
accuracy_history_AUC=zeros(1,kFold);
sv_history=zeros(1,kFold);

grid_errorA=zeros(21,21);
grid_errorN=zeros(21,21);
grid_AUC=zeros(21,21);
grid_accuracy=zeros(21,21);
grid_std=zeros(21,21);
grid_nrosvs=zeros(21,21);

%----------------------------------------------------------------
gammaFunction=@(x) 2.^x;
paramC=2.^(-5:15);
paramGamma=gammaFunction(5:-1:-15);

%--------------------------------------------------------------------------------------

i=1;
for C=2.^(-5:15)
    
    j=1;
    for gamma=gammaFunction(5:-1:-15)
        
        for k = 1:kFold
            
            [~,~,~, nsV, Err_Rate,  Err_RateA,Err_RateN,~,AUC,~]=getStatisticsSVM(Training(k,:), Validation(k,:),kernelOption, gamma,C);
            %[nsV, Err_Rate,  Err_RateA,Err_RateN,AUC]
            
            accuracy_history(k)=100-Err_Rate;
            accuracy_history_errorA(k)=Err_RateA;
            accuracy_history_errorN(k)=Err_RateN;
            accuracy_history_AUC(k)=AUC;
            sv_history(k)=nsV;
            
            
        end
        
        grid_accuracy(i,j)=mean(accuracy_history);
        grid_errorA(i,j)=mean(accuracy_history_errorA);
        grid_errorN(i,j)=mean(accuracy_history_errorN);        
        grid_AUC(i,j)=mean(accuracy_history_AUC);
        grid_std(i,j)=std(accuracy_history);
        grid_nrosvs(i,j)=mean(sv_history);
        j=j+1;
    end
    i=i+1;
end

%--------------------------------------------------------------------------------------

[index_C,index_gamma]=bestParammeters(grid_accuracy, grid_nrosvs);
bestC     = paramC(index_C);
bestGamma = paramGamma(index_gamma);
accuracy  = grid_accuracy(index_C,index_gamma);
standarDev    = grid_std(index_C,index_gamma);
nrosvs    = grid_nrosvs(index_C,index_gamma);
Err_RateA =grid_errorA(index_C,index_gamma);
Err_RateN= grid_errorN(index_C,index_gamma);
AUC=grid_AUC(index_C,index_gamma);
end



%----------------------Search the best parammeteres---return the indexes---------------------------------------
function [f,c]=bestParammeters(grid_acuracy, grid)
[row,cols]=size(grid_acuracy);
zeroMatrix=zeros(row,cols);
index=find(grid_acuracy==max(grid_acuracy(:)));
zeroMatrix(index)=1;

% uncommnet to solve ties with number of sv's
% % to solve ties, find minimum nSV position
% zeroMatrix=grid_nrosvs.*zeroMatrix;
% [fil,col]=find(zeroMatrix==min(zeroMatrix(index)));
% f=fil(1);
% c=col(1);

% to solve ties, find the maximum AUC position
zeroMatrix=grid.*zeroMatrix;
[fil,col]=find(zeroMatrix==max(zeroMatrix(index)));
f=fil(1);
c=col(1);

end
