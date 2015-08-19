function [b,alph1,pos, nSV, Err_Rate, Err_RateA,Err_RateN,ConfMat,AUC,stat]=getStatisticsSVM(Training, Validation,kernelOption, kernelParam,C)

% What kind of kernel is it?
% is not a kernel fuzzy
if (kernelOption == 1 ||kernelOption == 2 || kernelOption == 3)    
    X=(Training{1}+Training{2})/2;
    Z=(Validation{1}+Validation{2})/2;
else
    X=Training; Z=Validation;
end
%labels
y=Training{3}; yV=Validation{3};
% compute training kernel G  // getKernel(3,Xapp,Xapp,(1/gamma^2))
% is equal to svmkernel(Xapp,'gaussian',gamma);
G=getKernel(kernelOption,X,X,kernelParam);
% Train a SVM

[b,alph1,pos]=SVM_MONQP(G,y,C);
% compute validation kernel V
if (kernelOption == 1 ||kernelOption == 2 || kernelOption == 3)    
V=getKernel(kernelOption,X(pos,:),Z,kernelParam);   
else
   X={Training{1}(pos,:), Training{2}(pos,:)};
   V=getKernel(kernelOption,X,Z,kernelParam);
end

% svm test ...
ypred = V'*(y(pos).*alph1) + b;
% statistics = compare ypred agains yV
nSV=length(pos);
[Err_Rate, Err_RateA,Err_RateN,ConfMat,AUC,stat]=computeError(yV,ypred);


%--------------------------------
%TestCase
% Example
% n = 500; % upto n = 10000;
% sigma=1.4;
% [Xapp,yapp,Xtest,ytest]=dataset_KM('checkers',n,n^2,sigma);
% kernelOption=3; kernelParam=0.5;
% G=getKernel(kernelOption,X,X,kernelParam);
% [b,alph1,pos]=SVM_MONQP(G,y,C)
% V=getKernel(kernelOption,Z,X(pos,:),kernelParam);
% ypred = V*(y(pos).*alph1) + b;
% [Err_Rate, Err_RateA,Err_RateN,ConfMat,AUC,stat]=computeError(yV,ypred);



