%substitution pseudometric between fuzzy sets kernel
% June 2015
function G=kerD1(dataX,dataZ,gamma)
% input dataX = cell {rigth interval values, left interval values}, also  data Z
% output G = kernel matrix

%compute parameters for gaussian fuzzy sets
factor=1; % default
if (factor ==1)
    stdX=abs (dataX{1}-dataX{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))    
    stdZ=abs (dataZ{1}-dataZ{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))  
else
    stdX=abs (dataX{1}-dataX{2})*0.1667; %sigma=(l-r)*1/6
    stdZ=abs (dataZ{1}-dataZ{2})*0.1667; %sigma=(l-r)*1/6
end
stdX(stdX==0)=0.00000000001;
stdZ(stdZ==0)=0.00000000001;
X=(dataX{1}+dataX{2})./2;
Z=(dataZ{1}+dataZ{2})./2;

%------------
[m,~]=size(X);
[p,~]=size(Z);

G=zeros(m,p);


for i=1:m
    for j=1:p
        D1=prod(  ( X(i,:)- Z(j,:)).*( X(i,:)- Z(j,:)) )/( stdX(i,:)*stdX(i,:)' +  stdZ(j,:)*stdZ(j,:)');        
        G(i,j)=exp(-  gamma*D1^2 );
    end
end
