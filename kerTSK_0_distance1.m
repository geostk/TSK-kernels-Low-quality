% substitution fuzzy distance kernel with TSK_0 
% % Eq 6.16 (Thesis)                
% (it is not published yet)

function G=kerTSK_0_distance1(dataX,dataZ,gamma)
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
        diff_XZ=( X(i,:)- Z(j,:)).*( X(i,:)- Z(j,:));
        den_XZ=stdX(i,:).*stdX(i,:) +  stdZ(j,:).*stdZ(j,:);
        XZ=exp(- 0.5*sum(diff_XZ./den_XZ)  )*dot(X(i,:),Z(j,:));               
        D=2-2*XZ;        
        G(i,j)=exp(- 0.5*gamma*D^2  );
    end
end  		

