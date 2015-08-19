function [X1N,X2N,yN]=preprocessing(X1,X2,y,val)
j=1;
%duplicate data
for i=1:length(val)
    if strcmp(val(i),'1}')                
    else
        X1N(j,:)=X1(i,:);
        X2N(j,:)=X2(i,:);
        yN(j)=y(i);
        j=j+1;
    end
end
yN=grp2idx(yN);
%set labels to +1 , -1
yN(find(yN==2))=-1;