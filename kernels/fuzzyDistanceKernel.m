function k=fuzzyDistanceKernel(X,Y, param)
%   implements a kernel on fuzzy sets with substitution of distance  D=sum(abs(X-Y))/sum(abs(X+Y)) 
%   input:  X and Y are a Nxp and Mxp matrices, each row is a fuzzy set with p
%   elements, such X(i,j) is the degree membership of element j of the fuzzy set i%          
%
%    output kernel matriz k of size NxN with the similarity measures
%   
%    Example:
%    Asumme that X contain 10 fuzzy sets and Y contain 15 fuzzy sets:
%     %X=rand(10,100); Y=rand(15,100);
%      
%     k=fuzzyDistanceKernel(X,Y,1)

[m,n]=size(X);
[p,q]=size(Y);

    for i=1:m
        for j=1:p            
            D=sum(abs(X(i,:)-Y(j,:)))/sum(abs(X(i,:)+Y(j,:)));            
             k(i,j)=exp(-  param*D^2 );                                                
        end
    end


