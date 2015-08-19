function k=intersectionKernel(X,Y,TNorm)
%   Positive definite kernel functions on fuzzy sets. Guevara J, et al.
%   input:  X and Y are a Nxp matrices, each row is a fuzzy set with p
%   elements, such X(i,j) is the degree membership of element j of the fuzzy set i
%           TNorm:  1= minimum 2=produto 3=lukasiewics 4=drástico
%
%    output kernel matriz k of size NxN with the similarity measures
%
%    Example:
%    Asumme that X contain 10 fuzzy sets and Y contain 15 fuzzy sets:
%     %X=rand(10,100); Y=rand(15,100);
%      
%     k=intersectionKernel(X,Y,1);

[m,n]=size(X);
[p,q]=size(Y);

if TNorm==1
    for i=1:m
        for j=1:p
            k(i,j)=sum(min(X(i,:),Y(j,:))); %think to normalize
        end
    end
end
if TNorm==2
    for i=1:m
        for j=1:p
            k(i,j)=sum(X(i,:).*Y(j,:)); %think to normalize
        end
    end
end