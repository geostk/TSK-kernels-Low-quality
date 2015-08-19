function G = fuzzyCrossProductKernel(varargin)
%   The cross product kernel on fuzzy sets
%   fuzzyCrossProductKernel(X,Z,opt1,opt2,param1,param2)
%   inputs: 
%   X, Z are cells containing M and N fuzzy sets, respectively
%        Each fuzzy set is modeled by a matriz of size Ni*(D+1) (Ni
%        d-simensional points for the fuzzy set i).
%        The last column is the degree of membershiop of the point to the fuzzy
%        set i.
%   opt1 and opt2 are the kernel options opt1=opt2=[1=linear, 2=polynomial, 3=gaussian]
%   param1,param2 are the kernel parameters for both kernels respectively
X =varargin{1};
Z=varargin{2};
opt1=varargin{3};
opt2=varargin{4};
param1=cell2mat(varargin(3:end));
param2=cell2mat(varargin(3:end));

[m,~]=size(X);
[p,~]=size(Z);

if opt==1&& opt==1
    
    
    
end

if opt==1&& opt==3