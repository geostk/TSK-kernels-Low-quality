function G = getKernel(option, varargin)
%varargin = X,Z, parameters=[...]
X =varargin{1};
Z=varargin{2};
kernelParam=cell2mat(varargin(3:end));
beta_par=0.5;
switch lower(option)
    %--------------------------------------
    % p(x| mu,cov)=1/L, i.e uniform
    case 1% linear kernel
        %disp('linear kernel')
        G=X*Z';
    case 2 %polinomial kernel
        %G=((X*Z'+kernelParam(1)).^kernelParam(2));
        G=((X*Z').^kernelParam);
    case 3 % gaussian kernel
        G=exp(-0.5*kernelParam*sqdistAll(X,Z));
        %-------------------------------------------------------------------
        % kernels TSK nonsingleton (antecedent part of rules) gaussian MF
        %-------------------------------------------------------------------
    case 4 % gaussian TSK no parameter       
        G=kerTSK_0(X,Z);
    case 5 % gaussian TSK with parameter
        G=kerTSK_1(X,Z,kernelParam);
    case 6 % gaussian TSK with parameter (is kerTSK_0 with lambda parameter as in RBF)
        G=kerTSK_2(X,Z,kernelParam);
        %-------------------------------------------------------------------
        % kernels TSK nonsingleton (antecedent + consequent part of rules)
        %-------------------------------------------------------------------
    case 7 % Eq 6.16 (Thesis)
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=kerTSK_0(X,Z).*(XX*ZZ');        
    case 8 % same as Eq 6.16 (Thesis) but with kernel TSK1
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=kerTSK_1(X,Z,kernelParam).*(XX*ZZ');
    case 9 % same as Eq 6.16 (Thesis) but with kernel TSK2
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=kerTSK_2(X,Z,kernelParam).*(XX*ZZ');
    case 10 % same as Eq 6.17 (Thesis) but with kernel TSK0
        stdX=abs (X{1}-X{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))
        stdZ=abs (Z{1}-Z{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))
        stdX(stdX==0)=0.00000000001; %12 x 4
        stdZ(stdZ==0)=0.00000000001; %2 x 4
        mX=(X{1}+X{2})./2;           %12 x 4
        mZ=(Z{1}+Z{2})./2;           %2 x 4
        
        
        
        stdX_2=stdX.*stdX; 
        stdZ_2=stdZ.*stdZ;
        
        [m,~]=size(stdX);
        [n,~]=size(stdZ);
        GG=zeros(m,n);
        for i=1:m
            for j=1:n
                GG(i,j)=dot(mX(i), (mX(i).*stdZ_2(j)+mZ(j).*stdX_2(i))/(stdX_2(i)+stdZ_2(j)));
            end
        end
        
        % as theta paremeters are for fuzzy sets in the rule, then:
        %theta=(mX.*(stdZ.*stdZ) + mZ.*(stdZ.*stdZ))./((stdZ.*stdZ) +(stdZ.*stdZ));
        G=kerTSK_0(X,Z).*GG;
    case 11 % Eq 6.17 (Thesis)
        
        stdX=abs (X{1}-X{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))
        stdZ=abs (Z{1}-Z{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))
        stdX(stdX==0)=0.00000000001; %12 x 4
        stdZ(stdZ==0)=0.00000000001; %2 x 4
        mX=(X{1}+X{2})./2;           %12 x 4
        mZ=(Z{1}+Z{2})./2;           %2 x 4
        
        
        
        stdX_2=stdX.*stdX; 
        stdZ_2=stdZ.*stdZ;
        
        [m,~]=size(stdX);
        [n,~]=size(stdZ);
        GG=zeros(m,n);
        for i=1:m
            for j=1:n
                GG(i,j)=dot(mX(i), (mX(i).*stdZ_2(j)+mZ(j).*stdX_2(i))/(stdX_2(i)+stdZ_2(j)));              
            end
        end
        
              
        G=kerTSK_1(X,Z,kernelParam).*GG;
        
    case 12 % same as Eq 6.17 (Thesis) but with kernel TSK2
        stdX=abs (X{1}-X{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))
        stdZ=abs (Z{1}-Z{2})/2.3548; %sigma=(l-r)*/2*sqrt(2*log(2))
        stdX(stdX==0)=0.00000000001; %12 x 4
        stdZ(stdZ==0)=0.00000000001; %2 x 4
        mX=(X{1}+X{2})./2;           %12 x 4
        mZ=(Z{1}+Z{2})./2;           %2 x 4
        
        
        
        stdX_2=stdX.*stdX; 
        stdZ_2=stdZ.*stdZ;
        
        [m,~]=size(stdX);
        [n,~]=size(stdZ);
        GG=zeros(m,n);
        for i=1:m
            for j=1:n
                GG(i,j)=dot(mX(i), (mX(i).*stdZ_2(j)+mZ(j).*stdX_2(i))/(stdX_2(i)+stdZ_2(j)));           
            end
        end
        
        G=kerTSK_2(X,Z,kernelParam).*GG;
    case 13 % Eq 6.18 (Thesis)
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=kerTSK_0(X,Z).*exp(-kernelParam*sqdistAll(XX,ZZ));
        
    case 14 % same as Eq 6.18 (Thesis) but with kernel TSK1
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        
        %median heuristic
        if length(XX)+length(ZZ) > 1000
            kPer=1000;
            indX=randperm(length(XX),kPer) ;
            indZ=randperm(length(ZZ),kPer) ;
            MM=[XX(indX,:);ZZ(indZ,:)];
        else
            MM=[XX;ZZ];
        end
        M=sqdistAll(MM,MM);
        lambda=median(M(:));
        %-----------------------
        
        G=kerTSK_1(X,Z,kernelParam).*exp(-(1/lambda)*sqdistAll(XX,ZZ));
        
    case 15 % same as Eq 6.18 (Thesis) but with kernel TSK2
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        %median heuristic
        if length(XX)+length(ZZ) > 1000
            kPer=1000;
            indX=randperm(length(XX),kPer) ;
            indZ=randperm(length(ZZ),kPer) ;
            MM=[XX(indX,:);ZZ(indZ,:)];
        else
            MM=[XX;ZZ];
        end
        M=sqdistAll(MM,MM);
        lambda=median(M(:));
        %-----------------------
        G=kerTSK_2(X,Z,kernelParam).*exp(-(1/lambda)*sqdistAll(XX,ZZ));
        
        
        
        %--------------------------------------------------------------------
        % substitution fuzzy distance kernel
        %--------------------------------------------------------------------
    case 16  % kernel TSK_0 as distance substitution kernel
        
        G=kerTSK_0_distance(X,Z,kernelParam);
        
    case 17 % kernel TSK_0 (Eq 6.16 (Thesis)) as distance substitution kernel
        G=kerTSK_0_distance1(X,Z,kernelParam);
        
                        
    case 18 % lemma 17 thesis
        G=kerD1(X,Z,kernelParam);
    case 19 % lemma 17 thesis
        G=kerD2(X,Z,kernelParam);
        % multiple kernel fuzzy
        
    case 20 % lemma 17 thesis (variation of)
        G=kerD3(X,Z,kernelParam);
        
    case 21 % lemma 17 thesis (variation of)
        G=kerD4(X,Z,kernelParam);
        
        %----------------------------------------------------------------------
        % multiples kernels
        %----------------------------------------------------------------------
        
    case 22 % kernel id = 6
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=beta_par*kerTSK_0(X,Z)+(1-beta_par)*exp(-0.5*kernelParam*sqdistAll(XX,ZZ));
        
    case 23 % kernel id = 7
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=beta_par*exp(-0.5*kernelParam*sqdistAll(XX,ZZ)).*kerTSK_0(X,Z)+(1-beta_par)*exp(-0.5*kernelParam*sqdistAll(XX,ZZ));
        
    case 24 % kernel id = 8
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=beta_par*kerD1(X,Z,kernelParam)+(1-beta_par)*exp(-0.5*kernelParam*sqdistAll(XX,ZZ));
        
    case 25 % kernel id = 9
        XX=(X{1}+X{2})./2;
        ZZ=(Z{1}+Z{2})./2;
        G=beta_par*kerD2(X,Z,kernelParam)+(1-beta_par)*exp(-0.5*kernelParam*sqdistAll(XX,ZZ));
        
    otherwise
        disp('Unknown method.')
        
        % fuzzyDistanceKernel.m intersectionKernel.m, considerar cada
        % variable como por ejemplo X{i,4} = gaussmf(sample,[ (temp(2)-temp(1))/(2*sqrt(2*log(2))) ,mean(temp)]);
        
      %---------------------------------
      % fuzzy distance with metrics
      % intersection kernel
      % distance substitution with intersection kernel
      % cross product kernel
      % distance substitution with cross product kernels
end




