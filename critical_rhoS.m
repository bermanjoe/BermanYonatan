function [lim1r,medr,lim2r]=critical_rhoS(XS1,XS2,NN,nquantiles)
%
% Input
% XS1: Cross section 1 (initial)
% XS2: Cross section 2 (terminal)
% NN: Number of iterations (for bootsrapping)
% nquantiles: Resolution of analysis (number of quantiles)
%
% Output
% The analysis result is a distribution of critical rank correlation
% values, which is given in the variable 'tmp'. The output of this function
% is the 5th percentile of the distribution (lim1r), the median (medr) and
% the 95th percentile (lim2r)
%

tmp=zeros(NN,1);
pp=linspace(0,1,nquantiles);

% Range of rank correlations considered, with higher resolution above 0.5
rcs=[0.01:0.01:0.49 0.5:0.005:0.695 0.7:0.001:0.998];

for j=1:NN
    
    rc=rcs(1);count=1;
    flag=0;
    nn=min([length(XS1) length(XS2)]);
    
    % For bootstrapping we consider 80% of the length of the cross section
    nn=min([nn ceil(8*nn/10)]);

    tmp1=randperm(length(XS1));a=XS1(tmp1(1:nn));
    tmp2=randperm(length(XS2));b=XS2(tmp2(1:nn));
    
    while (flag==0)
        
        [w1,w2]=couple_vecs_plackett(a,b,plackett_theta(rc));
        absmobs=absmobg_quantiles(w1,w2,nquantiles);
        
        gim=zeros(nquantiles,1);
        gimu=zeros(nquantiles,1);
        for iii=2:nquantiles-1
            gim(iii)=(1/pp(iii))*trapz(pp(1:iii),absmobs(1:iii));
            gimu(iii)=(1/(1-pp(iii)))*trapz(pp(iii:end),absmobs(iii:end));
        end
        if (sum(gim>gimu))>=(nquantiles-2)
            rc=rcs(count+1);count=count+1;
            if rc>=0.998 %critical rho above 0.998 is defined as 1
                rc=1;
                flag=1;
            end
        else
            flag=1;
        end
    end
    tmp(j)=rc;
end

% Create output
lim1r=prctile(tmp,5);
lim2r=prctile(tmp,95);
medr=prctile(tmp,50);