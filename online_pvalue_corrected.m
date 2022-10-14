function pvalue = online_pvalue_corrected(stat, observed)
%function pvalue = online_pvalue_corrected(stat, observed)
%
% INPUT
% stat:     test statistic computed over transpositions. It is of size 
%           # of transpostions) x # of edges
% observed: it's the maximum observed test satistic over all edges 
%
% Computes the pvalue (corrected for multiple comparisions) based on the 
% collection of statistics value and observation sequentically. 
% It is needed for various applications
% where we need to know how p-values change. The seqencial p-value 
% computation algorithm is given in 
%
% Chung, M.K., Xie, L, Huang, S.-G., Wang, Y., Yan, J., Shen, L. 2019 Rapid 
% acceleration of the permutation test via transpositions, International 
% Workshop on Connectomics in NeuroImaging, in press. 
% http://www.stat.wisc.edu/~mchung/papers/chung.2019.CNI.pdf
%
% The code gives a single pvalue over the whole networks. 
% If you want to compute the pvalue at each edge, use online_pvalue.m
%
% This code is downloaded from
% https://github.com/laplcebeltrami/transpositions
%
% (C) 2022 Moo K. Chung  mkchung@wisc.edu
% University of Wisconsin-Madison
%
% Update histroy: 2019 May created
%                 2021 Oct. 11 Validation done 
%                 2022 Oct. 14 Multiple comparsions correction added
% 
%
% Validaiton against built-in Matlab function normrnd.m with 10000 N(0,1)
% random numbers. 
% norminv(0.95,0,1) = 1.6449 treshold corresponding to pvalue = 1-0.95 probability
% norminv(0.5,0,1) = 0
% norminv(0.05,0,1) = -1.6449
%
% Based on the above ground truth, we should get p-values
%
% stat=normrnd(0,1,10000,1);
% observed =-1.6449  %1.6449  
% pvalue = online_pvalue(stat, observed)
% pvalue(end) 
% 0.0482  0.4889   0.0489

% pvalue can be computed iteratively as
% pvalue(i+1) = (pvalue(i) * i  + ( t(i+1)>=observed ) / (i+1)

n=size(stat,1); % number of statistics.
l=size(stat,2); % numver of variables/connections

pvalue=zeros(n,1); %pvalue will be iteratively updated 

maxstat = max(stat, [], 2); %maximum statistics over whole brain network

if observed>=0
    
    pvalue(1)= (maxstat(1)>=observed); %initial p-value. It is either 0 or 1.
    for i=2:n
        pvalue(i) = (pvalue(i-1) * (i-1) + (maxstat(i) >= observed))/i;
    end
    
else %observed<0
    
    pvalue(1)= (maxstat(1)<=observed); %initial p-value. It is either 0 or 1.
    for i=2:n
        pvalue(i) = (pvalue(i-1) * (i-1) + (maxstat(i) <=observed))/i;
    end
    
end


