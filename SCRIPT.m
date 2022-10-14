%% The script peforms the two-sample t-test using the transposition test
%%  reported in
%
% Chung, M.K., Xie, L, Huang, S.-G., Wang, Y., Yan, J., Shen, L. 2019 Rapid 
% acceleration of the permutation test via transpositions, International 
% Workshop on Connectomics in NeuroImaging, in press. 
% http://www.stat.wisc.edu/~mchung/papers/chung.2019.CNI.pdf
%
% (C) 2022 Moo K. Chung
%  mkchung@wisc.edu
% University of Wisconsin-Madison
%
% The code is downloaded from https://github.com/laplcebeltrami/transpositions
%
% 2022 October 14

%The transposition test on the whole brain or network can be done by vecotrizing data.  
%Suppose we generated the following simulated vector data, where group 1 has m=10 subjects 
% and group 2 has n=12 subjects over l=1000 edges. 

m=10; n=12; l=1000; per_t=100000;
x= normrnd(0,1,m,l);   %x is size 10 (# of subjects) x 1000 (# of edges)
y= normrnd(1,1,n,l); %y is size 12 (# of subjects) x 1000 (# of edges)

%Then two-sample test via the transposition test is done as
[stat_t, time_t] = test_transpose(x,y,per_t); 
%stat_t is size 100000 (# of transpostions) x 1000 (# of edges)

%% Local inference
%We are interested in determining edges that shows the connectivity differences.
%Edges with smaller p-values are the edges where the connectivity difference occurs. 

%The p-value for t-test is then computed as
observed=(mean(x)-mean(y)).*sqrt(m*n*(m+n-2)./((m+n)*((m-1)*var(x)+(n-1)*var(y)))); 
          %This is the observed two-sample t-statistic at each edge
pvalue_t = online_pvalue(stat_t, observed);
%pvalue_t is size 100000 (# of transpostions) x 1000 (# of edges).
%The p-value at the end of 100000 transpositions is given as 
% pvalue_t(end,:), which gives the p-values at each edge.
%To find edges that gives p-value smaller than 0.001
find(pvalue_t(end,:)<0.001)

%The plot of how p-value converges at each edge. 
figure; plot(pvalue_t(:,1)) 

%% Global inference
%However, if we are interested in testing if two networks are different 
%using the above two sample t-test, the p-value formula has to be modified
%to account for multiple comparsisions using sup t-stat formula in page 48
%(equation right above the Validation section, where sup notation is used)
%in http://www.stat.wisc.edu/~mchung/papers/chung.2019.CNI.pdf.

observed=max(observed) %maximum observed t-statistic over the whole brain
pvalue = online_pvalue_corrected(stat_t, observed);
%The corrected p-value is of size 100000 (# of transpostions) x 1 
%At the end of transpostions, we have the corrected p-value 
pvalue(end)

%The convergence plot for corrected p-value is 
figure; plot(pvalue)




