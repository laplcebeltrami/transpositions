function [stat_t, time_t] = test_transpose(x,y,per_t)
%function [stat_t, time_t] = test_transpose(x,y,per_t)
% The function computes the two-sample t-statitic of vector data using the transposition
% test. If follows the method explained in 
%
% Chung, M.K., Xie, L, Huang, S.-G., Wang, Y., Yan, J., Shen, L. 2019 Rapid 
% acceleration of the permutation test via transpositions, International 
% Workshop on Connectomics in NeuroImaging, in press. 
% http://www.stat.wisc.edu/~mchung/papers/chung.2019.CNI.pdf
%
% INPUT
% x    : input data of size m x l (m= number of subjects, l=number of data per subject) 
% y    : input data of size n x l (n= number of subjects, l=number of data per subject)
% per_t: number of transpositions
%
% OUTPUT
% stat_t:  two-sample t-statistic of all transpositions
% time_t:  run time it took to compute the statistics
%
%
% This code is downloaded from
% http://www.stat.wisc.edu/~mchung/transpositions
%
% (C) 2019 Moo K. Chung, Yixian Wang  
%  mkchung@wisc.edu
% University of Wisconsin-Madison


tic

m=size(x,1);
n=size(y,1);

l=size(x,2); %dimension of vector
stat_t=zeros(per_t,l); %computed statistic over transpostions will be saved as stat_t.



pi1=randi(m,1,per_t); %generate random transposition indices in group 1
pi2=randi(n,1,per_t); %generate random transposition indices in group 2

z=[x;y];
z=z(randperm(m+n));
x=z(1:m);y=z(m+1:m+n); %initial random shuffle

%------------------
%initial values of sum and squared sum functions
mx1=sum(x);my1=sum(y);%initial summation (with no division)
mean_x=mean(x); mean_y=mean(y);
vx1=sum((x-mean_x).*(x-mean_x)); vy1=sum((y-mean_y).*(y-mean_y));%initial squared sum



%random sequantical transposition


for i=1:per_t
    
    a=x(pi1(i));x(pi1(i))=y(pi2(i));
    b=y(pi2(i));y(pi2(i))=a;%tranpose one element between x and y
    
    mx2=mx1+b-a;my2=my1+a-b;%update summation function
    
    vx2=vx1+(mx1.^2-mx2.^2)/m+b.^2-a.^2; %update squared sum function
    vy2=vy1+(my1.^2-my2.^2)/n+a.^2-b.^2; 
    
    stat_t(i,:)=(mx2/m-my2/n).*sqrt(m*n*(m+n-2)./((m+n)*(vx2+vy2)));%update t-statistic
    
    mx1=mx2;my1=my2;vx1=vx2;vy1=vy2;%prepare for the next transposition
    %   i=i+1;
    %   time=toc;
end

toc
time_t = toc;
