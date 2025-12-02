function [coeff,p_val]=autocor_r(X,m)
% AUTOCOR computes cross-correlations and autocorrelations
%         A=autocor_r(X,m) takes a T x N matrix X and constructs
%
%                       [ B[0] ]
%                       [ B[1] ]
%                     A=[ B[2] ]
%                       [ B[3] ]
%                       [  ... ]
%                       [ B[m] ]
%                
%                 where the (i,j) element of B[k] is corr(X[i,t],X[j,t-k]),
%                 with i,j in [1,...,N].
%                 The matrix A has dimensions N*(m+1) *N
%                 The matrix B contains the p-values associated with the
%                 cross-correlations coefficients.
%
%                 Ellen R. McGrattan, 2-1-99, modified by Romain Houssa to
%                 include p-values to reject the null hypothesis.
%
[T,N]=size(X);
if T-m<=0;
    str='The Number of observations in the time series must be larger than m';
    error(str);
end
b=zeros(N);
bp=zeros(N);

coeff=zeros(N*(m+1),N);
p_val=zeros(N*(m+1),N);

for k=1:m+1;
    for i=1:N;
        for j=1:N;
            [x,x1]=corrcoef([X(k:T,i),X(1:T-k+1,j)]);
            b(i,j)=x(2,1);
            bp(i,j)=x1(2,1);
        end
    end
    coeff((k-1)*N+1:k*N,:)=b
    p_val((k-1)*N+1:k*N,:)=bp
end