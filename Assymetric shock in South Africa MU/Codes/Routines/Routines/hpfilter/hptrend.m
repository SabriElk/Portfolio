function  [yt,yd]=hptrend(y,phi)
%HPTREND  Detrend time series with Hodrick-Prescott method.
%         [yt,yd]=hptrend(y,phi) finds the series yt(:,i) that minimizes:
%
%         sum  (y(t,i)-yt(t,i))^2 + phi sum [(yt(t+1,i)-yt(t,i))- ..
%        t=1:T                        t=2:T-1      -(yt(t,i)-yt(t-1,i))]^2
%
%         for each column i=1,...k in y.   A larger phi results in a 
%         smoother trend series.  For quarterly data Hodrick and Prescott 
%         (1980) use phi=1600.
%
%         Also returned are the series with no trend:  yd(:,i)=y(:,i)-
%         yt(:,i), i=1,...k
%
%          

%         Ellen R. McGrattan,  4-23-87
%         Revised, 3-22-92, ERM

%         References
%         ----------
%         [1] Hodrick, Robert J. and Edward C. Prescott, ``Post-War
%             U.S. Business Cycles: An Empirical Investigation,''
%             Working Paper, Carnegie-Mellon University, November, 1980.
%         [2] Prescott, Edward C., ``Theory Ahead of Business Cycle 
%             Measurement,'' QUARTERLY REVIEW, Federal Reserve Bank 
%             of Minneapolis, Fall 1986.
%      
[T,k]=size(y);
if T<301;
  A = [1+phi   -2*phi   phi                      zeros(1,T-3);
      -2*phi  1+5*phi  -4*phi  phi               zeros(1,T-4);
       phi    -4*phi                             zeros(1,T-2);
       0        phi                              zeros(1,T-2);
                        zeros(T-8,T);
      zeros(1,T-2)                                   phi        0;
      zeros(1,T-2)                                 -4*phi     phi;
      zeros(1,T-4)                    phi  -4*phi  1+5*phi -2*phi;     
      zeros(1,T-3)                           phi   -2*phi   1+phi];
  A(3:T-2,3:T-2) = phi*(diag(ones(T-6,1),2)+diag(ones(T-6,1),-2)) ... 
                -4*phi*(diag(ones(T-5,1),1)+diag(ones(T-5,1),-1)) ... 
             +(1+6*phi)*diag(ones(T-4,1));
else
  i   = [1,1,2,2,1:T-2,2:T-2,3:T-2,3:T-1,3:T  ,T-1,T-1,T  ,T]';
  j   = [1,2,1,2,3:T  ,3:T-1,3:T-2,2:T-2,1:T-2,T-1,T  ,T-1,T]';
  s   = [1+phi,-2*phi,-2*phi,1+5*phi,phi*ones(1,T-2),-4*phi*ones(1,T-3), ...
         (1+6*phi)*ones(1,T-4),-4*phi*ones(1,T-3),phi*ones(1,T-2), ...
         1+5*phi,-2*phi,-2*phi,1+phi]'; 
  A   = sparse(i,j,s);
end;
yt = A\y;
yd = y-yt;


