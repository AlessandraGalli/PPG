function [H]=Hk(x)
L=400;      
N=1000;    
K=N-L+1;    
xc=x(1:K); 
xr=x(K:N); 
H=hankel(xc,xr);
end 
