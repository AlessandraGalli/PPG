%-----------------------------%%-----------------------------%%
%  Authors: Alessandra Galli, Claudio Narduzzi, Giada Giorgi. %
%        Instrumentation and Measurement Research Group       %
%                    University of Padova                     %
%-----------------------------%%-----------------------------%%
function [H]=Hk(x)
L=400;      
N=1000;    
K=N-L+1;    
xc=x(1:K); 
xr=x(K:N); 
H=hankel(xc,xr);
end 
