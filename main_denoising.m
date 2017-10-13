%-----------------------------%%-----------------------------%%
%  Authors: Alessandra Galli, Claudio Narduzzi, Giada Giorgi. %
%        Instrumentation and Measurement Research Group       %
%                    University of Padova                     %
%-----------------------------%%-----------------------------%%
close all
clear all

filename = '.\DATABASE\Training_data\DATA_05_TYPE02.MAT';
matObj = matfile(filename);
Tracce=matObj.sig;

fn='.\DATABASE\Training_data\DATA_05_TYPE02_BPMtrace.MAT';
matObj=matfile(fn);
BPM0=matObj.BPM0;

if(min(size(Tracce))==6)    %training trace
    ECG=Tracce(1,:);
    PPG1=Tracce(2,:);
    PPG2=Tracce(3,:);
    AX=Tracce(4,:);
    AY=Tracce(5,:);
    AZ=Tracce(6,:);
else                       % competition trace
    PPG1=Tracce(1,:);
    PPG2=Tracce(2,:);
    AX=Tracce(3,:);
    AY=Tracce(4,:);
    AZ=Tracce(5,:);
end

N=length(PPG1); 
K=1000;         
DK=250;         
i=1; fin=0;
n=length(Tracce);
while fin < n - DK

    in=(i-1)*DK+1;
    fin=in+K-1;
    ppg1=PPG1(in:fin); 
    ppg2=PPG2(in:fin);
    accx=AX(in:fin);
    accy=AY(in:fin);
    accz=AZ(in:fin);
    [temp1,temp2]=signal_denoising(ppg1,ppg2,accx,accy,accz);
    [z_ppg1(i) z_ppg2(i)] = frequency_analysis(ppg1,ppg2,temp1,temp2);
    i=i+1;
end 
