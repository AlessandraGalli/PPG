%-----------------------------%%-----------------------------%%
%  Authors: Alessandra Galli, Claudio Narduzzi, Giada Giorgi. %
%        Instrumentation and Measurement Research Group       %
%                    University of Padova                     %
%-----------------------------%%-----------------------------%%
clc
clear all
close all

TOT_EST=[];
TOT_TRUE=[];
PPG1=[];
PPG2=[];

for k=1:12    
    
    eval(['load(''./DATI_INTERMEDI/TR',num2str(k),''');'])   
    [res]=kalman_filter(z_ppg1,z_ppg2);


    errore.v1(k)=mean(abs(BPM0(1:end)-res(1:end)'));
    errore.v2(k)=mean(abs(BPM0(1:end)-res(1:end)')./(BPM0(1:end)))*100;
    errore.v3(k)=max(abs(BPM0(1:end)-res(1:end)'));
    errore.v4(k)=std(abs(BPM0(1:end)-res(1:end)'));
 
    errore.r1_2s(k)=mean(abs(BPM0(1:end-1)-res(2:end)'));
    errore.r2_2s(k)=mean(abs(BPM0(1:end-1)-res(2:end)')./(BPM0(1:end-1)))*100;
    errore.r3_2s(k)=max(abs(BPM0(1:end-1)-res(2:end)'));
    errore.r4_2s(k)=std(abs(BPM0(1:end-1)-res(2:end)'));
   
    errore.r1_4s(k)=mean(abs(BPM0(1:end-2)-res(3:end)'));
    errore.r2_4s(k)=mean(abs(BPM0(1:end-2)-res(3:end)')./(BPM0(1:end-2)))*100;
    errore.r3_4s(k)=max(abs(BPM0(1:end-2)-res(3:end)'));
    errore.r4_4s(k)=std(abs(BPM0(1:end-2)-res(3:end)'));

    TOT_EST=[TOT_EST res(3:end)];
    TOT_TRUE=[TOT_TRUE, BPM0(1:end-2)'];

    PPG1=[PPG1 z_ppg1(3:end)];
    PPG2=[PPG2 z_ppg2(3:end)];


end 


