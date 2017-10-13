%-----------------------------%%-----------------------------%%
%  Authors: Alessandra Galli, Claudio Narduzzi, Giada Giorgi. %
%        Instrumentation and Measurement Research Group       %
%                    University of Padova                     %
%-----------------------------%%-----------------------------%%
function[p1,p2]=frequency_analysis(ppg1,ppg2,xr_1,xr_2)
   
    fS=125;                 
    Nfft=2^13;                                                             
    F=fS/Nfft;                                                              
    w=hann(length(xr_1))';                                                 
    w2=hann(length(xr_2))';                                                

    Xr_1=abs(fft(xr_1.*w,Nfft)/sum(w));                                    
    Xr_2=abs(fft(xr_2.*w2,Nfft)/sum(w2));                                  

    start=round(60/60/F);                                              
    fine=round(180/60/F);                                              
    [~, ind_pk_1]=findpeaks(Xr_1(start:fine),'sortstr','descend');      
    [~, ind_pk_2]=findpeaks(Xr_2(start:fine),'sortstr','descend');     

    w=hann(length(ppg1))';                                                 
    w2=hann(length(ppg2))';                                                
    X1=abs(fft(ppg1.*w,Nfft)/sum(w));                                       
    X2=abs(fft(ppg2.*w2,Nfft)/sum(w2));                                    
    [~, ind1]=findpeaks(X1(start:fine),'sortstr','descend');             
    [~, ind2]=findpeaks(X2(start:fine),'sortstr','descend');              


    if size(ind_pk_1)>0                                                     
        p1=(ind_pk_1(1)+start-2)*F*60;                                         
        picco1=(ind1(1)+start-2)*F*60;                                        
        if abs(2*picco1-p1)<10                                                
            p1=picco1;                                                         
        end 
    else
        p1=NaN;
    end

    if size(ind_pk_2)>0                                                     
        p2=(ind_pk_2(1)+start-2)*F*60;
        picco2=(ind2(1)+start-2)*F*60;
         if abs(2*picco2-p2)<10
             p2=picco2;
         end 
    else
        p2=NaN;
    end 
end
