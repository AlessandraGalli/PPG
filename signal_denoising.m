%-----------------------------%%-----------------------------%%
%  Authors: Alessandra Galli, Claudio Narduzzi, Giada Giorgi. %
%        Instrumentation and Measurement Research Group       %
%                    University of Padova                     %
%-----------------------------%%-----------------------------%%
function [ xr_1,xr_2 ] = signal_denoising(ppg1,ppg2,accx,accy,accz)
 %% Construction of Hankel matrices and SVD
    H1=Hk(ppg1);                 
    [U1,S1,V1]=svd(H1);             
    eig1=diag(S1);                   

    H2=Hk(ppg2);                 
    [U2,S2,V2]=svd(H2);          
    eig2=diag(S2);                    

    H_x=Hk(accx);               
    [U_x,S_x,~]=svd(H_x);          
    eig_x=diag(S_x);                  

    H_y=Hk(accy);                
    [U_y,S_y,~]=svd(H_y);        
    eig_y=diag(S_y);               

    H_z=Hk(accz);               
    [U_z,S_z,~]=svd(H_z);       
    eig_z=diag(S_z);              


    ppg_aut_1=auto_grandi(eig1);
    ppg_aut_2=auto_grandi(eig2);
    acc_aut_x=auto_grandi(eig_x);
    acc_aut_y=auto_grandi(eig_y);
    acc_aut_z=auto_grandi(eig_z);

    %% Correlation matrices
    Matrix_Cor_1x=(U1')*(U_x);   % ppg1 e accx
    Matrix_Cor_1y=(U1')*(U_y);   % ppg1 e accy
    Matrix_Cor_1z=(U1')*(U_z);   % ppg1 e accx
    Matrix_Cor_2x=(U2')*(U_x);   % ppg2 e accx
    Matrix_Cor_2y=(U2')*(U_y);   % ppg2 e accy
    Matrix_Cor_2z=(U2')*(U_z);   % ppg2 e accz 


    Matrix_Cor_1x=Matrix_Cor_1x(1:ppg_aut_1,1:acc_aut_x);
    Matrix_Cor_1y=Matrix_Cor_1y(1:ppg_aut_1,1:acc_aut_y);
    Matrix_Cor_1z=Matrix_Cor_1z(1:ppg_aut_1,1:acc_aut_z);
    Matrix_Cor_2x=Matrix_Cor_2x(1:ppg_aut_2,1:acc_aut_x);
    Matrix_Cor_2y=Matrix_Cor_2y(1:ppg_aut_2,1:acc_aut_y);
    Matrix_Cor_2z=Matrix_Cor_2z(1:ppg_aut_2,1:acc_aut_z);

    Matrix_Cor_1x=(Matrix_Cor_1x.^2)';
    Matrix_Cor_1y=(Matrix_Cor_1y.^2)';
    Matrix_Cor_1z=(Matrix_Cor_1z.^2)';
    Matrix_Cor_2x=(Matrix_Cor_2x.^2)';
    Matrix_Cor_2y=(Matrix_Cor_2y.^2)';
    Matrix_Cor_2z=(Matrix_Cor_2z.^2)';

    SUM_1x=max(Matrix_Cor_1x);
    SUM_1y=max(Matrix_Cor_1y);
    SUM_1z=max(Matrix_Cor_1z);
    SUM_2x=max(Matrix_Cor_2x);
    SUM_2y=max(Matrix_Cor_2y);
    SUM_2z=max(Matrix_Cor_2z); 
    
    SUM_1=max([SUM_1x; SUM_1y; SUM_1z]);
    SUM_2=max([SUM_2x; SUM_2y; SUM_2z]);

    %% Reconstruction
    
    SOGLIA=0.6;
    Sr_1=zeros(size(S1)); % matrix 701x300
    Sr_2=zeros(size(S2)); % matrix 701x300

    % ------------------ PPG1 --------------------
    cont=0;
    for i=1:length(SUM_1)
        if SUM_1(i)<SOGLIA     
            Sr_1(i,i)=eig1(i);
            cont=cont+1;
        end
    end

    if cont==0 
        for i=1:length(SUM_1)
           Sr_1(i,i)=eig1(i);
        end 
    end

    mH_1=U1*Sr_1*V1';

    for i=1:ppg_aut_1
           temp(i,:)=mH_1(:,i)/sqrt(eig1(i)); 
    end

    if size(temp,1)>1
        xr_1=sum(temp);
    else
        xr_1=temp;
    end

   % ------------------ PPG2 --------------------
    cont=0; 
    for j=1:length(SUM_2)
        if SUM_2(j)<SOGLIA
            Sr_2(j,j)=eig2(j);
            cont=cont+1;
        end
    end 
    if cont==0 
        for j=1:length(SUM_2)
           Sr_2(j,j)=eig2(j);
        end 
    end

    mH_2=U2*Sr_2*V2';

    for j=1:ppg_aut_2
           temp2(j,:)=mH_2(:,j)/sqrt(eig2(j));
    end
    if size(temp2,1)>1
        xr_2=sum(temp2);
    else
        xr_2=temp2;
    end
end

