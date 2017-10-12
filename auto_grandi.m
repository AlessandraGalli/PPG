function [num_aut]=auto_grandi(eigen)

eig1=eigen(1:10);       
dif=-diff(eig1);        
media=mean(dif);       
index=find(dif>media); 
ind=index(end);        
p=dif(ind);         
if p>2*(eigen(1)-eigen(2))
    num_aut=ind;
else
    num_aut=10;
end 
num_aut=max(4,num_aut);
end
    