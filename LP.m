function F = LP(interaction,SD,SM,nl,nd,alpha)


PL=interaction;			% interaction_ori 原始关联矩阵
PD=interaction';
P0=interaction;
PZ =interaction;

% PL=re_md_l;			% interaction_ori 原始关联矩阵
% PD=re_md_d';
% P0=re_md_l;
% PZ = re_md_d;

delta = 1;


LPM = zeros(nl,nd);
LPD = zeros(nl,nd);
num_l = 1;

while  (delta > 0.01 && num_l<1000)
   
    PL1 = alpha * SM * P0 + (1-alpha) * PL;
    LPM = LPM + PL1;
    delta =abs(sum(sum((abs(PL1)-abs(PL)))));
    PL = PL1;
    num_l = num_l+1;
    
end


delta = 1;

num_p = 1;

while  (delta > 0.01 && num_p<1000)
    
    PD1 = alpha * (SD * PD)' + (1-alpha) * PZ;
    LPD = LPD + PD1;
    delta =abs(sum(sum((abs(PD1)-abs(PZ)))));

    PZ =PD1;
    num_p= num_p+1;
    
end


% 归一化
for i=1:nd

     MAX_D = max(LPD(:,i));
    MIN_D = min(LPD(:,i));

        LPD(:,i) = (LPD(:,i) - MIN_D) / (MAX_D - MIN_D); 
    

end
 for i =1:nl
     MAX_M = max(LPM(i,:));
     MIN_M = min(LPM(i,:));
 
     LPM(i,:) = (LPM(i,:) - MIN_M) / (MAX_M - MIN_M);
  
 end

beta = 1/2;
F = beta * LPM + (1-beta) * LPD;