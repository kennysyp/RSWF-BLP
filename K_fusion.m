function [K_lnc,K_dis] = K_fusion(disSim,lncSim,kd,km,k,beta)
[~,n] = size(lncSim);
[~,m] = size(disSim);

 
 delta_lnc = 5.5;
 delta_rna = 1;
 delta_km = 5.5;
 
 delta_rna_km = 1;
 
 delta_dis = 8;
 
 delta_d = 1;
 delta_kd = 8;
 
 delta_d_kd = 1;
 delta_dis_p = 0;
 delta_lnc_p = 0;
 delta_km_p = 0;
 delta_kd_p = 0;
t = 1;
I_lnc = eye(n,n);

 P_lncSim = previous_l(lncSim);
 P_km = previous_l(km);
%计算P_km
% P_km = km ./ (repmat(sum(km,2),1,n));
% P_lncSim = (P_lncSim + P_lncSim') / 2;
%  P_km = (P_km + P_km') / 2;
 P_km_0 = P_km;
  P_lncSim_0 = P_lncSim;
   P_km_p = P_km;
  P_lncSim_p = P_lncSim;
  P_km_ini = P_km;
  P_lncSim_ini = P_lncSim;

%  [S_lnc,U_lnc] = matrix_dec(lncSim,40,10);
%  [S_km,U_km] = matrix_dec(km,40,10);
%  [n,s] = size(U_lnc);
%  Y_lnc = U_lnc ./ repmat(sum(U_lnc,2),1,s);
%  Y_km = U_km ./ repmat(sum(U_km,2),1,s);
% W_lnc = Y_lnc*Y_lnc';
% W_km = Y_km*Y_km';
% delta_lnc_p = 0;%TR(P_lncSim);
% delta_km_p = 0;%TR(P_km);
% 
% U_lnc_ini = U_lnc;
% U_km_ini = U_km;
% U_lnc_0 = U_lnc;
% U_km_0 = U_km;
%S_lnc_0 = S_lnc;
% S_km_0 = S_km;
 %E_lnc = eye(n,s);

% 初始化 L_lncSim、L_km
L_lncSim = zeros(n,n);
L_km = zeros(n,n);

L_lncSim = knn(P_lncSim,k);
L_km = knn(P_km,k);
L_lncSim_0 = L_lncSim;
L_km_0 = L_km;

% end
% delta_l_p = 0;
% delta_p = 0;
while  t <= 2 %&& (delta_lnc > 8 || delta_km > 8)  %((delta_lnc > 1e-3 || delta_km > 1e-3) && t <= 100) 
     P_lncSim_0 = 0.5 * L_lncSim* P_lncSim + 0.5 * P_lncSim_ini;
     P_km_0 = 0.5 * L_km* P_km + 0.5 * P_km_ini;
    P_lncSim = beta*L_lncSim_0 * P_km_0  * L_lncSim_0' + (1-beta)*P_km_ini+I_lnc;
    P_km = beta*L_km_0 * P_lncSim_0 * L_km_0' + (1-beta)*P_lncSim_ini+I_lnc;
%      U_lnc_0 = 0.5 * W_lnc * U_lnc_0 + 0.5 * U_lnc_ini;
%      U_km_0 = 0.5 * W_km * U_km_0 + 0.5 * U_km_ini;
%     U_lnc = Y_lnc * U_km_0' * Y_lnc+ E_lnc;
%     U_km = Y_km * U_lnc_0' * Y_km+ E_lnc;  
%     
%     P_lncSim_1 = P_lncSim ./ repmat(sum(P_lncSim,2),1,n);
%     P_lncSim = (P_lncSim_1 + P_lncSim_1') / 2;
%     
%      P_km_1 = P_km ./ repmat(sum(P_km,2),1,n);
%      P_km = (P_km_1 + P_km_1') / 2;
%      Y_lnc = U_lnc ./ repmat(sum(U_lnc,2),1,s);
%      Y_km = U_km ./ repmat(sum(U_km,2),1,s);
%     W_lnc = U_lnc * U_lnc';
%     %S_lnc = zzh(W_lnc);
%     W_km = U_km * U_km';
     %S_km = zzh(W_km);
    
    %lncSim相似性矩阵的迭代次数判断
    delta_lnc = sqrt(sum(sum((P_lncSim-P_lncSim_p).^2)));
%     P_lncSim = previous_l(P_lncSim);
%     P_km = previous_l(P_km);
    %delta_lnc = TR(P_lncSim-P_lncSim_0);
%       delta_rna = delta_lnc - delta_lnc_p;
%       delta_lnc_p = delta_lnc;
    
   delta_km = sqrt(sum(sum((P_km-P_km_p).^2)));
    %delta_km = sqrt(norm(P_km-P_km_0));
    %km矩阵的迭代次数判断
     %delta_km = TR(P_km - P_km_0);
%      delta_rna_km = delta_km - delta_km_p;
%       delta_km_p = delta_km;
%    
      %U_lnc_0 = U_lnc;
      
     % P_lncSim_0 = (P_lncSim_0 + P_lncSim_0') / 2;
%        P_lncSim_0 = P_lncSim;
%        P_km_0 = P_km;
      %U_km_0 = U_km;
%      S_lnc_0 = S_lnc;
%      S_km_0 = S_km;
        
%        P_km_0 = (P_km_0 + P_km_0') / 2; 
        L_lncSim = knn(P_lncSim,k);
        L_km = knn(P_km,k);

       
   % P_km_0 = (P_km_0 + P_km_0') / 2;
%     P_lncSim_p = P_lncSim;
%      P_km_p = P_km;
     t=t+1;
end

%  W_lnc = U_lnc * U_lnc';
% S_lnc = zzh(W_lnc);
%  W_km = U_km * U_km';
%  S_km = zzh(W_km);
%  K_lnc = (S_lnc + S_km) / 2;
 K_lnc = (P_lncSim + P_km);% ./ repmat(sum(P_lncSim + P_km,2),1,n);
 K_lnc = K_lnc ./ (repmat(sum(K_lnc,2),1,n));
% K_lnc = (K_lnc + K_lnc') / 2;
%K_lnc = GuiYi(K_lnc); 
K_lnc = (K_lnc + K_lnc'+I_lnc ) / 2;
%%
I_dis = eye(m,m);

% 计算 P_disSim  P_kd
P_disSim = previous_l(disSim);
P_kd = previous_l(kd);
% P_kd = kd ./ (repmat(sum(kd,2),1,m));
% P_disSim = (P_disSim + P_disSim') / 2;
%  P_kd = (P_kd + P_kd') / 2;
% 

 P_disSim_0 = P_disSim;
 P_kd_0 = P_kd;
 P_disSim_ini = P_disSim;
 P_kd_ini = P_kd;
% 
% L_disSim = zeros(m,m);
% L_kd = zeros(m,m);
t = 1;
%  [S_dis,U_dis] = matrix_dec(disSim,40,10);
%  [S_kd,U_kd] = matrix_dec(kd,40,10);
%  [n,s] = size(U_dis);
% Y_dis = U_dis ./ repmat(sum(U_dis,2),1,s);
% Y_kd = U_kd ./ repmat(sum(U_kd,2),1,s);
% W_dis = Y_dis*Y_dis';
% W_kd = Y_kd*Y_kd';
% 
% U_dis_ini = U_dis;
% U_kd_ini = U_kd;
% U_dis_0 = U_dis;
% U_kd_0 = U_kd;
% 
  %E_dis = eye(n,s);
L_disSim = knn(P_disSim,k);
L_kd = knn(P_kd,k);
L_disSim_0 = L_disSim;
L_kd_0 = L_kd;
% A = sum(P_disSim,2);
% B = sum(P_kd,2);



while  t <= 2 %&& (delta_dis > 7 || delta_kd > 8) %( (delta_dis > 1e-3 || delta_kd > 1e-3) && t <= 100  )
   P_disSim_0 = 0.5 * L_disSim * P_disSim + 0.5 * P_disSim_ini;
   P_kd_0 = 0.5 * L_kd * P_kd + 0.5 * P_kd_ini;
   P_kd = beta*L_kd_0 * P_disSim_0 * L_kd_0' + (1-beta)*P_disSim_ini + I_dis;         
   P_disSim =beta*L_disSim_0 * P_kd_0  * L_disSim_0' + (1-beta) * P_kd_ini + I_dis; 

    
%       U_dis_0 = 0.5 * W_dis * U_dis_0 + 0.5 * U_dis_ini;
%       U_kd_0 = 0.5 * W_kd * U_kd_0 + 0.5 * U_kd_ini; 
%         U_dis = Y_dis * U_kd_0' *  Y_dis + E_dis;
%         U_kd = Y_kd * U_dis_0' * Y_kd+ E_dis;
%     
%     P_disSim_1 = P_disSim ./ repmat(sum(P_disSim,2),1,m);
%     P_disSim = (P_disSim_1 + P_disSim_1') / 2;
    
%     P_kd_1 = P_kd ./ repmat(sum(P_kd,2),1,m);
%     P_kd = (P_kd_1 + P_kd_1') / 2;
%       Y_dis = U_dis ./ repmat(sum(U_dis,2),1,s);
%       Y_kd = U_kd ./ repmat(sum(U_kd,2),1,s);
%       W_dis = U_dis * U_dis';
%       W_kd = U_kd * U_kd';
   %disSim相似性矩阵的迭代次数判断
      delta_dis = sqrt(sum(sum((P_disSim-P_disSim_0).^2)));
      %delta_kd = sqrt(sum(sum((U_kd-U_kd_0).^2)));
%        delta_dis = TR(P_disSim - P_disSim_0);
        delta_d = delta_dis - delta_dis_p;
%        delta_dis_p = delta_dis;
%        P_disSim_0 = P_disSim;
   
%    
     delta_kd = sqrt(sum(sum((P_kd-P_kd_0).^2)));
   %kd矩阵的迭代次数判断
%       delta_kd = TR(P_kd - P_kd_0);
        delta_d_kd = delta_kd - delta_kd_p;
%        delta_kd_p = delta_kd;
%        P_kd_0 = P_kd;
% U_dis_0 = U_dis;
% U_kd_0 = U_kd; 
   
   
    
   
   
%     P_disSim_0 = previous_l(P_disSim);
%     P_disSim_0 = (P_disSim_0 + P_disSim_0') / 2;
%     
%     P_kd_0 = previous_l(P_kd);
%     P_kd_0 = (P_kd_0 + P_kd_0') / 2;
    
     L_disSim = knn(P_disSim,k);
     L_kd = knn(P_kd,k); 
%      P_disSim_0 =  P_disSim;
%      P_kd_0 = P_kd;
    t=t+1;
end
%    W_dis = U_dis * U_dis';
%   S_dis = zzh(W_dis);
%    W_kd = U_kd * U_kd';
%   S_kd = zzh(W_kd);
%    K_dis = (S_dis + S_kd) / 2;
   K_dis = (P_disSim + P_kd); %./ repmat(sum(P_disSim + P_kd,2),1,m);
   K_dis = K_dis ./ repmat(sum(K_dis,2),1,m);
 %K_dis = (K_dis + K_dis') / 2;
  %K_dis = K_dis ./ (repmat(sum(K_dis,2),1,m));
 K_dis = (K_dis + K_dis'+I_dis) / 2;
    function delta = TR(Sim)
        [~,S,~] = svd(Sim);
        delta = sum(sum(S));
    end
end






    