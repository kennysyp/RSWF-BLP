function [K_lnc,K_dis] = K_fusion(disSim,lncSim,kd,km,k,beta)
[~,n] = size(lncSim);
[~,m] = size(disSim);
t = 1;
I_lnc = eye(n,n);
P_lncSim = previous_l(lncSim);
P_km = previous_l(km);
%计算P_km

 P_km_0 = P_km;
 P_lncSim_0 = P_lncSim;
 P_km_p = P_km;
 P_km_ini = P_km;
 P_lncSim_ini = P_lncSim;



% 初始化 L_lncSim、L_km
L_lncSim = zeros(n,n);
L_km = zeros(n,n);

L_lncSim = knn(P_lncSim,k);
L_km = knn(P_km,k);
L_lncSim_0 = L_lncSim;
L_km_0 = L_km;


while  t <= 2  
    P_lncSim_0 = 0.5 * L_lncSim* P_lncSim + 0.5 * P_lncSim_ini;
    P_km_0 = 0.5 * L_km* P_km + 0.5 * P_km_ini;
    P_lncSim = beta*L_lncSim_0 * P_km_0  * L_lncSim_0' + (1-beta)*P_km_ini+I_lnc;
    P_km = beta*L_km_0 * P_lncSim_0 * L_km_0' + (1-beta)*P_lncSim_ini+I_lnc;
 

    L_lncSim = knn(P_lncSim,k);
    L_km = knn(P_km,k);
    t=t+1;
end


 K_lnc = (P_lncSim + P_km);
 K_lnc = K_lnc ./ (repmat(sum(K_lnc,2),1,n));

K_lnc = (K_lnc + K_lnc'+I_lnc ) / 2;
%%
I_dis = eye(m,m);

% 计算 P_disSim  P_kd
P_disSim = previous_l(disSim);
P_kd = previous_l(kd);

 P_disSim_0 = P_disSim;
 P_kd_0 = P_kd;
 P_disSim_ini = P_disSim;
 P_kd_ini = P_kd;

t = 1;

L_disSim = knn(P_disSim,k);
L_kd = knn(P_kd,k);
L_disSim_0 = L_disSim;
L_kd_0 = L_kd;




while  t <= 2 
   P_disSim_0 = 0.5 * L_disSim * P_disSim + 0.5 * P_disSim_ini;
   P_kd_0 = 0.5 * L_kd * P_kd + 0.5 * P_kd_ini;
   P_kd = beta*L_kd_0 * P_disSim_0 * L_kd_0' + (1-beta)*P_disSim_ini + I_dis;         
   P_disSim =beta*L_disSim_0 * P_kd_0  * L_disSim_0' + (1-beta) * P_kd_ini + I_dis; 
   
    L_disSim = knn(P_disSim,k);
    L_kd = knn(P_kd,k); 

    t=t+1;
end

   K_dis = (P_disSim + P_kd); 
   K_dis = K_dis ./ repmat(sum(K_dis,2),1,m);

   K_dis = (K_dis + K_dis'+I_dis) / 2;
 
end






    