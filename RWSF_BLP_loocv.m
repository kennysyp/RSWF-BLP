clc;clear;
% format long
%% 原始处理
lncSim = load('D:\MATLAB\Bidirectional_label_propagation\data\lncRNAsimilarity.txt');   %lncRNA 表达相似性
disSim = load('D:\MATLAB\Bidirectional_label_propagation\data\diseasesimilarity.txt');  % disease 语义相似性
interaction = load('D:\MATLAB\Bidirectional_label_propagation\data\known_lncRNA_disease_interaction.txt');

%% 数据处理方法

[nl,nd] = size(interaction);
interaction_ori = interaction; 
  n=1;
  m=1;

 [km,kd] = gaussiansimilarity(interaction_ori,nl,nd);
 [M_re,D_re] = K_fusion(disSim,lncSim,kd,km,7,0.4);

%% 主方法

F_ori = LP(interaction_ori,D_re,M_re,nl,nd,0.2);

index=find(interaction_ori==1); 
%%  留一交叉验证（LOOCV）
for u=1:length(index)
     u   
    interaction(index(u))=0;
   

%% 数据处理方法
%% computing Gaussian interaction profile kernel of lncRNAs

[km,kd] = gaussiansimilarity(interaction,nl,nd); 


[M_re,D_re] = K_fusion(disSim,lncSim,kd,km,7,0.4);

%% 主方法
  
F = LP(interaction,D_re,M_re,nl,nd,0.2);

%%

F_ori(index(u))=F(index(u));
interaction = interaction_ori;
    
end
 pre_label_score = F_ori(:);
 label_y = interaction_ori(:);
 auc=roc_1(pre_label_score,label_y,'red');

