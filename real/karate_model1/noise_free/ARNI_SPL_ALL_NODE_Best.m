clc;
clear;
close all;
MODEL={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
BASIS={'polynomial_diff','fourier_diff','power_series','RBF'};
addpath('Functions/')

%BASIS={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};
th=0.001;
disp('Initiating reconstruction...');
disp('Reading data...');

data = dlmread('Data/data.dat');
connectivity=dlmread('Data/connectivity.dat');
ts_param=dlmread('Data/ts_param.dat');
data=data';
para=load('Result_karate_model_1_undirected.txt');

tic;
adjacency=zeros(size(connectivity,1),size(connectivity,2));

for NODE=1:size(connectivity,1)
lambda=0.001;
v=ones(1,ts_param(1,1)*ts_param(1,2));
v1=zeros(1,ts_param(1,1)*ts_param(1,2));
Loss=zeros(1,ts_param(1,1)*ts_param(1,2));
result_auc=0;
t=0;
% b=1;
while (t<100) %&& (b==1)
    %%% Hard / Linear
    for i=1:length(v)
        if Loss(i)<lambda
            v(i)=1-Loss(i)/lambda;
        else
            v(i)=0;
        end
    end
%     %%% Logistic
%     for i=1:length(v)
%         v(i)=(1+exp(-lambda))/(1+exp(Loss(i)-lambda));
%     end
    if isequal(v1,v)||(result_auc==1)
%         b=0;
        break
    else
        [list,cost,Loss,vec,AUC]=reconstruct_spl_ALL_NODE_best(MODEL{1},NODE,BASIS{para(NODE,3)},para(NODE,2),connectivity,ts_param,data,th,v);
        t=t+1;
        %lambda=lambda+mean(Loss);
        lambda=1.25*lambda;
        %lambda=lambda*1.5;
        v1=v;
        if AUC > result_auc
            result_auc=AUC;
            result_vec=vec;
        end
    end
end
adjacency(NODE,:)=result_vec;
end
toc;

connectivity(connectivity~=0)=1;
[FPR,TPR,~,AUC]=perfcurve(reshape(connectivity,size(connectivity,1)*size(connectivity,2),1),reshape(adjacency,size(connectivity,1)*size(connectivity,2),1),1);
disp(AUC);

figure('Name','Reconstruction for all unit ');
plot(FPR,TPR,'LineWidth',2.5,'Color',[0,0.7,0.9])
title('Receiver-Operating-Characteristic Curve');
xlabel('False Positives Rate');
ylabel('True Positives Rate');
text(0.4,0.5,['AUC score=',num2str(AUC)])