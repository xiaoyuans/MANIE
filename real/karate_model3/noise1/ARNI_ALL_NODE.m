clc;
clear; 
close all;
addpath('Functions/')

MODEL={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
BASIS={'polynomial_diff','fourier_diff','power_series','RBF'};
%BASIS={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};
th=0.001;
disp('Initiating reconstruction...');
disp('Reading data...');

data = dlmread('Data/data.dat');
connectivity=dlmread('Data/connectivity.dat');
ts_param=dlmread('Data/ts_param.dat');
data=data';
para=load('Result_karate_model_3_undirected_uncircle.txt');

tic;
adjacency=zeros(size(connectivity,1),size(connectivity,2));

%%%%%%%%%%
for NODE=1:size(connectivity,1)
[list,cost,vec]=reconstruct_ALL_NODE(MODEL{3},NODE,BASIS{para(NODE,3)},para(NODE,2),ts_param,data,th);
adjacency(NODE,:)=vec;
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