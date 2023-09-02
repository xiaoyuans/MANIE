clc;
clear; 
MODEL={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
BASIS={'polynomial_diff','fourier_diff','power_series','RBF'};
addpath('Functions/')

% BASIS={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};
th=0.001;

% disp('Initiating reconstruction...');
% disp('Reading data...');

data = dlmread('Data/data.dat');
connectivity=dlmread('Data/connectivity.dat');
ts_param=dlmread('Data/ts_param.dat');
data=data';

fidout=fopen('Result_model_1_25_4_directed.txt','w');

for NODE=1:25

Max_Auc=0;
Corr_Order=0;
Corr_Basic=0;

for ORDER=4:10
% tic;
for i =1:4
if Max_Auc==1 
    break;
else
[list,cost,FPR,TPR,AUC]=reconstruct(MODEL{1},NODE,BASIS{i},ORDER,connectivity,ts_param,data,th);

if Max_Auc<AUC
Max_Auc=AUC;
Corr_Order=ORDER;
Corr_Basic=i;
disp(Corr_Order);
disp(Corr_Basic);
disp(AUC);
figure('Name',['Reconstruction for unit ',num2str(NODE)]);
ax1=subplot(2,1,1);
plot(cost,'-o','LineWidth',2.5,'Color',[0,0.7,0.9],'MarkerFaceColor',[0,0.7,0.9]);
title({'Evolution of Fitting Costs';'(Actual connections shown above)'});
xlabel('# Inferred Interactions');
ylabel('Cost');

for j=1:length(cost)
    text(1.01*j,1.03*cost(j),num2str(list(j)))
end

ax2=subplot(2,1,2);
plot(FPR,TPR,'LineWidth',2.5,'Color',[0,0.7,0.9])
title('Receiver-Operating-Characteristic Curve');
xlabel('False Positives Rate');
ylabel('True Positives Rate');
text(0.4,0.5,['AUC score=',num2str(AUC)])
end
end
end
% toc;
end

fprintf(fidout,'%2d %2d %2d %8f\n',NODE,Corr_Order,Corr_Basic,Max_Auc);

end

fclose(fidout);