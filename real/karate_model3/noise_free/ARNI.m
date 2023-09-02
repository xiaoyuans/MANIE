clear;
MODEL={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
BASIS={'polynomial_diff','fourier_diff','power_series','RBF'};
addpath('Functions/')

% BASIS={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};
ORDER=5;
NODE=1;
th=0.001;
disp('Initiating reconstruction...');
disp('Reading data...');
data = dlmread('Data/data.dat');
connectivity=dlmread('Data/connectivity.dat');
ts_param=dlmread('Data/ts_param.dat');
data=data';
tic;
[list,cost,FPR,TPR,AUC]=reconstruct(MODEL{3},NODE,BASIS{3},ORDER,connectivity,ts_param,data,th);
toc;
disp(AUC);
figure('Name',['Reconstruction for unit ',num2str(NODE)]);
ax1=subplot(2,1,1);
plot(cost,'-o','LineWidth',2.5,'Color',[0,0.7,0.9],'MarkerFaceColor',[0,0.7,0.9]);
title({'Evolution of Fitting Costs';'(Actual connections shown above)'});
xlabel('# Inferred Interactions');
ylabel('Cost');

for i=1:length(cost)
    text(1.01*i,1.03*cost(i),num2str(list(i)))
end

ax2=subplot(2,1,2);
plot(FPR,TPR,'LineWidth',2.5,'Color',[0,0.7,0.9])
title('Receiver-Operating-Characteristic Curve');
xlabel('False Positives Rate');
ylabel('True Positives Rate');
text(0.4,0.5,['AUC score=',num2str(AUC)])
