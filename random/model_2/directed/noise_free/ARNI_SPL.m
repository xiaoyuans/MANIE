clear;
MODEL={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
BASIS={'polynomial_diff','fourier_diff','power_series','RBF'};
addpath('Functions/')

% BASIS={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};
ORDER=4;
NODE=25;%%% 对节点NODE的连接进行重构
th=0.001;% Stopping criterium: decrease it to recover longer list of possible links
lambda=0.02;
disp('Initiating reconstruction...');
disp('Reading data...');

data = dlmread('Data/data.dat');
connectivity=dlmread('Data/connectivity.dat');
ts_param=dlmread('Data/ts_param.dat');
data=data';

v=ones(1,ts_param(1,1)*ts_param(1,2));
v1=zeros(1,ts_param(1,1)*ts_param(1,2));
Loss=zeros(1,ts_param(1,1)*ts_param(1,2));
t=0;
result_auc=0;
% b=1;
tic;
while (t<20) %&& (b==1)
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
    if v1==v
%         b=0;
        break
    else
        [list,cost,Loss,FPR,TPR,AUC]=reconstruct_spl(MODEL{1},NODE,BASIS{2},ORDER,connectivity,ts_param,data,th,v);
        t=t+1;
        lambda=lambda+mean(Loss);
        %lambda=1.25*lambda;
        %lambda=lambda*1.5;
        v1=v;
        if AUC > result_auc
            result_auc=AUC;
            result_list=list;
            result_cost=cost;
            result_FPR=FPR;
            result_TPR=TPR;
            result_epoch=t;
        end            
    end
end
toc;
disp(result_auc);
disp(result_epoch);
% link=[];
% for i=1:length(cost)
%     if cost(i)>0.4
%         link=[link,list(i)];
%     else
%         break
%     end
% end

figure('Name',['Reconstruction for unit ',num2str(NODE)]);
ax1=subplot(2,1,1);
plot(result_cost,'-o','LineWidth',2.5,'Color',[0,0.7,0.9],'MarkerFaceColor',[0,0.7,0.9]);
title({'Evolution of Fitting Costs';'(Actual connections shown above)'});
xlabel('# Inferred Interactions');
ylabel('Cost');

for i=1:length(result_cost)
    text(1.01*i,1.03*result_cost(i),num2str(result_list(i)))
end

ax2=subplot(2,1,2);
plot(result_FPR,result_TPR,'LineWidth',2.5,'Color',[0,0.7,0.9])
title('Receiver-Operating-Characteristic Curve');
xlabel('False Positives Rate');
ylabel('True Positives Rate');
text(0.4,0.5,['AUC score=',num2str(result_auc)])
