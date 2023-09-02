load('BA.txt');
load('X2.txt');
W=BA;
Wlearn=X2;
W1=[];
W2=[];
for i = 1:length(W)
    for j = 1:length(W)
        W1=[W1;W(i,j)];
        W2=[W2;Wlearn(i,j)];
    end
end
[FPR,TPR,~,AUC]=perfcurve(W1,W2,1);

figure();
plot(FPR,TPR,'LineWidth',2.5,'Color',[0,0.7,0.9])
title('Receiver-Operating-Characteristic Curve');
xlabel('False Positives Rate');
ylabel('True Positives Rate');
text(0.4,0.5,['AUC score=',num2str(AUC)])

[sre,srne]=estsrall(W,Wlearn,0);% mean(RSX{1},3) 对矩阵的第三维取平均
[sre,srne,AUC]