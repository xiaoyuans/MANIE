clear;
clc;
format short;
Ymatrix=dlmread('Ymatrix.txt');
%%% dlmread从由filename指定的具有分割标志的ASCII文件中读取数值数据，使用一个分隔符作为格式化的文件的暗示信息。默认分隔符使用逗号(,)。
%%% dlmread(‘filename’, delimiter) 指定了分隔符。使用\t代表制表符tab分隔。
Ymatrix = Ymatrix';
Xmatrix=load('Xmatrix.txt');
Amatrix=dlmread('Amatrix.txt');
%%
test_total=1;% 重构十次
eq_total=1;    %  # of attempts
[~,num_nodes]=size(Ymatrix);
for j=2:num_nodes
    for i=j:-1:2
        Xmatrix(i,j)=Xmatrix(i-1,j);%%%把邻接矩阵的对角线上方的值都往下挪一个，占用对角线的0的位置
    end
end
Xmatrix=Xmatrix(2:end,:);%%%忽略对角线上的自循环0值，即各点的邻接向量中省略掉其自身
%RS = cell(eq_total);
for i=1:eq_total
    RSX{i}=zeros(num_nodes-1,num_nodes,test_total);%%%用来存储十次的重构结果
end
data_x=1.0/eq_total:1.0/eq_total:1.0;
data_x=data_x*0.35;
% data_x=0.5:0.1/(eq_total-1):0.6;

detect_total=num_nodes;

testX = zeros(num_nodes-1,num_nodes);
avg_edge=zeros(eq_total,test_total);
avg_no_edge=zeros(eq_total,test_total);
data_total=num_nodes;
tic
t=0;
Loss=zeros(size(Ymatrix,1),num_nodes);
loss=zeros(1,size(Ymatrix,1));
v=ones(size(Ymatrix,1),1);
v1=zeros(size(Ymatrix,1),1);
lambda=1;
while(t<20)
    for i=1:length(v)
        if loss(i)<lambda
            v(i)=1-loss(i)/lambda;
        else
            v(i)=0;
        end
    end
    %%% Logistic
%     for i=1:length(v)
%         v(i)=(1+exp(-lambda))/(1+exp(loss(i)-lambda));
%     end
    if v1==v
        %b=0;
        break
    else
        for node=1:num_nodes
            % node
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  revised
            A=Amatrix((node-1)*data_total+1:(node)*data_total,:);
            A(:,node)=[];
            dif=Ymatrix(:,node);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% revised
            [temp,len]=size(A);%%% A: 34*33
            norm=zeros(1,len);
            for i=1:len
                norm(i)=sqrt(sum(A(:,i).*A(:,i)));
                A(:,i)=A(:,i)/norm(i);%%%归一化
            end
            AMA=A;
            B=dif;
            NOR=norm;
            for x_no=1:eq_total
                data_needed=data_x(x_no);
%                 K=len;
                K =round(data_needed*len);%%% round 四舍五入
                for test_no=1:test_total%%% 执行十次重构过程
                    j=1;
                    % K equations to get all parameters
                    selected_no=[];
                    i=0;
                    while i<K %%%选出K组不同的数据来实现重构，并把选择的组别放入selected_no里
                        test_selected_no=ceil(len*rand);%%% ceil 向上取整
                        if isempty(find(selected_no==test_selected_no, 1))
                            i=i+1;
                            selected_no(i)=test_selected_no;
                        end
                    end
                    A=AMA(selected_no,:);
                    y=B(selected_no,j);
                    x0=A'*y;
                    xp = L1_SPL_1(x0, A, [], y, v(selected_no,j));
                    xp =xp./NOR';             % xp is what we need, you can plot it to check the results
                    xp(find(isnan(xp)==1)) = 0;
                    Loss(:,node)=abs(B-AMA*xp);
                    RSX{x_no}(:,node,test_no)=xp;
                end
            end
        end
        loss=mean(Loss,2);
        t=t+1;
%         lambda=lambda+mean(loss);
        lambda=1.25*lambda;
%         lambda=lambda*1.5;
        v1=v;
    end
end
toc
%% evaluate the success rate
[sre,srne]=estsrall(Xmatrix,mean(RSX{1},3),0);% mean(RSX{1},3) 对矩阵的第三维取平均
W=Xmatrix;
Wlearn=mean(RSX{1},3);
W1=[];
W2=[];
for i = 1:size(W,1)
    for j = 1:size(W,2)
        W1=[W1;W(i,j)];
        W2=[W2;Wlearn(i,j)];
    end
end
[FPR,TPR,~,AUC]=perfcurve(W1,W2,1);
[sre,srne,AUC]