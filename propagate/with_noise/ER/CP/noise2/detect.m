clear;
clc;
format short;
Ymatrix=dlmread('Ymatrix.txt');
Ymatrix = Ymatrix';
Xmatrix=load('Xmatrix.txt');
Amatrix=dlmread('Amatrix.txt');
test_total=1;
eq_total=1;
[~,num_nodes]=size(Ymatrix);
for j=2:num_nodes
    for i=j:-1:2
        Xmatrix(i,j)=Xmatrix(i-1,j);
    end
end
Xmatrix=Xmatrix(2:end,:);
for i=1:eq_total
    RSX{i}=zeros(num_nodes-1,num_nodes,test_total);
end
data_x=1.0/eq_total:1.0/eq_total:1.0;
data_x=data_x*0.35;
%data_x=0.5:0.1/(eq_total-1):0.6;
detect_total=num_nodes;
testX = zeros(num_nodes-1,num_nodes);
avg_edge=zeros(eq_total,test_total);
avg_no_edge=zeros(eq_total,test_total);
data_total=num_nodes;
tic
for node=1:num_nodes
    A=Amatrix((node-1)*data_total+1:(node)*data_total,:);
    A(:,node)=[];
    dif=Ymatrix(:,node);
    [temp,len]=size(A);
    norm=zeros(1,len);
    for i=1:len
        norm(i)=sqrt(sum(A(:,i).*A(:,i)));
        A(:,i)=A(:,i)/norm(i);
    end
    AMA=A;
    B=dif;
    NOR=norm;
    for x_no=1:eq_total
        data_needed=data_x(x_no);
        K=len;
        % K =round(data_needed*len);
        for test_no=1:test_total
            j=1;
            selected_no=[];
            loss=zeros(1,K);
            i=0;
            while i<K
                test_selected_no=ceil(len*rand);              
                if isempty(find(selected_no==test_selected_no, 1))
                    i=i+1;
                    selected_no(i)=test_selected_no;
                end
            end
            selected_no1=selected_no;
            v=zeros(1,length(selected_no));
            t=0;
            A1=AMA(selected_no,:);
            y1=B(selected_no,j);                       
            while (t<500)
                A=AMA(selected_no1,:);
                y=B(selected_no1,j);
                x0=A'*y;
                xp = L1(x0, A, [], y);
                for i=1:length(selected_no)
                    loss(i)=abs(y1(i)-A1(i,:)*xp);
                end
                lambda=mean(loss);
%                 lambda=0.75*max(loss);
                for i=1:length(selected_no1)
                    if loss(i)<lambda
                        v(i)=1;%-loss(i)/lambda;
                    else
                        v(i)=0;
                    end
                end
                selected_no1=selected_no(v==1);% selected_no(v==0)=[];
                t=t+1;
            end           
            xp =xp./NOR';
            xp(find(isnan(xp)==1)) = 0;
            RSX{x_no}(:,node,test_no)=xp;
        end
    end
end
toc
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