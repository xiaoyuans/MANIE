function [list,cost,vec]=reconstruct_ALL_NODE(MODEL,NODE,BASIS,ORDER,ts_param,x,th)

models={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
bases={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};

if any(ismember(models,MODEL))==0
    disp('ERROR: MODEL must be a valid string: kuramoto1, kuramoto2, michaelis_menten, roessler');
elseif any(ismember(bases,BASIS))==0
    disp('ERROR: BASIS must be a valid string: polynomial, polynomial_diff, fouries, fourier_diff, power_series, RBF');
    
else
    S=ts_param(1,1);
    M=ts_param(1,2);
    [N,~]=size(x);    
    % Estimating time derivatives and constructing input matrices
%     disp('Estimating time derivatives and constructing input matrices...');
    Xtemp=[];
    DX=[];
    for s=1:S
        Ytemp=zeros(N,M);
        DY=zeros(N,M);
        for n=1:N
            for m=1:M
                Ytemp(n,m)=(x(n,m+(s-1)*(M+1))+x(n,m+1+(s-1)*(M+1)))*0.5;
                DY(n,m)=(-x(n,m+(s-1)*(M+1))+x(n,m+1+(s-1)*(M+1)))*1/(1.0);
            end
        end
        Xtemp=[Xtemp,Ytemp];
        DX=[DX,DY];
    end
    
    switch MODEL
        
        case 'roessler'            
            X=Xtemp;           
            % Beginning of reconstruction algorithm
%             disp('Performing ARNI...');
            Y=basis_expansion(X,ORDER,BASIS,NODE);
            nolist=1:N;
            list=[];
            cost=[];
            b=1;
            vec=zeros(1,N);
            while (~isempty(nolist)) && (b==1)
                % Composition of inferred subspaces
                Z=[];
                for n=1:length(list)
                    Z=[Z;Y(:,:,list(n))];
                end
                
                % Projection on remaining composite spaces
                P=zeros(length(nolist),2);
                cost_err=zeros(length(nolist),1);
                for n=1:length(nolist)
                    % Composition of a possible space
                    R=[Z;Y(:,:,nolist(n))];
                    % Error of projection on possible composite space
                    P(n,1)=std(DX(3*(NODE-1)+1,:)-DX(3*(NODE-1)+1,:)*pinv(R)*R);%%% pinv求其伪逆矩阵
                    P(n,2)=nolist(n);
                    % Fitting cost of possible composite space
                    cost_err(n,1)=1/M *norm(DX(3*(NODE-1)+1,:)-DX(3*(NODE-1)+1,:)*pinv(R)*R);
                    R=[];
                end
                
                if std(P(:,1))<th
                    b=0;
                    break
                    
                else
                    % Selection of composite space which minimizes
                    % projection error
                    [MIN,block]=min(P(:,1));
                    list=[list,P(block,2)];
                    nolist(nolist==P(block,2))=[];
                    vec(1,P(block,2))=MIN;
                    cost=[cost,cost_err(block,1)];
                end
                
            end
            % End of reconstruction algorithm                
        otherwise
            
            if (strcmp(MODEL,'kuramoto1')) || (strcmp(MODEL,'kuramoto2')) %%%判断两个字符串是否相同
                % Transforming data coming from phase oscillators
                X=mod(Xtemp,2*pi);%%% 取余数
            else
                X=Xtemp;
            end
            
            % Beginning of reconstruction algorithm
%             disp('Performing ARNI...');
            Y=basis_expansion(X,ORDER,BASIS,NODE);
            nolist=1:N;
            list=[];
            cost=[];
            b=1;
            vec=zeros(1,N);
            while (~isempty(nolist)) && (b==1)
                % Composition of inferred subspaces
                Z=[];
                for n=1:length(list)
                    Z=[Z;Y(:,:,list(n))];
                end
                
                % Projection on remaining composite spaces
                P=zeros(length(nolist),2);
                cost_err=zeros(length(nolist),1);
                for n=1:length(nolist)
                    % Composition of a possible space
                    R=[Z;Y(:,:,nolist(n))];
                    % Error of projection on possible composite space
                    P(n,1)=std(DX(NODE,:)-DX(NODE,:)*pinv(R)*R);
                    P(n,2)=nolist(n);
                    % Fitting cost of possible composite space
                    cost_err(n,1)=1/M *norm(DX(NODE,:)-DX(NODE,:)*pinv(R)*R);
                    R=[];
                end
                
                if std(P(:,1))<th
                    b=0;
                    break
                else
                    % Selection of composite space which minimizes
                    % projection error
                    [MIN,block]=min(P(:,1));
                    list=[list,P(block,2)];
                    nolist(nolist==P(block,2))=[];
                    vec(1,P(block,2))=MIN;
                    cost=[cost,cost_err(block,1)];
                end
                
            end
            % End of reconstruction algorithm
            
    end
end