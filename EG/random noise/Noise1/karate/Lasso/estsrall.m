function [sre,srne,rs] = estsrall( Xmatrix,detrs,flag)
    [m,n] = size(Xmatrix);
    sre=0;
    srne=0;
    rs = zeros(m,n);
    for i=1:n
		if flag==0
	        rs(:,i)=estsr_single(Xmatrix(:,i),detrs(:,i),m);%%%设置截断方式，将产生的邻接01向量保存在rs中
		elseif flag==1		
			rs(:,i)=estsr_single_less(Xmatrix(:,i),detrs(:,i),m);
		else		
			rs(:,i)=estsr_single_cut(Xmatrix(:,i),detrs(:,i),m);
		end
    end
    for i=1:n
        if sum(Xmatrix(:,i))==0
			sre = sre +1;
			srne = srne +1;
		else
			sre=sre+sum(rs(:,i).*Xmatrix(:,i))/sum(Xmatrix(:,i));
			srne=srne+sum((1-rs(:,i)).*(1-Xmatrix(:,i)))/sum(1-Xmatrix(:,i));
        end
    end
    sre = sre/n;
    srne = srne/n;
    %sre=sum(sum(rs.*Xmatrix))/sum(sum(Xmatrix));
    %srne=sum(sum((1-rs).*(1-Xmatrix)))/sum(sum(1-Xmatrix));
    
end

function [xadj] = estsr_single(xm,dtrs,m)%%%截断方式，产生最终重构出来的01邻接向量
	if sum(xm)==0
		xadj=zeros(m,1);
	else
    a = xm.*dtrs;
    a = unique(a);%%% unique函数返回的是和a中一样的值，但是没有重复元素。产生的结果向量按升序排序。   
    [k,~] = size(a);
    flag=0;
    xadj = zeros(m,1);
    for i=2:k
        adj = zeros(m,1);
        index = find(dtrs>=a(i,1));
        adj(index)=1;
%         sum(adj)
        sre = sum(adj.*xm)/sum(xm);
        srne=sum((1-adj).*(1-xm))/sum(1-xm);
        if sre*srne>flag
            flag = sre*srne;
            xadj = adj;
        end
    end
	end
end

function [xadj] = estsr_single_cut(~,dtrs,m)
	xadj=ones(m,1);
	sep=10;
	sepxy=[];
	y=dtrs;
	ymin=min(y);
	ymax=max(y);
	x=linspace(ymin,ymax,sep);
	yy=hist(y,x);
	for sep1=1:sep
		sepxy=[sepxy,length( find(( yy(max(sep1-2,1):min(sep1+2,sep))-yy(sep1) ) >= 0))];
	end
	[~,sepid]=max(sepxy);
	index=find(y<x(sepid));
	xadj(index)=0;
end


function [xadj] = estsr_single_less(~,dtrs,m)
	xadj=ones(m,1);
	%index=find(dtrs<0.9);
	index=find(dtrs<(0.3*max(dtrs)));
	xadj(index)=0;
end

