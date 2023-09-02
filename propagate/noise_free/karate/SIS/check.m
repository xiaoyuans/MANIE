load('karate.txt');
x1=karate;
load('Xmatrix.txt');
x2=Xmatrix;
k=0;
for i=1:size(x2,1)
    for j=1:size(x2,2)
        if x1(i,j)~=x2(i,j)
            k=k+1;
        end
    end
end
disp(k);
%%数据集都是对称的矩阵