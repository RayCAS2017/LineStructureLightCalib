function linePixel = steger(EinRDR)
%采用steger算法提取光条中心
 im = double(EinRDR);
 im = im/max(im(:)); %%归一化
 bw=im2bw(im, 0.4);  %%二值化
 bw = bwmorph(bw, 'clean', 1);    % 去孤立点
[x0,y0]=find(bw==1);  %%找到白点坐标
sigma=25/sqrt(3);%求Hessian矩阵之前需要对图像进行高斯滤波，高斯滤波时，设置高斯方差σ<（w/√3），w是光条宽度
[Dx,Dy,Dxx,Dxy,Dyy] = Hessian2D(im,sigma);
[eigenvalue1, eigenvalue2, eigenvectorx, eigenvectory]=eig2image(Dxx, Dxy, Dyy); 
%eigenvectorx,eigenvectory:Hessian矩阵最大特征值对应的特征向量对应于光条的法线方向 
t = -(Dx.*eigenvectorx + Dy .* eigenvectory) ./...  
    (Dxx .* eigenvectorx.^2 + 2*Dxy.*eigenvectorx.*eigenvectory + Dyy.*eigenvectory.^2 );   
px = t.*eigenvectorx;  
py = t.*eigenvectory;
[candidateX1, candidateY1] = find(px >= -0.5 & px <= 0.5 & py >= -0.5 & py <= 0.5 & bw==1);  
%判断：如果(px,py)∈[?0.5,0.5]×[?0.5,0.5]
%即一阶导数为零的点位于当前像素内，
%且(nx,ny)方向的二阶导数大于指定的阈值，
%则该点(candidataX1,candidataY1）为光条的中心点，candidataX1+px,candidataY1+py 则为亚像素坐标。
linePixel_t = [candidateX1, candidateY1];
for i=1:size(candidateX1,1)
    m1=candidateX1(i,1);
    n1=candidateY1(i,1);
    px1(i,1)=px(m1,n1);
    py1(i,1)=py(m1,n1);
    x1(i,1)=m1+px(m1,n1);
    x2(i,1)=py(m1,n1)+n1;
end
linePixel = [x1,x2];
%获得激光线像素数点坐标
end

function [Dx,Dy,Dxx,Dxy,Dyy] = Hessian2D(I,Sigma)
%构造高斯模板
if nargin < 2, Sigma = 1; end
[X,Y]   = ndgrid(-round(3*Sigma):round(3*Sigma));
DGaussx  = 1/(2*pi*Sigma^4)*(-X).* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussy  = 1/(2*pi*Sigma^4)*(-Y).* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussxx = 1/(2*pi*Sigma^4) * (X.^2/Sigma^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussxy = 1/(2*pi*Sigma^6) * (X .* Y)           .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussyy = DGaussxx';
%卷积
Dx  = imfilter(I,DGaussx,'conv');
Dy  = imfilter(I,DGaussy,'conv');
Dxx = imfilter(I,DGaussxx,'conv');
Dxy = imfilter(I,DGaussxy,'conv');
Dyy = imfilter(I,DGaussyy,'conv');
end

function [Lambda1,Lambda2,Ix,Iy]=eig2image(Dxx,Dxy,Dyy)
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v2x = 2*Dxy; v2y = Dyy - Dxx + tmp;
% Normalize
%标准化
mag = sqrt(v2x.^2 + v2y.^2); i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);
% The eigenvectors are orthogonal
%实对称矩阵性质：不同特征值对应的特征向量是正交的
v1x = -v2y; 
v1y = v2x;
% Compute the eigenvalues
%计算特征值
mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);
% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
%按绝对值abs(Lambda1)<abs(Lambda2)对特征值排序
check=abs(mu1)>abs(mu2);
Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);
Ix=v1x; Ix(check)=v2x(check);
Iy=v1y; Iy(check)=v2y(check);
end

