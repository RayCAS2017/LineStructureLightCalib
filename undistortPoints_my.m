function points_undistorted = undistortPoints_my(points, Kc, D, max_iters, error_thresh)
% 坐标[u,v],points 是一个n*2的矩阵
% o--->u,x
% |
% v,y
if nargin < 5
    max_iters=1000;
    error_thresh = 1e-20;
end
fx = Kc(1,1);
fy = Kc(2,2);
cx = Kc(1,3);
cy = Kc(2,3);
%坐标系
%--->x
%|
%y
% 图像其次坐标[u,v,1]
uv = cat(2, points, ones(size(points,1),1));
% 归一化平面畸变坐标[xnd,ynd]
xnyn_d = Kc\uv';
X0 = xnyn_d(1:2,:)';
Xd = xnyn_d(1:2,:)';
for i = 1:max_iters
    if i == 1
        Xk = X0;
    else
        Xk = Xk1;
    end
    [Xk1,error] = single_iter(Xk,Xd,D);
    if error < error_thresh
        break
    end

end
points_undistorted = [Xk1(:,1)*fx+cx, Xk1(:,2)*fy+cy];

end

function [Xk1,error] = single_iter(Xk,Xd,D)
r2 = Xk(:,1).^2 + Xk(:,2).^2;
FXk = 1 + D(1)*r2+D(2)*r2.^2;
GXk = [2*D(4)*Xk(:,1).*Xk(:,2)+D(5)*(r2+2*Xk(:,1).^2), D(4)*(r2+2*Xk(:,2).^2)+2*D(5)*Xk(:,1).*Xk(:,2)];
Temp = Xd-GXk;
Xk1 = [Temp(:,1)./FXk, Temp(:,2)./FXk];
error = sqrt(sum((Xk1-Xk).* (Xk1-Xk),2));
error = sum(error);
end






 