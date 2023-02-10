function [coors_camera] = get_points_camera_coord(coors_image, Kc, nl, d)
% Zc = -d/(a1*xn+a2*yn+a3)
% Xc = Zc*xn
% Yc = Zc*yn

num_points = size(coors_image,1);
coors_camera = zeros(num_points,3);

for i = 1:num_points
    % 1、获取归一化平面坐标
    xi = coors_image(i,1);
    yi = coors_image(i,2);
    pi = [xi, yi, 1]';
    pn = Kc\pi;
    % 2、获取相机坐标
    Zc = -d/(nl(1)*pn(1)+nl(2)*pn(2)+nl(3));
    Xc = Zc*pn(1);
    Yc = Zc*pn(2);
    coors_camera(i,:) = [Xc,Yc,Zc];
 
end
