clear all
close all
clc

load('intrinsic.mat')
% imageFileName='D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_56_45_515.bmp';
% imageFileName='D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_55_34_404.bmp';
% imageFileName='D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_55_10_099.bmp';
imageFileNames={'D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_56_45_515.bmp',...
    'D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_55_34_404.bmp',...
    'D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_55_10_099.bmp',...
    };
squareSize = 20; % mm
num_image = length(imageFileNames);
light_vps=zeros(num_image,2);
light_cb_intersects_cords = cell(num_image,1);
light_cb_intersects_points = cell(num_image,1);
for i = 1:num_image
    [light_vp,light_cb_intersects_cord, light_insects] = get_info_in_single_frame(imageFileNames{i}, 0);
    light_vps(i,:) = light_vp([2,1],:)'; % vp空间坐标[x,y]
    light_cb_intersects_cords{i} = light_cb_intersects_cord; % 激光线与棋盘格交点在棋盘格局部坐标
    light_cb_intersects_points{i} = light_insects(:,[2,1]); % 交点图像坐标[x,y]
end

%matlab 图像空间坐标系
% o--->x
% |
% y
% 拟合激光平面消隐线
light_vl = polyfit(light_vps(:,1),light_vps(:,2),1);
f_lvl = figure('Name', 'Light vanish line');
plot(light_vps(:,2),light_vps(:,1),'r*');
hold on
xx = round(min(light_vps(:,1))-100):round(max(light_vps(:,1))+100);
yy =  light_vl(1)*xx + light_vl(2);
plot(yy,xx,'Color','black');

%由消隐线获得激光平面法向量（a1,a2,a3）
Kc = cameraParams4.Intrinsics.IntrinsicMatrix';
lp = [light_vl(1),-1,light_vl(2)]';
nl0 = Kc'*lp;
nl0 = nl0/nl0(3);

[nl,d, d0] = get_light_plane(nl0, Kc, light_cb_intersects_cords, light_cb_intersects_points, squareSize);

save light_plane nl d nl0 d0;

