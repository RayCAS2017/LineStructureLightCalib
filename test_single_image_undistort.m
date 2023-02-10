clear all
close all
clc

load('intrinsic.mat')

% imageFileName='D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_56_45_515.bmp';
imageFileName='D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_55_34_404.bmp';
% imageFileName='D:\Projects\LineStructureLightCalib\data\checkerboard_light_selected\2023-01-17_13_55_10_099.bmp';
[light_vp,light_cb_intersects_cord, light_insects] = get_info_in_single_frame_undistort(imageFileName,cameraParams4, 1);