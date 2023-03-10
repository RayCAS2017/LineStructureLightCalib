#### 简介
实现了文献[1]的算法，基于棋盘格，利用消隐线和交比不变性，实现3D线结构光相机的标定。

#### 运行环境
windows10, matlab2022

#### 标定步骤
1、利用matlab标定相机内参和畸变系数。

**Notes**:

可以多拍些图片，matlab会自动剔除一些棋盘格角点检测效果不好的图片，同时可以手动一步一步剔除投影误差比较大的图像，不断优化相机标定，将标定结果保存为intrinsic.mat;

2、利用带激光线条的棋盘格图片（data/checkerboard_light_selected）标定激光平面,分为两个版本:1)main_func_calibrate.m和2）main_func_calibrate_undistort.m, 后者做了畸变矫正。

**Notes**：

i)理论上只要两幅就可以标定出激光平面，但实际中可以多采一些图，采集的过程中，棋盘格不要与相机正视，要有一定角度，形成透视变换，这样好检测消隐点，同时激光线条最好完整穿过棋盘格。一副好的激光-棋盘格图片，程序应该能够检测出尽量多的棋盘格的角点，使得棋盘格的水平和垂直线簇能够拟合，从而得到所有的角点。一个简单判断这副图片是否满足要求的方法是，利用matlab自带的相机标定工具，加载这些图片，看角点检测情况，从而决定取哪些图。

![相机标定](https://github.com/RayCAS2017/LineStructureLightCalib/raw/main/assets/matlab_camera_cali.jpg)

ii) 棋盘格的水平线簇和垂直线簇不可能都交于一点，可以通过看投影误差看相交的怎么样。程序代码中，剔除了一些投影误差比较大的线；

![消隐点投影误差](https://github.com/RayCAS2017/LineStructureLightCalib/raw/main/assets/vp_proj_error.jpg)

iii) 激光线条中心像素点的提取。为了不受周边高亮区域的影响，程序自动提取了一个大致靶面。这个时候，可能棋盘格边缘白色高亮的地方也检测为激光像素点。所以，在拟合激光线条的时候，又做了一层处理，只提取了棋盘格角点最小外接四边形内的像素点。如下图所示，绿色为区域内激光像素点，蓝色为区域外激光像素点。

![激光线提取](https://github.com/RayCAS2017/LineStructureLightCalib/raw/main/assets/detect_laser_line.jpg)

iv)采用文献[1]中式（12），做Levenberg-Marguqrdt进行非线性优化，联合优化线结构光的4个参数（a1,a2,a3，d）,确实提升了精度。

#### 验证步骤
文件main_verify.mlx实现了对标定结果的验证。采用针形标靶数据（data/needlepoint）,其中有不同高度和不同间距的细针。采用畸变矫正标定和验证的结果，比不用的结果，精度提升比较多。如下图所示，两个光点之间的距离是40mm,上面是不用畸变矫正的结果，下面是使用畸变矫正的结果。

![验证结果](https://github.com/RayCAS2017/LineStructureLightCalib/raw/main/assets/verify_results.jpg)

#### 改进项

1、可以利用针形、锯齿等标靶，再采用文献[1]中式（12），对激光平面进行一轮优化，精度应该能够变得更高。

2、算法程序标定的是相机世界坐标系下的线结构光平面。在实际使用过程中，可能还需要标定个相机世界坐标系到一个测量世界坐标系（例如，线结构光平面）的外参。


#### 参考文献：
1、基于共面靶标的线结构光传感器标定新方法

2、[最小二乘求多条线交点wiki](https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection)

