function [vl_light_inters, light_cb_intersects_cord, light_insects] = get_info_in_single_frame(imageFileName,is_show)

I = imread(imageFileName);

% Detect the checkerboard points
[imagePoints, boardSize] = detectCheckerboardPoints(I);% 对应空间坐标[y,x]
% 剔除NaN点
j = 0;
for i = 1:size(imagePoints, 1)
    if ~isnan(imagePoints(i,1)) && ~isnan(imagePoints(i,2))
        j = j+1;
        imagePoints1(j,:) = imagePoints(i,:);
    end
end

% Display detected points
if is_show
    J = insertText(I, imagePoints1,1:size(imagePoints1, 1)); 
    J = insertMarker(J, imagePoints1, 'o', 'Color','red', 'Size', 5); 
    J = insertMarker(J, [500, 1000], 'o', 'Color','green', 'Size', 10); 
    f1 = figure('Name','Detected checkboard points');
    imshow(J); 
    title(sprintf('Detected a %d x %dCheckerboard', boardSize));
end


if length(size(I))==3 
    im = rgb2gray(I);
    im = im2single(im);
else
    im = im2single(I);
end
im0=im;

%对应空间坐标[y, x]
[m, n] = size(im);


%Detect the checkboard horizontal and vertical lines
%1、线聚类的阈值统计
points_min_xy = zeros(size(imagePoints1, 1),2);
for i = 1 : size(imagePoints1, 1)
    min_dx = 10000;
    min_dy = 10000;
    for j = 1 : size(imagePoints1, 1)
        if j == i
            continue
        else
            dx = abs(imagePoints1(i,1) - imagePoints1(j,1));
            if dx < min_dx
                min_dx = dx;
            end
            dy = abs(imagePoints1(i,2) - imagePoints1(j,2));
            if dy < min_dy 
                min_dy = dy;
            end
        end
    end
    points_min_xy(i,1) = min_dx;
    points_min_xy(i,2) = min_dy;
end
thesh_dxy = 5*mean(points_min_xy,1);

squre_dis = zeros(size(imagePoints1, 1),1);
for i = 1 : size(imagePoints1, 1)
    min_squre_dis = 10000;
    for j = 1 : size(imagePoints1, 1)
        if j == i
            continue
        else
            dx = abs(imagePoints1(i,1) - imagePoints1(j,1));            
            dy = abs(imagePoints1(i,2) - imagePoints1(j,2));
            dis = sqrt(dy^2+dx^2);
            if dis < min_squre_dis
                min_squre_dis = dis;
            end
            
        end
    end
    squre_dis(i) = min_squre_dis;
end
thresh_dis = 2.5*mean(squre_dis);
%2、棋盘格水平方向线
line_horizon = zeros(boardSize(1)-1, boardSize(2)-1, 2);
line_horizon_index = zeros(boardSize(1)-1, boardSize(2)-1);
used = zeros(1,size(imagePoints1, 1));
for l = 1:size(line_horizon,1)
    for i = 1:size(imagePoints1,1)
        if line_horizon(l,1,1) == 0 && used(i) ==0
            line_horizon(l,1,:) = imagePoints1(i,:);
            used(i) = 1;
            line_horizon_index(l,1) = i;
        end
    end
    last_pt = line_horizon(l,1,:);
    for p = 2:size(line_horizon,2)
        for i = 1:size(imagePoints1,1)
            if used(i) == 1
                continue
            else
                dx = abs(imagePoints1(i,1) - last_pt(1));
                dy = abs(imagePoints1(i,2) - last_pt(2));
                dis = sqrt(dy^2+dx^2);
                if dy < thesh_dxy(2) && dis < thresh_dis
                    line_horizon(l,p,:) = imagePoints1(i,:);
                    line_horizon_index(l,p) = i;
                    last_pt = imagePoints1(i,:);
                    used(i) = 1;
                    break
                end
            end
        end
    end
end

line_horizon_paras = zeros(size(line_horizon_index,1),2);
for i = 1:size(line_horizon_index)    
    line_horizon_cur_index = line_horizon_index(i,:,:);
    line_cur_points_index = line_horizon_cur_index(line_horizon_cur_index>0);
    line_cur_points = imagePoints1(line_cur_points_index,:);
    x = line_cur_points(:,1);
    y = line_cur_points(:,2);
    line_cur_para = polyfit(y, x, 1);
    line_horizon_paras(i,:) = line_cur_para;
end
%3、棋盘格垂直方向线
line_vertical = zeros(boardSize(2)-1, boardSize(1)-1, 2);
line_vertical_index = zeros(boardSize(2)-1, boardSize(1)-1);
used = zeros(1,size(imagePoints1, 1));
for l = 1:size(line_vertical,1)
    for i = 1:size(imagePoints1,1)
        if line_vertical(l,1,1) == 0 && used(i) ==0
            line_vertical(l,1,:) = imagePoints1(i,:);
            used(i) = 1;
            line_vertical_index(l,1) = i;
        end
    end
    last_pt = line_vertical(l,1,:);
    for p = 2:size(line_vertical,2)
        for i = 1:size(imagePoints1,1)
            if used(i) == 1
                continue
            else
                dx = abs(imagePoints1(i,1) - last_pt(1));
                dy = abs(imagePoints1(i,2) - last_pt(2));
                dis = sqrt(dy^2+dx^2);
                if dx < thesh_dxy(1) && dis < thresh_dis
                    line_vertical(l,p,:) = imagePoints1(i,:);
                    line_vertical_index(l,p) = i;
                    last_pt = imagePoints1(i,:);
                    used(i) = 1;
                    break
                end
            end
        end
    end
end
line_vertical_paras = zeros(size(line_vertical_index,1),2);
for i = 1:size(line_vertical_index)    
    line_vertical_cur_index = line_vertical_index(i,:,:);
    line_cur_points_index = line_vertical_cur_index(line_vertical_cur_index>0);
    line_cur_points = imagePoints1(line_cur_points_index,:);
    x = line_cur_points(:,1);
    y = line_cur_points(:,2);
    line_cur_para = polyfit(y, x, 1);
    line_vertical_paras(i,:) = line_cur_para;
end

if is_show
    f3 = figure('Name','Detect checkboard lines and vanish points');
    imshow(im0);
    hold on 
    xx=1:1:n;
    for i = 1:size(line_horizon_index)
        yy = line_horizon_paras(i,1)*xx + line_horizon_paras(i,2);
        plot(yy, xx, 'blue')
        hold on
    end
    for i = 1:size(line_vertical_index)
        yy = line_vertical_paras(i,1)*xx + line_vertical_paras(i,2);
        plot(yy, xx, 'green')
        hold on
    end
end



%Detect the checkboard vanish points
%1、寻找水平线消失点
num_line_horizon = size(line_horizon_paras,1);
start_points = ones(num_line_horizon)*round(n/3);
end_points = ones(num_line_horizon)*round(2*n/3);
A_horizon = zeros(num_line_horizon,2);
B_horizon = zeros(num_line_horizon,2);
for i = 1:num_line_horizon
    A_horizon(i,1) = start_points(i);
    A_horizon(i,2) = line_horizon_paras(i,1)*start_points(i) + line_horizon_paras(i,2);
    B_horizon(i,1) = end_points(i);
    B_horizon(i,2) = line_horizon_paras(i,1)*end_points(i) + line_horizon_paras(i,2);
end
[X_h,P_h,R_h,x_h,p_h,l_h] = lineXline(A_horizon,B_horizon);
y_interset_h = x_h{2};
x_interset_h= x_h{1};
%2、寻找垂直线消失点
num_line_vertical = size(line_vertical_paras,1);
start_points = ones(num_line_vertical)*round(n/3);
end_points = ones(num_line_vertical)*round(2*n/3);
A_vertical = zeros(num_line_vertical,2);
B_vertical = zeros(num_line_vertical,2);
for i = 1:num_line_vertical
    A_vertical(i,1) = start_points(i);
    A_vertical(i,2) = line_vertical_paras(i,1)*start_points(i) + line_vertical_paras(i,2);
    B_vertical(i,1) = end_points(i);
    B_vertical(i,2) = line_vertical_paras(i,1)*end_points(i) + line_vertical_paras(i,2);
end
[X_v,P_v,R_v,x_v,p_v,l_v] = lineXline(A_vertical,B_vertical);
y_interset_v = x_v{2};
x_interset_v = x_v{1}; 
%3、水平线消失点投影误差
proj_error_h = zeros(num_line_horizon,1);
for i = 1:num_line_horizon
    pt1 = zeros(2,1);
    pt2 = zeros(2,1);
    pt1(1) = x_interset_h - 5*i;
    pt1(2) = line_horizon_paras(i,1)*pt1(1) + line_horizon_paras(i,2);
    pt2(1) = x_interset_h + 5*i;
    pt2(2) = line_horizon_paras(i,1)*pt2(1) + line_horizon_paras(i,2);
    tmpt = [0,-1;1,0];
    nl = tmpt*(pt2-pt1)/sqrt((pt2-pt1)'*(pt2-pt1));
    inter_pt = [x_interset_h;y_interset_h];
    dis_cur = sqrt((inter_pt-pt1)'*nl*nl'*(inter_pt-pt1));
    proj_error_h(i) = dis_cur;
end
%4、垂直线消失点投影误差
proj_error_v = zeros(num_line_vertical,1);
for i = 1:num_line_vertical
    pt1 = zeros(2,1);
    pt2 = zeros(2,1);
    pt1(1) = x_interset_v - 5*i;
    pt1(2) = line_vertical_paras(i,1)*pt1(1) + line_vertical_paras(i,2);
    pt2(1) = x_interset_v + 5*i;
    pt2(2) = line_vertical_paras(i,1)*pt2(1) + line_vertical_paras(i,2);
    tmpt = [0,-1;1,0];
    nl = tmpt*(pt2-pt1)/sqrt((pt2-pt1)'*(pt2-pt1));
    inter_pt = [x_interset_v;y_interset_v];
    dis_cur = sqrt((inter_pt-pt1)'*nl*nl'*(inter_pt-pt1));
    proj_error_v(i) = dis_cur;
end
% 5、根据投影误差，剔除一些误差比较大的线，重新计算消失点
% 5.1、水平线消失点
proj_error_h_std = std(proj_error_h);
line_h_selected_index = proj_error_h < proj_error_h_std;
line_h_selected = line_horizon_paras(line_h_selected_index,:);
if length(line_h_selected) < 3
    [B,I] = mink(proj_error_h,3);
    line_h_selected = line_horizon_paras(I,:);
end
num_line_select_h = size(line_h_selected,1);
start_points = ones(num_line_select_h)*round(n/3);
end_points = ones(num_line_select_h)*round(2*n/3);
A_horizon = zeros(num_line_select_h,2);
B_horizon = zeros(num_line_select_h,2);
for i = 1:num_line_select_h
    A_horizon(i,1) = start_points(i);
    A_horizon(i,2) = line_h_selected(i,1)*start_points(i) + line_h_selected(i,2);
    B_horizon(i,1) = end_points(i);
    B_horizon(i,2) = line_h_selected(i,1)*end_points(i) + line_h_selected(i,2);
end
[X_h,P_h,R_h,x_h,p_h,l_h] = lineXline(A_horizon,B_horizon);
y_interset_h_2 = x_h{2};
x_interset_h_2 = x_h{1};
%5.2、垂直线消失点
proj_error_v_std = std(proj_error_v);
line_v_selected_index = proj_error_v < proj_error_v_std;
line_v_selected = line_vertical_paras(line_v_selected_index,:);
if length(line_v_selected) < 3
    [B,I] = mink(proj_error_v,3);
    line_v_selected = line_vertical_paras(I,:);
end
num_line_select_v = size(line_v_selected,1);
start_points = ones(num_line_select_v)*round(n/3);
end_points = ones(num_line_select_v)*round(2*n/3);
A_vertical = zeros(num_line_select_v,2);
B_vertical = zeros(num_line_select_v,2);
for i = 1:num_line_select_v
    A_vertical(i,1) = start_points(i);
    A_vertical(i,2) = line_v_selected(i,1)*start_points(i) + line_v_selected(i,2);
    B_vertical(i,1) = end_points(i);
    B_vertical(i,2) = line_v_selected(i,1)*end_points(i) + line_v_selected(i,2);
end

[X_v,P_v,R_v,x_v,p_v,l_v] = lineXline(A_vertical,B_vertical);
y_interset_v_2 = x_v{2};
x_interset_v_2 = x_v{1};

if is_show
    f4=figure('Name','Vanish point project error');
    h_x_plot = round(x_interset_h_2-100):round(x_interset_h_2+100);
    proj_error_h_2 = zeros(num_line_select_h,1);
    for i = 1:num_line_select_h
        h_y_plot = line_h_selected(i,1)*h_x_plot +  line_h_selected(i,2);
        subplot(1,2,1)
        plot(h_y_plot,h_x_plot)
        hold on
        pt1 = zeros(2,1);
        pt2 = zeros(2,1);
        pt1(1) = x_interset_h_2 - 5*i;
        pt1(2) = line_h_selected(i,1)*pt1(1) + line_h_selected(i,2);
        pt2(1) = x_interset_h_2 + 5*i;
        pt2(2) = line_h_selected(i,1)*pt2(1) + line_h_selected(i,2);
        tmpt = [0,-1;1,0];
        nl = tmpt*(pt2-pt1)/sqrt((pt2-pt1)'*(pt2-pt1));
        inter_pt = [x_interset_h_2;y_interset_h_2];
        dis_cur = sqrt((inter_pt-pt1)'*nl*nl'*(inter_pt-pt1));
        proj_error_h_2(i) = dis_cur;
        text(pt1(2), pt1(1), strcat(num2str(dis_cur),'\rightarrow'), 'HorizontalAlignment', 'right', 'FontSize',10);
        hold on
    end
    plot(y_interset_h_2,x_interset_h_2,'*','Color','green', 'MarkerSize',5)
    % 垂直线消失点投影误差
    v_x_plot = round(x_interset_v_2-100):round(x_interset_v_2+100);
    proj_error_v_2 = zeros(num_line_select_v,1);
    for i = 1:num_line_select_v
        v_y_plot = line_v_selected(i,1)*v_x_plot +  line_v_selected(i,2);
        subplot(1,2,2)
        plot(v_y_plot,v_x_plot)
        hold on
        pt1 = zeros(2,1);
        pt2 = zeros(2,1);
        pt1(1) = x_interset_v_2 - 5*i;
        pt1(2) = line_v_selected(i,1)*pt1(1) + line_v_selected(i,2);
        pt2(1) = x_interset_v_2 + 5*i;
        pt2(2) = line_v_selected(i,1)*pt2(1) + line_v_selected(i,2);
        tmpt = [0,-1;1,0];
        nl = tmpt*(pt2-pt1)/sqrt((pt2-pt1)'*(pt2-pt1));
        inter_pt = [x_interset_v_2;y_interset_v_2];
        dis_cur = sqrt((inter_pt-pt1)'*nl*nl'*(inter_pt-pt1));
        proj_error_v_2(i) = dis_cur;

        text(pt1(2), pt1(1), strcat(num2str(dis_cur),'\rightarrow'), 'HorizontalAlignment', 'right', 'FontSize',10);

        hold on
    end
    plot(y_interset_v_2,x_interset_v_2,'*','Color','green', 'MarkerSize',5)
end


%Detect light line
%1、截取棋盘格区域roi,排除其他高亮区域对于激光线的提取
topleft=[round(min(imagePoints1(:,2))), round(min(imagePoints1(:,1)))];
bottomright = [round(max(imagePoints1(:,2))), round(max(imagePoints1(:,1)))];
im(1:topleft(1),:) = 0;
im(:,1:topleft(2)) = 0;
im(:,bottomright(2):n) = 0;
im(bottomright(1):m,:) = 0;
%2、利用steger提取激光线中心点
level = 0.9;
BW = imbinarize(im,level);
linepixel = steger(BW);
linepixel(:,[1,2]) = linepixel(:,[2,1]);
% 进一步约束激光点在角点区域内
% 1）提取棋盘格角点的4个端点
horizon_line1 = line_horizon_paras(1,:);
horizon_line2 = line_horizon_paras(end,:);
vertical_line1 = line_vertical_paras(1,:);
vertical_line2 = line_vertical_paras(end,:);
vertices = zeros(5,2);
% vertice1 = two_line_intersect(horizon_line1,vertical_line1);
% vertice2 = two_line_intersect(horizon_line1,vertical_line2);
% vertice3 = two_line_intersect(horizon_line2,vertical_line1);
% vertice4 = two_line_intersect(horizon_line2,vertical_line2);
vertices(1,:) = two_line_intersect(horizon_line1,vertical_line1);
vertices(2,:) = two_line_intersect(horizon_line1,vertical_line2);
vertices(3,:) = two_line_intersect(horizon_line2,vertical_line2);
vertices(4,:) = two_line_intersect(horizon_line2,vertical_line1);
vertices(5,:) = vertices(1,:) ;
% 2)判断激光像素点是否在这4个端点围成的4变形内部
x = linepixel(:,1);
y = linepixel(:,2);
in = inpolygon(y,x,vertices(:,1),vertices(:,2));
x1 = x(in);
y1 = y(in);
%3、拟合激光线
% line_light=polyfit(y, x, 1);
line_light=polyfit(y1, x1, 1);

if is_show
    f2=figure('Name','Detect light line');
    subplot(121)
    imshow(im);
    hold on 
%     for i = 1:4
%         plot(vertices(i,2),vertices(i,1),'r*', 'MarkerSize',10);
%         hold on
%     end
    plot(vertices(:,2),vertices(:,1),'r','LineWidth',1.5);
    plot(x(in),y(in),'g+');
    hold on;
    plot(x(~in),y(~in),'bo');
    hold on
    xx=1:1:n;
    yy = line_light(1)*xx + line_light(2);
    plot(yy, xx,'red')
    subplot(122)
    imshow(BW);
end


% Detect the intersection of the vanish line and light line
vp_h = [x_interset_h_2; y_interset_h_2];
vp_v = [x_interset_v_2;y_interset_v_2];
vl_k = (vp_h(2) - vp_v(2))/(vp_h(1)-vp_v(1));
vl_b = (vp_h(2)*vp_v(1)-vp_v(2)*vp_h(1))/(vp_v(1)-vp_h(1));
vl = [vl_k;vl_b];

vl_light_inters_x = (line_light(2)-vl(2))/(vl(1)-line_light(1));
vl_light_inters_y = vl(1)*vl_light_inters_x+vl(2);
vl_light_inters = [vl_light_inters_x;vl_light_inters_y];
if is_show
    figure(f3);
    hold on
    plot(vp_h(2),vp_h(1),'*','Color','black', 'MarkerSize', 10)
    hold on
    plot(vp_v(2),vp_v(1),'*','Color','black', 'MarkerSize', 10)
    hold on
    plot([vp_h(2), vp_v(2)],[vp_h(1), vp_v(1)],'Color','black')
    hold on
    plot(vl_light_inters(2),vl_light_inters(1),'*','Color','red', 'MarkerSize', 10)
    hold on
    light_x = n;
    light_y = line_light(1)*light_x + line_light(2);
    plot_light_point = [light_x, light_y];
    plot([vl_light_inters(2), plot_light_point(2)],[vl_light_inters(1), plot_light_point(1)],'Color','black')
end


% Detect the checkerboard corners
cb_corners = zeros(size(line_vertical_paras,1),size(line_horizon_paras,1),2);
light_insects = zeros(size(line_vertical_paras,1),2);
for i=1:size(line_vertical_paras,1)
    light_insects(i,:) = two_line_intersect(line_vertical_paras(i,:), line_light);
    for j = 1:size(line_horizon_paras,1)
        cb_corners(i,j,:) = two_line_intersect(line_vertical_paras(i,:),line_horizon_paras(j,:));
    end
end

if is_show
    f5 = figure('Name', 'Checkerboard corners and light intersects');
    imshow(im0);
    hold on
    for i=1:size(line_vertical_paras,1)
        plot(light_insects(i,2),light_insects(i,1),'r*', 'MarkerSize', 10);
        hold on;
        text(light_insects(i,2), light_insects(i,1), strcat(num2str(i),'\rightarrow'), 'HorizontalAlignment', 'right', 'FontSize',15, 'Color', 'red');
        hold on;
        for j = 1:size(line_horizon_paras,1)
            plot(cb_corners(i,j,2),cb_corners(i,j,1),'g+', 'MarkerSize', 10);
            hold on;
            text(cb_corners(i,j,2), cb_corners(i,j,1), strcat('(',num2str(i),',',num2str(j),')','\rightarrow'), 'HorizontalAlignment', 'right', 'FontSize',15, 'Color', 'green');
            hold on;
        end
    end
    
end


%根据交比不变性求激光和棋盘格交点的局部物理坐标
cb_h_n = size(cb_corners,2); %棋盘格水平线数
light_cb_intersects_cord = zeros(size(light_insects,1),2);
for li = 1:size(light_insects,1)
    light_cb_intersects_cord(li,1)=li-1;
    % 当前激光与棋盘格交点在图像中的坐标
    m= light_insects(li,:);
    m = squeeze(m);
    % 当前激光与棋盘格交点在棋盘格中的坐标
    M = [];
    % 取A、B、C三个点
    count_n = 0;
    for i = 1:cb_h_n-2
        a = cb_corners(li,i,:);
        a = squeeze(a)';
        A = [li-1, i-1];
        for j = i+1:cb_h_n-1
            b = cb_corners(li,j,:);
            b = squeeze(b)';
            B = [li-1,j-1];
            for k=j+1:cb_h_n
                c = cb_corners(li,k,:);
                c = squeeze(c)';
                C = [li-1,k-1];
                ac = sqrt((a-c)*(a-c)');
                bc = sqrt((b-c)*(b-c)');
                am = sqrt((a-m)*(a-m)');
                bm = sqrt((b-m)*(b-m)');
                CrossRatio = (ac/bc)/(am/bm);
                AC = A(2)-C(2);
                BC = B(2)-C(2);
                AM_divide_BM = (AC/BC)/CrossRatio;
                M_cur = (A(2)-AM_divide_BM*B(2))/(1-AM_divide_BM);
                count_n = count_n+1;
                M(count_n)=M_cur;
            end
        end
    end
    for tt = 1:3
        M_mean = mean(M);
        M_std = std(M);
        M = M(abs(M-M_mean)<0.6*M_std);
    end    
    light_cb_intersects_cord(li,2) = mean(M);
end





