function [nl,d,d0,ds] = get_light_plane(nl0, Kc, light_cb_intersects_cords,light_cb_intersects_points, squareSize)
% nl0：初始结构光平面法向量
% Kc: 摄像机内参矩阵
% light_cb_intersects_cord: 激光和棋盘格交点局部坐标
% light_cb_intersects_points: 激光和棋盘格交点图像坐标点，[x,y] 
% squareSize:棋盘格方块尺寸
%matlab 图像空间坐标系
% o--->x
% |
% y
num_image = length(light_cb_intersects_cords);
ds = zeros(num_image,1);
d_MNs = zeros(500,1);
M_n_s = zeros(500,3);
N_n_s = zeros(500,3);
a1 = nl0(1);
a2 = nl0(2);
a3 = nl0(3);

count_d = 0;
for i = 1:num_image
    cur_intersects_cord = light_cb_intersects_cords{i}*squareSize;
    cur_intersects_points = light_cb_intersects_points{i};
    cur_d_numerator = 0;
    cur_d_denominator = 0;
    for j = 1:size(cur_intersects_cord,1)-1
        % M点局部坐标
        M_cb = cur_intersects_cord(j,:);
        %M点图像坐标点
        m = cur_intersects_points(j,:);
        %其次话
        m(3)=1;
        % M点归一化平面坐标
        if size(m,1) == 1
            m = m';
        end
        M_n = Kc\m;
        f_M = M_n/(a1*M_n(1)+a2*M_n(2)+a3);
        for k = j+1:size(cur_intersects_cord,1)
            % N点局部坐标
            N_cb = cur_intersects_cord(k,:);
            %M点图像坐标点
            n = cur_intersects_points(k,:);
            %其次话
            n(3)=1;
            % M点归一化平面坐标
            if size(n,1) == 1
                n = n';
            end
            N_n = Kc\n;
            f_N = N_n/(a1*N_n(1)+a2*N_n(2)+a3);
            d_MN = norm(M_cb-N_cb);
            
            cur_d_numerator = cur_d_numerator + d_MN^2;
            cur_d_denominator = cur_d_denominator + norm(f_M-f_N)^2;
            count_d = count_d + 1;
            d_MNs(count_d) = d_MN;
            M_n_s(count_d,:) = M_n';
            N_n_s(count_d,:) = N_n';
%             d_cur = d_MN/(norm(f_M-f_N));
%             
%             ds(count_d) = d_cur;
        end

    end
    cur_d = sqrt(cur_d_numerator/cur_d_denominator);
    ds(i) = cur_d;

end
d_MNs = d_MNs(1:count_d);
M_n_s = M_n_s(1:count_d,:);
N_n_s = N_n_s(1:count_d,:);
index_select = find(~isnan(d_MNs));
d_MNs = d_MNs(index_select);
M_n_s = M_n_s(index_select,:);
N_n_s = N_n_s(index_select,:);

d0 = mean(ds(find(~isnan(ds))));
x0 = [a1, a2, a3, d0];
% 利用Levenberg-Marquardt算法，做分线性优化a1,a2,a3,d

% f_Ms = 1./(x(1)*M_n_s(:,1)+x(2)*M_n_s(:,2)+x(3)).*M_n_s;
% f_Ns = 1./(x(1)*N_n_s(:,1)+x(2)*N_n_s(:,2)+x(3)).*N_n_s;
% sum((sum((f_Ms - f_Ns).^2,2)*x(4)^2-d_MNs.^2).^2)
gfun = @(x)sum((sum((1./(x(1)*M_n_s(:,1)+x(2)*M_n_s(:,2)+x(3)).*M_n_s - 1./(x(1)*N_n_s(:,1)+x(2)*N_n_s(:,2)+x(3)).*N_n_s).^2,2)*x(4)^2-d_MNs.^2).^2);
options = optimoptions('lsqnonlin','Display','iter');
options.Algorithm = 'levenberg-marquardt';
[x,resnorm,residual,exitflag,output] = lsqnonlin(gfun,x0,[],[],options);
output
x = x/x(3);
nl = x(1:3);
d = x(4);
end









