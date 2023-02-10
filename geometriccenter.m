% 输入：光斑二值蒙版
% 输出：中心坐标
function coor = geometriccenter(mask)

% 连通域
[label,num] = bwlabel(mask);

% 计算
coor = zeros(num,2);
for n = 1:num
    [x,y] = find(label==n);
    coor(n,:) = [mean(x),mean(y)];
end

end
