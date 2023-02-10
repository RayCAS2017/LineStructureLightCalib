% 输入：原始灰度图像，光斑二值蒙版
% 输出：中心坐标
function coor = gravitycenter(im,mask)

% 连通域
[M,~] = size(im);
[label,num] = bwlabel(mask);
coor = zeros(num,2);

for n = 1:num
    [x,y] = find(label==n);
    % 取出对应点灰度
    idx = (y-1)*M+x; % matlab是列优先
    imtmp = im(idx);
    % 计算权值
    w = (imtmp/sum(imtmp))';   
    coor(n,:) = [w*x,w*y];
end

end
