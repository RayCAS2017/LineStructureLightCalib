% 输入：原始灰度图像，光斑二值蒙版
% 输出：中心坐标
function coor = gausscenter(im,mask)

% 连通域
[label,num] = bwlabel(mask);

% 计算
coor = zeros(num,2);
for n = 1:num
    [x,y] = find(label==n);
    % 生成计算矩阵
    m_iN = length(x);
    tmp_A = zeros(m_iN,1);
    tmp_B = zeros(m_iN,5);    
    for k = 1:m_iN        
        pSrc = im(x(k),y(k));
        if pSrc>0
            tmp_A(k) = pSrc*log(pSrc);
        end
        tmp_B(k,1) = pSrc ;
        tmp_B(k,2) = pSrc*x(k);
        tmp_B(k,3) = pSrc*y(k);
        tmp_B(k,4) = pSrc*x(k)*x(k);
        tmp_B(k,5) = pSrc*y(k)*y(k);
    end
    
    % QR分解
    Vector_A = tmp_A;
    matrix_B = tmp_B;
    [Q,R] = qr(matrix_B);
    
    % 求解中心
    S = Q'*Vector_A;
    S = S(1:5);
    R1 = R(1:5,1:5);
    C = R1\S;   
    coor(n,:) = [-0.5*C(2)/C(4),-0.5*C(3)/C(5)];
end

end
