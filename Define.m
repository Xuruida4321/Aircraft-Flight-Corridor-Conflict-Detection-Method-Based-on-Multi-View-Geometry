tic; % 开始计时
% 定义点A和B的坐标
xA = 3; yA = 3; zA = 8; % A的坐标
xB = 10; yB = 2; zB = 6; % B的坐
% xB = 15; yB = 2; zB = 7; % B的坐标

% 定义椭圆柱体的尺寸参数
a1 = 5; b1 = 3; h1 = 10; % A的长轴、短轴和高度
a2 = 5; b2 = 3; h2 = 10; % B的长轴、短轴和高度

% 定义椭圆柱体的方向
D_A = 45; % A的方向角（度）
P_A = 30; % A的俯仰角（度）
D_B = 190; % B的方向角（度）
P_B = 50; % B的俯仰角（度）

% 计算旋转矩阵
R_A = makeRotationMatrix(D_A, P_A);
R_B = makeRotationMatrix(D_B, P_B);

% 定义一个辅助函数来创建旋转矩阵

% 创建椭圆柱体A
[theta, z] = meshgrid(linspace(0, 2*pi, 50), linspace(0, h1, 50));
X_A = a1*cos(theta);
Y_A = b1*sin(theta);
Z_A = z;
for i = 1:size(Z_A, 1)
    for j = 1:size(Z_A, 2)
        point = R_A * [X_A(i,j); Y_A(i,j); Z_A(i,j)];
        X_A(i,j) = point(1);
        Y_A(i,j) = point(2);
        Z_A(i,j) = point(3);
    end
end
X_A = X_A + xA;
Y_A = Y_A + yA;
Z_A = Z_A + zA;

% 创建椭圆柱体B
[theta, z] = meshgrid(linspace(0, 2*pi, 50), linspace(0, h2, 50));
X_B = a2*cos(theta);
Y_B = b2*sin(theta);
Z_B = z;
for i = 1:size(Z_B, 1)
    for j = 1:size(Z_B, 2)
        point = R_B * [X_B(i,j); Y_B(i,j); Z_B(i,j)];
        X_B(i,j) = point(1);
        Y_B(i,j) = point(2);
        Z_B(i,j) = point(3);
    end
end
X_B = X_B + xB;
Y_B = Y_B + yB;
Z_B = Z_B + zB;

% 绘制椭圆柱体A和B
figure;
surf(X_A, Y_A, Z_A);
hold on;
surf(X_B, Y_B, Z_B);
hold off;
xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');
title('两个椭圆柱体');
rotate3d on; % 激活旋转模式




% 调用检测函数
isIntersecting = checkIntersection(xA, yA, zA, a1, b1, h1, R_A, xB, yB, zB, a2, b2, h2, R_B);
if isIntersecting
    disp('两个椭圆柱体相交。');
else
    disp('两个椭圆柱体不相交。');
end
elapsedTime = toc; % 停止计时并获取经过的时间

fprintf('计算用时：%.2f 秒。\n', elapsedTime);

% 检测椭圆柱体是否相交
function isIntersect = checkIntersection(xA, yA, zA, a1, b1, h1, R_A, xB, yB, zB, a2, b2, h2, R_B)
    % 创建检测用的空间网格
    stepSize = 0.1;
    [xGrid, yGrid, zGrid] = meshgrid(min(xA, xB):stepSize:max(xA+a1, xB+a2),...
                                     min(yA, yB):stepSize:max(yA+b1, yB+b2),...
                                     min(zA, zB):stepSize:max(zA+h1, zB+h2));

    % 初始化标记
    isIntersect = false;

    % 检测网格中的点是否同时在两个椭圆柱体内部
    for i = 1:numel(xGrid)
        if isInsideEllipticalCylinder([xGrid(i), yGrid(i), zGrid(i)], xA, yA, zA, a1, b1, h1, R_A) && ...
           isInsideEllipticalCylinder([xGrid(i), yGrid(i), zGrid(i)], xB, yB, zB, a2, b2, h2, R_B)
            isIntersect = true;
            break;
        end
    end
end

% 检测点是否在椭圆柱体内部的函数
function isInside = isInsideEllipticalCylinder(point, xC, yC, zC, a, b, h, R)
    % 将点转换到椭圆柱体的局部坐标系
    localPoint = R' * (point' - [xC; yC; zC]);

    % 检测点是否在局部坐标系定义的椭圆柱体内
    isInside = (localPoint(1)^2 / a^2 + localPoint(2)^2 / b^2) <= 1 && ...
               localPoint(3) >= 0 && localPoint(3) <= h;
end

function R = makeRotationMatrix(D, P)
    Rx1 = [1 0 0; 0 0 -1; 0 1 0]; % 绕X轴旋转90度    
    % 绕Z轴旋转，调整方向角
    Rz = [cosd(D) -sind(D) 0; sind(D) cosd(D) 0; 0 0 1];

    % 绕X轴旋转，调整俯仰角
    Rx = [1 0 0; 0 cosd(P) -sind(P); 0 sind(P) cosd(P)];

    % 最终旋转矩阵：先绕Z轴旋转，再绕X轴旋转
    R = Rx * Rz * Rx1
end