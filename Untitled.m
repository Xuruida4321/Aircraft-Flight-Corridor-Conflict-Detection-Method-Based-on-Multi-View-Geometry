% ��ȡExcel�ļ�
[inputParams, ~, ~] = xlsread('input_parameters.xlsx');
% ����ȡ�����ݴ�С��ȷ�������Ԥ�ڣ����磬100��x16�У�
disp(size(inputParams));
results = zeros(30, 2); % ����ÿ����������һ�������һ��ʱ��ֵ
% results = zeros(size(inputParams, 1), 2);
% ѭ������ÿ���������
% for idx = 1:size(inputParams, 1)
for idx = 1:30
    % ��ȡ��ǰ�еĲ���
    xA = inputParams(idx, 1); yA = inputParams(idx, 2); zA = inputParams(idx, 3);
    xB = inputParams(idx, 4); yB = inputParams(idx, 5); zB = inputParams(idx, 6);
    a1 = inputParams(idx, 7); b1 = inputParams(idx, 8); h1 = inputParams(idx, 9);
    a2 = inputParams(idx, 10); b2 = inputParams(idx, 11); h2 = inputParams(idx, 12);
    D_A = inputParams(idx, 13); P_A = inputParams(idx, 14);
    D_B = inputParams(idx, 15); P_B = inputParams(idx, 16);
tic; % ��ʼ��ʱ
% �����A��B������
% xA = 3; yA = 3; zA = 8; % A������
% %xB = 10; yB = 2; zB = 6; % B����
% xB = 10; yB = 2; zB = 6; % B������

% % ������Բ����ĳߴ����
% a1 = 5; b1 = 3; h1 = 10; % A�ĳ��ᡢ����͸߶�
% a2 = 5; b2 = 3; h2 = 10; % B�ĳ��ᡢ����͸߶�

% % ������Բ����ķ���
% D_A = 45; % A�ķ���ǣ��ȣ�
% P_A = 30; % A�ĸ����ǣ��ȣ�
% D_B = 190; % B�ķ���ǣ��ȣ�
% P_B = 50; % B�ĸ����ǣ��ȣ�

% ������ת����
R_A = makeRotationMatrix(D_A, P_A);
R_B = makeRotationMatrix(D_B, P_B);

% ����һ������������������ת����

% ������Բ����A
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

% ������Բ����B
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

createInvisibleFigure(X_A, Y_A, Z_A, X_B, Y_B, Z_B);


% ���ü�⺯��
isIntersecting = checkIntersection(xA, yA, zA, a1, b1, h1, R_A, xB, yB, zB, a2, b2, h2, R_B);
if isIntersecting
    disp('������Բ�����ཻ��');
else
    disp('������Բ���岻�ཻ��');
end
elapsedTime = toc; % ֹͣ��ʱ����ȡ������ʱ��

fprintf('������ʱ��%.2f �롣\n', elapsedTime);
% ������ʱ����¼���
    elapsedTime = toc;
%     results(i, 1) = isIntersecting; % �����������ļ�����
%     results(i, 2) = elapsedTime;
    results(idx, 1) = isIntersecting; % �������ļ�����
    results(idx, 2) = elapsedTime; % ���Ǽ���ʱ��
    fprintf('��ǰ�����%d, ʱ�䣺%.4f\n', isIntersecting, elapsedTime);
    
end
% �����д��Excel�ļ�
disp('���ս����');
disp(results);
% �����Բ�����Ƿ��ཻ
function isIntersect = checkIntersection(xA, yA, zA, a1, b1, h1, R_A, xB, yB, zB, a2, b2, h2, R_B)
    % ��������õĿռ�����
    stepSize = 0.1;
    [xGrid, yGrid, zGrid] = meshgrid(min(xA, xB):stepSize:max(xA+a1, xB+a2),...
                                     min(yA, yB):stepSize:max(yA+b1, yB+b2),...
                                     min(zA, zB):stepSize:max(zA+h1, zB+h2));

    % ��ʼ�����
    isIntersect = false;

    % ��������еĵ��Ƿ�ͬʱ��������Բ�����ڲ�
    for i = 1:numel(xGrid)
        if isInsideEllipticalCylinder([xGrid(i), yGrid(i), zGrid(i)], xA, yA, zA, a1, b1, h1, R_A) && ...
           isInsideEllipticalCylinder([xGrid(i), yGrid(i), zGrid(i)], xB, yB, zB, a2, b2, h2, R_B)
            isIntersect = true;
            break;
        end
    end
end

% �����Ƿ�����Բ�����ڲ��ĺ���
function isInside = isInsideEllipticalCylinder(point, xC, yC, zC, a, b, h, R)
    % ����ת������Բ����ľֲ�����ϵ
    localPoint = R' * (point' - [xC; yC; zC]);

    % �����Ƿ��ھֲ�����ϵ�������Բ������
    isInside = (localPoint(1)^2 / a^2 + localPoint(2)^2 / b^2) <= 1 && ...
               localPoint(3) >= 0 && localPoint(3) <= h;
end

function R = makeRotationMatrix(D, P)
    Rx1 = [1 0 0; 0 0 -1; 0 1 0]; % ��X����ת90��    
    % ��Z����ת�����������
    Rz = [cosd(D) -sind(D) 0; sind(D) cosd(D) 0; 0 0 1];

    % ��X����ת������������
    Rx = [1 0 0; 0 cosd(P) -sind(P); 0 sind(P) cosd(P)];

    % ������ת��������Z����ת������X����ת
    R = Rx * Rz * Rx1
end
function createInvisibleFigure(X_A, Y_A, Z_A, X_B, Y_B, Z_B)
    % ����һ�����ɼ���ͼ�δ���
    fig = figure('Visible', 'off');

    % ������Բ����A��B
    surf(X_A, Y_A, Z_A);
    hold on;
    surf(X_B, Y_B, Z_B);
    hold off;

    % �����������ǩ�ͱ���
    xlabel('X��');
    ylabel('Y��');
    zlabel('Z��');
    title('������Բ����');

    % ������תģʽ
    rotate3d on;

    % �����Ҫ����ͼ�Σ�����������ʹ�� saveas ����
    % saveas(fig, 'myFigure.png');

    % �ر�ͼ�δ���
    close(fig);
end