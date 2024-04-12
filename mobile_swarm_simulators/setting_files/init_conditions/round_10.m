
% 周回用 1０台

% 台数の定義
posset_Na = 10;

% 配置の定義．
posset_x = zeros(2,5);
posset_y = zeros(2,5);

%posset_x(1,:) = -4;
%posset_x(2,:) = -5;
posset_x(1,:) = -6;
posset_x(2,:) = -7;
posset_y(1,:) = [6.4 5.2 4.0 2.8 1.6];
posset_y(2,:) = [7.1 5.9 4.7 3.5 2.3];
%posset_y(3,:) = [6.4 5.2 4.0 2.8 1.6];
%posset_y(4,:) = [7.1 5.9 4.7 3.5 2.3];

posset_x = reshape(posset_x.',[posset_Na,1]);
posset_y = reshape(posset_y.',[posset_Na,1]);