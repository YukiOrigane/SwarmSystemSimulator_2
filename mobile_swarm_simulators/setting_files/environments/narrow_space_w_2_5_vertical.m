% 壁あり

% 環境サイズ決定
envset_xmin = -8;
envset_xmax = 8;
envset_ymin = -8;
envset_ymax = 8;

% 障害物
L1_ = 5;    % 狭くなりはじめ-広くなり終わりまでの距離
L2_ = 5;    % 狭所長さ
L3_ = 2.5;  % 狭所幅, w
% L3_ = obj.param.w;

% [[x1,y1];[x2,y2]]の形で線分を指定していく
envset_wall_segments(:,:,1) = [[envset_xmin,envset_ymin];[envset_xmin,envset_ymax]];
envset_wall_segments(:,:,2) = [[envset_xmin,envset_ymax];[-L1_/2,envset_ymax]];
envset_wall_segments(:,:,3) = [[-L1_/2,envset_ymax];[-L2_/2,L3_/2]];
envset_wall_segments(:,:,4) = [[-L2_/2,L3_/2];[+L2_/2,L3_/2]];
envset_wall_segments(:,:,5) = [[+L2_/2,L3_/2];[+L1_/2,envset_ymax]];
envset_wall_segments(:,:,6) = [[+L1_/2,envset_ymax];[envset_xmax,envset_ymax]];
envset_wall_segments(:,:,7) = [[envset_xmax,envset_ymax];[envset_xmax,envset_ymin]];
envset_wall_segments(:,:,8) = [[envset_xmax,envset_ymin];[+L1_/2,envset_ymin]];
envset_wall_segments(:,:,9) = [[+L1_/2,envset_ymin];[+L2_/2,-L3_/2]];
envset_wall_segments(:,:,10) = [[+L2_/2,-L3_/2];[-L2_/2,-L3_/2]];
envset_wall_segments(:,:,11) = [[-L2_/2,-L3_/2];[-L1_/2,envset_ymin]];
envset_wall_segments(:,:,12) = [[-L1_/2,envset_ymin];[envset_xmin,envset_ymin]];

envset_upwall = [3,4,5];    % 上限を決める壁（枠除く）