
classdef MobileRobots2dSimulator < Simulator
    % 時空間のシミュレータ
    properties
        % システムの変数を記載
        t_vec   % 固有時刻ベクトル
        x       % ロボット位置 [台数,空間次元,時刻]
        dxdt    % ロボット速さ [台数,空間次元,時刻]
        u       % 入力
        G       % グラフオブジェクト．MATLABのgraph参照
        wall     % 壁セグメントの集合 [開始/終了,空間次元,セグメント]
        variances   % 位置の分散の時刻履歴 [台数,1,時刻]
        is_stop  % 停止状態か？ [台数,1,時刻]
    end

    methods
        function obj = MobileRobots2dSimulator()
            % コンストラクタ（宣言時に呼ばれる）
            obj@Simulator();    % 親クラスのコンストラクタも呼び出す
            obj = obj.setDefaultParameters();       % パラメタのデフォルト値を設定
        end

    %%%%%%%%%%%%%% 初期設定まわり %%%%%%%%%%%%%

        function obj = setDefaultParameters(obj)
            % パラメータとデフォルト値を設定
            obj = obj.setDefaultParameters@Simulator();   % スーパークラス側の読み出し
            %%%%%%%% シミュレータの基本情報 %%%%%%%
            obj.param.dt = 0.05;    % 刻み時間
            obj.param.Nt = 400;    % 計算するカウント数
            obj.param.Na = 20;       % エージェント数
            obj.param.space_x = [-8 8]; % 空間のサイズを定義
            obj.param.space_y = [-8 8];
            %%%%%%%% システムパラメータ %%%%%%%%
            % obj.param.K = 1;       % ゲイン
            obj.param.rv = 1.6;      % 観測範囲
            obj.param.stop_timehistry = 256;    % 停止検知用時間幅
            obj.param.stop_threshold = 0.01;    % 停止検知用閾値
            obj.param.adjacency_method = "distance";    % 隣接行列の作り方
            %%%%%%%% 読み込みファイル名 %%%%%%%%
            obj.param.environment_file = "setting_files/environments/narrow_space.m";  % 環境ファイル
            obj.param.placement_file = "setting_files/init_conditions/narrow_20.m";    % 初期位置ファイル
            %%%%%%%%%%%%%% 初期値 %%%%%%%%%%%%%
            obj.param.x_0 = rand(obj.param.Na, 2);
            obj.param.dxdt_0 = 0;%zeros(obj.param.Na, 2);
            obj.param.initial_pos_variance = 0;  % 初期位置の分散
        end
        
        function obj = initializeVariables(obj)
            % 各種変数を初期化．シミュレーションをやり直す度に必ず呼ぶこと
            % 状態変数の定義と初期値の代入を行うこと
            obj.t_vec = 0:obj.param.dt:obj.param.dt*(obj.param.Nt-1); % 時刻ベクトルの定義
            obj.x(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);    % 状態変数の定義
            obj.param.x_0 = obj.param.x_0 + (2*rand(obj.param.Na, 2)-1)*obj.param.initial_pos_variance;
            obj.x(:,:,1) = obj.param.x_0;   % 初期値の代入
            obj.dxdt(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);    % 状態変数の定義
            obj.dxdt(:,:,1) = obj.param.dxdt_0;   % 初期値の代入
            obj.u(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);    % 入力の履歴
            obj.variances(:,:,:) = ones(obj.param.Na, 1, obj.param.Nt);  % 位置の分散
            obj.is_stop(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);  % 停止判定
        end
        
        function obj = readSettingFiles(obj)
            % 環境や初期配置の初期設定ファイルを読み込む
            % initializeVariablesの前に読み込む！！
            %%%%%% 初期配置の読み込み %%%%%%
            run(obj.param.placement_file);
            obj.param.Na = posset_Na;               % posset_Naに台数
            obj.param.x_0 = [posset_x, posset_y];   % posset_x[台数,1],posset_y[台数,1]に初期位置
            %%%%%% 環境情報の読み込み %%%%%%
            run(obj.param.environment_file);
            obj.param.space_x = [envset_xmin, envset_xmax]; % envset_xmin,maxにx方向の端の値を指定
            obj.param.space_y = [envset_ymin, envset_ymax];
            obj.wall = permute(envset_wall_segments,[1,2,3]); % 障害物セグメント情報を読み取り
        end

        function obj = defineSystem(obj)
            % システムの定義が必要な場合はここ．シミュレーションをやり直すたびに呼ぶこと
            
        end


    %%%%%%%%%%%%%%%%%%%% 時間更新まわり %%%%%%%%%%%%%%%%%%

        function obj = simulate(obj)
            % シミュレーション本体
            disp("シミュレーションを開始します...")
            tic
            for t = 1:obj.param.Nt-1
                % ループ毎の更新をここに
                obj.G = obj.calcGraph(t);
                obj = obj.calcControlInput(t);  % 入力の計算
                obj.x(:,:,t+1) = obj.x(:,:,t) + obj.param.dt*obj.dxdt(:,:,t);   % オイラー法による更新
                obj.dxdt(:,:,t+1) = obj.dxdt(:,:,t) + obj.param.dt*obj.u(:,:,t);
            end % for
            toc
        end % simulate
        
        function G_ = calcGraph(obj,t)
            % 所定時刻におけるグラフを更新
            arguments
                obj
                t   % 時刻
            end
            X = repmat(obj.x(:,1,t),1,obj.param.Na);    % x座標を並べた行列
            Y = repmat(obj.x(:,2,t),1,obj.param.Na);    % y座標を並べた行列
            distances_ = (X-X.').^2 + (Y-Y.').^2;  % ユークリッド距離の２乗．X-X.'でx座標の差分が得られる
            if obj.param.adjacency_method == "delaunay"
                DT_ = delaunayTriangulation(obj.x(:,:,t));   % ドロネー分割の実施
                E_ = edges(DT_);                             % エッジを取得
                dist_ = distances_(sub2ind(size(distances_),E_(:,1),E_(:,2)));        % 各エッジのユークリッド距離
                E_enable_ = dist_<obj.param.rv^2;           % 観測距離よりも短いエッジを選択
                G_ = graph;
                G_ = addnode(G_,obj.param.Na);
                G_ = addedge(G_, E_(E_enable_,1),E_(E_enable_,2)); % エッジからグラフを生成
            else
                % 隣接行列はロボット間距離が観測範囲rvよりも小さいかどうか．対角要素は無視してグラフを作成
                G_ = graph(distances_<obj.param.rv^2, 'omitselfloops');
            end
        end

        function obj = calcControlInput(obj,t)
            % 入力の生成．継承して使う
            arguments
                obj
                t    % 時刻
            end
            disp("WARNIG : 継承前クラスのメソッドcalcControlInputが呼び出されている可能性があります")
            u_t = zeros(obj.param.Na, 2);   % 時刻tにおける入力
            %u_t(12,:) = [0.1,0.1];
            obj.u(:,:,t) = u_t;
        end

        function relative_vectors = calcVectorToWalls(obj,t)
            % 壁セグメントまでの相対位置
            % @return relative_vectors [ロボット数, 空間次元, 壁セグメント数]
            arguments
                obj
                t   % 時刻
            end
            Nwall_ = length(obj.wall); % 壁セグメントの数
            relative_vectors = zeros(obj.param.Na, 2, Nwall_);
            for p_ = 1:Nwall_   % セグメント毎に計算
                arg_ = atan2(obj.wall(2,2,p_)-obj.wall(1,2,p_),obj.wall(2,1,p_)-obj.wall(1,1,p_)); % セグメントのx軸からの偏角
                Rot_ = [cos(arg_), sin(arg_); -sin(arg_), cos(arg_)];   % 回転行列
                rotated_x_ = Rot_ * obj.x(:,:,t).'; % 回転後のエージェント座標．[空間次元, ロボット数]なので注意
                rotated_segment_ = Rot_ * obj.wall(:,:,p_).';  % 回転後のセグメント端点. y座標は同じになるはず．[空間次元, [開始/終了]]なので注意
                rv__ = zeros(2,obj.param.Na);   % 相対ベクトルを保持する中間変数．[空間次元,ロボット数]
                % 開始点のx座標よりもロボットのx座標が小さければ，開始点との相対位置
                rv__(:,:) = ( rotated_x_(1,:)<rotated_segment_(1,1) ).* (rotated_segment_(:,1)-rotated_x_(:,:));
                % 終了点のx座標よりもロボットのx座標が大きければ，終了点との相対位置
                rv__(:,:) = rv__(:,:) + ( rotated_x_(1,:)>rotated_segment_(1,2) ).* (rotated_segment_(:,2)-rotated_x_(:,:));
                % 終了点と開始点の間にロボットのx座標があるならば，y座標の差
                rv__(:,:) = rv__(:,:) + ( (rotated_x_(1,:)>=rotated_segment_(1,1)).*(rotated_x_(1,:)<=rotated_segment_(1,2)) ).*([rotated_x_(1,:);rotated_segment_(2,2)*ones(1,obj.param.Na)]-rotated_x_(:,:));
                
                relative_vectors(:,:,p_) = permute(Rot_.'*rv__,[2,1]);
            end
        end

        function obj = stopDetect(obj,t,is_debug)
            % 停止状態を検知
            arguments
                obj
                t                   % 時刻
                is_debug = false    % デバッグか？
            end
            if t < obj.param.stop_timehistry
                return  % データがたまっていなかったらリターン 迷いどころ
            end
            obj.variances(:,:,t) = 1/2*sum(var(obj.x(:,:,t-obj.param.stop_timehistry+1:t),0,3), 2) ; % 位置の分散を取得．x方向とy方向のやつを足し合わせ
            obj.is_stop(:,:,t) = obj.variances(:,:,t) < obj.param.stop_threshold;    % 位置の分散が閾値より低かったら停止状態と判断
        end

    %%%%%%%%%%%%%%%%%%%%% 解析まわり %%%%%%%%%%%%%%%%%%
    function obj = minimumDistanceCheck(obj, stepwidth, num)
        % 時刻歴の中での最小距離を出す．cbfのバリデーション用に
        % @brief 座標xは計算済みであることを前提
        arguments
            obj
            stepwidth = 1:obj.param.Nt  % 対象となる時刻歴
            num = 1:obj.param.Na        % 対象となるエージェントの集合
        end
        minimum_distance_agents = zeros(1,obj.param.Nt);    % エージェント間の最小距離
        minimum_distance_wall = zeros(1,obj.param.Nt);      % エージェント-壁間の最小距離
        for t = stepwidth
            Xa_ij(:,:,1) = repmat(obj.x(num,1,t),1,obj.param.Na)-repmat(obj.x(num,1,t),1,obj.param.Na).';   % 各方向相対位置ベクトル
            Xa_ij(:,:,2) = repmat(obj.x(num,2,t),1,obj.param.Na)-repmat(obj.x(num,2,t),1,obj.param.Na).';
            minimum_distance_agents(1,t) = min(vecnorm(Xa_ij,2,3)+100*eye(length(num)),[],'all');           % 対角に値を足した上で，最小距離を出す
            Xw_ij = obj.calcVectorToWalls(t);   % 壁までの相対位置ベクトル
            minimum_distance_wall(1,t) = min(vecnorm(Xw_ij(num,:,:),2,2),[],'all'); % 最小距離
        end
        figure
        subplot(2,1,1)
        plot(stepwidth,minimum_distance_agents);
        ylabel("between Robots")
        subplot(2,1,2)
        plot(stepwidth,minimum_distance_wall);
        ylabel("between Robot and Walls")
        xlabel("Timestep")
    end

    %%%%%%%%%%%%%%%%%%%%% 描画まわり %%%%%%%%%%%%%%%%%%

        function obj = placePlot(obj, t, view_edge, val)
            % ロボットの位置プロット
            arguments
                obj
                t               % 時刻
                view_edge = false                    % エッジ表示するか？
                val = zeros(obj.param.Na,1)         % ロボットに特徴づける値．これに応じてロボットの色が変わる
            end
            if view_edge    % エッジ表示するなら
                obj.showEdges(t);
            end
            obj.showWalls();    % 壁表示
            hold on
            scatter(obj.x(:,1,t),obj.x(:,2,t),120,val,'filled','MarkerEdgeColor','k'); % 散布図表示
            xlim(obj.param.space_x);    % 描画範囲を決定
            ylim(obj.param.space_y);
            %pbaspect([1 1 1])             % 縦横のアスペクト比を合わせる
            daspect([1 1 1])
            colormap(gca,"cool")
            hold on
        end

        function obj = numberPlacePlot(obj, t, view_edge)
            % ロボットの位置とロボットナンバーを描画
            arguments
                obj
                t               % 時刻
                view_edge = false                    % エッジ表示するか？
            end
            obj.placePlot(t,view_edge);
            text(obj.x(:,1,t)-0,obj.x(:,2,t)-0.5,string(1:obj.param.Na),'FontSize',12);
        end

        function obj = showEdges(obj,t)
            G_ = obj.calcGraph(t);  % グラフ計算
            e_ = table2array(G_.Edges);          % エッジ取得
            x_ = obj.x(:,1,t); 
            y_ = obj.x(:,2,t); 
            line(x_(e_).', y_(e_).','Color',"#0072BD",'LineWidth',1); % エッジの描画
        end
        
        function obj = showWalls(obj)
            % 壁の描画
            if isempty(obj.wall)    % 壁無かったら描画しない
                return
            end
            line(permute(obj.wall(:,1,:),[1,3,2]), permute(obj.wall(:,2,:),[1,3,2]), 'Color',"k",'LineWidth',1)
        end
        
        function obj = variancePlot(obj,num)
            % 位置の分散の時刻プロット
            arguments
                obj
                num = 8 % 表示対象のエージェント
            end
            
            figure
            plot(1:obj.param.Nt, permute(obj.variances(num,1,:),[3,1,2]))
            l = legend(string(num));
            l.NumColumns = 4;
            %ylim([-100,100])
            xlim([0,obj.param.Nt])
            ylabel("Variance of Position m^2]")
            xlabel("TIme Step")
        end

        function obj = moviePlot(obj,t)
            % 動画用プロット
            delete(gca)
            obj.placePlot(t,true);
            text(obj.param.space_x(2)*0.7, obj.param.space_y(2)*0.8, "t = "+string(t), 'FontSize',12);
            hold off
        end
        
        function obj = subplot(obj)
            % 一括描画用の最低限プロット
            % plot(obj.t_vec, obj.x(:,:));
            % xlabel("時刻 t [s]")
            % legend(["1","2","3"])
        end
        
        function obj = generateMovie(obj, filename, speed)
            arguments
                obj
                filename string = "movie.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
            obj.makeMovie(@obj.moviePlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end

    end % methods
end