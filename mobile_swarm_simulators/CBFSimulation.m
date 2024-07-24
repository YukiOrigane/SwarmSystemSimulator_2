
classdef CBFSimulation < MobileRobots2dSimulator
    
    properties
        cbf             % CBFのインスタンス
        lambda_history_lower  % ラグランジュ乗数履歴 [台数,次元,時刻]
        lambda_history_upper  % ラグランジュ乗数履歴 [台数,次元,時刻]
    end

    methods
        %%%%%%% 初期設定まわり %%%%%%%
        function obj = CBFSimulation()
            % コンストラクタ（宣言時に呼ばれる）
            obj@MobileRobots2dSimulator();    % 親クラスのコンストラクタも呼び出す
            obj = obj.setDefaultParameters();       % パラメタのデフォルト値を設定
            obj.cbf = CollisionAvoidanceCBF();              % CBF
        end

        function obj = setDefaultParameters(obj)
            obj = obj.setDefaultParameters@MobileRobots2dSimulator();   % スーパークラス側の読み出し
            % CBF %
            obj.param.cbf_rs = 0.8;         % 安全距離
            obj.param.cbf_gamma = 5;        % ナイーブパラメタ
            obj.param.cbf_lb = [-2; -2];    % 入力下限      
            obj.param.cbf_ub = [2; 2];      % 入力上限
        end
        
        function obj = initializeVariables(obj)
            % 各種変数を初期化．シミュレーションをやり直す度に必ず呼ぶこと
            % 状態変数の定義と初期値の代入を行うこと
            obj = obj.initializeVariables@MobileRobots2dSimulator();   % スーパークラス側の読み出し
            obj.lambda_history_lower = zeros(obj.param.Na,2,obj.param.Nt);
            obj.lambda_history_upper = zeros(obj.param.Na,2,obj.param.Nt);
        end
        
        function obj = defineSystem(obj)

        end

        %%%%%%%% 時間更新 %%%%%%%%%
        function obj = calcControlInput(obj,t)
            % 入力の生成．継承して使う
            arguments
                obj
                t    % 時刻
            end
            obj.showSimulationTime(t);
            u_t = zeros(obj.param.Na, 2);   % 時刻tにおける入力
            u_nominal = zeros(obj.param.Na, 2); % CBFをかける前のノミナル入力
            Adj = full(adjacency(obj.G));   % 隣接行列

            %
            X = repmat(obj.x(:,1,t),1,obj.param.Na);  % 位置x
            Y = repmat(obj.x(:,2,t),1,obj.param.Na);  % 位置y
            X_ij = X.'-X;   % 相対位置 X(i,j) = x(j)-x(i)
            Y_ij = Y.'-Y;

            % ノミナル入力の決定
            u_nominal = zeros(obj.param.Na,2);%-obj.param.k*obj.x(:,:,t);

            %%%% CBF %%%%
            % 詳細はCollisionAvoidanceCBF.mを参照

            cbf_(obj.param.Na) = obj.cbf;
            lambda_lower_t = obj.lambda_history_lower(:,:,t);   % lambda_lower_tを初期化し，要素数を確定
            lambda_upper_t = obj.lambda_history_upper(:,:,t);
            param_ = obj.param;
            for i = 1:obj.param.Na
            %parfor i = 1:obj.param.Na
                % CBFの適用 %
                cbf_(i) = cbf_(i).setParameters(1,obj.param.cbf_rs,obj.param.dt,obj.param.cbf_gamma,false);
                cbf_(i) = cbf_(i).addConstraints(-obj.x(:,:,t), -obj.dxdt(:,:,t));
                cbf_(i) = cbf_(i).addInputMinMaxConstraint(param_.cbf_lb,param_.cbf_ub);
                [u_t(i,:),lambda_] = cbf_(i).apply(u_nominal(i,:));
                lambda_lower_t(i,:,1) = (lambda_.lower).';
                lambda_upper_t(i,:,1) = (lambda_.upper).';
                cbf_(i) = cbf_(i).clearConstraints();
            end
            obj.lambda_history_lower(:,:,t) = lambda_lower_t;
            obj.lambda_history_upper(:,:,t) = lambda_upper_t;
            obj.u(:,:,t) = u_t;
        end
        
        %%%%%% 描画 %%%%%%
        function obj = movementPlot(obj, t)
            % ロボットの位置プロット
            arguments
                obj
                t               % 時刻
            end
            hold on
            scatter(obj.x(:,1,t),obj.x(:,2,t),120,0,'filled','MarkerEdgeColor','k'); % 散布図表示
            plot(permute(obj.x(:,1,:),[3,2,1]),permute(obj.x(:,2,:),[3,2,1]),'LineWidth',0.8)
            rectangle(...
                'Position', [-1, -1, 2, 2],...
                'Curvature',[1 1],...
                'FaceColor','k')
            xline(1.2,'--')
            xlim(obj.param.space_x);    % 描画範囲を決定
            ylim(obj.param.space_y);
            daspect([1 1 1])
            colormap(gca,"cool")
            hold on
        end

    end
end