classdef CollisionAvoidanceCBF
    %COLLISIONAVOIDANCECBF 衝突回避用CBFによる入力制約を行うクラス
    %   2次2d系を想定したCBFの実施
    
    properties
        m               % 質量
        rs             % 安全距離（定数．１次のクラスK関数）
        Tc             % 制御周期 s
        gamma           % ナイーブ定数（2次のクラスK関数）
        devide = false  % 分割化されているCBFか？ true なら Lghの前の係数が2に
        A               % 非線形制約の左辺の行列(Au \leq b) [制約数, 空間次元]
        b               % 非線形制約の右辺のベクトル [制約数, 1]
        options         % 最適化オプション
        enable_cbf = true   % CBFを使うか？（falseにするのはデバッグのとき）
    end
    
    methods
        function obj = CollisionAvoidanceCBF()
            %COLLISIONAVOIDANCECBF このクラスのインスタンスを作成
            obj.m = 1;
            obj.rs = 1;
            obj.Tc = 1;
            obj.gamma = 1;
            obj.enable_cbf = true;
            obj.devide = false;
            obj = obj.clearConstraints();   % Aとbをクリア
            obj.options = optimoptions("quadprog","Display","off"); % ２次計画法のコマンドウィンドウ出力をなしに
        end

        function obj = setParameters(obj,m_,rs_,Tc_,gamma_,devide_)
            %   各種パラメータの値をセット
            arguments
                obj
                m_          % 質量
                rs_         % 安全距離
                Tc_         % 制御周期
                gamma_      % ナイーブパラメタ
                devide_ = false    % 分割CBFか？
            end
            obj.m = m_;
            obj.rs = rs_;
            obj.Tc = Tc_;
            obj.gamma = gamma_;
            obj.devide = devide_;
        end
        
        function obj = enable(obj)
            % CBF制約の適用を有効化
            obj.enable_cbf = true;
        end

        function obj = disable(obj)
            % CBF制約の適用を無向化
            obj.enable_cbf = false;
        end

        function obj = addConstraints(obj,x_,dotx_)
            arguments
                obj
                x_      % 相対位置 [制約数, 空間次元]
                dotx_   % 相対速度 [制約数, 空間次元]
            end
            xnorm = vecnorm(x_,2,2);
            Lfh = sum(dotx_.*x_,2)./xnorm + (dotx_(:,1).*x_(:,2)-dotx_(:,2).*x_(:,1)).^2./(xnorm.^3)*obj.Tc;    % [制約数,1]ベクトル
            Lgh = x_./xnorm*obj.Tc; % [制約数,2]ベクトル
            gammah = obj.gamma*(xnorm + sum(dotx_.*x_,2)./xnorm*obj.Tc - obj.rs);   % [制約数,1]ベクトル
            if obj.devide == true
                obj.A = [obj.A; 2*Lgh];  % 分割してるなら２倍．A,bのサイズがループ毎に変わるので計算時間的には改善の余地あり
            else
                obj.A = [obj.A; Lgh];    % 分割してなければそのまま
            end
            obj.b = [obj.b; Lfh + gammah];
        end

        function obj = clearConstraints(obj)
            % 制約不等式をクリア
            obj.A = [];
            obj.b = [];
        end

        function u = apply(obj,u_nominal)
            % ノミナル入力に対しCBF制約を適用
            % @return 制約後の入力 [1,空間次元]のベクトル
            arguments
                obj
                u_nominal   % ノミナル入力．[1,空間次元]のベクトル
            end
            if obj.enable_cbf == false
                u = u_nominal;          % CBF使わないとき
                return
            end
            dim = length(u_nominal);
            delta_u = zeros(dim,1); % u-\bar{u}ノミナル入力との差 [空間次元,1]ベクトル
            [delta_u,~,flag] = quadprog(eye(dim), zeros(dim,1), obj.A, obj.b-obj.A*u_nominal.',[],[],[],[],[],obj.options);    % A\delta u\leq b-A\hat{u}
            u = (u_nominal.' + delta_u).';
            %{  %　デバッグ
            if flag ~= 1
                 disp("CBF missied "+string(flag))
                u = 0*u_nominal;    % だめなら入力なしに
            else
                u = (u_nominal.' + delta_u).';  % CBF解あり
            end
            %}
        end
    end
end

