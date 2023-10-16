
classdef SwarmWithWaveInteractionSimple < SwarmWithWaveInteractionSimulation
    
    properties

    end

    methods
        %%%%%%% 初期設定まわり %%%%%%%
        function obj = SwarmWithWaveInteractionSimple()
            % コンストラクタ（宣言時に呼ばれる）
            obj@SwarmWithWaveInteractionSimulation();    % 親クラスのコンストラクタも呼び出す
            obj = obj.setDefaultParameters();       % パラメタのデフォルト値を設定
        end

        function obj = setDefaultParameters(obj)
            obj = obj.setDefaultParameters@SwarmWithWaveInteractionSimulation();   % スーパークラス側の読み出し
            obj.param.Nt = 100;
        end
        
        function obj = initializeVariables(obj)
            % 各種変数を初期化．シミュレーションをやり直す度に必ず呼ぶこと
            % 状態変数の定義と初期値の代入を行うこと
            obj = obj.initializeVariables@SwarmWithWaveInteractionSimulation();   % スーパークラス側の読み出し
        end
        
        function obj = defineSystem(obj)
            % システムの定義が必要な場合はここ．シミュレーションをやり直すたびに呼ぶこと
            obj = obj.defineSystem@SwarmWithWaveInteractionSimulation();   % スーパークラス側の読み出し
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

            %%%% COSの更新 %%%%
            obj.cos = obj.cos.setGraph(obj.G);  % グラフを渡す
            obj.cos = obj.cos.setPosition(obj.x(:,:,t),t);  % 位置を渡す
            obj.cos = obj.cos.stepSimulate(t);  % COS側の更新
            
            %%%% デッドロック判定とその利用 %%%%
            obj = obj.stopDetect(t);    % 停止検知
            if obj.param.deadlock_source == "cos"
                obj.is_deadlock(:,1,t) = obj.is_stop(:,1,t).*obj.cos.is_deadlock(:,1,t);    % COSによるデッドロック判定+停止条件
            elseif obj.param.deadlock_source == "stop"
                obj.is_deadlock(:,1,t) = obj.is_stop(:,1,t);            % 停止検知をデッドロック判定として利用
            end
            if obj.param.do_kp_adjust % デッドロックに基づくkp調整
                obj = obj.kpAdjust(t);
            end
            obj = obj.tripUpdate(t);    % trip状態の更新

            %%%% 勾配追従力 %%%%
            u_p = zeros(obj.param.Na, 2);
            pos_index = round( (obj.x(:,:,t)-repmat([obj.param.space_x(1) obj.param.space_y(1)],obj.param.Na,1))./obj.attract_field.param.dx ) + 1;
            % ロボットのポジションインデックスの計算 round( (x-x_min)/dx )+1 結果は[エージェント数,空間次元]
            for i = 1:obj.param.Na
                %u_p(i,:) = obj.param.kp*( obj.attract_field.cx(pos_index(i,1), pos_index(i,2))*[1 0] + obj.attract_field.cy(pos_index(i,1), pos_index(i,2))*[0 1] );
                cx = obj.attract_field.cx(pos_index(i,1), pos_index(i,2));
                cy = obj.attract_field.cy(pos_index(i,1), pos_index(i,2));
                if obj.param.attract_force_type == "field_xonly"
                    u_p(i,:) = obj.param.kp*obj.kp_adjust(i,:,t)*( cx*[1 0] )/norm([cx,cy]);   % 誘導場をx方向のみ利用
                elseif obj.param.attract_force_type == "field_xy"
                    u_p(i,:) = obj.param.kp*obj.kp_adjust(i,:,t)*( cx*[1 0] + cy*[0 1] )/norm([cx,cy]);  % 誘導場を利用
                elseif obj.param.attract_force_type == "linear_fbx"
                    u_p(i,:) = obj.param.kp*obj.kp_adjust(i,:,t)*( [1 0] );  % ただのx方向への単位ベクトル
                elseif obj.param.attract_force_type == "trip"
                    u_p(i,:) = obj.param.kp*( (obj.trip_state(i,1,t)==0)*[1 0] +...
                                                (obj.trip_state(i,1,t)==1)*[0 -1] +...
                                                (obj.trip_state(i,1,t)==2)*[-1 0] +...
                                                (obj.trip_state(i,1,t)==3)*[0 1]);  % [0,1,2,3] : [right, down, left, up]
                end
                % u_p = kp( c_x ex + c_y ey )
            end
            %u_p = zeros(obj.param.Na, 2);   %%%% 注意！！入力削除
            %%%% 群形成力 %%%%
            u_f = zeros(obj.param.Na, 2);   % 群形成力
            X = repmat(obj.x(:,1,t),1,obj.param.Na);  % 位置x
            Y = repmat(obj.x(:,2,t),1,obj.param.Na);  % 位置y
            X_ij = X.'-X;   % 相対位置 X(i,j) = x(j)-x(i)
            Y_ij = Y.'-Y;
            D_ij = sqrt(X_ij.^2+Y_ij.^2)/obj.param.rc + eye(obj.param.Na);   % 正規化相対距離 (零割防止のため対角に1がならぶ)
            u_f(:,1) = obj.param.kf*sum( full(adjacency(obj.G)).*(-X_ij/obj.param.rc).*(D_ij.^(-3)-D_ij.^(-2)).*exp(-D_ij) ,2);
            u_f(:,2) = obj.param.kf*sum( full(adjacency(obj.G)).*(-Y_ij/obj.param.rc).*(D_ij.^(-3)-D_ij.^(-2)).*exp(-D_ij) ,2);

            % ノミナル入力の決定
            u_nominal = u_p + u_f;

            %%%% CBF %%%%
            % 詳細はCollisionAvoidanceCBF.mを参照
            x_io = obj.calcVectorToWalls(t);    % 壁との相対位置ベクトル
            for i = 1:obj.param.Na
                % ロボット間衝突回避CBF %
                obj.cbf = obj.cbf.setParameters(1,obj.param.cbf_rs,obj.param.dt,obj.param.cbf_gamma,true);
                x_ij = obj.x(:,:,t) - obj.x(i,:,t);          % 相対位置ベクトル
                dxdt_ij = obj.dxdt(:,:,t) - obj.dxdt(i,:,t); % 相対速度ベクトル
                obj.cbf = obj.cbf.addConstraints([x_ij(Adj(:,i)==1,1), x_ij(Adj(:,i)==1,2)], [dxdt_ij(Adj(:,i)==1,1), dxdt_ij(Adj(:,i)==1,2)]);
                % 隣接ロボットとの相対ベクトルのみCBF制約として利用
                % 壁との衝突回避CBF %
                obj.cbf = obj.cbf.setParameters(1,obj.param.cbf_rs,obj.param.dt,obj.param.cbf_gamma,false);
                obj.cbf = obj.cbf.addConstraints(permute(x_io(i,:,:),[3,2,1]), -repmat(obj.dxdt(i,:,t),length(x_io(i,:,:)),1));
                % 壁との相対位置ベクトルと，自身の速度ベクトル(壁との相対速度ベクトル)をCBFに入れる
                % CBFの適用 %
                u_t(i,:) = obj.cbf.apply(u_nominal(i,:));
                obj.cbf = obj.cbf.clearConstraints();
            end

            %%%% 最終的な入力の生成 %%%%
            obj.u(:,:,t) = u_t - obj.param.kd*obj.dxdt(:,:,t);  % CBF後に粘性が入っている…

            %%%% デバッグ %%%%
            if obj.param.is_debug_view
                figure
                obj.placePlot(t,true);
                hold on
                quiver(obj.x(:,1,t),obj.x(:,2,t),u_p(:,1),u_p(:,2));    % 勾配追従力プロット
            end
            %%%% 連結性の判定 %%%%
            Lap_ = full(laplacian(obj.G));  % グラフラプラシアン
            [~,Sigma] = eig(Lap_);  % 固有値展開
            sigma_ = sort(diag(Sigma));
            %if (sigma_(2) <= 10^-5)        % 0303時点で使っていないので一度外す
            %    obj.is_connected = false;
            %end
        end

            %%%%%%%%%%%%% 描画まわり %%%%%%%%%%%%%%%
    end

end % clasdef