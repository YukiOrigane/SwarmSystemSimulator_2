
classdef SwarmWithWaveInteractionSimulation < MobileRobots2dSimulator
    
    properties
        cbf             % CBFのインスタンス
        attract_field   % 誘導場生成クラスのインスタンス
        cos             % 結合振動子のインスタンス
        is_connected = true     % グラフが連結かどうか
        kp_adjust       % 勾配追従力の調整係数[台数,1,時刻]
        is_deadlock     % デッドロック状態か？ COSからもらったり停止検出を使ったり [台数,1,時刻]
    end

    methods
        %%%%%%% 初期設定まわり %%%%%%%
        function obj = SwarmWithWaveInteractionSimulation()
            % コンストラクタ（宣言時に呼ばれる）
            obj@MobileRobots2dSimulator();    % 親クラスのコンストラクタも呼び出す
            obj = obj.setDefaultParameters();       % パラメタのデフォルト値を設定
            obj.cbf = CollisionAvoidanceCBF();              % CBF
            obj.attract_field = AttractantFieldGenerator(); % 誘導場生成クラス
            obj.cos = WaveInteractionSimulator();   % COS
        end

        function obj = setDefaultParameters(obj)
            obj = obj.setDefaultParameters@MobileRobots2dSimulator();   % スーパークラス側の読み出し
            obj.param.kp = 1;   % 勾配追従力ゲイン
            obj.param.kf = 1;   % 群形成力ゲイン
            obj.param.rc = 1;   % 群形成力の平衡距離
            obj.param.kd = 1;   % 粘性項
            obj.param.attract_force_type = "field_xy";% x方向のみ誘導力の形式
            obj.param.is_debug_view = false;    % デバッグ表示をするか？
            % kp調整 %
            obj.param.deadlock_source = "cos";   % デッドロック判定のソースは？
            obj.param.do_kp_adjust = false; % デッドロック時のkp調整を実施？
            obj.param.kp_adjust_out = -1;  % デッドロック時外側
            obj.param.kp_adjust_in = 2.0;   % デッドロック時内側
            obj.param.adjust_stepwith = 100;% 最後のデッドロック発生から何ステップ調整を発動するか
            % CBF %
            obj.param.cbf_rs = 0.8;         % 安全距離
            obj.param.cbf_gamma = 5;        % ナイーブパラメタ
        end
        
        function obj = initializeVariables(obj)
            % 各種変数を初期化．シミュレーションをやり直す度に必ず呼ぶこと
            % 状態変数の定義と初期値の代入を行うこと
            obj = obj.initializeVariables@MobileRobots2dSimulator();   % スーパークラス側の読み出し
            obj.cos = obj.cos.setParam("phi_0",rand(obj.param.Na,1));
            obj.kp_adjust = ones(obj.param.Na,1,obj.param.Nt);          % kp補正係数．標準は1
            obj.is_deadlock = zeros(obj.param.Na,1,obj.param.Nt);
        end
        
        function obj = defineSystem(obj)
            % システムの定義が必要な場合はここ．シミュレーションをやり直すたびに呼ぶこと
            %%% 誘導場をつくる %%%
            obj.attract_field = obj.attract_field.setParam("environment_file",obj.param.environment_file);  % 環境情報を誘導場クラスに渡す
            obj.attract_field = obj.attract_field.readSettingFiles();   % 設定ファイル読み出し
            obj.attract_field = obj.attract_field.initializeVariables();
            obj.attract_field = obj.attract_field.defineSystem();
            obj.attract_field = obj.attract_field.simulate();           % simulateまで回す．以降，obj.attract_field.c(:,:,end)に誘導場情報がある
            %%% COSの各種設定 %%%
            obj.cos = obj.cos.setParam("dt",obj.param.dt);  % 基本設定を共有しておく
            obj.cos = obj.cos.setParam("Nt",obj.param.Nt);
            obj.cos = obj.cos.setParam("Na",obj.param.Na);
            obj.cos = obj.cos.setParam("interaction_type","wave");
            obj.cos = obj.cos.initializeVariables();
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
                obj.is_deadlock(:,1,t) = obj.cos.is_deadlock(:,1,t);    % COSによるデッドロック判定を利用
            elseif obj.param.deadlock_source == "stop"
                obj.is_deadlock(:,1,t) = obj.is_stop(:,1,t);            % 停止検知をデッドロック判定として利用
            end
            if obj.param.do_kp_adjust % デッドロックに基づくkp調整
                obj = obj.kpAdjust(t);
            end

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
                elseif obj.param.attract_force_type == "linear_fbxy"
                    u_p(i,:) = obj.param.kp*obj.kp_adjust(i,:,t)*( [1 0] - [0 1]*obj.x(i,2,t) );  % x方向への単位ベクトルと，y=0への線形FB
                end
                % u_p = kp( c_x ex + c_y ey )
            end
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
                quiver(obj.x(:,1,t),obj.x(:,2,t),u_p(:,1),u_p(:,2),1.5);    % 勾配追従力プロット
            end
            %%%% 連結性の判定 %%%%
            Lap_ = full(laplacian(obj.G));  % グラフラプラシアン
            [~,Sigma] = eig(Lap_);  % 固有値展開
            sigma_ = sort(diag(Sigma));
            if (sigma_(2) <= 10^-5)
                obj.is_connected = false;
            end
        end

        function obj = kpAdjust(obj,t)
            % デッドロック情報に基づくkp調整の実施
            arguments
                obj
                t
            end
            if t<obj.param.adjust_stepwith+1
                return  % データの蓄積不十分ならなにもしない
            end
            adjust_do_ = sum(obj.is_deadlock(:,1,t-obj.param.adjust_stepwith:t),3); % 時刻幅内の履歴にデッドロック状態があるか？(論理和のかわりにsum)
            for i = 1:obj.param.Na
                if adjust_do_(i)>0  % 調整する場合
                    if obj.is_deadlock(i,1,t)
                        % 今デッドロック状態なら，デッドロック状態に応じて調整
                        obj.kp_adjust(i,1,t) = obj.cos.is_edge(i,2,t)*obj.param.kp_adjust_out + ~obj.cos.is_edge(i,2,t)*obj.param.kp_adjust_in;
                    else
                        % 今はデッドロックでないのなら，前の調整値を利用
                        obj.kp_adjust(i,1,t) = obj.kp_adjust(i,1,t-1);
                    end
                else
                    obj.kp_adjust(i,1,t) = 1;   % 調整しないなら1
                end
            end
            % 調整するなら，out:kp_adjust=kp_adjust_out, in:kp_adjust=kp_adjust_in
            % 否定~の方が積.*より演算が先
        end
            
            %%%%%%%%%%%%%%%% 解析 %%%%%%%%%%%%%%%
            function n = obtainNumberOfPassedRobots(obj,narrow_end,t)
                % 通り抜けに成功したエージェント数
                arguments
                    obj
                    narrow_end = 0      % 狭所の終わり
                    t = obj.param.Nt    % 判定時刻
                end
                n = sum(obj.x(:,1,t)>narrow_end);
            end

            %%%%%%%%%%%%% 描画まわり %%%%%%%%%%%%%%%
        function obj = phasePlacePlot(obj, t)
            % ロボットの位置プロット
            arguments
                obj
                t               % 時刻
            end
            delete(gca)
            obj = obj.placePlot(t,true, sin( obj.cos.phi(:,1,t) ) );
            clim([-1,1])
            colorbar
            text(obj.param.space_x(2)*0.7, obj.param.space_y(2)*0.8, "t = "+string(t), 'FontSize',12);
            hold off
        end

        function obj = edgeDeadlockPlot(obj,t,dim,is_delete)
            % ロボットの相対位置推定+デッドロック判定プロット
            arguments
                obj
                t               % 時刻
                dim = 2         % 次元
                is_delete = true
            end
            if (is_delete)
                delete(gca)
            end
            obj = obj.placePlot(t,false, (-1+2*obj.cos.is_edge(:,dim,t)).*obj.cos.is_deadlock(:,1,t));
            clim([-1,1])
            colorbar
            text(obj.param.space_x(2)*0.65, obj.param.space_y(2)*0.8, "t = "+string(t), 'FontSize',12);
            hold off
        end

        function obj = kpAdjustPlot(obj,num)
            % kp_adjustの値をプロット
            arguments
                obj
                num = 8 % 表示対象のエージェント
            end
            figure
            plot(1:obj.param.Nt, permute(obj.kp_adjust(num,1,:),[3,1,2]))
            l = legend(string(num));
            l.NumColumns = 2;
            ylim([-1 2.1])
            xlim([0,1000])
            ylabel("k_p adjust")
            xlabel("TIme Step")
        end

        function obj = edgeJudgePlot(obj, t, dim)
            % ロボットの位置プロット
            arguments
                obj
                t               % 時刻
                dim = 2         % 次元
            end
            delete(gca)
            obj = obj.placePlot(t,true, obj.cos.is_edge(:,dim,t));
            clim([0,1])
            colorbar
            text(obj.param.space_x(2)*0.7, obj.param.space_y(2)*0.8, "t = "+string(t), 'FontSize',12);
            hold off
        end
        
        function obj = trajectryJudgePlot(obj, step_width, num)
            % 軌跡を指定した時間幅でプロット
            % kp_adjustを利用して色を変える
            arguments
                obj
                step_width            % 軌跡の描画範囲
                num = 1:obj.param.Na  % 描画対象のエージェント
            end
            %delete(gca)
            for t = step_width  % 軌跡の描画．色の都合でやむなくfoe文
                for i = num
                    l = line([obj.x(i,1,t-1) obj.x(i,1,t)],[obj.x(i,2,t-1) obj.x(i,2,t)]);
                    l.LineWidth = 0.7;
                    if obj.kp_adjust(i,1,t) == obj.param.kp_adjust_out
                        l.Color = "#FF0000";    % 外なら赤
                        l.LineWidth = 1.5;
                    elseif obj.kp_adjust(i,1,t) == obj.param.kp_adjust_in
                        l.Color = "#0000FF";    % 内なら青
                        l.LineWidth = 1.5;
                    else
                        l.Color = "#111111";    % それ以外なら黒
                    end
                end
            end
            hold on
            obj = obj.placePlot(step_width(end),false);  % ステップ幅の最終時刻におけるエージェント位置をプロット
            title("t=["+string(step_width(1)-1)+", "+string(step_width(end))+"] step")
        end
        
        function obj = deadlockDetectionPlot(obj,source,stepwidth,num)
            % 各エージェントのデッドロック検出時刻を提示
            arguments
                obj
                source {mustBeMember(source,["result","cos","stop"])} = "result"  % なにを表示？
                stepwidth = 1:obj.param.Nt  % 時刻幅
                num = 1:obj.param.Na        % エージェント集合
            end
            if source == "result"
                judge_ = obj.is_deadlock;
            elseif source == "cos"
                judge_ = obj.cos.is_deadlock;
            elseif source == "stop"
                judge_ = obj.is_stop;
            end
            figure
            %{
            for n = num
                l = line([stepwidth(1) stepwidth(end)],[n n]);
                l.LineWidth = 0.7;
                l.Color = "#CCCCCC";
                hold on
                t_deadlock = stepwidth(permute(judge_(n,1,stepwidth),[3,1,2])==1); % デッドロックと判定されていた時刻
                plot(t_deadlock, repmat(n,length(t_deadlock)),'o','Color',"#D95319",'MarkerFaceColor',"#D95319");
            end
            %}
            
            imagesc(permute(1-judge_(num,1,stepwidth),[1,3,2]));
            colormap("gray")
            
            hold off
            ylabel("Robots Number")
            xlabel("Timestep")
           
        end

        function obj = generateMovieEstimate(obj, filename, speed)
            % 相対位置判定のムービーを生成
            arguments
                obj
                filename string = "movie.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
             %obj.makeMovie(@obj.edgeJudgePlot, obj.param.dt, obj.param.Nt, filename, speed, true);
            obj.makeMovie(@obj.edgeDeadlockPlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end

        function obj = generateMovie(obj, filename, speed)
            arguments
                obj
                filename string = "movie.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
            obj.makeMovie(@obj.phasePlacePlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end
        
    end

end % clasdef