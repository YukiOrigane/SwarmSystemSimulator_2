classdef WaveInteractionSimulator < Simulator
    %WAVEINTERACTIONSIMULATOR このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        % システムの変数を記載
        t_vec     % 固有時刻ベクトル
        phi       % 位相 [台数,1,時刻]
        xi        % ラプラシアン固有モード量 [台数,1,時刻]
        sigma     % 固有値 [台数,1,時刻]
        %dphidt   % ロボット速さ [台数,1,時刻]
        %u         % 入力
        G         % グラフオブジェクト．MATLABのgraph参照
        dphi_touch  % ソフトタッチ用．最後にエッジが追加された時刻の位相値
        dphi_release% ソフトリリース用．最後にエッジが削除された時刻の位相差（ソフトタッチ込）
        L_pre       % 一時刻前のラプラシアン行列
        x        % エージェント座標
        phi_x     % 位相の方向微分 [台数,空間次元,時刻]
        is_edge  % 自身が端っこか？ [台数,空間次元,時刻]
        peaks       % ピークの大きさ [台数,モード数,時刻]
        peak_freqs  % ピークの位置 [台数,モード数,時刻]
        is_deadlock % 自身がデッドロック状態か判定 [台数,1,時刻]
        peak_variances_db  % ピークの分散 [台数,モード数,時刻]
        freq_variances  % ピーク周波数の分散 [台数,モード数,時刻]
        potential_energy    % 全ポテンシャルエネルギー [時刻]
        phase_variances % 位相値の分散（拡散相互作用の場合）[台数,1,時刻]
    end
    
    methods
        function obj = WaveInteractionSimulator()
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
            %%%%%%%% システムパラメータ %%%%%%%%
            % 振動子系そのもの %
            % obj.param.K = 1;       % ゲイン
            obj.param.kappa = 10;      % 結合強度
            rng(1)
            obj.param.omega_0 = 5*rand(obj.param.Na,1);      % 固有角速度
            obj.param.gamma = 0;        % 粘性
            obj.param.interaction_type = "wave";    % 相互作用の形
            obj.param.is_normalize = false;     % 相互作用の正規化をするか？
            obj.param.use_softtouch = true;     % ソフトタッチを使うか？
            obj.param.use_softrelease = false;     % ソフトリリースを使うか？
            obj.param.delta_release = 0.1;     % ソフトリリースの位相差下限値
            % 各種推定 %
            obj.param.do_estimate = false;
            obj.param.is_judge_continuous = false;  % 内外判定結果を連続量にするか？
            obj.param.time_histry = 2048;     % パワースペクトラムで，どれくらい前の時刻情報まで使うか？
            obj.param.minimum_store = 64;     % ここまでデータたまるまではスタートしない
            obj.param.power_threshold_dB = 10^-10;
            obj.param.prominence_threshold_dB = 1;
            obj.param.peak_memory_num = 2;   % ピーク情報を何次まで記録するか
            obj.param.power_variance_db = 10^-3; % デッドロック判定時のパワー分散閾値
            obj.param.freq_variance_hz = 10^-5;  % デッドロック判定時の周波数分散閾値
            obj.param.deadlock_stepwith = 100;  % デッドロック判定．何ステップ分の定常状態を要請するか？
            obj.param.deadlock_usepower = true; % 判定にパワーも使うか？
            obj.param.use_peak_matching = false;    % 推定にピークマッチングを導入するか？
            obj.param.number_of_matching_peaks = 3; % ピークマッチングを何点でとるか？
            obj.param.diffusion_phase_var_threshold = 1e-6;    % 拡散相互作用の場合の，デッドロック検出閾値
            %%%%%%%% 読み込みファイル名 %%%%%%%%
            %obj.param.environment_file = "setting_files/environments/narrow_space.m";  % 環境ファイル
            %obj.param.placement_file = "setting_files/init_conditions/narrow_20.m";    % 初期位置ファイル
            %%%%%%%%%%%%%% 初期値 %%%%%%%%%%%%%
            obj.param.x_0 = zeros(obj.param.Na, 2);
            obj.param.phi_0 = zeros(obj.param.Na, 1);
            obj.param.dphidt_0 = zeros(obj.param.Na, 1);
            %obj.param.dxdt_0 = zeros(obj.param.Na, 2);
        end
        
        function obj = initializeVariables(obj)
            % 各種変数を初期化．シミュレーションをやり直す度に必ず呼ぶこと
            % 状態変数の定義と初期値の代入を行うこと
            if obj.param.use_peak_matching
                obj.param.peak_memory_num = obj.param.number_of_matching_peaks;  % ピークマッチングを使用する場合，3つモードを記録しておく
                disp("ピークマッチングを使うため，ピークは"+string(obj.param.number_of_matching_peaks)+"つ保存されます")
            end
            obj.t_vec = 0:obj.param.dt:obj.param.dt*(obj.param.Nt-1); % 時刻ベクトルの定義
            obj.phi(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);    % 状態変数の定義
            obj.xi(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);    % 状態変数の定義
            obj.sigma(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);    % 状態変数の定義
            obj.phi_x(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);    % 状態変数の定義
            obj.phi(:,:,1) = obj.param.phi_0;   % 初期値の代入
            obj.x(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);    % 状態変数の定義
            obj.is_edge(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);  % 内外変数
            obj.is_deadlock(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);  % デッドロック判定
            obj.peaks(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % ピークの大きさ
            obj.peak_freqs(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % ピークの大きさ
            obj.peak_variances_db(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % ピークの分散
            obj.freq_variances(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % ピーク位置の分散
            obj.potential_energy(:,1) = zeros(obj.param.Nt,1);  % ポテンシャルエネルギー
            obj.dphi_touch(:,:) = zeros(obj.param.Na,obj.param.Na); % 最終接続時の位相差
            obj.dphi_release(:,:) = zeros(obj.param.Na,obj.param.Na); % 最終切断時の位相差
            obj.L_pre(:,:) = zeros(obj.param.Na,obj.param.Na);
            obj.phase_variances(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);    % 各エージェントの位相の分散
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
                obj.x(:,:,t+1) = obj.x(:,:,t);  % 単発で動かすときは基本位置変えない
            end % for
            toc
        end % simulate

        function obj = stepSimulate(obj,t)
            % 他のシミュレーションに組み込む用．1ステップ分だけ更新を行う
            % @brief グラフラプラシアンを使うので，事前にsetGraph等でグラフ構造を与えておくこと．
            arguments
                obj
                t   % 時刻
            end
            L = full(laplacian(obj.G));
            if (obj.param.interaction_type == "wave")
                %%% 振動的相互作用 %%%
                if(t>2)
                    if obj.param.is_normalize   % 相互作用を正規化するか？
                        obj.phi(:,:,t+1) = 1/(1+obj.param.gamma/2*obj.param.dt)*(2*obj.phi(:,1,t) ...
                        +(obj.param.gamma/2*obj.param.dt-1)*obj.phi(:,1,t-1) ...
                        -obj.param.kappa./(degree(obj.G)+(degree(obj.G)==0)).*obj.param.dt^2.*L*obj.phi(:,1,t)); % 陽解法によるステップ更新
                    else    % しない
                        obj.phi(:,:,t+1) = 1/(1+obj.param.gamma/2*obj.param.dt)*(2*obj.phi(:,1,t) ...
                        +(obj.param.gamma/2*obj.param.dt-1)*obj.phi(:,1,t-1) ...
                        -obj.param.kappa*obj.param.dt^2*L*obj.phi(:,1,t)); % 陽解法によるステップ更新
                    end
                    Phi_i = repmat(obj.phi(:,:,t),1,obj.param.Na);  % 位相ベクトルを横に並べたもの．行列表現のphi_i相当
                    Phi_j = Phi_i.';    % 位相ベクトルを縦に並べたもの．行列表現のphi_j相当
                    if obj.param.use_softtouch  % ソフトタッチ
                        deltaL = L-obj.L_pre;   % Lの差分
                        deltaL_mask_all = ( deltaL - diag(diag(deltaL)) ) ~= 0; % deltaLのエッジ変更部分をマスク
                        deltaL_mask_add = ( deltaL - diag(diag(deltaL)) ) < 0;  % deltaLのエッジ追加部分をマスク
                        obj.dphi_touch = obj.dphi_touch-obj.dphi_touch.*deltaL_mask_all;    % エッジが変化した部分のdphi_touchをゼロに
                        obj.dphi_touch = obj.dphi_touch + (Phi_j-Phi_i).*deltaL_mask_add;   % エッジ追加時の位相差を記録
                        dphi_touch_memory = obj.dphi_touch;
                        obj.phi(:,:,t+1) = obj.phi(:,:,t+1) - obj.param.kappa*obj.param.dt^2*sum(obj.dphi_touch,2);
                        % エネルギー計算
                        phi_diff = (Phi_j-Phi_i);
                        obj.potential_energy(t) = 1/2*obj.param.kappa*sum((full(adjacency(obj.G)).*(phi_diff-obj.dphi_touch)).^2,'all')/2;
                    else
                        % エネルギー計算
                        obj.potential_energy(t) = 1/2*obj.param.kappa*obj.phi(:,1,t).'*L*obj.phi(:,1,t);
                    end
                    if obj.param.use_softrelease    % ソフトリリース
                        deltaL = L-obj.L_pre;   % Lの差分
                        deltaL_mask_cut = ( deltaL - diag(diag(deltaL)) ) > 0;% deltaLのエッジ削除部分をマスク
                        releasing_pair = (abs(obj.dphi_release-Phi_i).*(obj.dphi_release~=0)) > obj.param.delta_release;
                        obj.dphi_release = obj.dphi_release.*releasing_pair;    % 位相差が運動エネルギーになったらゼロに
                        releasing_pair = releasing_pair | deltaL_mask_cut;
                        obj.dphi_release = obj.dphi_release + (Phi_j-dphi_touch_memory).*deltaL_mask_cut;
                        obj.phi(:,:,t+1) = obj.phi(:,:,t+1) + obj.param.kappa*obj.param.dt^2*sum((obj.dphi_release-Phi_i).*releasing_pair,2);
                        obj.potential_energy(t) = obj.potential_energy(t) + obj.param.kappa*sum(((obj.dphi_release-Phi_i).*releasing_pair).^2,'all')/2;
                    end
                    if obj.param.do_estimate == true % 推定の実施
                        % xがsetされていることを要確認
                        obj = obj.calcPartialDerivative(t); % 位相の空間微分の計算
                        obj = obj.relativePositionEstimate(t);  % 相対位置推定
                        obj.showSimulationTime(t);
                    end
                else
                    obj.phi(:,:,t+1) = obj.phi(:,:,t);
                end
                obj.L_pre = L;
            elseif obj.param.interaction_type == "diffusion"
                %%% 拡散相互作用 %%%
                obj.phi(:,:,t+1) = obj.phi(:,:,t) + obj.param.dt*(obj.param.omega_0 ...
                    -obj.param.kappa*full(laplacian(obj.G))*obj.phi(:,:,t));
                obj = obj.judgeDeadlockDiffusion(t);
            end
            
            %%% 固有モード関連の計算 %%%
            if obj.param.is_normalize
                D = diag(diag(L));
                L = pinv(D)*L;
                [P_,D_] = eig(L);
            else
                [P_,D_] = eig(full(laplacian(obj.G)));
            end
            obj.xi(:,1,t) = P_.'*obj.phi(:,1,t);
            obj.sigma(:,1,t) = diag(D_);
        end
        
        function obj = calcPartialDerivative(obj,t)
            % 位相変数の空間微分の計算
            % 座標をsetPosition関数などで与えておくこと
            X_ij = repmat(obj.x(:,1,t).',obj.param.Na,1) - repmat(obj.x(:,1,t),1,obj.param.Na); % x方向相対ベクトル
            Y_ij = repmat(obj.x(:,2,t).',obj.param.Na,1) - repmat(obj.x(:,2,t),1,obj.param.Na); % y方向相対ベクトル
            phi_ij = repmat(obj.phi(:,1,t).',obj.param.Na,1) - repmat(obj.phi(:,1,t),1,obj.param.Na);
            R_ij = sqrt(X_ij.^2+Y_ij.^2) + eye(obj.param.Na);   % 相対距離ベクトル
            invR_ij = full(adjacency(obj.G)).*(1./R_ij);
            obj.phi_x(:,1,t) = sum(X_ij.*invR_ij.*phi_ij,2)./sum(abs(X_ij).*invR_ij,2);
            obj.phi_x(:,2,t) = sum(Y_ij.*invR_ij.*phi_ij,2)./sum(abs(Y_ij).*invR_ij,2);
            % 以下３行でphi_xのNaN要素を0に変換
            phi_x_ = obj.phi_x(:,:,t);
            phi_x_(isnan(phi_x_)) = 0;
            obj.phi_x(:,:,t) = phi_x_;
        end

        function obj = setGraph(obj, G_)
            % 外部で作ったグラフを与える
            arguments
                obj
                G_  % グラフオブジェクト
            end
            obj.G = G_;
        end

        function obj = setPosition(obj,x_,t)
            % 外部で計算した座標を渡す
            obj.x(:,:,t) = x_;
        end

        function G_ = calcGraph(obj,t)
            % 所定時刻におけるグラフを更新
            arguments
                obj
                t   % 時刻
            end
            X = repmat(obj.x(:,1,t),1,obj.param.Na);    % x座標を並べた行列
            Y = repmat(obj.x(:,2,t),1,obj.param.Na);    % y座標を並べた行列
            distances = (X-X.').^2 + (Y-Y.').^2;  % ユークリッド距離の２乗．X-X.'でx座標の差分が得られる
            % 隣接行列はロボット間距離が観測範囲rvよりも小さいかどうか．対角要素は無視してグラフを作成
            G_ = graph(distances<obj.param.rv^2, 'omitselfloops');
        end

        %%%%%%%%%%%%%%%%%%%%% 解析まわり %%%%%%%%%%%%%%%%%%
        function obj = relativePositionEstimate(obj,t,debug_agents)
            % 相対位置推定を行う
            arguments
                obj
                t                   % 時刻
                debug_agents = [];  % デバッグ用の描画を行うエージェント集合．空ならデバッグ描画なし
            end
            
            if t<obj.param.minimum_store    % 蓄積データ少ない間は推定しない
                return
            end
            if t>obj.param.time_histry
                % 時刻が推定に使うデータ点数より多いかどうかで，使う時刻幅を変える
                t_start_ = t-obj.param.time_histry;
            else
                t_start_ = 1;
            end
            % 各位相情報に関するパワースペクトラム p_は [周波数,チャンネル]となっているので注意
            [p_,f_] = pspectrum(permute(obj.phi(:,1,t_start_:t),[3,1,2]), obj.t_vec(t_start_:t));
            [px_,fx_] = pspectrum(permute(obj.phi_x(:,1,t_start_:t),[3,1,2]), obj.t_vec(t_start_:t));
            [py_,fy_] = pspectrum(permute(obj.phi_x(:,2,t_start_:t),[3,1,2]), obj.t_vec(t_start_:t));
            for i = 1:obj.param.Na  % エージェント毎回し
                [~,peak_index] = findpeaks(10*log10(p_(:,i)),"MinPeakHeight",obj.param.power_threshold_dB, "MinPeakProminence",obj.param.prominence_threshold_dB);    % ピーク検出
                [~,peakx_index] = findpeaks(10*log10(px_(:,i)),"MinPeakHeight",obj.param.power_threshold_dB, "MinPeakProminence",obj.param.prominence_threshold_dB);    % ピーク検出
                [~,peaky_index] = findpeaks(10*log10(py_(:,i)),"MinPeakHeight",obj.param.power_threshold_dB, "MinPeakProminence",obj.param.prominence_threshold_dB);    % ピーク検出
                p_i_ = p_(:,i);
                peak = p_i_(peak_index);  % ピーク値取得
                % TODO : 自励の処理を書く
                if length(peak)<obj.param.peak_memory_num
                    n_ = length(peak);
                else
                    n_ = obj.param.peak_memory_num;
                end
                obj.peaks(i,1:n_,t) = peak(1:n_);
                obj.peak_freqs(i,1:n_,t) = f_(peak_index(1:n_));
                if isempty(peakx_index)  % ピークがemptyの場合は最低周波数でピーク0に
                    pksx_ = 0;
                    peakx_index = 1;
                end
                if isempty(peaky_index)
                    pksy_ = 0;
                    peaky_index = 1;
                end
                pxi_ = px_(:,i);    % 論理取り出しをするためにベクトルに
                pyi_ = py_(:,i);
                %fxi_ = fx_(:,i);
                %fyi_ = fy_(:,i);
                index_xwin = pxi_(peakx_index)>pyi_(peakx_index);   % xのピークの内，yの値より高かったもの
                index_ywin = pyi_(peaky_index)>pxi_(peaky_index);
                if(sum(index_xwin)==0) % x側のピークが勝てる位置がなかった
                    maxindexx_ = peakx_index(1);
                else
                    win_indexx_ = peakx_index(index_xwin);
                    maxindexx_ = win_indexx_(1);
                end
                if(sum(index_ywin)==0) % y側のピークが勝てる位置がなかった
                    maxindexy_ = peaky_index(1);
                else
                    win_indexy_ = peaky_index(index_ywin);
                    maxindexy_ = win_indexy_(1);
                end
                %[~, maxindexx_] = max(pksx_);    % 偏微分側の最大ピークインデックスを得る
                %[~, maxindexy_] = max(pksy_);    % 注意) 2次モードや3次モードが最大になってしまう場合，修正の必要がある
                obj.is_edge(i,1,t) = 0;
                obj.is_edge(i,2,t) = 0;
                if (obj.param.is_judge_continuous)  % 判定結果は連続？論理値？
                    obj.is_edge(i,1,t) = p_(maxindexx_,i)-px_(maxindexx_,i);
                    obj.is_edge(i,2,t) = p_(maxindexy_,i)-py_(maxindexy_,i);
                else
                    obj.is_edge(i,1,t) = p_(maxindexx_,i)>px_(maxindexx_,i);
                    obj.is_edge(i,2,t) = p_(maxindexy_,i)>py_(maxindexy_,i);
                end
                %obj.is_edge(i,1,t) = p_(peakx_index(maxindexx_),i)>px_(peakx_index(maxindexx_),i)*sqrt(obj.param.kappa)/(2*pi*fx_(peakx_index(maxindexx_)));   % 最大
                %obj.is_edge(i,2,t) = p_(peakx_index(maxindexy_),i)>py_(peakx_index(maxindexy_),i)*sqrt(obj.param.kappa)/(2*pi*fy_(peaky_index(maxindexy_)));
                % 補正項の詳細
                % ピーク周波数f[Hz]としてエージェント長l. \mu次モードについて l = \mu\sqrt{\kappa}/{2f}
                % 微分時に\pi/l倍されているはずなので，l/\pi = \um\sqrt{\kappa}/{2\pi
                % f}を描ければいいのではと．一旦\mu = 1

                if ismember(i,debug_agents) %デバッグ用描画
                    figure
                    plot(f_,10*log(p_(:,i)));
                    hold on
                    plot(fx_,10*log(px_(:,i)));
                    plot(fy_,10*log(py_(:,i)));
                    plot(fx_(maxindexx_)*ones(2,1),10*log([p_(maxindexx_,i); px_(maxindexx_,i)]),'o');
                    plot(fy_(maxindexy_)*ones(2,1),10*log([p_(maxindexy_,i); py_(maxindexy_,i)]),'o');
                    xlabel("周波数 Hz")
                    xlim([0,5])
                    ylabel("パワー dB")
                    legend("\phi","\phi_x","\phi_y","x方向判定位置","y方向判定位置")
                    title("i = "+string(i)+", l_x = "+string(sqrt(obj.param.kappa)/2/fx_(maxindexx_))+", l_y = " + string(sqrt(obj.param.kappa)/2/fy_(maxindexy_)));
                end
            end
            obj = obj.judgeDeadlock(t); % デッドロック判定
        end

        function obj = judgeDeadlock(obj,t)
            % deadlock判定
            % @brief is_deadlock変数に1か0を返す
            % @brief 時刻tにおけるpeakの計算後に呼び出すこと
            if t < obj.param.minimum_store+obj.param.deadlock_stepwith
                return  % データがたまっていなかったらリターン
            end
            %if t>700
            %    disp("debug")
            %end

            %%% 振動的相互作用 %%%
            peak_variances_ = zeros(obj.param.Na,obj.param.peak_memory_num);    % ピークの大きさの分散
            freq_variances_ = zeros(obj.param.Na,obj.param.peak_memory_num);    % ピークの位置の分散
            
            if obj.param.use_peak_matching  % 余計にピーク点を記録しておき，周波数の近いピークを同一ピークとみなす
                % 1個前と今の差分行列
                freq_Diff_ = repmat(permute((obj.peak_freqs(:,:,t)+(obj.peak_freqs(:,:,t)==0).*100),[2,3,1]),1,obj.param.number_of_matching_peaks,1)-repmat(permute(obj.peak_freqs(:,:,t-1),[3,2,1]),obj.param.number_of_matching_peaks,1,1);
                peak_Diff_ = repmat(permute(10*log10(obj.peaks(:,:,t)+(obj.peaks(:,:,t)==0).*10^-10),[2,3,1]),1,obj.param.number_of_matching_peaks,1)-repmat(permute(10*log10(obj.peaks(:,:,t-1)),[3,2,1]),obj.param.number_of_matching_peaks,1,1);
                [MF_,col_] = mink(freq_Diff_,1,2,"ComparisonMethod","abs");   % 各ピーク点について，1つ前と最も近い点を探す．列インデックスも取る
                row_ = repmat((1:obj.param.number_of_matching_peaks).', 1,1,obj.param.Na);
                oku_ = repmat(permute(1:obj.param.Na,[1,3,2]),obj.param.number_of_matching_peaks,1,1);
                linear_Index = sub2ind(size(freq_Diff_),row_,col_,oku_);
                MP_ = peak_Diff_(linear_Index);              % 該当するパワー差も抽出
                [MinkFreq,row_] = mink(MF_,obj.param.number_of_matching_peaks-1,1,"ComparisonMethod","abs");    % 近い順にピーク点のペアいくつかの周波数差分を抽出
                col_ = repmat(ones(obj.param.number_of_matching_peaks-1,1), 1,1,obj.param.Na);
                oku_ = repmat(permute(1:obj.param.Na,[1,3,2]),obj.param.number_of_matching_peaks-1,1,1);
                linear_Index = sub2ind(size(MP_),row_,col_,oku_);
                MinkPower = MP_(linear_Index);
                % ピーク変動の2乗和を計算
                peak_variances_ = permute(sum(MinkPower.^2,1)/(obj.param.number_of_matching_peaks-1),[3,1,2]);
                freq_variances_ = permute(sum(MinkFreq.^2,1)/(obj.param.number_of_matching_peaks-1),[3,1,2]);
                obj.peak_variances_db(:,1,t) = peak_variances_;
                obj.freq_variances(:,1,t) = freq_variances_;    % 分散とは書いてあるが，この場合変化量の二乗和で分散ではない
                pv_mean_ = mean(obj.peak_variances_db(:,:,t-obj.param.deadlock_stepwith+1:t),3);
                fv_mean_ = mean(obj.freq_variances(:,:,t-obj.param.deadlock_stepwith+1:t),3);
                if obj.param.deadlock_usepower == true
                    obj.is_deadlock(:,:,t) = (pv_mean_(:,1)<obj.param.power_variance_db).*(fv_mean_(:,1)<obj.param.freq_variance_hz);
                else
                    obj.is_deadlock(:,:,t) = fv_mean_(:,1)<obj.param.freq_variance_hz;
                end
            else
                peak_variances_ = var(10*log10(obj.peaks(:,:,t-obj.param.deadlock_stepwith+1:t)),0,3);   % 時刻に沿った分散を計算．N-1で正規化
                freq_variances_ = var(obj.peak_freqs(:,:,t-obj.param.deadlock_stepwith+1:t),0,3);
                if obj.param.deadlock_usepower == true
                    obj.is_deadlock(:,:,t) = prod(peak_variances_<obj.param.power_variance_db,2).*prod(freq_variances_<obj.param.freq_variance_hz,2);
                else
                    obj.is_deadlock(:,:,t) = prod(freq_variances_<obj.param.freq_variance_hz,2);
                end
                obj.peak_variances_db(:,:,t) = peak_variances_;
                obj.freq_variances(:,:,t) = freq_variances_;
            end
            %obj.is_deadlock(:,:,t) = 1; %% DEBUG CODE HERE !!!!
            % 各モードの大きさ，周波数について全ての分散が閾値を下回っていたら，デッドロックと判定
        end

        function obj = judgeDeadlockDiffusion(obj,t)
            if t < obj.param.minimum_store+obj.param.deadlock_stepwith
                return  % データがたまっていなかったらリターン
            end
                %%% 拡散的相互作用の場合 %%%
                %phase_variances_ = var(obj.phi(:,1,t-obj.param.minimum_store:t)-phi_dc_,0,3);
                phase_variances_ = var(diff(obj.phi(:,1,t-obj.param.minimum_store:t),1,3),0,3);
                obj.is_deadlock(:,:,t) = phase_variances_<obj.param.diffusion_phase_var_threshold;
                obj.phase_variances(:,:,t) = phase_variances_;
        end

        %%%%%%%%%%%%%%%%%%%%% 描画まわり %%%%%%%%%%%%%%%%%%

        function obj = plot(obj,is_power)
            % ロボットの位置プロット
            arguments
                obj
                is_power = false    % power形式のプロット
            end
            figure
            %plot(obj.t_vec, permute(obj.phi(:,1,:),[1,3,2]))
            plot(1:obj.param.Nt, permute(obj.phi(:,1,:),[1,3,2]))
            if is_power == true
             plot(1:obj.param.Nt, permute(10*log10(abs(obj.phi(:,1,:))),[1,3,2]))
            end
            obj.saveFigure(gcf, "phase");
        end

        function obj = partPlot(obj,t)
            % 動画用の，指定時刻までの位相をプロットする関数
            arguments
                obj
                t      % 時刻
            end
            if t == 1
                return
            end
            plot(1:t, permute(obj.phi(:,1,1:t),[1,3,2]))
            hold on
            xline(t,'r')
            hold off
            xlim([1,obj.param.Nt]);
        end
        
        function obj = phaseTimeVariancePlot(obj)
            arguments
                obj
            end
            
            plot(1:obj.param.Nt, permute(obj.phase_variances(:,1,:),[1,3,2]))
            xlim([1,obj.param.Nt]);
        end

        function obj = energyPlot(obj,t)
            % 動画用の，指定時刻までのエネルギーをプロットする関数
            arguments
                obj
                t = obj.param.Nt      % 時刻
            end
            if t == 1
                return
            end
            phi_dot = gradient(permute(obj.phi(:,:,1:t),[1,3,2]))/obj.param.dt;     % 中心差分．両端だけ片側差分
            kinetic_energy = 1/2*phi_dot.^2;
            plot(1:t,sum(kinetic_energy))
            hold on
            plot(1:t,obj.potential_energy(1:t))
            plot(1:t,obj.potential_energy(1:t).'+sum(kinetic_energy))
            hold on
            xline(t,'r')
            hold off
            xlim([1,obj.param.Nt]);
            legend(["$K$","$U$","$K+U$"],'Interpreter','latex')
            xlabel("timestep")
            ylim([0 300])
            ylabel("Energy")
        end

        function obj = spectrumPlot(obj,t,view_eigen,num)
            % 指定エージェントのスペクトラムを描画
            arguments
                obj
                t       % 時刻
                view_eigen = true   % 固有値ベースの真値を描画するか？
                num = [9,10]    % エージェント番号
            end
            if t<obj.param.minimum_store    % 蓄積データ少ない間は推定しない
                return
            end
            if t>obj.param.time_histry
                % 時刻が推定に使うデータ点数より多いかどうかで，使う時刻幅を変える
                t_start_ = t-obj.param.time_histry;
            else
                t_start_ = 1;
            end
            % 各位相情報に関するパワースペクトラム p_は [周波数,チャンネル]となっているので注意
            [p,f] = pspectrum(permute(obj.phi(num,1,t_start_:t),[3,1,2]), obj.t_vec(t_start_:t));
            plot(f,10*log10(p));
            hold on
            if view_eigen   % 固有値に基づく真値の描画
                xline(sqrt(abs(obj.param.kappa*permute(obj.sigma(2:3,1,t),[3,1,2])))/2/pi,'--k',"$f_"+string((2:3)-1)+"$",'Interpreter','latex','LineWidth',0.5,'FontSize',14)
            end
            for mu = 1:obj.param.peak_memory_num
                plot(obj.peak_freqs(num,mu,t),10*log10(obj.peaks(num,mu,t)),'o');
            end
            hold off
            text(max(f)*0.7, 0, "t = "+string(t), 'FontSize',12);
            ylim([-100,20])
            xlim([0,10])
            legend(string(num))
        end

        function phaseSpacePlot(obj,t_end,mu_list,tile_layout,basis,G_)
            % 1ステップ分の相平面プロット
            arguments
                obj
                t_end = obj.param.Nt     % 時刻
                mu_list = 0:3             % モード
                tile_layout = [2,2]     % 描画タイルの配置
                basis {mustBeMember(basis,["eigen","static"])} = "eigen"    % 基底の取り方
                G_ = graph()          % 使うグラフ
            end
            t_vec_ = 1:t_end;
            E = permute(obj.xi(:,1,:),[3,1,2]); %．[時刻,モード,1]なので注意
            K = zeros(obj.param.Nt, obj.param.Na);
            K(2:end,:) = permute(((obj.xi(:,1,2:end)-obj.xi(:,1,1:end-1))./obj.param.dt),[3,1,2]);

            mu_cnt = 1; % 描画用インデックス
            for mu = mu_list
                subplot(tile_layout(1),tile_layout(2),mu_cnt)
                plot(E(t_vec_,mu+1), K(t_vec_,mu+1))
                hold on
                plot(E(t_vec_(end),mu+1), K(t_vec_(end),mu+1),'*')
                hold off
                mu_cnt = mu_cnt + 1;
                title("$\mu = "+string(mu)+", t = "+string(t_end)+"$",'Interpreter','latex')
                xlabel("$\xi_"+string(mu)+"(t)$",'Interpreter','latex')
                ylabel("$\dot{\xi}_"+string(mu)+"(t)$",'Interpreter','latex')
            end
            %obj.saveFigure(gcf, "phase_space");
        end

        function peakAndFreqPlot(obj,num)
            % 特定エージェントのピーク及びピーク周波数の時刻履歴を，【エージェント毎に】プロット
            arguments
                obj
                num = 8 % 表示対象のエージェント
            end
            figure
            for i = 1:length(num)
                subplot(length(num),1,i)
                plot(1:obj.param.Nt, permute(10*log10(obj.peaks(num(i),:,:)),[2,3,1]))
                hold on
                legend(string(1:obj.param.peak_memory_num))
                ylim([-100,100])
                %xlim([0,1000])
                ylabel("Power of Peaks [dB]")
                xlabel("TIme Step")
                title("i="+string(num(i)))
            end
            figure
            for i = 1:length(num)
                subplot(length(num),1,i)
                plot(1:obj.param.Nt, permute(obj.peak_freqs(num(i),:,:),[2,3,1]))
                legend(string(1:obj.param.peak_memory_num))
                %ylim([-100,100])
                %xlim([0,1000])
                ylabel("Frequency of Peaks [Hz]")
                xlabel("TIme Step")
                title("i="+string(num(i)))
            end
        end

        function peakAndFreqPlot2(obj,num)
            % 特定エージェントのピーク及びピーク周波数の時刻履歴を，【ピーク毎に】プロット
            arguments
                obj
                num = 8 % 表示対象のエージェント
            end
            
            for mu = 1:obj.param.peak_memory_num
                figure
                plot(1:obj.param.Nt, permute(10*log10(obj.peaks(num,mu,:)),[3,1,2]))
                l = legend(string(num));
                l.NumColumns = 4;
                ylim([-100,100])
                %xlim([0,1000])
                ylabel("Power of Peaks [dB]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
            
            for mu = 1:obj.param.peak_memory_num
                figure
                plot(1:obj.param.Nt, permute(obj.peak_freqs(num,mu,:),[3,1,2]))
                l = legend(string(num));
                l.NumColumns = 4;
                %ylim([-100,100])
                %xlim([0,1000])
                ylabel("Frequency of Peaks [Hz]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
        end
        
        function obj = compareEstimateAndEigenFreq(obj,num,modes)
            % 特定エージェントのピーク周波数と，固有値から求まる真の空間周波数を比較
            arguments
                obj
                num = 8 % 表示対象のエージェント
                modes = 2:3 % 表示する固有モード
            end
            figure
            for mu = 1:obj.param.peak_memory_num
                figure
                plot(1:obj.param.Nt, permute(obj.peak_freqs(num,mu,:),[3,1,2]))
                hold on
                plot(1:obj.param.Nt, sqrt(obj.param.kappa*permute(obj.sigma(modes,1,:),[3,1,2]))/2/pi,'--')
                hold off
                leg_str = ["$f^*_{"+string(num)+","+string(mu)+"}$", "$\sqrt{\kappa\sigma_"+string(modes-1)+"}/2\pi$"];
                l = legend(leg_str,'Interpreter','latex');
                l.NumColumns = 4;
                %ylim([-100,100])
                %xlim([0,1000])
                ylabel("Frequency[Hz]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
            obj.saveFigure(gcf, "estimate_eigen");
        end

        function obj = deadlockPlot(obj,num)
            % デッドロック判定の時系列結果を表示
            arguments
                obj
                num = 8 % 表示対象のエージェント
            end
            figure
            plot(1:obj.param.Nt, permute(obj.is_deadlock(num,1,:),[3,1,2]))
            l = legend(string(num));
            l.NumColumns = 2;
            ylim([-0.1 1.1])
            xlim([0,obj.param.Nt])
            ylabel("is deadlock")
            xlabel("TIme Step")
            obj.saveFigure(gcf, "deadlockPlot");
        end

        function obj = variancePlot(obj,num)
            % ピークの大きさ及び分散の時刻プロット
            arguments
                obj
                num = 8 % 表示対象のエージェント
            end
            
            for mu = 1:obj.param.peak_memory_num
                figure
                plot(1:obj.param.Nt, permute(obj.peak_variances_db(num,mu,:),[3,1,2]))
                l = legend(string(num));
                l.NumColumns = 4;
                %ylim([-100,100])
                xlim([0,obj.param.Nt])
                ylabel("Variance of Peak Power [dB^2]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
            
            for mu = 1:obj.param.peak_memory_num
                figure
                plot(1:obj.param.Nt, permute(obj.freq_variances(num,mu,:),[3,1,2]))
                l = legend(string(num));
                l.NumColumns = 4;
                %ylim([-100,100])
                xlim([0,obj.param.Nt])
                ylabel("Variance of Peak Frequency [Hz^2]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
            obj.saveFigure(gcf, "variance");
        end

        function obj = meanVariancePlot(obj,num)
            arguments
                obj
                num = 8 % 表示対象のエージェント
            end
            pv_mean = zeros(obj.param.Na, obj.param.peak_memory_num-1, obj.param.Nt);
            fv_mean = zeros(obj.param.Na, obj.param.peak_memory_num-1, obj.param.Nt);
            for t = 1:obj.param.Nt
                if (t>obj.param.deadlock_stepwith)
                    pv_mean_(:,:,t) = mean(obj.peak_variances_db(:,:,t-obj.param.deadlock_stepwith+1:t),3);
                    fv_mean_(:,:,t) = mean(obj.freq_variances(:,:,t-obj.param.deadlock_stepwith+1:t),3);
                end
            end
            mu = 1;
            figure
            plot(1:obj.param.Nt, permute(pv_mean_(num,mu,:),[3,1,2]))
            %l = legend(string(num));
            xlim([0,obj.param.Nt])
            ylabel("Variance of Peak Power [dB^2]")
            xlabel("TIme Step")
            yline(obj.param.power_variance_db,'--r')
            l = legend([string(num),"threshold"]);
            title("mode "+string(mu))

            figure
            plot(1:obj.param.Nt, permute(fv_mean_(num,mu,:),[3,1,2]))
            xlim([0,obj.param.Nt])
            ylabel("Variance of Peak Frequency [Hz^2]")
            xlabel("TIme Step")
            yline(obj.param.freq_variance_hz,'--r')
            l = legend([string(num),"threshold"]);
            title("mode "+string(mu))
        end

        function obj = phaseGapPlot(obj)
            % ロボットの位置プロット
            arguments
                obj
            end
            figure
            plot(obj.t_vec, permute(obj.phi(:,1,:),[1,3,2])-mean( permute(obj.phi(:,1,:),[1,3,2]), 1 ))
        end

        function obj = phaseACPlot(obj)
            % 固有角速度分をカットした位相時間履歴
            arguments
                obj
            end
            t_vec_ = 0:1:obj.param.Nt-1;
            phi_dc_ = obj.param.omega_0.*permute(obj.t_vec,[1,3,2]);    % 固有角速度による上昇分をカット  
            figure
            %plot(t_vec_, permute(obj.phi(:,1,:)-phi_dc_,[1,3,2]));
            plot(t_vec_, [zeros(obj.param.Na,1),permute(diff(obj.phi(:,1,:),1,3),[1,3,2])]);
        end

        function obj = generatePhaseSpaceMovie(obj,filename, speed)
            arguments
                obj
                filename string = "phase_space.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
            obj.makeMovie(@obj.phaseSpacePlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end

        function obj = generateEnergyMovie(obj,filename, speed)
            arguments
                obj
                filename string = "energy.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
            obj.makeMovie(@obj.energyPlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end

        function obj = generateSpectrumMovie(obj,filename, speed)
            arguments
                obj
                filename string = "movie.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
            obj.makeMovie(@obj.spectrumPlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end

        function obj = generatePhaseMovie(obj,filename, speed)
            arguments
                obj
                filename string = "phase.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
            obj.makeMovie(@obj.partPlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end
    end
end

