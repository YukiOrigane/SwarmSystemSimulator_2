classdef WaveInteractionSimulator < Simulator
    %WAVEINTERACTIONSIMULATOR このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        % システムの変数を記載
        t_vec     % 固有時刻ベクトル
        phi       % 位相 [台数,1,時刻]
        sigma     % 固有値 [台数,1,時刻]
        %dphidt   % ロボット速さ [台数,1,時刻]
        %u         % 入力
        G         % グラフオブジェクト．MATLABのgraph参照
        x        % エージェント座標
        phi_x     % 位相の方向微分 [台数,空間次元,時刻]
        is_edge  % 自身が端っこか？ [台数,空間次元,時刻]
        peaks       % ピークの大きさ [台数,モード数,時刻]
        peaks_x       % ピークの大きさ [台数,モード数,時刻]
        peaks_y       % ピークの大きさ [台数,モード数,時刻]
        peak_freqs  % ピークの位置 [台数,モード数,時刻]
        peak_x_freqs  % ピークの位置 [台数,モード数,時刻]
        peak_y_freqs  % ピークの位置 [台数,モード数,時刻]
        is_deadlock % 自身がデッドロック状態か判定 [台数,1,時刻]
        is_deadlock_variance % 分散による判定 [台数,3(phi,phi_x,phi_y),時刻]
        is_deadlock_periodic % 周期性による判定 [台数,3(phi,phi_x,phi_y),時刻]
        peak_variances_db  % ピークの分散 [台数,モード数,時刻]
        freq_variances  % ピーク周波数の分散 [台数,モード数,時刻]
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
            %%%%%%%% シミュレータの基本情報 %%%%%%%
            obj.param.dt = 0.05;    % 刻み時間
            obj.param.Nt = 400;    % 計算するカウント数
            obj.param.Na = 20;       % エージェント数
            %%%%%%%% システムパラメータ %%%%%%%%
            % 振動子系そのもの %
            % obj.param.K = 1;       % ゲイン
            obj.param.kappa = 10;      % 結合強度
            obj.param.omega_0 = [5; 5];      % 固有角速度
            obj.param.gamma = 0;        % 粘性
            obj.param.interaction_type = "wave";    % 相互作用の形
            % 各種推定 %
            obj.param.do_estimate = false;
            obj.param.is_judge_continuous = false;  % 内外判定結果を連続量にするか？
            obj.param.time_histry = 2048;     % パワースペクトラムで，どれくらい前の時刻情報まで使うか？
            obj.param.minimum_store = 64;     % ここまでデータたまるまではスタートしない
            obj.param.power_threshold = 10^-10;
            obj.param.peak_memory_num = 2;   % ピーク情報を何次まで記録するか
            obj.param.power_variance_db = 10^-3; % デッドロック判定時のパワー分散閾値
            obj.param.freq_variance_hz = 10^-5;  % デッドロック判定時の周波数分散閾値
            obj.param.deadlock_stepwith = 100;  % デッドロック判定．何ステップ分の定常状態を要請するか？
            obj.param.deadlock_stepwith_periodic = 512;  % デッドロック判定．周期性検出窓
            obj.param.periodic_coeff_threshold = 0.5;   % 周期性検出時の自己相関閾値
            obj.param.periodic_minimum_shift = 10;      % 周期性検出時に，これ以下のシフトは除外する
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
            obj.t_vec = 0:obj.param.dt:obj.param.dt*(obj.param.Nt-1); % 時刻ベクトルの定義
            obj.phi(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);    % 状態変数の定義
            obj.sigma(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);    % 状態変数の定義
            obj.phi_x(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);    % 状態変数の定義
            obj.phi(:,:,1) = obj.param.phi_0;   % 初期値の代入
            obj.x(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);    % 状態変数の定義
            obj.is_edge(:,:,:) = zeros(obj.param.Na, 2, obj.param.Nt);  % 内外変数
            obj.is_deadlock(:,:,:) = zeros(obj.param.Na, 1, obj.param.Nt);  % デッドロック判定
            obj.is_deadlock_variance(:,:,:) = zeros(obj.param.Na, 3, obj.param.Nt);
            obj.is_deadlock_periodic(:,:,:) = zeros(obj.param.Na, 3, obj.param.Nt);
            obj.peaks(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % ピークの大きさ 
            obj.peaks_x(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % x微分
            obj.peaks_y(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % y微分
            obj.peak_freqs(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    %
            obj.peak_x_freqs(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % 
            obj.peak_y_freqs(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % 
            obj.peak_variances_db(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % ピークの分散
            obj.freq_variances(:,:,:) = zeros(obj.param.Na, obj.param.peak_memory_num, obj.param.Nt);    % ピーク位置の分散
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
            if (obj.param.interaction_type == "wave")
                %%% 振動的相互作用 %%%
                if(t>2)
                    obj.phi(:,:,t+1) = 1/(1+obj.param.gamma/2*obj.param.dt)*(2*obj.phi(:,1,t) ...
                        +(obj.param.gamma/2*obj.param.dt-1)*obj.phi(:,1,t-1) ...
                        -obj.param.kappa*obj.param.dt^2*full(laplacian(obj.G))*obj.phi(:,1,t)); % 陽解法によるステップ更新
                    if obj.param.do_estimate == true % 推定の実施
                        % xがsetされていることを要確認
                        obj = obj.calcPartialDerivative(t); % 位相の空間微分の計算
                        obj = obj.relativePositionEstimate(t);  % 相対位置推定
                        obj.showSimulationTime(t);
                    end
                else
                    obj.phi(:,:,t+1) = obj.phi(:,:,t);
                end
            elseif obj.param.interaction_type == "diffusion"
                %%% 拡散相互作用 %%%
                obj.phi(:,:,t+1) = obj.phi(:,:,t) + obj.param.dt*(obj.param.omega_0 ...
                    -obj.param.kappa*full(laplacian(obj.G))*obj.phi(:,:,t));
            end
            [~,D_] = eig(full(laplacian(obj.G)));
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
                [peak,peak_index,n_] = obj.findFillPeaks(p_(:,i),obj.param);
                [peak_x,peak_x_index,nx_] = obj.findFillPeaks(px_(:,i),obj.param);
                [peak_y,peak_y_index,ny_] = obj.findFillPeaks(py_(:,i),obj.param);
                obj.peaks(i,1:n_,t) = peak(1:n_);
                obj.peaks_x(i,1:nx_,t) = peak_x(1:nx_);
                obj.peaks_y(i,1:ny_,t) = peak_y(1:ny_);
                obj.peak_freqs(i,1:n_,t) = f_(peak_index(1:n_));
                obj.peak_x_freqs(i,1:nx_,t) = fx_(peak_x_index(1:nx_));
                obj.peak_y_freqs(i,1:ny_,t) = fy_(peak_y_index(1:ny_));

                if isempty(peak_x_index)  % ピークがemptyの場合は最低周波数でピーク0に
                    peak_x_index = 1;
                end
                if isempty(peak_y_index)
                    peak_y_index = 1;
                end
                pxi_ = px_(:,i);    % 論理取り出しをするためにベクトルに
                pyi_ = py_(:,i);
                %fxi_ = fx_(:,i);
                %fyi_ = fy_(:,i);
                index_xwin = pxi_(peak_x_index)>pyi_(peak_x_index);   % xのピークの内，yの値より高かったもの
                index_ywin = pyi_(peak_y_index)>pxi_(peak_y_index);
                if(sum(index_xwin)==0) % x側のピークが勝てる位置がなかった
                    maxindexx_ = peak_x_index(1);
                else
                    win_indexx_ = peak_x_index(index_xwin);
                    maxindexx_ = win_indexx_(1);
                end
                if(sum(index_ywin)==0) % y側のピークが勝てる位置がなかった
                    maxindexy_ = peak_y_index(1);
                else
                    win_indexy_ = peak_y_index(index_ywin);
                    maxindexy_ = win_indexy_(1);
                end
                obj.is_edge(i,1,t) = 0;
                obj.is_edge(i,2,t) = 0;
                if (obj.param.is_judge_continuous)  % 判定結果は連続？論理値？
                    obj.is_edge(i,1,t) = p_(maxindexx_,i)-px_(maxindexx_,i);
                    obj.is_edge(i,2,t) = p_(maxindexy_,i)-py_(maxindexy_,i);
                else
                    obj.is_edge(i,1,t) = p_(maxindexx_,i)>px_(maxindexx_,i);
                    obj.is_edge(i,2,t) = p_(maxindexy_,i)>py_(maxindexy_,i);
                end
                %obj.is_edge(i,1,t) = p_(peak_x_index(maxindexx_),i)>px_(peak_x_index(maxindexx_),i)*sqrt(obj.param.kappa)/(2*pi*fx_(peak_x_index(maxindexx_)));   % 最大
                %obj.is_edge(i,2,t) = p_(peak_x_index(maxindexy_),i)>py_(peak_x_index(maxindexy_),i)*sqrt(obj.param.kappa)/(2*pi*fy_(peak_y_index(maxindexy_)));
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
            %obj = obj.judgeDeadlock(t); % デッドロック判定
            obj = obj.judgeDeadlockWithPeriodic(t); % デッドロック判定
        end

        function [peaks_,indeces_,n_] = findFillPeaks(~,p_,param_)
            % FFT結果に対してピーク検出を行う．十分な数のピークがなかったら0で埋める
            arguments
                ~
                p_  % FFT結果（パワー）
                param_
            end
            [peaks_,indeces_] = findpeaks(p_(:),"MinPeakHeight",param_.power_threshold);    % ピーク検出
            if length(peaks_)<param_.peak_memory_num
                n_ = length(peaks_);
            else
                n_ = param_.peak_memory_num;
            end
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
            peak_variances_ = var(10*log10(obj.peaks(:,:,t-obj.param.deadlock_stepwith+1:t)),0,3);   % 時刻に沿った分散を計算．N-1で正規化
            freq_variances_ = var(obj.peak_freqs(:,:,t-obj.param.deadlock_stepwith+1:t),0,3);
            obj.is_deadlock(:,:,t) = prod(peak_variances_<obj.param.power_variance_db,2).*prod(freq_variances_<obj.param.freq_variance_hz,2);
            obj.peak_variances_db(:,:,t) = peak_variances_;
            obj.freq_variances(:,:,t) = freq_variances_;
            % 各モードの大きさ，周波数について全ての分散が閾値を下回っていたら，デッドロックと判定
        end

        function obj = judgeDeadlockWithPeriodic(obj,t)
            % 周期性を加味したデッドロック判定
            % @brief is_deadlock変数に1か0を返す
            % @brief 時刻tにおけるpeakの計算後に呼び出すこと
            if t > obj.param.minimum_store+obj.param.deadlock_stepwith
                freq_variances_ = obj.calcPeakVariance(obj.peak_freqs,t,obj.param);  % 時刻に沿った分散を計算．N-1で正規化
                freq_x_variances_ = obj.calcPeakVariance(obj.peak_x_freqs,t,obj.param);
                freq_y_variances_ = obj.calcPeakVariance(obj.peak_y_freqs,t,obj.param);
                obj.is_deadlock_variance(:,1,t) = prod(freq_variances_<obj.param.freq_variance_hz,2);
                obj.is_deadlock_variance(:,2,t) = prod(freq_x_variances_<obj.param.freq_variance_hz,2);
                obj.is_deadlock_variance(:,3,t) = prod(freq_y_variances_<obj.param.freq_variance_hz,2);
            end
            if t > obj.param.minimum_store+obj.param.deadlock_stepwith_periodic
                for i = 1:obj.param.Na
%                     f_ = permute(obj.peak_freqs(i,1,t-obj.param.deadlock_stepwith_periodic+1:t),[3,1,2]);
%                     fx_ = permute(obj.peak_x_freqs(i,1,t-obj.param.deadlock_stepwith_periodic+1:t),[3,1,2]);
%                     fy_ = permute(obj.peak_y_freqs(i,1,t-obj.param.deadlock_stepwith_periodic+1:t),[3,1,2]);
%                     [c_,lags_] = xcorr(f_-mean(f_),'normalized');    % １次ピークに限定することに注意
%                     [cx_,lags_x_] = xcorr(fx_-mean(fx_),'normalized');
%                     [cy_,lags_y_] = xcorr(fy_-mean(fy_),'normalized');
%                     [corr_peak_,corr_loc_] = findpeaks(c_,lags_,'MinPeakProminence',obj.param.periodic_coeff_threshold/2);
%                     [corr_peak_x_,corr_loc_x_] = findpeaks(cx_,lags_x_,'MinPeakProminence',obj.param.periodic_coeff_threshold/2);
%                     [corr_peak_y_,corr_loc_y_] = findpeaks(cy_,lags_y_,'MinPeakProminence',obj.param.periodic_coeff_threshold/2);
%                     obj.is_deadlock_periodic(i,1,t) = max([corr_peak_(corr_loc_>obj.param.periodic_minimum_shift);0]) > obj.param.periodic_coeff_threshold;     % ピークがない場合の[]>0 = [] を避けるために，[[],0]>0 = 0 とした
%                     obj.is_deadlock_periodic(i,2,t) = max([corr_peak_x_(corr_loc_x_>obj.param.periodic_minimum_shift);0]) > obj.param.periodic_coeff_threshold;
%                     obj.is_deadlock_periodic(i,3,t) = max([corr_peak_y_(corr_loc_y_>obj.param.periodic_minimum_shift);0]) > obj.param.periodic_coeff_threshold;
                    obj.is_deadlock_periodic(i,1,t) = obj.calcMaxPeriodic(obj.peak_freqs, i, t, obj.param) > obj.param.periodic_coeff_threshold;
                    obj.is_deadlock_periodic(i,2,t) = obj.calcMaxPeriodic(obj.peak_x_freqs, i, t, obj.param) > obj.param.periodic_coeff_threshold;
                    obj.is_deadlock_periodic(i,3,t) = obj.calcMaxPeriodic(obj.peak_y_freqs, i, t, obj.param) > obj.param.periodic_coeff_threshold;
                end
            end
            obj.is_deadlock(:,:,t) = prod(obj.is_deadlock_variance(:,:,t),2) + sum(obj.is_deadlock_periodic(:,:,t),2);
        end

        function var_ = calcPeakVariance(~,freq_,t_,param_)
            % 1行だが，解析でも使うので関数化．
            var_ = var(freq_(:,:,t_-param_.deadlock_stepwith+1:t_),0,3);
        end

        function max_peak_ = calcMaxPeriodic(~,freq_,i_,t_,param_)
            arguments
                ~
                freq_
                i_
                t_
                param_
            end
            if t_ == 1200
                disp("1200")
            end
            f_ = permute(freq_(i_,1,t_-param_.deadlock_stepwith_periodic+1:t_),[3,1,2]);   % １次ピークに限定することに注意
            [c_,lags_] = xcorr(f_-mean(f_),'normalized');
            [corr_peak_,corr_loc_] = findpeaks(c_,lags_,'MinPeakProminence',param_.periodic_coeff_threshold/2);
            if (length(corr_peak_)>3)
                max_peak_ = 0;      % ピーク数多すぎたら外す
            else
                max_peak_ = max([corr_peak_(corr_loc_>param_.periodic_minimum_shift);0]);
            end
        end

        %%%%%%%%%%%%%%%%%%%%% 描画まわり %%%%%%%%%%%%%%%%%%

        function obj = plot(obj)
            % ロボットの位置プロット
            arguments
                obj
            end
            figure
            plot(obj.t_vec, permute(obj.phi(:,1,:),[1,3,2]))
        end

        function obj = spectrumPlot(obj,t,view_eigen,num)
            % 指定エージェントのスペクトラムを描画
            arguments
                obj
                t       % 時刻
                view_eigen = true; % 固有値に基づく真値をプロットするか？
                num = [24,32,40]    % エージェント番号
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
                xline(sqrt(abs(obj.param.kappa*permute(obj.sigma(2:5,1,t),[3,1,2])))/2/pi,'--k',"$f_"+string((2:5)-1)+"$",'Interpreter','latex','LineWidth',0.5,'FontSize',14)
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

        function spectrumPlotDiff(obj,t,view_eigen,num)
            % 指定エージェントのスペクトラムを，空間微分含めて描画
            arguments
                obj
                t       % 時刻
                view_eigen = true; % 固有値に基づく真値をプロットするか？
                num {mustBeNumeric} = 32;    % エージェント番号
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
            [px,fx] = pspectrum(permute(obj.phi_x(num,1,t_start_:t),[3,1,2]), obj.t_vec(t_start_:t));
            [py,fy] = pspectrum(permute(obj.phi_x(num,2,t_start_:t),[3,1,2]), obj.t_vec(t_start_:t));
            plot(f,10*log10(p));
            hold on
            plot(fx,10*log10(px));
            plot(fy,10*log10(py));
            if view_eigen   % 固有値に基づく真値の描画
                xline(sqrt(abs(obj.param.kappa*permute(obj.sigma(2:5,1,t),[3,1,2])))/2/pi,'--k',"$f_"+string((2:5)-1)+"$",'Interpreter','latex','LineWidth',0.5,'FontSize',14)
            end
            for mu = 1:obj.param.peak_memory_num
                plot(obj.peak_freqs(num,mu,t),10*log10(obj.peaks(num,mu,t)),'o');
                plot(obj.peak_x_freqs(num,mu,t),10*log10(obj.peaks_x(num,mu,t)),'o');
                plot(obj.peak_y_freqs(num,mu,t),10*log10(obj.peaks_y(num,mu,t)),'o');
            end
            hold off
            text(max(f)*0.7, 0, "t = "+string(t), 'FontSize',12);
            ylim([-100,20])
            xlim([0,10])
            legend(["\phi","\partial_x\phi","\partial_y\phi"])
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
                xlim([0,1000])
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
                xlim([0,1000])
                ylabel("Frequency of Peaks [Hz]")
                xlabel("TIme Step")
                title("i="+string(num(i)))
            end
        end

        function peakAndFreqPlot2(obj,num,diff_)
            % 特定エージェントのピーク及びピーク周波数の時刻履歴を，【ピーク毎に】プロット
            arguments
                obj
                num = 8 % 表示対象のエージェント
                diff_ {mustBeMember(diff_,["","x","y"])} = ""
            end
            
            if diff_ == "x"
                p_ = obj.peaks_x;
                f_ = obj.peak_x_freqs;
            elseif diff_ == "y"
                p_ = obj.peaks_y;
                f_ = obj.peak_y_freqs;
            else
                p_ = obj.peaks;
                f_ = obj.peak_freqs;
            end

            for mu = 1:obj.param.peak_memory_num
                figure
                plot(1:obj.param.Nt, permute(10*log10(p_(num,mu,:)),[3,1,2]))
                l = legend(string(num));
                l.NumColumns = 4;
                ylim([-100,100])
                xlim([0,1500])
                ylabel("Power of Peaks [dB]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
            
            for mu = 1:obj.param.peak_memory_num
                figure
                plot(1:obj.param.Nt, permute(f_(num,mu,:),[3,1,2]))
                l = legend(string(num));
                l.NumColumns = 4;
                %ylim([-100,100])
                xlim([0,1500])
                ylabel("Frequency of Peaks [Hz]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
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
            xlim([0,1000])
            ylabel("is deadlock")
            xlabel("TIme Step")
        end

        function obj = peakVariancePlot(obj,num,dim)
            arguments
                obj
                num         % 表示対象のエージェント
                dim = 1     % 表示対象のモード
            end
            figure
            freq_variances_ = zeros(length(num),length(dim),obj.param.Nt);
            freq_x_variances_ = zeros(length(num),length(dim),obj.param.Nt);
            freq_y_variances_ = zeros(length(num),length(dim),obj.param.Nt);
            for t = obj.param.minimum_store+obj.param.deadlock_stepwith:obj.param.Nt
%                 freq_variances_(:,:,t) = var(obj.peak_freqs(num,dim,t-obj.param.deadlock_stepwith+1:t),0,3);   % 時刻に沿った分散を計算．N-1で正規化
%                 freq_x_variances_(:,:,t) = var(obj.peak_x_freqs(num,dim,t-obj.param.deadlock_stepwith+1:t),0,3);
%                 freq_y_variances_(:,:,t) = var(obj.peak_y_freqs(num,dim,t-obj.param.deadlock_stepwith+1:t),0,3);
                freq_variances_(:,:,t) = obj.calcPeakVariance(obj.peak_freqs(num,dim,:),t,obj.param);
                freq_x_variances_(:,:,t) = obj.calcPeakVariance(obj.peak_x_freqs(num,dim,:),t,obj.param);
                freq_y_variances_(:,:,t) = obj.calcPeakVariance(obj.peak_y_freqs(num,dim,:),t,obj.param);
            end
            for i = 1:length(dim)
                ylabel_str = "variance of \phi"+["","_x","_y"]+" (Hz^2)";
                fv_(:,:,1) = permute(freq_variances_(:,i,:),[1,3,2]);
                fv_(:,:,2) = permute(freq_x_variances_(:,i,:),[1,3,2]);
                fv_(:,:,3) = permute(freq_y_variances_(:,i,:),[1,3,2]);
                for j = 1:3
                    subplot(3,length(dim),j+3*(i-1))
                    semilogy(1:obj.param.Nt,fv_(:,:,j))
                    hold on
                    legend(string(num))
                    xlabel("timestep")
                    ylabel(ylabel_str(j))
                    yline(obj.param.freq_variance_hz)
                    hold off
                end
            end
        end

        function obj = peakPeriodicPlot(obj,num)
            arguments
                obj
                num         % 表示対象のエージェント
            end
            figure
            periodic_ = zeros(length(num),3,obj.param.Nt);
            for t = obj.param.minimum_store+obj.param.deadlock_stepwith_periodic:obj.param.Nt
                for i = 1:length(num)
                    periodic_(i,1,t) = obj.calcMaxPeriodic(obj.peak_freqs, num(i), t, obj.param);
                    periodic_(i,2,t) = obj.calcMaxPeriodic(obj.peak_x_freqs, num(i), t, obj.param);
                    periodic_(i,3,t) = obj.calcMaxPeriodic(obj.peak_y_freqs, num(i), t, obj.param);
%                     obj.is_deadlock_periodic(i,2,t) = obj.calcMaxPeriodic(obj.peak_x_freqs, i, t, obj.param) > obj.param.periodic_coeff_threshold;
%                     obj.is_deadlock_periodic(i,3,t) = obj.calcMaxPeriodic(obj.peak_y_freqs, i, t, obj.param) > obj.param.periodic_coeff_threshold;
                end
            end
            ylabel_str = "periodic of \phi"+["","_x","_y"]+" (Hz^2)";
            for j = 1:3
                subplot(3,1,j)
                plot(1:obj.param.Nt,permute(periodic_(:,j,:),[1,3,2]))
                hold on
                legend(string(num))
                xlabel("timestep")
                ylabel(ylabel_str(j))
                yline(obj.param.periodic_coeff_threshold)
                hold off
            end
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
                xlim([0,1000])
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
                xlim([0,1000])
                ylabel("Variance of Peak Frequency [Hz^2]")
                xlabel("TIme Step")
                title("mode "+string(mu))
            end
        end

        function obj = phaseGapPlot(obj)
            % ロボットの位置プロット
            arguments
                obj
            end
            figure
            plot(obj.t_vec, permute(obj.phi(:,1,:),[1,3,2])-mean( permute(obj.phi(:,1,:),[1,3,2]), 1 ))
        end

        function obj = generateSpectrumMovie(obj,filename, speed)
            arguments
                obj
                filename string = "movie.mp4" % 保存するファイル名
                speed = 1       % 動画の再生速度
            end
            obj.makeMovie(@obj.spectrumPlot, obj.param.dt, obj.param.Nt, filename, speed, true);
        end
    end
end

