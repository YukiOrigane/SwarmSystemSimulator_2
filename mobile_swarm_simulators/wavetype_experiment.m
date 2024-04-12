

% こちらのファイルを実行する
% 一発実施する

%% 最初に1回だけ回すもの．おまじない
addpath(genpath("../../../SwarmSystemSimulator_2/"))    % パスを通す
simulation = SwarmWithWaveInteractionSimple();                 % オブジェクトの定義
simulation.setFigureProperty("large");                  % 描画の基本設定を変更

%% 誘導場の生成


%% シミュレーションの実施 : 単発
simulation = simulation.setParam("environment_file","setting_files/environments/square.m");   % パラメタ変更
simulation = simulation.setParam("placement_file","setting_files/init_conditions/round_20.m");   % パラメタ変更
%simulation = simulation.setParam("placement_file","setting_files/init_conditions/round_10.m");   % パラメタ変更
%simulation = simulation.setParam("placement_file","setting_files/init_conditions/free_Na_3.m");   % パラメタ変更
%simulation = simulation.setParam("placement_file","setting_files/init_conditions/missed_Na_10.m");   % パラメタ変更
%simulation = simulation.setParam("placement_file","setting_files/init_conditions/round_1.m");   % パラメタ変更
%simulation = simulation.setParam("environment_file","setting_files/environments/narrow_space_hosome_w_4.m");   % パラメタ変更
%simulation = simulation.setParam("placement_file","setting_files/init_conditions/narrow_20.m");   % パラメタ変更
%simulation = simulation.setParam("environment_file","setting_files/environments/free.m");   % パラメタ変更
%simulation = simulation.setParam("placement_file","setting_files/init_conditions/narrow_20.m");   % パラメタ変更
% COS %
simulation.cos = simulation.cos.setParam("interaction_type","diffusion");
%simulation.cos = simulation.cos.setParam("interaction_type","wave");
simulation.cos = simulation.cos.setParam("kappa",100);
%simulation.cos = simulation.cos.setParam("dt",0.01);
simulation.cos = simulation.cos.setParam("gamma",0);
simulation.cos = simulation.cos.setParam("self_exitation_gain", 0);
%simulation.cos = simulation.cos.setParam("self_exitation_threshold", -60);
simulation.cos = simulation.cos.setParam("do_estimate",true);
simulation.cos = simulation.cos.setParam("time_histry",256);
simulation.cos = simulation.cos.setParam("power_threshold_dB",-60); % ピークの高さ閾値
simulation.cos = simulation.cos.setParam("prominence_threshold_dB",5); % ピークのプロミネンス閾値
simulation.cos = simulation.cos.setParam("use_peak_matching", true);   % ピークマッチングを利用
simulation.cos = simulation.cos.setParam("number_of_matching_peaks", 3);    % ピークマッチングを何点でとるか？
simulation.cos = simulation.cos.setParam("peak_memory_num",1);
simulation.cos = simulation.cos.setParam("deadlock_usepower",true);% true
simulation.cos = simulation.cos.setParam("is_normalize",false);  % 相互作用の正規化を行う
simulation.cos = simulation.cos.setParam("deadlock_stepwith",100);  % デッドロック判定期間を延長
simulation.cos = simulation.cos.setParam("use_softtouch",true);    % ソフトタッチ使うか？
simulation.cos = simulation.cos.setParam("use_softrelease",true);  % ソフトリリース使うか？
simulation.cos = simulation.cos.setParam("delta_release",0.1);      % ソフトリリースの下限値
simulation.cos = simulation.cos.setParam("power_variance_db",2*10^-3);
%simulation.cbf = simulation.cbf.disable();
% 停止検知 %
simulation = simulation.setParam("stop_timehistry",256);
%simulation = simulation.setParam("stop_threshold",10^-3);
%simulation = simulation.setParam("stop_threshold",10^-2);
simulation = simulation.setParam("stop_threshold",0.5);
% Swarm %
simulation = simulation.setParam("kp",8);   % Swarm : 勾配追従力ゲイン
%simulation = simulation.setParam("kp",0);   % ！！！停止注意！！！
simulation = simulation.setParam("kf",0);  % Swarm : 群形成力ゲイン
simulation = simulation.setParam("kd",10);   % Swarm : 粘性ゲイン
simulation = simulation.setParam("Nt",3000);
simulation = simulation.setParam("dt",0.05);
simulation = simulation.setParam("is_debug_view",false);
simulation = simulation.setParam("initial_pos_variance", 0);
%simulation = simulation.setParam("attract_force_type", "linear_fbx");
simulation = simulation.setParam("attract_force_type", "trip");
simulation = simulation.setParam("adjacency_method", "distance");   % 隣接行列の作り方
% CBF %
simulation = simulation.setParam("cbf_rs", 0.8);  % 安全距離
simulation = simulation.setParam("cbf_gamma", 5); % ナイーブパラメタ
% kp調整 %
%simulation = simulation.setParam("deadlock_source","cos");
simulation = simulation.setParam("deadlock_source","stop");
%simulation = simulation.setParam("do_kp_adjust",false);  % kp調整を実施？
simulation = simulation.setParam("do_kp_adjust",true);  % kp調整を実施？
simulation = simulation.setParam("kp_adjust_out",-0.3);
%simulation = simulation.setParam("kp_adjust_in",-0.3);
simulation = simulation.setParam("kp_adjust_in",1.2);
simulation = simulation.setParam("adjust_stepwith",80);
%simulation = simulation.setParam("dxdt_0",[[0 0];[0 0]]);   % パラメタ変更
% trip %
%simulation = simulation.setParam("trip_mode","straight");
simulation = simulation.setParam("trip_mode","round");
% 本番 %
simulation = simulation.readSettingFiles(); % 設定ファイルの読み込み
simulation = simulation.initializeVariables();  % 初期値の計算
simulation = simulation.defineSystem();  % システム設定（誘導場の生成）
simulation = simulation.simulate(); % シミュレーションの実施
% データ保存先フォルダの作成と指定 %
folder_name = "data"+string(datetime('now','Format','yyyyMMdd/HH_mm_ss'));
% folder_name = "Nt3000_stop";
simulation = simulation.setParam("folder_name",folder_name);
simulation.cos = simulation.cos.setParam("folder_name",folder_name);
mkdir(folder_name);
save(folder_name+"/simulation.mat")

%% 描画とか
figure
%simulation.edgeDeadlockPlot(1,2);
%simulation.numberPlacePlot(3121,true);
%simulation.placePlot(1,false);
simulation.cos = simulation.cos.phaseGapPlot();
% simulation.cos = simulation.cos.partPlot(3500);
% simulation = simulation.generateMovieEstimate([],10);
%simulation = simulation.generateMovieTrip("movie_fast.mp4",10);
% simulation = simulation.generateMovieDetection("detection.mp4",20);
%simulation.cos = simulation.cos.generatePhaseMovie("phase.mp4",10);
%simulation = simulation.generateMovieEstimate();
simulation = simulation.setParam("is_debug_view",true);
%simulation = simulation.calcControlInput(8000);
%simulation.cos.relativePositionEstimate(150,[8,9,10]);  % 推定デバッグ表示
%simulation.cos.peakAndFreqPlot([1,6,8]);   % エージェント毎ピーク履歴
%simulation.cos.peakAndFreqPlot2([1,2]);   % モード毎ピーク履歴
% simulation.cos.spectrumPlot(7227,true,[1,2]);   % 特定時刻スペクトラムプロット
%simulation.cos.generateSpectrumMovie("spectrum.mp4",5);
% simulation.cos.deadlockPlot([1,5:20]);
% simulation.cos.variancePlot([9,10]);
%simulation.cos.meanVariancePlot([9,10]);
% simulation.kpAdjustPlot([1,5:20]);
% simulation.minimumDistanceCheck();
% simulation.deadlockDetectionPlot("stop");
%simulation.deadlockDetectionPlotColor();
%simulation.cos.phaseTimeVariancePlot();
% simulation.stopDetect(600);
% simulation.variancePlot([1:10]);
% simulation.plotPositionVariance();
%simulation.obtainNumberOfPassedRobots();
%simulation.cos.compareEstimateAndEigenFreq([1],2);
simulation.saveFigure(gcf, "estimate_eigen");

%simulation.cos.phaseSpacePlot(2500,0:9,[4,3]);
simulation.saveFigure(gcf, "phase_space");

figure
simulation.cos.energyPlot();
simulation.saveFigure(gcf, "energy");
%simulation.cos.generateEnergyMovie("energy.mp4", 10);
%simulation.cos.generatePhaseSpaceMovie("phase_space.mp4",1);
%{
t = 7217;
i = 1;
t_start_ = t-simulation.cos.param.time_histry;
[p_,f_] = pspectrum(permute(simulation.cos.phi(i,1,t_start_:t),[3,1,2]), simulation.cos.t_vec(t_start_:t));
figure
plot(f_,10*log10(p_));
text(max(f_)*0.7, 0, "t = "+string(t), 'FontSize',12);
ylim([-100,20])
xlim([0,10])
%}
%{
t_list = [1 600 1200 1800];
for t = t_list
    figure
    simulation.tripPlot(t);
    xticks(-8:4:8)
    yticks(-8:4:8)
end
figure
subplot(3,2,1)
simulation.tripPlot(1);
subplot(3,2,2)
simulation.tripPlot(600);
subplot(3,2,3)
simulation.tripPlot(1200);
subplot(3,2,4)
simulation.tripPlot(1800);
subplot(3,2,5)
simulation.tripPlot(2400);
subplot(3,2,6)
simulation.tripPlot(3000);
%}
%{ 
figure
subplot(2,2,1)
simulation.trajectryJudgePlot([2:400]);
subplot(2,2,2)
simulation.trajectryJudgePlot([401:800]);
subplot(2,2,3)
simulation.trajectryJudgePlot([801:1200]);
subplot(2,2,4)
simulation.trajectryJudgePlot([1201:1600]);
figure
subplot(1,2,1)
simulation.trajectryJudgePlot([1601:2000]);
subplot(1,2,2)
simulation.trajectryJudgePlot([2001:3000]);
%}
%% シミュレーションの実施 : 回す
%{
addpath(genpath("../../../SwarmSystemSimulator_2/"))    % パスを通す
simulation = SwarmWithWaveInteractionSimulation();                 % オブジェクトの定義
simulation.setFigureProperty("large");                  % 描画の基本設定を変更


env_list = ["setting_files/environments/narrow_space_hosome_w_4.m"];
source_list = ["cos","stop"];
number_of_sets = length(env_list)*length(source_list);
sim_per_sets = 3;
results = table('size',[number_of_sets,3],'VariableTypes',["string","string",'cell']);
results.Properties.VariableNames = ["env","source","passing robots"];
sets_count = 0;
for env = env_list
    for source = source_list
        sets_count = sets_count+1;
        passing = zeros(1,sim_per_sets);
        for n = 1:sim_per_sets  % 1セット
            clc
            disp(string(sets_count)+"/"+string(number_of_sets)+"セット，"+string(n)+"/"+string(sim_per_sets)+"試行");

            %%% シミュレーション %%%
            simulation = simulation.setParam("environment_file",env);   % パラメタ変更
            simulation = simulation.setParam("placement_file","setting_files/init_conditions/narrow_20.m");   % パラメタ変更
            % COS %
            simulation.cos = simulation.cos.setParam("kappa",80);
            simulation.cos = simulation.cos.setParam("do_estimate",true);
            simulation.cos = simulation.cos.setParam("time_histry",256);
            simulation.cos = simulation.cos.setParam("power_threshold",10^-7);
            %simulation.cbf = simulation.cbf.disable();
            % 停止検知 %
            simulation = simulation.setParam("stop_timehistry",256);
            simulation = simulation.setParam("stop_threshold",10^-3);
            % Swarm %
            simulation = simulation.setParam("kp",8);   % Swarm : 勾配追従力ゲイン
            simulation = simulation.setParam("kf",0);  % Swarm : 群形成力ゲイン
            simulation = simulation.setParam("kd",10);   % Swarm : 粘性ゲイン
            simulation = simulation.setParam("Nt",3000);
            simulation = simulation.setParam("is_debug_view",false);
            simulation = simulation.setParam("initial_pos_variance", 0.2);
            simulation = simulation.setParam("attract_force_type", "trip");
            simulation = simulation.setParam("trip_mode","straight");
            % CBF %
            simulation = simulation.setParam("cbf_rs", 0.8);  % 安全距離
            simulation = simulation.setParam("cbf_gamma", 5); % ナイーブパラメタ
            % kp調整 %
            if source == "stop"
                simulation = simulation.setParam("deadlock_source","stop");
                simulation = simulation.setParam("kp_adjust_in",-0.3);
            elseif source == "cos"
                simulation = simulation.setParam("deadlock_source","cos");
                simulation = simulation.setParam("kp_adjust_in",-0.3);
            elseif source == "cos_inout"
                simulation = simulation.setParam("deadlock_source","cos");
                simulation = simulation.setParam("kp_adjust_in",1.2);
            end
            %simulation = simulation.setParam("deadlock_source",source);
            simulation = simulation.setParam("do_kp_adjust",true);  % kp調整を実施？
            simulation = simulation.setParam("kp_adjust_out",-0.3);
            %simulation = simulation.setParam("kp_adjust_in",-0.3);
            %simulation = simulation.setParam("kp_adjust_in",1.2);
            simulation = simulation.setParam("adjust_stepwith",80);
            %simulation = simulation.setParam("dxdt_0",[[0 0];[0 0]]);   % パラメタ変更
            % trip %
            simulation = simulation.setParam("trip_mode","straight");
            % 本番 %
            simulation = simulation.readSettingFiles(); % 設定ファイルの読み込み
            simulation = simulation.initializeVariables();  % 初期値の計算
            simulation = simulation.defineSystem();  % システム設定（誘導場の生成）
            simulation = simulation.simulate(); % シミュレーションの実施
            
            save("data_0204/"+string(sets_count)+"_"+string(n)+".mat","simulation");
            passing(1,n) = simulation.obtainNumberOfPassedRobots();
        end
        results(sets_count,:) = {env,source,num2cell(passing,[1 2])};
    end
end

%% 複数回解析
cond(2) = RunningCondition();
cond(1).setParameter("source")
results_new = table('size',[number_of_sets*sim_per_sets,3],'VariableTypes',["string","string",'double']);
results_new.Properties.VariableNames = ["env","source","passing robots"];
i = 0;
for sets = 1:3%number_of_sets
    dat = cell2mat(table2array(results(sets,"passing robots")));
    if results.env(sets) == "setting_files/environments/narrow_space_hosome_w_2_1.m"
        env_str = "w=2.1";
    elseif results.env(sets) == "setting_files/environments/narrow_space_w_2.m"
        env_str = "w=2.0";
    elseif results.env(sets) == "setting_files/environments/narrow_space_w_2_5.m"
        env_str = "w=2.5";
    end
    for n = 1:sim_per_sets
        i = i+1;
        results_new(i,:) = {env_str, results.source(sets), dat(n)};
    end
end

% 箱ひげ図
boxchart(categorical(results_new.env),results_new.("passing robots"),'GroupByColor',results_new.source);
legend
ylabel("Number of Passing Robots")
xlabel("Width of Aisle")

% 検定
rows = (results_new.env == "w=2.1") .* (results_new.source == "stop");
dat1 = table2array(results_new(rows==1,"passing robots"));
rows = (results_new.env == "w=2.1") .* (results_new.source == "cos_inout");
dat2 = table2array(results_new(rows==1,"passing robots"));
[h,p] = ttest2(dat1,dat2)
%}