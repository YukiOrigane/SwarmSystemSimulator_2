

% こちらのファイルを実行する
% 一発実施する

%% 最初に1回だけ回すもの．おまじない
addpath(genpath("../../../SwarmSystemSimulator_2/"))    % パスを通す
simulation = SwarmWithWaveInteractionSimulation();                 % オブジェクトの定義
simulation.setFigureProperty("large");                  % 描画の基本設定を変更

%% 誘導場の生成


%% シミュレーションの実施 : 単発
simulation = simulation.setParam("environment_file","setting_files/environments/narrow_space_w_2_5_vertical_short.m");   % パラメタ変更
simulation = simulation.setParam("placement_file","setting_files/init_conditions/narrow_40.m");   % パラメタ変更
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
simulation = simulation.setParam("Nt",2000);
simulation = simulation.setParam("is_debug_view",false);
simulation = simulation.setParam("initial_pos_variance", 0);
%simulation = simulation.setParam("attract_force_type", "linear_fbx");
simulation = simulation.setParam("attract_force_type", "linear_fbxy");
% CBF %
simulation = simulation.setParam("cbf_rs", 0.8);  % 安全距離
simulation = simulation.setParam("cbf_gamma", 5); % ナイーブパラメタ
% kp調整 %
simulation = simulation.setParam("deadlock_source","cos");
simulation = simulation.setParam("do_kp_adjust",true);  % kp調整を実施？
simulation = simulation.setParam("kp_adjust_out",-0.3);
%simulation = simulation.setParam("kp_adjust_in",-0.3);
simulation = simulation.setParam("kp_adjust_in",1.2);
simulation = simulation.setParam("adjust_stepwith",80);
%simulation = simulation.setParam("dxdt_0",[[0 0];[0 0]]);   % パラメタ変更
% 本番 %
simulation = simulation.readSettingFiles(); % 設定ファイルの読み込み
simulation = simulation.initializeVariables();  % 初期値の計算
simulation = simulation.defineSystem();  % システム設定（誘導場の生成）
simulation = simulation.simulate(); % シミュレーションの実施
%% 描画とか
figure
simulation.edgeDeadlockPlot(2000,2);
simulation.placePlot(1);
% simulation.cos = simulation.cos.plot();
% simulation = simulation.generateMovieEstimate();
simulation = simulation.generateMovieEstimate("0528_movie_2times.mp4",2);
simulation = simulation.setParam("is_debug_view",true);
simulation = simulation.calcControlInput(10);
simulation.cos.relativePositionEstimate(750,[8,9,10]);  % 推定デバッグ表示
simulation.cos.peakAndFreqPlot([8,9,10]);   % エージェント毎ピーク履歴
simulation.cos.peakAndFreqPlot2([1,5:20]);   % モード毎ピーク履歴
% simulation.cos.spectrumPlot(750,8);   % 特定時刻スペクトラムプロット
% simulation.cos.generateSpectrumMovie();
% simulation.cos.deadlockPlot([1,5:20]);
% simulation.cos.variancePlot([1,5:20]);
% simulation.kpAdjustPlot([1,5:20]);
% simulation.minimumDistanceCheck();
% simulation.deadlockDetectionPlot("result");
% simulation.stopDetect(600);
% simulation.variancePlot([1,5:20]);
simulation.obtainNumberOfPassedRobots();
 
figure
subplot(2,2,1)
simulation.edgeDeadlockPlot(1,2,false);
subplot(2,2,2)
simulation.edgeDeadlockPlot(225,2,false);
subplot(2,2,3)
simulation.edgeDeadlockPlot(450,2,false);
subplot(2,2,4)
simulation.edgeDeadlockPlot(675,2,false);
figure
subplot(1,2,1)
simulation.trajectryJudgePlot([1601:2000]);
subplot(1,2,2)
simulation.trajectryJudgePlot([2001:3000]);

%% シミュレーションの実施 : 回す
addpath(genpath("../../../SwarmSystemSimulator_2/"))    % パスを通す
simulation = SwarmWithWaveInteractionSimulation();                 % オブジェクトの定義
simulation.setFigureProperty("large");                  % 描画の基本設定を変更


env_list = ["setting_files/environments/narrow_space_w_2.m",...
    "setting_files/environments/narrow_space_w_1_5.m",...
    "setting_files/environments/narrow_space_w_2_5.m"];
source_list = ["cos","stop"];
number_of_sets = length(env_list)*length(source_list);
sim_per_sets = 25;
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
            simulation.cos = simulation.cos.setParam("kappa",100);
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
            simulation = simulation.setParam("Nt",2000);
            simulation = simulation.setParam("is_debug_view",false);
            simulation = simulation.setParam("initial_pos_variance", 0.2);
            simulation = simulation.setParam("attract_force_type", "linear_fbx");
            % CBF %
            simulation = simulation.setParam("cbf_rs", 0.8);  % 安全距離
            simulation = simulation.setParam("cbf_gamma", 5); % ナイーブパラメタ
            % kp調整 %
            simulation = simulation.setParam("deadlock_source",source);
            simulation = simulation.setParam("do_kp_adjust",true);  % kp調整を実施？
            simulation = simulation.setParam("kp_adjust_out",-0.3);
            simulation = simulation.setParam("kp_adjust_in",-0.3);
            %simulation = simulation.setParam("kp_adjust_in",1.2);
            simulation = simulation.setParam("adjust_stepwith",80);
            %simulation = simulation.setParam("dxdt_0",[[0 0];[0 0]]);   % パラメタ変更
            % 本番 %
            simulation = simulation.readSettingFiles(); % 設定ファイルの読み込み
            simulation = simulation.initializeVariables();  % 初期値の計算
            simulation = simulation.defineSystem();  % システム設定（誘導場の生成）
            simulation = simulation.simulate(); % シミュレーションの実施
            
            save("data_0202/"+string(sets_count)+"_"+string(n)+".mat","simulation");
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
for sets = 1:number_of_sets
    dat = cell2mat(table2array(results(sets,"passing robots")));
    if results.env(sets) == "setting_files/environments/narrow_space_w_1_5.m"
        env_str = "w=1.5";
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
rows = (results_new.env == "w=2.0") .* (results_new.source == "cos");
dat1 = table2array(results_new(rows==1,"passing robots"));
rows = (results_new.env == "w=2.0") .* (results_new.source == "stop");
dat2 = table2array(results_new(rows==1,"passing robots"));
[h,p] = ttest2(dat1,dat2)
