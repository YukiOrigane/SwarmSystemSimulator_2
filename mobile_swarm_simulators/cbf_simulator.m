
%% 最初に1回だけ回すもの．おまじない
addpath(genpath("../../../SwarmSystemSimulator_2/"))    % パスを通す

gamma_list = [0.1,0.2,0.5,1.0,2.0,5.0,10];
for gamma = gamma_list

%% シミュレーションの実施 : 単発
simulation = CBFSimulation();                 % オブジェクトの定義
simulation.setFigureProperty("large");                  % 描画の基本設定を変更

simulation = simulation.setParam("environment_file","setting_files/environments/free.m");   % パラメタ変更
simulation = simulation.setParam("placement_file","setting_files/init_conditions/free_Na_1.m");   % パラメタ変更

simulation = simulation.setParam("Nt",200);

% CBF %
simulation = simulation.setParam("cbf_rs", 1.0);  % 安全距離
simulation = simulation.setParam("cbf_gamma", gamma); % ナイーブパラメタ
simulation = simulation.setParam("cbf_lb", []); % 入力下限 ex) [-10; -10]
simulation = simulation.setParam("cbf_ub", []); % 入力上限 ex) [10; 10]

% 初期条件 %
%simulation = simulation.setParam("dxdt_0", 2*[2;-5]/vecnorm([2;-5]));
simulation = simulation.setParam("dxdt_0", 2*[0;-5]/vecnorm([0;-5]));

% 本番 %
simulation = simulation.readSettingFiles(); % 設定ファイルの読み込み
rng(5);     % 乱数の固定
simulation = simulation.initializeVariables();  % 初期値の計算
simulation = simulation.defineSystem();  % システム設定（誘導場の生成）
simulation = simulation.simulate(); % シミュレーションの実施
%% 描画とか
f = figure('Position',[100 100 950 400]);
subplot(1,2,1)
simulation.movementPlot(200);
%simulation = simulation.generateMovie("0716_Na01_motion.mp4",8);
xlabel("x (m)")
ylabel("y (m)")

subplot(1,2,2)
plot(simulation.t_vec,permute(simulation.dxdt(1,:,:),[2,3,1]));
hold on
plot(simulation.t_vec,permute(vecnorm(simulation.dxdt(1,:,:),2,2),[2,3,1]),'Color',"#7E2F8E");
legend(["$\dot x$","$\dot y$","speed"],'Interpreter','latex')
xlabel("time (s)")
ylim([-2 2])
ylabel("Velocity (m/s)")
saveas(f,"data_0716_cbf/gamma_"+string(gamma)+".png")

end
