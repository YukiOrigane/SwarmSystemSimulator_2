% mobileRobotsCOSSimulatorのテスタ
% 相互作用の変更と，内外推定のテストに使うことができる．
% こちらのファイルを実行する
% 一発実施する

%% 最初に1回だけ回すもの．おまじない
addpath(genpath("../../../SwarmSystemSimulator_2/"))    % パスを通す
simulation = mobileRobotCOSSimulator();                 % オブジェクトの定義
simulation.setFigureProperty("large");                  % 描画の基本設定を変更

%% シミュレーションの実施
rng(5); % 乱数シード固定
%simulation = simulation.setParam("K",2);   % パラメタ変更
%simulation.cos = simulation.cos.setParam("phi_0",rand(20,1));   % パラメタ変更
simulation.cos = simulation.cos.setParam("kappa",100);
simulation.cos = simulation.cos.setParam("interaction_type","wave");
simulation.cos = simulation.cos.setParam("do_estimate",true);
simulation.cos = simulation.cos.setParam("is_judge_continuous",false);
simulation = simulation.setParam("environment_file","setting_files/environments/free.m");   % パラメタ変更
%simulation = simulation.setParam("placement_file","setting_files/init_conditions/narrow_20.m");   % パラメタ変更
simulation = simulation.setParam("placement_file","setting_files/init_conditions/random_20.m");   % パラメタ変更
simulation = simulation.setParam("dxdt_0",zeros(20,2));   % パラメタ変更
simulation = simulation.setParam("rv",1.5);
simulation = simulation.setParam("Nt",3000);
simulation = simulation.readSettingFiles(); % 設定ファイルの読み込み
simulation = simulation.initializeVariables();  % 初期値の計算
simulation = simulation.defineSystem();  % システム設定（COSの初期化）

simulation = simulation.simulate(); % シミュレーションの実施
%simulation.phasePlacePlot(1);
%simulation.edgeJudgePlot(3999,1);
%simulation.cos.relativePositionEstimate(3999,[1,5,10]);
%simulation.generateMovie();
simulation.cos.plot();

% TODO : 格子状の初期配置でチェック