classdef Simulator
    %SIMULATOR シミュレータでよく使う機能をまとめたクラス
    %       パラメータの変更機能，動画の撮影機能などを入れる予定．継承して色々なシミュレータに利用してほしい
    
    properties
        param   % パラメータの入った構造体オブジェクト
    end
    
    methods (Access=public)
        function obj = Simulator()
            % コンストラクタ
            obj.param = struct;                     % paramオブジェクトの生成
            obj = obj.setFigureProperty('default'); % 図の設定を一旦デフォルトに戻す．bodeなどを扱うときにこっちで設定したプロパティが残っていると危険なので．
            obj = obj.setDefaultParameters();       % パラメタのデフォルト値を設定
        end
        
        function obj = setDefaultParameters(obj)
            % パラメータとデフォルト値を設定．継承で好きにいじってほしい
            obj.param.dt = 0.01;
        end

        function obj = setParam(obj, field, value)
            % 1つのパラメターを変更．フィールド名，ないものは警告出す．外から使う
            arguments
                obj             % Simulatorオブジェクト
                field string      % パラメータ名（文字ベクトル）
                value = 0       % パラメータ初期値．省略なら0が入る
            end
            if isfield(obj.param, field)
                obj.param.(field) = value;
            else
                disp("WARNIG：パラメータに "+field+" は存在しません");
            end
        end

        function obj = setMultiParam(obj, list)
            % 複数パラメータを同時に変更
            arguments
                obj         % Simulatorオブジェクト
                list        % フィールド名と値が並んだ配列．["Nt", 10, "use_CBF", true] みたいに使う
            end
            if mod(length(list),2) == 0
                paramN = length(list)/2;
                for n = 1:paramN
                    obj = obj.setParam(list(2*n-1), list(2*n));
                end
            else
                disp("WARNIG : パラメータリストの長さは，2の整数倍である必要があります");
            end
        end

        function obj = showSimulationTime(obj,t)
            % シミュレーション進行中，現在時刻をウィンドウに表示
            clc
            disp("t = "+string(t)+"/"+string(obj.param.Nt));
        end

        function obj = setFigureProperty(obj, preset)
            % 図のプロパティを変更．MATLABのオリジナルはちょっと文字が小さいので
            % 継承先では自分オリジナルのプリセットを作っても良いと思う
            arguments
                obj
                preset {mustBeMember(preset, ["default", "large"])} = "default" % 標準でdefault設定に
            end
            if preset == "default"
                set(0,"DefaultAxesFontSize",'default');
                set(0,"DefaultLineLineWidth",'default');
                set(0,"DefaultAxesXGrid",'default');
                set(0,"DefaultAxesYGrid",'default');
            end
            if preset == "large"
                % 論文掲載できる位の文字サイズ，線太さに．グリッドも標準でつける．
                set(0,"DefaultAxesFontSize",13);    % フォントサイズ13
                set(0,"DefaultLineLineWidth",2);    % 線の太さ2
                set(0,"DefaultAxesXGrid",'on');     % X軸方向のグリッドON
                set(0,"DefaultAxesYGrid",'on');     % Y軸方向のグリッドON
            end
        end
        
        function obj = makeMovie(obj, func, dt, Nt, filename, speed, capture_all, padding_frame)
            % 動画の作成
            arguments
                obj
                func                            % フレームごとの描画関数．シミュレーションオブジェクトと現在のカウントを渡す
                dt                              % シミュレーションの刻み時間
                Nt                              % シミュレーションカウント数
                filename string = "movie.mp4"   % 保存するファイル名
                speed = 1                       % 動画の再生速度
                capture_all = false             % 画面全体をキャプチャするか？
                padding_frame = round(1/dt/speed)   % 動画の開始終了前後に静止フレームを追加．使わない場合は0に
            end
            f = figure;
            disp("アニメーション描画を開始します")
            Nk = Nt/speed + 2*padding_frame;  % 前後に静止フレームを追加
            t_list = [1*ones(1,padding_frame),1:speed:Nt,Nt*ones(1,padding_frame)];   % 1とNtでパディング
            F(Nk) = struct('cdata',[],'colormap',[]);
            for k = 1:Nk
                func(t_list(k));
                drawnow;
                if capture_all
                    F(k) = getframe(f); % gcfをキャプチャ
                else
                    F(k) = getframe;    % gcaをキャプチャ
                end
                if strcmp(get(gcf,'currentcharacter'),'q')  % key stop
                    break; % 中止したいときはqを押す
                end
            end
            disp("アニメーション描画を終了します")
            v = VideoWriter(filename,'MPEG-4');
            v.FrameRate = round(1/dt);
            open(v);
            writeVideo(v,F(1:end));
            close(v);
            disp("Animation : 動画が保存されました")
        end

    end
end

