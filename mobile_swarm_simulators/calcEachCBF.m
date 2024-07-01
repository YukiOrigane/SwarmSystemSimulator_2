function [ui_t, lambda_] = calcEachCBF(param_,x_ij,dxdt_ij,x_io,dxdt_io,Adj_coli,ui_nominal,cbf_)
    % エージェント毎CBF計算．parfor計算のために分離
    % ロボット間衝突回避CBF %
    cbf_ = cbf_.setParameters(1,param_.cbf_rs,param_.dt,param_.cbf_gamma,true);
    cbf_ = cbf_.addConstraints([x_ij(Adj_coli==1,1), x_ij(Adj_coli==1,2)], [dxdt_ij(Adj_coli==1,1), dxdt_ij(Adj_coli==1,2)]);
    % 隣接ロボットとの相対ベクトルのみCBF制約として利用
    % 壁との衝突回避CBF %
    cbf_ = cbf_.setParameters(1,param_.cbf_rs,param_.dt,param_.cbf_gamma,false);
    cbf_ = cbf_.addConstraints(permute(x_io,[3,2,1]), -repmat(dxdt_io,length(x_io(:,:,:)),1));
    % 壁との相対位置ベクトルと，自身の速度ベクトル(壁との相対速度ベクトル)をCBFに入れる
    % 入力範囲の制限 %
    cbf_ = cbf_.addInputMinMaxConstraint(param_.cbf_lb,param_.cbf_ub);
    [ui_t,lambda_] = cbf_.apply(ui_nominal(1,:));
end