
beta_list = [0.2 0.5 0.8];
T_fft_list = [128 256 512 1024];
start_list = [1 60 120];
figure('Position', [100 100 500 400]);
for i = 1:length(beta_list)
    for j = 1:length(T_fft_list)
        for k = 1:length(start_list)
            N = T_fft_list(j);
            t_start = t_start_list(k);
            beta = beta_list(i);
            pspectrum(permute(simulation.cos.phi(48,1,t_start:t_start+N-1),[3,1,2]),0:simulation.param.dt:simulation.param.dt*(N-1),'Leakage',beta)
            xline(f(1:3),'--')
            xlim([0 5])
            saveas(gcf,"pspectrum_beta_"+string(beta)+"_T_fft_"+string(N)+"_start_"+string(t_start)+".png")
        end
    end
end