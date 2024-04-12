
set(0,"DefaultAxesFontSize",13);    % フォントサイズ13
set(0,"DefaultLineLineWidth",2);    % 線の太さ2
set(0,"DefaultAxesXGrid",'on');     % X軸方向のグリッドON
set(0,"DefaultAxesYGrid",'on');     % Y軸方向のグリッドON

f = [f1,f2];
p = [p1,p2];

figure
plot(f,10*log10(p));
hold on

ylim([-100,20])
xlim([0,10])

power_threshold_dB = -60;
prominence_threshold_dB = 10;
get_peak_number = 3;

[~,peak_i,~,prom] = findpeaks(10*log10(p(:,1)),"MinPeakHeight",power_threshold_dB, "MinPeakProminence",prominence_threshold_dB); 
peak_freqs(:,1) = f(peak_i(1:get_peak_number),1);
peak_powers(:,1) = p(peak_i(1:get_peak_number),1);
peak_prominences(:,1) = prom(1:get_peak_number);
[~,peak_i,~,prom] = findpeaks(10*log10(p(:,2)),"MinPeakHeight",power_threshold_dB, "MinPeakProminence",prominence_threshold_dB);
peak_freqs(:,2) = f(peak_i(1:get_peak_number),2);
peak_powers(:,2) = p(peak_i(1:get_peak_number),2);
peak_prominences(:,2) = prom(1:get_peak_number);

plot(peak_freqs(:,1), 10*log10(peak_powers(:,1)),'o','MarkerSize',15);
plot(peak_freqs(:,2), 10*log10(peak_powers(:,2)),'+','MarkerSize',15);

legend(["t = 7227","t=7217","t = 7227","t=7217"])