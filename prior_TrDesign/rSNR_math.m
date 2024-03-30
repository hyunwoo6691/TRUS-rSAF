clc;

num_syn = [3 5 7 9 11 13 15 17 19 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65];

depth_ = 1:7;
root_n = sqrt(num_syn);

% c = linspace(1,30,numel(num_syn));
c = num_syn;
%%
for d_idx = depth_
    fele5mm_noise_con = mean(fele5mm_noise(:,d_idx,1));
    fele5mm_signal_con = mean(fele5mm_signal(:,d_idx,1));
    
    fele5mm_noise_rsaf = mean(squeeze(fele5mm_noise(:,d_idx,2:end)));
    fele5mm_signal_rsaf = mean(squeeze(fele5mm_signal(:,d_idx,2:end)));
    
    signal_ratio = fele5mm_signal_rsaf/fele5mm_signal_con;
    noise_ratio = fele5mm_noise_con./fele5mm_noise_rsaf;
    
    snr_ = signal_ratio.*noise_ratio;
    
    c_sq = signal_ratio.*num_syn.^2;
    figure(1);
    subplot(1,4,1);
    plot(num_syn, signal_ratio,'LineWidth',2); hold on;
    subplot(1,4,2);
    plot(num_syn, sqrt(c_sq),'LineWidth',2); hold on;
    subplot(1,4,3);
    plot(num_syn, noise_ratio,'LineWidth',2); hold on;
    subplot(1,4,4);
    plot(num_syn, snr_,'LineWidth',2); hold on;
    plot(num_syn, c_sq./num_syn,'LineWidth',2);
    
end

subplot(1,4,1);
plot(num_syn, 1./(num_syn.^2),'color','r','LineWidth',2); hold off;
title('S_{rSAF}/S_{CON}');
legend('1cm','2cm','3cm', '4cm', '5cm', '6cm', '7cm', '1/N_{syn}^{2}','location','northeast');
set(gca,'FontSize',14,'FontWeight','bold');

subplot(1,4,2);
plot(num_syn, num_syn, 'color','k','LineWidth',2); hold off;
title('c');
legend('1cm','2cm','3cm', '4cm', '5cm', '6cm', '7cm', 'ref','location','northwest');
set(gca,'FontSize',14,'FontWeight','bold');

subplot(1,4,3);plot(num_syn, num_syn,'color','k','LineWidth',2); hold off;
title('N_{CON}/N_{rSAF}');
legend('1cm','2cm','3cm', '4cm', '5cm', '6cm', '7cm', 'ref','location','northwest');
set(gca,'FontSize',14,'FontWeight','bold');

subplot(1,4,4); hold off;
title('rSNR');
legend('1cm', 'ref-1cm', '2cm', 'ref-2cm', '3cm', 'ref-3cm', '4cm', 'ref-4cm', '5cm', 'ref-5cm', '6cm', 'ref-6cm', '7cm', 'ref-7cm', 'location','northwest');

set(gcf,'Position',[241 78 2630 882]);
set(gca,'FontSize',14,'FontWeight','bold');