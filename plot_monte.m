figure
subplot(2,3,1)

% i_snr = 5;
i_snr = 6;
% range_monte = 1:11;
range_monte = 1:99;

% figure
for ty = 1:5
    plot(1 - squeeze(mean(SPCidics{ty}(range_monte,i_snr,:), 1)),...
             squeeze(mean(TPRidics{ty}(range_monte,i_snr,:), 1)));%, col{ty})
    hold on;
end;
title('iDICS');
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204')
xlabel('1 - specificity');
ylabel('sensitivity');
xlim([0, 0.05]);
ylim([0, 1]);

% figure
subplot(2,3,2);
for ty = 1:5
    plot(1 - squeeze(mean(SPCdics{ty}(range_monte,i_snr,:), 1)),...
             squeeze(mean(TPRdics{ty}(range_monte,i_snr,:), 1)));%, col{ty + 3})
    hold on;
end
title('DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204')
xlabel('1 - specificity');
ylabel('sensitivity');
xlim([0, 0.05]);
ylim([0, 1]);

% figure
subplot(2,3,3);
for ty = 1:5
    plot(1 - squeeze(mean(SPCgcs_dics{ty}(range_monte,i_snr,:), 1)),...
             squeeze(mean(TPRgcs_dics{ty}(range_monte,i_snr,:), 1)));%, col{ty + 6})
    hold on;
end
title('GCS-DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204')
xlabel('1 - specificity');
ylabel('sensitivity');
xlim([0, 0.05]);
ylim([0, 1]);

% --------------------------------------------------------------------------- %

% figure
subplot(2,3,4);
for ty = 1:5
    plot(squeeze(mean(TPRidics{ty}(range_monte,i_snr,:), 1)),...
         squeeze(mean(PPVidics{ty}(range_monte,i_snr,:), 1)));%, col{ty})
    hold on;
end
title('iDICS');
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204')
xlabel('recall');
ylabel('precision');

% figure
subplot(2,3,5);
for ty=1:5
    plot( squeeze(mean(TPRdics{ty}(range_monte,i_snr,:), 1)),...
          squeeze(mean(PPVdics{ty}(range_monte,i_snr,:), 1)));%, col{ty + 3})
    hold on;
end
title('DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204')
xlabel('recall');
ylabel('precision');

% figure
subplot(2,3,6);
for ty=1:5
    plot( squeeze(mean(TPRgcs_dics{ty}(range_monte,i_snr,:), 1)),...
          squeeze(mean(PPVgcs_dics{ty}(range_monte,i_snr,:), 1)));%, col{ty + 6})
    hold on;
end;
title('GCS-DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204')
xlabel('recall');
ylabel('precision');
