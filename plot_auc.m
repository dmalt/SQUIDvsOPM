
figure
subplot(2,3,1)

i_snr = 1;
range_monte = 1:9;

% figure
for ty = 1:5
    auc_roc{ty} = calc_auc(1 - SPCidics{ty}(range_monte,:,:), TPRidics{ty}(range_monte,:,:));
    errorbar([InducedScale{:}], mean(auc_roc{ty},1), std(auc_roc{ty}, 1));
    hold on;
end;
title('iDICS');
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
xlabel('SNR');
ylabel('AUC ROC');
% xlim([0, 0.05]);
% ylim([0, 1]);

% figure
subplot(2,3,2);
for ty = 1:5
    auc_roc{ty} = calc_auc(1 - SPCdics{ty}(range_monte,:,:), TPRdics{ty}(range_monte,:,:));
    errorbar([InducedScale{:}], mean(auc_roc{ty},1), std(auc_roc{ty}, 1));
    hold on;
end
title('DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
xlabel('SNR');
ylabel('AUC ROC');
% xlim([0, 0.05]);
% ylim([0, 1]);

% figure
subplot(2,3,3);
for ty = 1:5
    auc_roc{ty} = calc_auc(1 - SPCgcs_dics{ty}(range_monte,:,:), TPRgcs_dics{ty}(range_monte,:,:));
    errorbar([InducedScale{:}], mean(auc_roc{ty},1), std(auc_roc{ty}, 1));
    hold on;
end
title('GCS-DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
xlabel('SNR');
ylabel('AUC ROC');
% xlim([0, 0.05]);
% ylim([0, 1]);

% --------------------------------------------------------------------------- %

% figure
subplot(2,3,4);
for ty = 1:5
    auc_prec{ty} = calc_auc(TPRidics{ty}(range_monte,:,:), PPVidics{ty}(range_monte,:,:));
    errorbar([InducedScale{:}], mean(auc_prec{ty},1), std(auc_prec{ty}, 1));
    hold on;
end
title('iDICS');
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
xlabel('SNR');
ylabel('AUC Precision-Recall');

% figure
subplot(2,3,5);
for ty=1:5
    auc_prec{ty} = calc_auc(TPRdics{ty}(range_monte,:,:), PPVdics{ty}(range_monte,:,:));
    errorbar([InducedScale{:}], mean(auc_prec{ty},1), std(auc_prec{ty}, 1));
    hold on;
end
title('DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
xlabel('SNR');
ylabel('AUC Precision-Recall');

% figure
subplot(2,3,6);
for ty=1:5
    auc_prec{ty} = calc_auc(TPRgcs_dics{ty}(range_monte,:,:), PPVgcs_dics{ty}(range_monte,:,:));
    errorbar([InducedScale{:}], mean(auc_prec{ty},1), std(auc_prec{ty}, 1));
    hold on;
end;
title('GCS-DICS')
legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
xlabel('SNR');
ylabel('AUC Precision-Recall');
