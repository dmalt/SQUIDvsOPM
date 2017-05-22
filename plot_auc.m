% Plot ROC-AUC vs SNR and pre-rec-AUC vs SNR
% for monte-carlo simulations of networks detection
% INPUT:
%   SPC<method_name>
%   TPR<method_name>
%   PPV<method_name>
%__________________________________________________

figure
% subplot(2,3,1)
% subplot(2,2,1)
% subplot(2,1,1)

range_snr = 1:size(InducedScale,2); 
% n_monte = 48;
n_monte = 48;
range_monte = 1:n_monte;

% figure
% Divide std by sqrt(n)
% for ty = 1:5
%     auc_roc{ty} = calc_auc(1 - SPCidics{ty}(range_monte,:,:), TPRidics{ty}(range_monte,:,:));
%     errorbar(range_snr, (mean(auc_roc{ty},1)), (std(auc_roc{ty}, 1)) / sqrt(n_monte));
%     hold on;
% end;
% title('iDICS');
% legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
% xlabel('SNR');
% ylabel('AUC ROC');
% xlim([0, 0.05]);
% ylim([0, 1]);

% figure
% subplot(2,3,2);
% for ty = 1:5
%     auc_roc{ty} = calc_auc(1 - SPCdics{ty}(range_monte,:,:), TPRdics{ty}(range_monte,:,:));
%     errorbar([InducedScale{:}], mean(auc_roc{ty},1), std(auc_roc{ty}, 1));
%     hold on;
% end
% title('DICS')
% legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
% xlabel('SNR');
% ylabel('AUC ROC');
% xlim([0, 0.05]);
% ylim([0, 1]);

% figure
% subplot(2,3,3);
% subplot(2,2,2);
% for ty = 1:5
%     auc_roc{ty} = calc_auc(1 - SPCgcs_dics{ty}(range_monte,:,:), TPRgcs_dics{ty}(range_monte,:,:));
%     errorbar(range_snr, (mean(auc_roc{ty},1)), (std(auc_roc{ty}, 1)));
%     hold on;
% end
% title('GCS-DICS')
% legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
% xlabel('SNR');
% ylabel('AUC ROC');
% xlim([0, 0.05]);
% ylim([0, 1]);

% --------------------------------------------------------------------------- %

% figure
% subplot(2,3,4);
% subplot(2,2,3);
% subplot(2,1,2);
for ty = 1:5
    auc_prec{ty} = calc_auc(TPRidics{ty}(range_monte,:,:), PPVidics{ty}(range_monte,:,:));
    errorbar(range_snr, (mean(auc_prec{ty},1)), (std(auc_prec{ty}, 1) / sqrt(n_monte)));
    hold on;
end
title('iDICS');
legend(sprintf('%s\n%s','SQUID', '(gradiometers)'), 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
xlabel('SNR');
ylabel('AUC Precision-Recall');

% figure
% subplot(2,3,5);
% for ty=1:5
%     auc_prec{ty} = calc_auc(TPRdics{ty}(range_monte,:,:), PPVdics{ty}(range_monte,:,:));
%     errorbar([InducedScale{:}], mean(auc_prec{ty},1), std(auc_prec{ty}, 1));
%     hold on;
% end
% title('DICS')
% legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
% xlabel('SNR');
% ylabel('AUC Precision-Recall');

% figure
% subplot(2,3,6);
% subplot(2,2,4);
% for ty=1:5
%     auc_prec{ty} = calc_auc(TPRgcs_dics{ty}(range_monte,:,:), PPVgcs_dics{ty}(range_monte,:,:));
%     errorbar(range_snr, (mean(auc_prec{ty},1)), (std(auc_prec{ty}, 1)));
%     hold on;
% end;
% title('GCS-DICS')
% legend('SQUID (gradiometers)', 'mSQUID', 'nOPM102', 'tOPM', 'nOPM204', 'location', 'best')
% xlabel('SNR');
% ylabel('AUC Precision-Recall');
