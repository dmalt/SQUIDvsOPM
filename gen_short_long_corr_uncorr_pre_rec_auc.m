% Generate figure with 4 plots of AUC vs SNR 
% for network detection with iDICS for following cases:
% - long networks with uncorrelated noise (dist > 3 cm)
% - short networks with uncorrelated noise (dist < 3cm)
% - long networks with spatially correlated noise
% - short networks with spatially correlated noise
%
% Images are generated in grayscale
% _____________________________________________________
% Author: dmalt, date: May 22 11:27 MSK 2017



data_path = '../data/white_noise/';
fnames = {'short_ntw_correlated_noise.mat',...
          'long_ntw_50_mc_trials_correlated_noise.mat',...
          'short_ntw_48_mc_trials.mat',...
          'long_ntw_50_mc_trials.mat'};

n_monte = 48;
range_monte = 1:n_monte;

l_style = {'-.k', '--k', ':k', '-k'};

figure;
for i_file = 1:length(fnames)
    scores = load([data_path, fnames{i_file}],...
                  'TPRidics', 'PPVidics', 'InducedScale');
    range_snr = 1:size(InducedScale,2); 

    subplot(2,2,i_file);

    for ty = 1:4
        auc_prec{ty} = calc_auc(scores.TPRidics{ty}(range_monte,:,:),...
                                scores.PPVidics{ty}(range_monte,:,:));

        errorbar(range_snr,...
                 (mean(auc_prec{ty},1)), (std(auc_prec{ty}, 1) / sqrt(n_monte)),...
                 l_style{ty}, 'linewidth', 2);
        % title(fnames(i_file));
        hold on;
    end
    % title('iDICS');
    legend('gSQUID204', 'mSQUID102', 'nOPM102', 'tOPM204', 'location', 'best')
    xlabel('SNR index', 'FontSize', 10);
    ylabel('AUC pre-rec', 'FontSize', 10);
end

% Print column and row titles. To put titles in the right place
% rescale the image
t1 = text(-12., 0.9, sprintf('%s\n%s','Correlated', '   noise'), 'FontSize', 14);
set(t1, 'rotation', 90);
t2 = text(-12., 0, sprintf('%s\n%s','Uncorrelated', '     noise'), 'FontSize', 14);
set(t2, 'rotation', 90);
t3 = text(-8.7, 1.6, 'Short range (< 3cm apart)', 'FontSize', 14);
t4 = text(1.9, 1.6, 'Long range (> 3cm apart)', 'FontSize', 14);
% t1 = text(-12, 1, 'Correlated noise', 'FontSize', 14);
% set(t1, 'rotation', 90);
