function [Induced_src, BrainNoise_src,...
          Fs, Ntr] = SimulateSrc(dPhi, alpha_in)


    if(nargin < 7)
        alpha = 0.05;
    else
        alpha = alpha_in;
    end;

    % if(nargin == 0)
    %     bNewBrainNoise = true; % whether or not to generate new brain noise
    % end;

    N_brain = 1000;

    % create normalized forward matrix, leaving only two components in the tangential plane 
    % [Ns, Nsites] = size(G2d);

    % [Ggen, XYZGenAct] = GetGeneratingFwd(XYZGen, G2d, R);

    nw = [1, 2]; % If you want to simulate all four networks use

    sp = zeros(2, 500);
    T = 500; % number of timeslices per trial
    t = 1:T;

    % synchrony profiles of networks, one of each  as specified in lines 34-37
    sp(1,:) = exp(-1e-8 * (abs(t - 150)) .^ 4);
    sp(2,:) = exp(-1e-8 * (abs(t - 300)) .^ 4);
    sp(3,:) = exp(-0.2e-8 * (abs(t - 225)) .^ 4);
    sp(4,:) = 0.5 * (sin(10 * t / 500) + 1);

    Ntr = 100; %  number of trials
    range = 1:T;
    clear Data;
    Fs = 500; % sampling rate
    BrainNoise_src   = zeros(N_brain, Ntr * T);
    % SensorNoise  = zeros(Ns, Ntr * T);
    Induced_src      = zeros(2, Ntr * T);

    F1 = 10; % Hz;
    t = linspace(0, 1, Fs);
    clear s;

    % if(~bNewBrainNoise)
    %     % BN = load('BrainNoise2017.mat');
    %     % BN = load(['../data/BrainNoise2017_', num2str(ty), '.mat']);
    %     BN = load(['../data/BrainNoise_src', '.mat']);
    % end;

    fprintf('Simulating trial data ...\n');
    fprintf('Current trial number (Max %d):', Ntr);
    PhaseShiftsOut = zeros(Ntr, 8);
    for tr = 1:Ntr
        
        phi1 = 2 * (rand - 0.5) * pi;
        phi_alpha = alpha * (rand - 0.5) * pi;

        PhaseShiftsOut(tr, 1:2) = [phi1, phi_alpha];
        rnd_phi12 = phi_alpha + dPhi;
        s{1}(1,:) = sin(2 * pi * F1 * t + phi1) .* sp(1,:);
        s{1}(2,:) = sin(2 * pi * F1 * t + phi1 + rnd_phi12) .* sp(1,:);

         % collect activity from the selected networks

        % induced = Ggen * s{1};
        induced = s{1};
        % induced = induced / sqrt(sum((induced(:) .^ 2)));
        Induced_src(:,range) = induced;
        
        
        % Change this
        % if(bNewBrainNoise)
            % brainnoise = GenerateBrainNoise_SQUIDvsOPM(G2d, T, 500, 1000, Fs);
        brainnoise = GenerateBrainNoise_src(T, 500, N_brain, Fs);
        % else
        %     brainnoise = BN.BrainNoise_src(:, range);
        % end;
        
        % brainnoise = brainnoise / sqrt(sum((brainnoise(:) .^ 2)));
        BrainNoise_src(:, range) = brainnoise;
        
        % sensornoise = randn(Ns, T);
        % sensornoise = sensornoise / sqrt(sum((sensornoise(:) .^ 2)));
        % SensorNoise(:, range) = sensornoise;
        range = range + T;

        if tr > 1
          for j = 0 : log10(tr - 1)
              fprintf('\b'); % delete previous counter display
          end
         end
         fprintf('%d', tr);
    end

    % if(bNewBrainNoise)
    %    save(['./BrainNoise_src', '.mat'], 'BrainNoise_src');
    % end;

    fprintf('\nDone.\n');
end




% function [ BrainNoise ] = GenerateBrainNoise_src(G, T, Tw, N, Fs)
function [SourceNoise] = GenerateBrainNoise_src(T, Tw, N, Fs)
    % N - nuber of noise sources on cortex
    % Tw - time window for better filtering

    % Nsrc = size(G, 2);
    % Generate random indices of noise sources.
    % SrcIndex = fix(rand(1, N) * Nsrc + 1);
    % Nsrc = fix(size(G, 2));


    q = randn(N, T + 2 * Tw);

    alpha_band  = [8,  12] / (Fs / 2);
    beta_band   = [15, 30] / (Fs / 2);
    gamma1_band = [30, 50] / (Fs / 2);
    gamma2_band = [50, 70] / (Fs / 2);


    [b_alpha, a_alpha] = butter(4, alpha_band);
    [b_beta, a_beta] = butter(4, beta_band);
    [b_gamma1, a_gamma1] = butter(4, gamma1_band);
    [b_gamma2, a_gamma2] = butter(4, gamma2_band);

    SourceNoise = 1 / mean(alpha_band)  * filtfilt(b_alpha,  a_alpha,  q')' +...
                  1 / mean(beta_band)   * filtfilt(b_beta,   a_beta,   q')' +...
                  1 / mean(gamma1_band) * filtfilt(b_gamma1, a_gamma1, q')' +...
                  1 / mean(gamma2_band) * filtfilt(b_gamma2, a_gamma2, q')';

    % BrainNoise = G(:,SrcIndex) * SourceNoise(:, Tw + 1 : Tw + T) / N;
    SourceNoise = SourceNoise(:, Tw + 1 : Tw + T) / N;
end
