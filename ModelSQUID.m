% clear all;
% close all;


import ups.spm_svd
import ups.CrossSpectralTimeseries
import ups.GenerateROC
import ups.ReduceToTangentSpace


% cd('/home/asus/MyProjects/SQUIDvsAM_MEG/Data');
InducedScale = {0.25, 0.5, 0.75, 1.0, 1.25}; 

data = matfile('../data/data_2D.mat') ;
data = data.data(1,1);
data = data{1};

GainSVDTh = 0.0001;

G_OPM = data.nOPM.L;
G_OPM = ReduceToTangentSpace(G_OPM, 'all');

G_mSQUID = data.mSQUID.L;
G_mSQUID = ReduceToTangentSpace(G_mSQUID, 'all');

G_SQUID = data.SQUID.L;
G_SQUID = ReduceToTangentSpace(G_SQUID, 'grad');

RL = data.sources{1}.p;
RR = data.sources{2}.p;
R = [RL; RR];

clear data;

data_big = matfile('../data/data_tOPM.mat');
data_big = data_big.data_tOPM(1,1);
data_big = data_big{1};

G_tOPM = data_big.tOPM.L;
G_tOPM = ReduceToTangentSpace(G_tOPM, 'all');

clear data_big;


data_nopm = matfile('../data/data_nOPM204.mat');
data_nopm = data_nopm.data_nOPM204(1, 1);
data_nopm = data_nopm{1};

G_nOPM204 = data_nopm.nOPM204.L;
G_nOPM204 = ReduceToTangentSpace(G_nOPM204, 'all');

clear data_nopm;

% Do rank reduction, use the same threshold for the two sensor types 
[ug_OPM, ~, ~] = spm_svd(G_OPM * G_OPM', GainSVDTh);
UP_OPM = ug_OPM';

[ug_mSQUID, ~, ~] = spm_svd(G_mSQUID * G_mSQUID', GainSVDTh);
UP_mSQUID = ug_mSQUID';

[ug_SQUID, ~, ~] = spm_svd(G_SQUID * G_SQUID', GainSVDTh);
UP_SQUID = ug_SQUID';

[ug_tOPM, ~, ~] = spm_svd(G_tOPM * G_tOPM', GainSVDTh);
UP_tOPM = ug_tOPM';

[ug_nOPM204, ~, ~] = spm_svd(G_nOPM204 * G_nOPM204', GainSVDTh);
UP_nOPM204 = ug_nOPM204';

G = {G_SQUID, G_mSQUID, G_OPM, G_tOPM, G_nOPM204};
UP = {UP_SQUID, UP_mSQUID, UP_OPM, UP_tOPM, UP_nOPM204};

[N_ch_squid, N_src_2] = size(G_SQUID);
[N_ch_msquid, ~] = size(G_mSQUID);
[N_ch_opm, ~] = size(G_OPM);
[N_ch_topm, ~] = size(G_tOPM);
[N_ch_nopm204, ~] = size(G_nOPM204);
N_ch = {N_ch_squid, N_ch_msquid, N_ch_opm, N_ch_topm, N_ch_nopm204};
N_src = N_src_2 / 2;

N_ch_p = {size(UP_SQUID,1), size(UP_mSQUID, 1),...
          size(UP_OPM,1), size(UP_tOPM, 1), size(UP_nOPM204,1)};

dec = 10;
N_src_dec = ceil(N_src / dec);

% % fix the phase lag 
phi = pi / 2 - pi / 20;
%generate random phase jitters and trials of brain noise (to save time)
PhaseShiftsIn = [];


Rdec = R(1:dec:end,:);

% for mc = 1:100

for mc = 1:100
for i_snr = 1:length(InducedScale)
    disp('MC -------------> ')
    disp(mc)
    dst = 0; 
    while(dst < 0.03)
        ind = fix(1 + rand(2,1) * (size(Rdec,1) - 2));
        dst = norm(Rdec(ind(1),:) - Rdec(ind(2),:))
    end

    for ty = 1:4

        range = 1:2;
        range_dec = 1:2;
        G_dec = zeros(N_ch{ty}, N_src_dec * 2);

        for i_src = 1:N_src_dec
            G_dec(:,range) = G{ty}(:, range_dec);
            range_dec = range_dec + 2 * dec;
            range = range + 2;
        end

        Gp_dec{ty} = UP{ty} * G_dec;
        Gp{ty} = UP{ty} * G{ty};
        % generate data from a network with randomly chosen nodes located @ R(ind,:); 
        %
        [Induced{ty},...
         BrainNoise{ty},...
         SensorNoise{ty},...
         Fs,...
         Ntr,...
         XYZGenOut,...
         Ggen{ty},...
         PhaseShiftsOut{ty}] = SimulateDataPhase_SQUIDvsOPM(Rdec(ind,:),...
                                                            phi,...
                                                            false,...
                                                            PhaseShiftsIn,...
                                                            Rdec, G_dec, 0.25, ty);

         % [C, ~, XYZGenOut] = ups.SimulateData(phi, 100, InducedScale, 0, false, G{ty}, R, UP{ty});

        Data0{ty} = BrainNoise{ty} + InducedScale{i_snr} * Induced{ty};
        % Data0{ty} = ups.ShufflePhases(Data0{ty});
        [bf, af] = butter(5, [8 12] / (0.5 * Fs), 'bandpass');
        % Filter in the band of interest
        Data{ty} = filtfilt(bf, af, Data0{ty}')';
        clear Data0;
        % Reshape the data in a 3D structure(Space x Time x Epochs)
        [Nch{ty}, Tcnt] = size(Data{ty});
        T = fix(Tcnt / Ntr);
        Nch{ty} = size(UP{ty}, 1);
        %reshape Data and store in a 3D array X
        X1{ty} = zeros(Nch{ty}, T, Ntr);
        range = 1:T;
        for i = 1:Ntr
            X1{ty}(:,:,i) = UP{ty} * Data{ty}(:,range);
            range = range + T;
        end;
        %% Calculate band cross-spectral matrix 
        CrossSpecTime{ty} = CrossSpectralTimeseries(X1{ty});
        % CrossSpecTime{ty} = C;
        C = reshape(mean(CrossSpecTime{ty}, 2), N_ch_p{ty}, N_ch_p{ty});


        % [Qidics{ty}, Psidics{ty}, IND{ty}] = iDICS_1D(C, Gp{ty}(:, 1:dec:end));
        % [~, ~, IND{ty}] = iDICS_1D(C, Gp_dec{ty});
        [A, Psidics{ty}, Qidics{ty}, IND{ty}] = ups.DICS((C), Gp_dec{ty}, 1000, true);

        [SPCidics{ty}(mc, i_snr, :),...
         TPRidics{ty}(mc, i_snr, :),...
         PPVidics{ty}(mc, i_snr, :)] = GenerateROC(Qidics{ty}', 0.015,...
                                                   R(1:dec:end,:),...
                                                   IND{ty}, 200, XYZGenOut, 1);


        [A, Psdics{ty}, Qdics{ty}, IND{ty}] = ups.DICS((C), Gp_dec{ty}, 1000, false);

        [SPCdics{ty}(mc, i_snr, :),...
         TPRdics{ty}(mc, i_snr, :),...
         PPVdics{ty}(mc, i_snr, :)] = GenerateROC(Qdics{ty}', 0.015, R(1:dec:end,:),...
                                           IND{ty}, 200, XYZGenOut, 1);


        [Qgcs_dics{ty}, IND{ty}] = ups.GCS_DICS((C), Gp_dec{ty}, 1000);

        [SPCgcs_dics{ty}(mc, i_snr, :),...
         TPRgcs_dics{ty}(mc, i_snr, :),...
         PPVgcs_dics{ty}(mc, i_snr, :)] = GenerateROC(Qgcs_dics{ty}', 0.015,...
                                                      R(1:dec:end,:),...
                                                      IND{ty}, 200, XYZGenOut, 1);

        % [~, ~, ~, Qpsiicos{ty}] = ps.T_PSIICOS(imag(C(:)), Gp_dec{ty}, 0.9, 350, 0, []);

        % [SPCpsiicos{ty}(mc, i_snr, :),...
        %  TPRpsiicos{ty}(mc, i_snr, :),...
        %  PPVpsiicos{ty}(mc, i_snr, :)] = GenerateROC(Qpsiicos{ty}', 0.015, R(1:dec:end,:),...
        %                                      IND{ty}, 200, XYZGenOut, 1);
        end
end

end
