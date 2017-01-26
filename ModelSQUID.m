% clear all;
% close all;


import ups.spm_svd
import ups.CrossSpectralTimeseries
import ups.GenerateROC
import ups.ReduceToTangentSpace


InducedScale = 0.3;
% cd('/home/asus/MyProjects/SQUIDvsAM_MEG/Data');

data = load('data_2D.mat') ;

GainSVDTh = 0.0001;

G_OPM = data.data{1}.nOPM.L;
G_mSQUID = data.data{1}.mSQUID.L;
G_SQUID = data.data{1}.SQUID.L;

RL = data.data{1}.sources{1}.p;
RR = data.data{1}.sources{2}.p;
R = [RL; RR];

clear data;

G_OPM = ReduceToTangentSpace(G_OPM, 'all');
G_SQUID = ReduceToTangentSpace(G_SQUID, 'grad');
G_mSQUID = ReduceToTangentSpace(G_mSQUID, 'all');

% Do rank reduction, use the same threshold for the two sensor types 
[ug_OPM, sg_OPM, vg_OPM] = spm_svd(G_OPM * G_OPM', GainSVDTh);
UP_OPM = ug_OPM';

[ug_mSQUID, sg_mSQUID, vg_mSQUID] = spm_svd(G_mSQUID * G_mSQUID', GainSVDTh);
UP_mSQUID = ug_mSQUID';

[ug_SQUID, sg_SQUID, vg_SQUID] = spm_svd(G_SQUID * G_SQUID', GainSVDTh);
UP_SQUID = ug_SQUID';

G = {G_SQUID, G_mSQUID, G_OPM};
UP = {UP_SQUID, UP_mSQUID, UP_OPM};

[N_ch_squid, N_src_2] = size(G_SQUID);
[N_ch_msquid, ~] = size(G_mSQUID);
[N_ch_opm, ~] = size(G_OPM);
N_ch = {N_ch_squid, N_ch_msquid, N_ch_opm};
N_src = N_src_2 / 2;

N_ch_p = {size(UP_SQUID,1), size(UP_mSQUID, 1), size(UP_OPM,1)};

dec = 10;
N_src_dec = ceil(N_src / dec);

% % fix the phase lag 
phi = pi / 2 - pi / 20;
%generate random phase jitters and trials of brain noise (to save time)
PhaseShiftsIn = [];


Rdec = R(1:dec:end,:);

% for mc = 1:100
for mc = 1:100
    dst = 0; 
    while(dst < 0.03)
        ind = fix(1 + rand(2,1) * (size(Rdec,1) - 2));
        dst = norm(Rdec(ind(1),:) - Rdec(ind(2),:))
    end

    for ty = 1:3

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
        % generate data from a network wit randomly chosen nodes located @ R(ind,:); 
        [Induced{ty},...
         BrainNoise{ty},...
         SensorNoise{ty},...
         Fs,...
         Ntr,...
         XYZGenOut,...
         Ggen{ty},...
         PhaseShiftsOut] = SimulateDataPhase_SQUIDvsOPM(Rdec(ind,:),...
                                                        phi,...
                                                        true,...
                                                        PhaseShiftsIn,...
                                                        Rdec, G_dec, 0.25);

         % [C, ~, XYZGenOut] = ups.SimulateData(phi, 100, InducedScale, 0, false, G{ty}, R, UP{ty});

        Data0{ty} = BrainNoise{ty} + InducedScale * Induced{ty};
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
        [A, Psidics{ty}, Qidics{ty}, IND{ty}] = ups.DICS((C), Gp_dec{ty}, 0.1, true);
        % [~, ~, ~, Qidics{ty}] = ps.T_PSIICOS(imag(C(:)), Gp_dec{ty}, 0.9, 350, 0, []);

        [SPCidics{ty}(mc,:),...
         TPRidics{ty}(mc,:),...
         PPVidics{ty}(mc,:)] = GenerateROC(Qidics{ty}', 0.015, R(1:dec:end,:),...
                                           IND{ty}, 200, XYZGenOut, 1);
        end
end

figure
col = {'c', 'm', 'y'};
for ty = 1:3
    plot(1 - mean(SPCidics{ty}, 1), mean(TPRidics{ty}, 1), col{ty})
    hold on;
end;
legend('SQUID (gradiometers)', 'mSQUID', 'OPM')
xlabel('1 - specificity');
ylabel('sensitivity');

figure
col = {'c', 'm', 'y'};
for ty = 1:3
    plot(mean(TPRidics{ty}, 1), mean(PPVidics{ty}, 1), col{ty})
    hold on;
end;
legend('SQUID (gradiometers)', 'mSQUID', 'OPM')
xlabel('recall');
ylabel('precision');
