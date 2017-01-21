clear all;
close all;


import ups.spm_svd
import ups.CrossSpectralTimeseries
import ups.GenerateROC


InducedScale = 1;
% cd('/home/asus/MyProjects/SQUIDvsAM_MEG/Data');

data = load('data.mat') ;

GainSVDTh = 0.0001;
G_OPM = data.data{1}.nOPM.L;
G_SQUID = data.data{1}.mSQUID.L;
RL = data.data{2}.sources{1}.p;
RR = data.data{2}.sources{2}.p;
R = [RL; RR];
% Do rank reduction, use the same threshold for the two sensor types 
%
[ug_OPM, sg_OPM, vg_OPM] = spm_svd(G_OPM*G_OPM', GainSVDTh);
UP_OPM = ug_OPM';
[ug_SQUID, sg_SQUID, vg_SQUID] = spm_svd(G_SQUID*G_SQUID', GainSVDTh);
UP_SQUID = ug_SQUID';
G = {G_SQUID, G_OPM};
UP = {UP_SQUID, UP_OPM};

% fix the phase lag 
phi = pi / 2 - pi / 20;
%generate random phase jitters and trials of brain noise (to save time)
PhaseShiftsIn = [];
dec = 10;

% [~, ~, ~, ~, ~,...
%  XYZGenOut,...
%  ~,...
%  PhaseShiftsOut] = SimulateDataPhase_SQUIDvsOPM(R([1,2],:),...
%                                                 phi,...
%                                                 true,...
%                                                 PhaseShiftsIn,...
%                                                 R,...
%                                                 G_SQUID);

Rdec = R(1:dec:end,:);

%for mc = 1:100
    
dst = 0; 
while(dst < 0.1)
    ind = fix(1 + rand(2,1) * (size(Rdec,1) - 2));
    dst = norm(Rdec(ind(1),:) - Rdec(ind(2),:));
end

for ty = 1:2
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
                                                    R, G{ty});

    Data0{ty} = BrainNoise{ty} + InducedScale * Induced{ty};
    [bf af] = butter(5,[8 12] / (0.5 * Fs), 'bandpass');
    % Filter in the band of interest
    Data{ty} = filtfilt(bf, af, Data0{ty}')';
    % Reshape the data in a 3D structure(Space x Time x Epochs)
    [Nch{ty}, Tcnt] = size(Data{ty});
    T = fix(Tcnt / Ntr);
    Nch{ty} = size(UP{ty}, 1);
    %reshape Data and store in a 3D array X
    X1{ty} = zeros(Nch{ty}, T, Ntr);
    range = 1:T;
    for i=1:Ntr
        X1{ty}(:,:,i) = UP{ty} * Data{ty}(:,range);
        range = range + T;
    end;
    %% Calculate band cross-spectral matrix 
    CrossSpecTime{ty} = CrossSpectralTimeseries(X1{ty});
    C = reshape(mean(CrossSpecTime{ty}, 2), Nch{ty}, Nch{ty});
    Gp{ty} = UP{ty} * G{ty};

    [Qidics{ty}, Psidics{ty}, IND{ty}] = iDICS_1D(C, Gp{ty}(:, 1:dec:end));

    [SPCidics{ty}(1,:),...
     TPRidics{ty}(1,:),...
     PPVidics{ty}(1,:)] = GenerateROC(Qidics{ty}, 0.01, R(1:dec:end,:),...
                                      IND{ty}, 100, XYZGenOut, 1);

end

figure
col = {'r','b'}
for ty = 1:2
    plot(1 - SPCidics{ty}, TPRidics{ty}, col{ty})
    hold on;
end;
legend('SQUID','OPM')
