clear all;
close all;
InducedScale = 10;
cd('/home/asus/MyProjects/SQUIDvsAM_MEG/Data');
data = load('data.mat') ;
GainSVDTh = 0.0001;
G_OPM = data.data{1}.nOPM.L;
G_SQUID = data.data{1}.mSQUID.L;
RL = data.data{2}.sources{1}.p;
RR = data.data{2}.sources{2}.p;
R = [RL; RR];
% Do rank reduction, use the same threshold for the two sensor types 
[ug_OPM sg_OPM vg_OPM] = spm_svd(G_OPM*G_OPM',GainSVDTh);
UP_OPM = ug_OPM';
[ug_SQUID sg_SQUID vg_SQUID] = spm_svd(G_SQUID*G_SQUID',GainSVDTh);
UP_SQUID = ug_SQUID';
% fix the phase lag 
phi = pi/2-pi/20;
%generate random phase jitters and trials of brain noise (to save time)
PhaseShiftsIn = [];
dec  = 10;

[Induced, BrainNoise, SensorNoise, Fs,Ntr,XYZGenOut,Ggen, PhaseShiftsOut]  = SimulateDataPhase_SQUIDvsOPM(R([1,2],:), phi, false, PhaseShiftsIn,R, G_OPM);

Rdec = R(1:dec:end,:);

%for mc = 1:100
    
    dst = 0; 
    while(dst<0.12)
        ind = fix(1+rand(2,1)*(size(Rdec,1)-2));
        dst = norm(Rdec(ind(1),:)-Rdec(ind(2),:));
    end
    % generate data from a network wit randomly chosen nodes located @ R(ind,:); 
    [Induced, BrainNoise, SensorNoise, Fs,Ntr,XYZGenOut,Ggen, PhaseShiftsOut]  = SimulateDataPhase_SQUIDvsOPM(Rdec(ind,:), phi, false, PhaseShiftsIn,R, G_OPM);
    
    Data0 = BrainNoise + InducedScale*Induced;
    [bf af] = butter(5,[8 12]/(0.5*Fs),'bandpass');
    % Filter in the band of interest
    Data = filtfilt(bf,af,Data0')';
    % Reshape the data in a 3D structure(Space x Time x Epochs)
    [Nch Tcnt] = size(Data);
    T = fix(Tcnt/Ntr);
    Nch = size(UP_OPM,1);
    %reshape Data and store in a 3D array X
    X1 = zeros(Nch,T,Ntr);
    range = 1:T;
    for i=1:Ntr
        X1(:,:,i) = UP_OPM*Data(:,range);
        range = range+T;
    end;
    %% Calculate band cross-spectral matrix 
    CrossSpecTime = CrossSpectralTimeseries(X1);
    C = reshape(mean(CrossSpecTime,2),Nch,Nch);
    Gp_OPM = UP_OPM*G_OPM;
    [Qidics, Psidics, IND] = iDICS_1D(C,Gp_OPM(:,1:dec:end));
    
    [SPCidics(1,:), TPRidics(1,:),PPVidics(1,:),NWCidics(1,:)] = GenerateROC(Qidics,0.001,R(1:dec:end,:), IND, 100,XYZGenOut,1);
    figure,
    plot(1-SPCidics, TPRidics)
    
    
    
end

