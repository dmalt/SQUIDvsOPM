function [Induced, BrainNoise, SensorNoise, Fs,Ntr,XYZGenOut,Ggen, PhaseShiftsOut] = SimulateDataPhase_SQUIDvsOPM(XYZGen,dPhi, bNewBrainNoise,PhaseShiftsIn,R, G, alpha_in)

bUsePhases = ~isempty(PhaseShiftsIn);

if(nargin < 7)
    alpha = 0.05;
else
    alpha = alpha_in;
end;

SSPVE = 0.98;
RAP = 1;% number of RAP-MUSIC iterations

if(nargin == 0)
    bNewBrainNoise = true; % whether or not to generate new brain noise
end;

% create normalized forward matrix, leaving only two components in the tangential plane 
[Ns, Nsites] = size(G);

for i=1:size(XYZGen, 1)
    d = repmat(XYZGen(i,:), size(R, 1), 1) - R;
    d = sum(d .* d, 2);
    [val ind] = min(d);
    XYZGenAct(i,:) = R(ind,:);
    Ggen(:,i) = G(:,ind);
    GenInd(i) = ind; 
end;

nw = [1, 2]; % If you want to simulate all four networks use

Ngen = 2; % total number of generators
sp = zeros(2, 500);
T = 500;   % number of timeslices per trial
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
BrainNoise   = zeros(Ns, Ntr * T);
SensorNoise  = zeros(Ns, Ntr * T);
Induced      = zeros(Ns, Ntr * T);
Evoked       = zeros(Ns, Ntr * T);
F1 = 10; % Hz;
t = linspace(0, 1, Fs);
clear s;

if(~bNewBrainNoise)
    BN = load('/home/asus/MyProjects/SQUIDvsAM_MEG/Data/BrainNoiseBiomag2017.mat');
end;

fprintf('Simulating trial data ...\n');
fprintf('Current trial number (Max %d):', Ntr);
PhaseShiftsOut = zeros(Ntr, 8);
for tr = 1:Ntr
    
    if(bUsePhases)
        phi1 = PhaseShiftsIn(tr, 1);
        phi_alpha = PhaseShiftsIn(tr, 2);
    else
        phi1= 2 * (rand - 0.5) * pi;
        phi_alpha = alpha * (rand - 0.5) * pi;
    end;
    PhaseShiftsOut(tr,1:2) = [phi1 phi_alpha];
    rnd_phi12 = phi_alpha + dPhi;
    s{1}(1,:) =     sin(2 * pi * F1 * t + phi1) .* sp(1,:);
    s{1}(2,:) =     sin(2 * pi * F1 * t + phi1 + rnd_phi12) .* sp(1,:);

     % collect activity from the selected networks
    induced = Ggen * s{1};
    induced = induced / sqrt(sum((induced(:) .^ 2)));
    Induced(:,range) = induced;
    
    
    if(bNewBrainNoise)
        brainnoise = GenerateBrainNoise_SQUIDvsOPM(G, T, 500, 1000, Fs);
    else
        brainnoise = BN.BrainNoise(:, range);
    end;
    
    brainnoise = brainnoise / sqrt(sum((brainnoise(:) .^ 2)));
    BrainNoise(:, range) = brainnoise;
    
    sensornoise = randn(size(brainnoise));
    sensornoise = sensornoise / sqrt(sum((sensornoise(:) .^ 2)));
    SensorNoise(:, range) = sensornoise;
    range = range + T;

    if tr > 1
      for j=0:log10(tr - 1)
          fprintf('\b'); % delete previous counter display
      end
     end
     fprintf('%d', tr);
end

if(bNewBrainNoise)
   save('./BrainNoise2017.mat','BrainNoise');
end;

fprintf('\nDone.\n');

XYZGenOut = XYZGenAct(nw,:);
