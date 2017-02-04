import ups.spm_svd
data = load('data_2D.mat') ;

GainSVDTh = 0.0001;

G_OPM = data.data{1}.nOPM.L;
G_mSQUID = data.data{1}.mSQUID.L;
G_SQUID = data.data{1}.SQUID.L;

[ug_OPM, sg_OPM, vg_OPM] = spm_svd(G_OPM * G_OPM', GainSVDTh);

[ug_mSQUID, sg_mSQUID, vg_mSQUID] = spm_svd(G_mSQUID * G_mSQUID', GainSVDTh);

[ug_SQUID, sg_SQUID, vg_SQUID] = spm_svd(G_SQUID * G_SQUID', GainSVDTh);

plot(sg_OPM);
