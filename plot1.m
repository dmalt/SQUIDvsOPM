m_opm = max(sg_OPM);
m_squid = max(sg_SQUID);
m_msquid = max(sg_mSQUID);

ss_opm = diag(sg_OPM / m_opm(1));
ss_squid = diag(sg_SQUID / m_squid(1));
ss_msquid = diag(sg_mSQUID / m_msquid(1));


T = 34

plot(1:T, ss_opm(1:T))
hold on;
plot(1:T, ss_squid(1:T));
hold on
plot(1:T, ss_msquid(1:T));

legend('OPM', 'SQUID (gradiometers)', 'mSQUID');
