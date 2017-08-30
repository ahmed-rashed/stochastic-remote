function NoisySignal=addNoise(signal,SNR_db)

noise=std(signal)*sqrt(db2mag(SNR_db))*randn(size(signal)); %randn produces noise whose std = 1
NoisySignal=signal+noise;