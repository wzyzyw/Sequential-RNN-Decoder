function rx_bits=demodulation(mod_type,data_mod,use_turbo,noiseVar)
if strcmp(mod_type,'BPSK')
    if use_turbo
        hDemod=comm.BPSKDemodulator('DecisionMethod','Log-likelihood ratio', 'Variance',noiseVar);
    else
        hDemod = comm.BPSKDemodulator();
    end
    bps=1; M=2^bps; % Number of bits per (modulated) symbol
elseif strcmp(mod_type,'QPSK')
    bps=2; M=2^bps; % Number of bits per (modulated) symbol
    if use_turbo
        hDemod=comm.QPSKDemodulator('DecisionMethod','Log-likelihood ratio', 'Variance',noiseVar,'BitOutput',true)
    else
        hDemod = comm.QPSKDemodulator('BitOutput',true);
    end
elseif strcmp(mod_type,'16QAM')
    bps=4;
    M=2^bps;
    if use_turbo
        hDemod=comm.RectangularQAMDemodulator('ModulationOrder',M, 'BitOutput',true,'NormalizationMethod','Average power', 'DecisionMethod','Log-likelihood ratio', 'Variance',noiseVar);
    else
        hDemod=comm.RectangularQAMDemodulator('ModulationOrder',M, 'BitOutput',true,'NormalizationMethod','Average power');
    end
end
rx_bits=step(hDemod,data_mod);
end