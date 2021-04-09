function data_mod=modulation(mod_type,msgint)
    if strcmp(mod_type,'BPSK')  
        hMod = comm.BPSKModulator();%default:pi/4,gray map
        bps=1; M=2^bps; % Number of bits per (modulated) symbol
    elseif strcmp(mod_type,'QPSK')
        bps=2; M=2^bps; % Number of bits per (modulated) symbol
        hMod = comm.QPSKModulator('BitInput',true);%default:pi/4,gray map
    elseif strcmp(mod_type,'16QAM')
        bps=4;
        M=2^bps;
        hMod=comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average power');
    end
    [row,col]=size(msgint)
    if col~=1
        msgint=msgint';
    end
    data_mod=step(hMod,msgint);
end