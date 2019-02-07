function f = fundFreq_ceps(s,~)

    fs = 16000;

    C=fft(log(abs(fft(s,fs))+eps));
    [~,q] = max(abs(C(floor(fs/1500):ceil(fs/50))));
    
    f = fs/(q+fs/1500);

end

