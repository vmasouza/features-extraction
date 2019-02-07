function [f, sp] = fundFreq_spec(s,~)

    fs = 16000;
    sp = abs(fft(s,fs));
    [~,f] = max(sp(1:floor(fs/2)));

end

