function sp = spectrum(s,~)
    fs = 16000;
    sp = abs(fft(s, fs));
    sp = sp(1:fs/2)';
end
