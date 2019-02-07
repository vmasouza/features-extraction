function Y = Dan_AutoCor(d,naccf)
    
    sr = 16000;
    xctime = 0.010;
    xclen = round(xctime * sr);
    xc = xcorr(d,xclen);
    xc = xc(xclen + [1:xclen]);
    Y = (xc(1+[1:naccf])/xc(1))';
    
end
