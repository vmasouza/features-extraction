function v = spectral(s, ~)

    fs = 16000;

    spec = abs(fft(s, fs));    
    spec = (spec(1:fs/2)) + eps;
    
    len = length(spec);
    
    % f0 = Max peak
    f0 = aux_fundFreq(spec,fs);
    %f0 = fundFreq_spec(s, fs);
    %f0 = fundFreq_ceps(s, fs);

    % Sum (Auxiliary)
    sum = 0;
    
    % Sum of squares (Auxiliary)
    sqsum = 0;
    
    % Sum of ln (Auxiliary)
    lnsum = 0;
    
    % Short-time energy
    %ste = 0;
    
    % Magnitude
    %mag = 0;
    
    % Spectral centroid
    sc = 0;
    
    % Spectral irregularity
    si = 0;
    % SI Modified
    %sim = 0;
    
    % Spectral flatness measure
    %sfm = 0;
    
    % Spectral flux
    flux = 0;
    
    for i = 1 : len
      
        sum = sum + spec(i);
        sqsum = sqsum + spec(i)^2;
        lnsum = lnsum + log(spec(i));

        sc = sc + i * spec(i);

        if (i>1 && i<len)
            si = si + abs(20*log10(spec(i)) - (20*log10(spec(i-1))+20*log10(spec(i)+20*log10(spec(i+1))))/3);
        end

        if (i>1)
            flux = flux + (spec(i)-spec(i-1))^2;
        end;
      
    end

    ste = sqsum/len;
    mag = sum/len;

    flux = sqrt(flux);
    
    sc = sc/sum;    
    
    %sim = flux/sqsum;
    sim = log(flux/sqsum); % Log to avoid too small variance
    
    sfm = (exp(lnsum/len) / (sum/len));
    
    % -------------------------------------------------------------------------------------------------------
    
    % Spectral roll-off
    sro_sum = 0.85*sum;
    sro = 0;
    
    % Variance
    var = 0;
    
    % Standard deviation
    %std = 0;
    
    % Skewness
    skw = 0;
    
    % Kurtosis
    kts = 0;
    
    for i = 1 : len
        
        if (sro_sum > 0)
            sro_sum = sro_sum - spec(i);
            sro = i;
        end
        
        aux = (spec(i)-mag).^2;
        var = var + aux;
        
        % (spec(i)-mea).^3
        aux = aux * (spec(i)-mag);
        skw = skw + aux;
        
        % (spec(i)-mea).^4
        aux = aux * (spec(i)-mag);
        kts = kts + aux;
        
    end
    
    var = var / (len-1);
    std = sqrt(var);
    
    aux = std.^3;
    skw = skw / ((len-1) * aux);
    
    % std^4
    aux = aux*std;
    kts = kts / ((len-1) * aux);
    
    [inharm, ts1, ts2, ts3] = aux_harmonics(spec, f0);
    
    v = [f0, ste, mag, sc, si, sim, sfm, flux, sro, var, std, skw, kts, inharm, ts1, ts2, ts3];
    
end


function f0 = aux_fundFreq(spec,fs,freqBeg,freqEnd)

    if (nargin < 3)
        freqEnd = 1500;
    end
    
    if (nargin < 4)
        freqBeg = 50;
    end 
    
    cepsBeg=floor(fs/freqEnd);
    cepsEnd=ceil(fs/freqBeg);
    
    C=fft(log(abs(spec)));
    [~, i] = max(abs(C(cepsBeg:cepsEnd)));
    f0 = fs/(cepsBeg+i);

end


function [inharm, ts1, ts2, ts3, obs, the] = aux_harmonics(spec, f0, nharm)

    spec = abs(spec);
    
    window = floor(f0/2);

    if (nargin < 3)
        nharm = ceil(length(spec)/f0) - 2;
    end
    
    % Theoric harmonics
    the = zeros(1,nharm);
    % Observed harmonics
    obs = zeros(1,nharm);
    
    for i = 1 : nharm
        the(i) = round(f0*i);
        [~, f] = max(spec(floor(f0*i-window):ceil(f0*i+window)));
        obs(i) = round(f+the(i));
    end

    % Sum of energy
    sum = 0;
    
    % Inharmonicity
    inharm = 0;
    
    % Tristimulus
    ts1 = spec(obs(1));
    ts2 = spec(obs(2)) + spec(obs(3)) + spec(obs(4));
    %ts3 = 0;

    for i = 1 : nharm
        sum = sum + spec(obs(i));
        inharm = inharm + (abs(obs(i)-the(i))/the(i));
    end
    
    ts3 = (sum-ts1-ts2)/sum;
    ts1 = ts1/sum;
    ts2 = ts2/sum;

end


