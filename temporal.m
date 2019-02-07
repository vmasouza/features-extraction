function v = temporal(signal,~)

    fs = 16000;
    len = length(signal);
    
    % Short-time Energy
    ste = 0;
    
    % Magnitude average
    mag = 0;
        
    % Root mean square
    %rms = 0;
    
    % Mean
    mean = 0;
    
    % Temporal centroid
    tc = 0;
    
    % Zero-crossing rate
    zcr = 0;
    
    % Interval
    itv_max = 0;
    itv_min = 0;
    %itv = 0;
    
    % Complexity estimate
    ce = 0;
    
    for i = 1 : len
        
      ste = ste + signal(i).^2;
      mag = mag + abs(signal(i));
      mean = mean + signal(i);
      
      tc = tc + i * abs(signal(i));
      
      if (i>1)
          zcr = zcr + abs(sign(signal(i)) - sign(signal(i-1)));
          ce = ce + (signal(i) - signal(i-1)).^2;
      end
      
      if (signal(i)>itv_max)
          itv_max = signal(i);
      end
      if (signal(i)<itv_min)
          itv_min = signal(i);
      end
      
    end

    ste = ste/len;
    rms = sqrt(ste);
    mag = mag/len;
    mean = mean/len;
    
    tc = (tc/mag);
    
    zcr = zcr / (len-1);

    itv = itv_max - itv_min;
    
    ce = sqrt(ce);
    
    % -------------------------------------------------------------------------------------------------------
    
    % Variance
    var = 0;
    
    % Standard deviation
    %std = 0;
    
    % Skewness
    skw = 0;
    
    % Kurtosis
    kts = 0;
    
    for i = 1 : len
        
        aux = (signal(i)-mean)^2;
        var = var + aux;
        
        % (signal(i)-mea).^3
        aux = aux * (signal(i)-mean);
        skw = skw + aux;
        
        % (signal(i)-mea).^4
        aux = aux * (signal(i)-mean);
        kts = kts + aux;        
        
    end
    
    var = var / (len-1);
    std = sqrt(var);
    
    aux = std.^3;
    skw = skw / ((len-1) * aux);
    
    % std^4
    aux = aux*std;
    kts = kts / ((len-1) * aux);
    
    v = [ste,mag,rms,mean,tc,zcr,itv,var,std,skw,kts,ce,len/fs]; % last att = duration

end

function p = sign(value)
    if value < 0
        p=0;
    else
        p=1;
    end
end
