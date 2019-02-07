function MC = Dan_MFCC(d,nmelc)

    sr = 16000;

    % Overlap in specgram
    fft_olap_fct = 0.75; %% .75

    ext_win = 0.1; %% .280

    % only use columns within this proportion of highest E
    magpropthresh = 0.0;

    % Mel cepstrum
    nmelf = 45;
    %nmelf = 34;
    %nmelc = 30;
    %mx = fft2melmx(fftlen,sr,nmelf);
    wmel = 1.0;
    fminmel = 0;
    % Upper limit (as proportion of 8000) of mel filters
    fmaxmel = 8000; %% 3000
    % full-range of Mel scale
    srscale = 1.0; %% 2.0;

    % Find central part
    win_e = round( 0.100 * sr );
    hop_e = round( 0.010 * sr );
    E = env(d, win_e, hop_e);
    [~,xx] = max(E);

    %  d = d/vv;

    midix = 1+hop_e*(xx-1);

    % Segment out 200ms
    %ext_win = 0.200;

    ext_l = round( ext_win * sr );

    % Make sure we don't overrun ends of file
    midix = max(round(ext_l/2),midix);
    midix = min(length(d)-(ext_l-round(ext_l/2)),midix);

    dx = d(-round(ext_l/2) + midix + [1:ext_l]);


    % Calc specgram
    %fftwintime = 0.100;
    fftwintime = 0.050;
    fftlen = 2^round(log(fftwintime*sr)/log(2));

    % Actual specgram
    DX = abs(specgram(dx,fftlen,sr,fftlen,round(fftlen*fft_olap_fct)));
    %DX = abs(specgram(dx,fftlen,sr,fftlen,fftlen*3/4));
    % Average along time
    DS = mean(DX(:,sum(DX)>magpropthresh*max(sum(DX))), 2);

    %DS = DX;
    
    persistent melmx
    if isempty(melmx)
        melmx = fft2melmx(fftlen,sr*srscale,nmelf,wmel,fminmel*srscale,fmaxmel*srscale);
    end
    MS = (melmx * (DS));
    % calcualte cepstrum
    MC = (real(fft(log([MS;MS(end-1:-1:2,:)]))));
    MC = MC(1:nmelc,:)';
    
end

function [wts,binfrqs] = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
    % [wts,frqs] = fft2melmx(nfft, sr, nfilts, width, minfrq, maxfrq, htkmel, constamp)
    %      Generate a matrix of weights to combine FFT bins into Mel
    %      bins.  nfft defines the source FFT size at sampling rate sr.
    %      Optional nfilts specifies the number of output bands required 
    %      (else one per "mel/width"), and width is the constant width of each 
    %      band relative to standard Mel (default 1).
    %      While wts has nfft columns, the second half are all zero. 
    %      Hence, Mel spectrum is fft2melmx(nfft,sr)*abs(fft(xincols,nfft));
    %      minfrq is the frequency (in Hz) of the lowest band edge;
    %      default is 0, but 133.33 is a common standard (to skip LF).
    %      maxfrq is frequency in Hz of upper edge; default sr/2.
    %      You can exactly duplicate the mel matrix in Slaney's mfcc.m
    %      as fft2melmx(512, 8000, 40, 1, 133.33, 6855.5, 0);
    %      htkmel=1 means use HTK's version of the mel curve, not Slaney's.
    %      constamp=1 means make integration windows peak at 1, not sum to 1.
    %      frqs returns bin center frqs.
    % 2004-09-05  dpwe@ee.columbia.edu  based on fft2barkmx

    if nargin < 2;     sr = 8000;      end
    if nargin < 3;     nfilts = 0;     end
    if nargin < 4;     width = 1.0;    end
    if nargin < 5;     minfrq = 0;     end  % default bottom edge at 0
    if nargin < 6;     maxfrq = sr/2;  end  % default top edge at nyquist
    if nargin < 7;     htkmel = 0;     end
    if nargin < 8;     constamp = 0;   end

    if nfilts == 0
      nfilts = ceil(hz2mel(maxfrq, htkmel)/2);
    end

    wts = zeros(nfilts, nfft/2+1);

    % Center freqs of each FFT bin
    fftfrqs = [0:(nfft/2)]/nfft*sr;

    % 'Center freqs' of mel bands - uniformly spaced between limits
    minmel = hz2mel(minfrq, htkmel);
    maxmel = hz2mel(maxfrq, htkmel);
    binfrqs = mel2hz(minmel+[0:(nfilts+1)]/(nfilts+1)*(maxmel-minmel), htkmel);

    for i = 1:nfilts
    %  fs = mel2hz(i + [-1 0 1], htkmel);
      fs = binfrqs(i+[0 1 2]);
      % scale by width
      fs = fs(2)+width*(fs - fs(2));
      % lower and upper slopes for all bins
      loslope = (fftfrqs - fs(1))/(fs(2) - fs(1));
      hislope = (fs(3) - fftfrqs)/(fs(3) - fs(2));
      % .. then intersect them with each other and zero
    %  wts(i,:) = 2/(fs(3)-fs(1))*max(0,min(loslope, hislope));
      wts(i,:) = max(0,min(loslope, hislope));

      % actual algo and weighting in feacalc (more or less)
    %  wts(i,:) = 0;
    %  ww = binbin(i+2)-binbin(i);
    %  usl = binbin(i+1)-binbin(i);
    %  wts(i,1+binbin(i)+[1:usl]) = 2/ww * [1:usl]/usl;
    %  dsl = binbin(i+2)-binbin(i+1);
    %  wts(i,1+binbin(i+1)+[1:(dsl-1)]) = 2/ww * [(dsl-1):-1:1]/dsl;
    % need to disable weighting below if you use this one

    end

    if (constamp == 0)
      % Slaney-style mel is scaled to be approx constant E per channel
      %wts = diag(2./(binfrqs(2+[1:nfilts])-binfrqs(1:nfilts)))*wts;
      wts = repmat(2./(binfrqs(2+[1:nfilts])-binfrqs(1:nfilts))',1, ...
                   nfft/2+1).*wts;
      % 0.14s vs 0.37s
    end

    % Make sure 2nd half of FFT is zero
    %wts(:,(nfft/2+2):nfft) = 0;
    % seems like a good idea to avoid aliasing
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = mel2hz(z, htk)
    %   f = mel2hz(z, htk)
    %   Convert 'mel scale' frequencies into Hz
    %   Optional htk = 1 means use the HTK formula
    %   else use the formula from Slaney's mfcc.m
    % 2005-04-19 dpwe@ee.columbia.edu

    if nargin < 2
      htk = 0;
    end

    if htk == 1
      f = 700*(10.^(z/2595)-1);
    else

      f_0 = 0; % 133.33333;
      f_sp = 200/3; % 66.66667;
      brkfrq = 1000;
      brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
      logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

      linpts = (z < brkpt);

      f = 0*z;

      % fill in parts separately
      f(linpts) = f_0 + f_sp*z(linpts);
      f(~linpts) = brkfrq*exp(log(logstep)*(z(~linpts)-brkpt));

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = hz2mel(f,htk)
    %  z = hz2mel(f,htk)
    %  Convert frequencies f (in Hz) to mel 'scale'.
    %  Optional htk = 1 uses the mel axis defined in the HTKBook
    %  otherwise use Slaney's formula
    % 2005-04-19 dpwe@ee.columbia.edu

    if nargin < 2
      htk = 0;
    end

    if htk == 1
      z = 2595 * log10(1+f/700);
    else
      % Mel fn to match Slaney's Auditory Toolbox mfcc.m

      f_0 = 0; % 133.33333;
      f_sp = 200/3; % 66.66667;
      brkfrq = 1000;
      brkpt  = (brkfrq - f_0)/f_sp;  % starting mel value for log region
      logstep = exp(log(6.4)/27); % the magic 1.0711703 which is the ratio needed to get from 1000 Hz to 6400 Hz in 27 steps, and is *almost* the ratio between 1000 Hz and the preceding linear filter center at 933.33333 Hz (actually 1000/933.33333 = 1.07142857142857 and  exp(log(6.4)/27) = 1.07117028749447)

      linpts = (f < brkfrq);

      z = 0*f;

      % fill in parts separately
      z(linpts) = (f(linpts) - f_0)/f_sp;
      z(~linpts) = brkpt+(log(f(~linpts)/brkfrq))./log(logstep);

    end
    
end

function E = env(X,W,H)
    % E = env(X,W,H)
    %   Find the envelope of a signal as the local RMS.
    %   X is the waveform.  W is the length of the (hann) window; H is
    %   the downsampling hop.
    % 2012-06-30 Dan Ellis dpwe@ee.columbia.edu

    if nargin < 2; W = 256; end
    if nargin < 3; H = W/2; end

    persistent swin
    if isempty(swin)
      swin = hann(W)/sum(hann(W));
    end

    %OLDWAY = 0;
    %if OLDWAY
    %  E = sqrt(conv(swin, X.^2));
    %  % downsample
    %  E = E(floor(W/2)+[1:H:(end-W)]);
    %else
      X2 = [zeros(floor(W/2)-1,1);X.^2;zeros(W-floor(W/2),1)];
      winbase = [0:H:length(X2)-W];
      ixs = repmat(winbase,length(swin),1) ...
            + repmat([1:length(swin)]',1,length(winbase));
      E = sqrt(swin'*X2(ixs));
    %end

    % new way gives identical result down to rounding error (1e-16)

end