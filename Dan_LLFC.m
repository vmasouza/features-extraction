function LLFC = Dan_LLFC(d,nlfcb)

    sr = 16000;

    % Overlap in specgram
    fft_olap_fct = 0.75; %% .75

    ext_win = 0.1; %% .280

    % only use columns within this proportion of highest E
    magpropthresh = 0.0;

    % Pitch projection
    lfcbpo = 6; % was 6
    lfcwdt = 1.2; % was 1.2
    lminf = 30; % was 30

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

    % lin-freq cepstrum (to find pitch)
    LFC = (real(fft(log([DS;DS(end-1:-1:2,:)]))));

    % Compress into log-spaced bins
    % Use fake settings on fft2logfmx to map from bin 7/100 up into
    % log-spaced bins
    persistent logfmx
    if isempty(logfmx)
        logfmx = fft2logfmx(size(LFC,1),1000,nlfcb,lminf,lfcbpo,lfcwdt);
    end
    LLFC = (logfmx * (LFC(1:size(LFC,1)/2+1,:)))';

end

function [mx,logffrqs] = fft2logfmx(nfft, sr, nfilts, fmin, bpo, width)
    % [wts,binfrqs] = fft2logfmx(nfft, sr, nfilts, fmin, bpo, width)
    %   Construct a weight mapping matrix to convert a linear frequency
    %   spectrogram into a log-frequency spectrogram.  The underlying
    %   FFT has nfft points (512), at sampling rate sr (16000), the
    %   output will have nfilts bins (60) with the lowest at fmin
    %   (100 Hz), and bpo bins per octave (12).
    %   width is a fact by which to "widen" filters (1.0).
    % 2012-07-02 Dan Ellis dpwe@ee.columbia.edu

    if nargin < 1; nfft = 512; end
    if nargin < 2; sr = 16000; end
    if nargin < 3; nfilts = 72; end
    if nargin < 4; fmin = 100; end
    if nargin < 5; bpo = 12; end
    if nargin < 6; width = 1.0; end

    % Ratio between adjacent frequencies in log-f axis
    fratio = 2^(1/bpo);

    % Freqs corresponding to each bin in FFT
    fftfrqs = [0:(nfft/2)]*(sr/nfft);
    nfftbins = nfft/2+1;

    % Freqs corresponding to each bin in log F output
    logffrqs = fmin * exp(log(2)*[0:(nfilts-1)]/bpo);

    % Bandwidths of each bin in log F
    logfbws = logffrqs * (fratio - 1) * width;

    % .. but bandwidth cannot be less than FFT binwidth
    logfbws = max(logfbws, sr/nfft);

    % Controls how much overlap there is between adjacent bands
    ovfctr = 0.5475;   % Adjusted by hand to make sum(mx'*mx) close to 1.0

    % Weighting matrix mapping energy in FFT bins to logF bins
    % is a set of Gaussian profiles depending on the difference in 
    % frequencies, scaled by the bandwidth of that bin
    freqdiff = ( repmat(logffrqs',1,nfftbins) - repmat(fftfrqs,nfilts,1) )./repmat(ovfctr*logfbws',1,nfftbins);
    mx = exp( -0.5*freqdiff.^2 );
    % Normalize rows by sqrt(E), so multiplying by mx' gets approx orig spec back
    mx = mx ./ repmat(sqrt(2*sum(mx.^2,2)), 1, nfftbins);
    
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