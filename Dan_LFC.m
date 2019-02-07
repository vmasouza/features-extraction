function LFC = Dan_LFC(d,nlfc)

    sr = 16000;

    % Overlap in specgram
    fft_olap_fct = 0.75; %% .75

    ext_win = 0.1; %% .280

    % only use columns within this proportion of highest E
    magpropthresh = 0.0;

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
    LFC = LFC(1:nlfc,:)';
    
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