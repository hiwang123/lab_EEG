function [TF] = tfa_morlet(td, fs, fmin, fmax, fstep)
%tfa_morlet: time frequency analysis using Morlet Wavelet.
%   [TF] = tfa_morlet(td, fs, fmin, fmax, fstep) returns the time-frequency table
%   calculated by using the Morlet Wavelet, where td is the input temporal data,
%   fs is the sampling rate, [fmin:fstep:fmax] constitute the resulted frequency
%   resolution.
%   
%
%   Yong-Sheng Chen

TF = [];
for fc=fmin:fstep:fmax
    MW = MorletWavelet(fc/fs); % calculate the Morlet Wavelet by giving the central freqency

    MWHL = (length(MW)-1)/2;    % half length of the Morlet Wavelet
    cr = conv(td, MW);
    
    cr = cr(MWHL+1:length(cr)-MWHL);
    
    TF = [TF; abs(cr)];
end

% % e.g.,
% % temp= tfa_morlet(td, 1000, 1, 55, 1);