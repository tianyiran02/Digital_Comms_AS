%   Noisy OFDM channel with QPSK and AWGN
%
%   Copyright (C) 2014 Tim Drysdale <tim.drysdale@glasgow.ac.uk>

%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU Lesser General Public License as 
%   published by the Free Software Foundation; either version 2 of the
%   License, or (at your option) any later version.

%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.

%	You should have received a copy of the GNU Lesser General Public License
%	along with this program; if not, write to the Free Software
%	Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

%   OFDM transmission and reception modified from Bernd's code:-
%	-- OFDM transmission and reception with energy dispersal
%   -- Copyright (C) 2013 Bernd Porr <mail@berndporr.me.uk>
%   -- OFDM transmission features
%   -- complex symbols
%   -- modulated and demodulated
%   -- frequencies in the FFT from neg to positive ones so that the 
%      spectrum is at at the upsampled fs/4
%   -- this is the version for just sending it and receiving it in OCTAVE
%   -- energy dispersal via xor with pseudorandom bitstream

% There are three main functions
% transmitOFDM(data,wavfile)
% output_wavfile = AWGNchannel(wavfile)
% rx_data = receiveOFDM(output_wavfile)

%Example of their use

% a = imread('greytee.png');
% wavfile('testtxfunc');
% transmitOFDM(data,wavfile);
% output_wavfile = AWGNchannel(wavfile);
% rx_data = receiveOFDM(output_wavfile);
% figure();
% imagesc(rx_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel modelling with AWGN (this is not a complete model of a real
% comms channel, but AWGN is often considered a good start, especially
% for a general case where specific fading mechanisms are not specified)
% Since the signal occupies the full band obtained by FFT, then
% there is no need to explicitly bandlimit or filter the noise.
% the noise is applied to meet a signal to noise ratio (SNR) criteria

function output_wavfile = AWGNchannel(wavfile,SNR)

    clean_signal = wavread(wavfile);
    desiredSNR = SNR; %dB
    plot_stuff = 1;

    % append '_dirtySNRXX' to the clean filename where XX... is the SNR
    dotspot = max(strfind(wavfile,'.'));
    if length(dotspot) == 0
        dotspot = length(wavfile);
    end
    dirty_wavfile = strcat(wavfile(1:dotspot-1),'_dirtySNR', ...
        num2str(desiredSNR),wavfile(dotspot:end));

    sigp = 10*log10(norm(clean_signal,2)^2/numel(clean_signal));
    snr = sigp-desiredSNR;
    noisep = 10^(snr/10);
    scaled_noise = sqrt(noisep)*randn(size(clean_signal));
    noiseFFT = fft(scaled_noise);

    dirty_signal = clean_signal + scaled_noise;

    wavwrite(dirty_signal,dirty_wavfile);

    output_wavfile = dirty_wavfile;

    if plot_stuff > 0
        nIFFT = 2048; %assume we still use the same channel
        offset = nIFFT * 3;
        Fs = 1; %normalise to Fs
        f = Fs/2*linspace(-1,1,length(clean_signal));       
        h=figure;
        subplot(3,2,1);
        plot(clean_signal(offset:offset+nIFFT./8));
        YL1 = ylim;
        ylabel('V(t)')
        xlabel('Time/dT')
        title('clean')
        subplot(3,2,5);
        plot(dirty_signal(offset:offset+nIFFT./8),'m');
        ylabel('V(t)')
        xlabel('Time/dT')
        title('dirty = signal + noise')
        subplot(3,2,2);
        plot(f,abs(fft(clean_signal)));
        ylabel('|V(f)|')
        xlabel('frequency/Fs')
        title('clean')
        subplot(3,2,6);
        plot(f,abs(fft(dirty_signal)),'m'); 
        YL2 = ylim;
        ylabel('|V(f)|')
        xlabel('frequency/Fs')
        title('dirty = signal + noise')
        subplot(3,2,3)
        plot(scaled_noise(offset:offset+nIFFT./8),'r');
        ylim(YL1);
        ylabel('V(t)')
        xlabel('Time/dT')
        title('noise')
        subplot(3,2,4)
        plot(f,abs(fft(scaled_noise)),'r')
        ylabel('|V(f)|')
        xlabel('frequency/Fs')
        title('noise')
        ylim(YL2);
    end

    saveas(h,'DigitalComm_AS_AWGNResult','eps');
    % Other noise aspects to consider
    % Burst noise - implentation thoughts :-
    % Add high SNR noise for about the length of a 1/256 - 1/128 of a 
    % symbol thus writing off a whole byte or two at a time, requiring
    % interleaving choose timing of bursts randomly, but seed so as to 
    % avoid the need for a full statistical analysis over multiple
    % transmission attempts, for 
    % a given burst noise occurence rate and severity level.

    % mimic the real channel from the assignment - 
    % add multipath, and frequency selective fading, but not Doppler

end % AWGNchannel

