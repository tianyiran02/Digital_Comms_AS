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

function transmitOFDM(data,wavfile)
% transmitOFDM
% inputs:
% data: 2D matrix of uint8 with 256 columns and any number of rows.
% wavfile: file to write the transmitted real signal
%
% Note we don't care what meaning the source data has - we assume it is 
% somehow already encoded into a binary stream, and we will use 256x8bits
% for each symbol we send. Hence the restriction to 2D arrays of uint8 with
% 256 columns is actually not restricting on what data you send.
% we are using energy dispersal so you can safely zero-pad the bytes of the 
% last symbol if your data does not encode to a multiple of 256 bytes.

    % modulation specific information - these values are interdependent so 
    % don't bother changing them unless you REALLY want to have some fun.
    nIFFT = 2048;
    nCyclic = nIFFT/4;
    nTotal = (nIFFT+nCyclic)*320; %NOT USED ELSEWHERE - ELSE ADJUST TO ARB LEN
    pilot_amplitude = 2;
    npilots = 16;
    pilotdistance = 256 / npilots;
    offset = nIFFT; % to have some delay before it starts
    bytes_per_symbol = 257;  
    
    % check that the input data is ok....
    % needs to be 256 columns of uint8, else throws an unrecoverable error
   
    validateattributes(data,{'uint8'},{'2d','ncols',bytes_per_symbol});

    % find out how many symbols we will send (this is the number of rows)
    nSymbols = size(data,1);

    % This offset is only used once per transmission
    signal = zeros(offset,1);
    offsetsignal = zeros(2*offset,1);
    longoffsetsignal = zeros(20*offset,1);
    
    % loop through all our symbols
    for y = 1:nSymbols

        row = data(y,:);

        spectrum = zeros(nIFFT,1);
        % we start with a negative frequency and then work
        % ourselves up to positive ones
        % we place them symmetrical around DC
        k = 1024+512-npilots/2;

        % we make sure that the random number generator always
        % generates the same random numbers for every line
        % this is only guaranteed to work on a single machine
        % i.e. you would need to save the PRBS to be able to
        % read the signal on a different machine
        rand ('seed',1);  

        % we loop through one line in the data
        for x = 1:bytes_per_symbol

            % let's convert one byte into a bitstream
            bitstream = zeros(8,1);
            greyvalue = row(x);

            r = floor( rand(1) * 255 );
            greyvalue = bitxor(greyvalue,r);  %energy dispersal

            for bit = 7:-1:0
                m = bitshift(1,bit);
                testbit = bitand(m,greyvalue);
                if ( testbit > 0 )
                    bitstream(8-bit) = 1;
                end %if
            end %for

            % now we have 8 bits which we distibute over 4 Fourier 
            % coefficients with the QPSK coding

            for cnum = 1:4

                % we add a pilot every 16 steps
                if (mod(k,pilotdistance)==(npilots/2))
                    % pilot signal, make it stronger than the payload data
                    % so that it's easy to recognise
                    spectrum(k) = pilot_amplitude;
                    k = k + 1;
                    if (k >= nIFFT)
                        k = 1;
                    end
                end

                spectrum(k) = bitstream(cnum*2-1)*2-1 + 1i*bitstream(cnum*2)*2-1i;
                k = k + 1;
                if (k >= nIFFT)
                    k = 1;
                end
            end

        end % for x

        % pilot signal, make it stronger than the payload data
        % so that it's easy to recognise
        spectrum(k) = pilot_amplitude;

        % create one symbol
        complex_symbol_dft = ifft(spectrum);
        % re-shuffle time so that negative time comes first, then pos time
        % negative time
        complex_symbol(1:nIFFT/2) = complex_symbol_dft(nIFFT/2+1:nIFFT);
        % postive time
        complex_symbol(nIFFT/2+1:nIFFT) = complex_symbol_dft(1:nIFFT/2);

        tx_symbol = zeros(length(complex_symbol)*2,1);

        % now we upsample at factor 2 and interleave
        % the I and Q signals
        s = 1;
        txindex = 1;
        for symbidx = 1:length(complex_symbol)
            tx_symbol(txindex) = s * real(complex_symbol(symbidx));
            txindex = txindex + 1;
            tx_symbol(txindex) = s * imag(complex_symbol(symbidx));
            txindex = txindex + 1;
            s = s * -1;
        end

        %cyclic prefix taken from the end of the signal
        cyclicPrefix = tx_symbol(nIFFT*2-nCyclic*2+1:nIFFT*2);

        %generate the complete symbol with cyclic prefix
        complete_tx_symbol = cat(1,cyclicPrefix,tx_symbol);

        %merge with the previously generated symbols
        if y == 140
            signal = cat(1,signal,complete_tx_symbol,longoffsetsignal);
        else
            signal = cat(1,signal,complete_tx_symbol,offsetsignal);
        end

    end % for y

    wavwrite(signal,8000,wavfile);
    h=figure;
    plot(linspace(0,1000*length(signal)/6000,length(signal)),signal);
    axis([0 1000*length(signal)/6000 min(signal) max(signal)]);
    title('Transmitted signal');
    xlabel('time/ms');
    ylabel('amplitude');
    saveas(h,'DigitalComm_AS_TRANSsignal','eps');
    
    
end %%% End of transmitOFDM

