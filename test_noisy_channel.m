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
close all
clear all
if 1
    data = imread('greytee.png');
    wavfile='testtxfunc.wav';
    transmitOFDM(data,wavfile);

    %clean ideal channel (never happens in real life!!)
    rx_data_clean = receiveOFDM(wavfile);

    SNR = 5; %dB
    output_wavfile = AWGNchannel(wavfile,5);
    rx_data = receiveOFDM(output_wavfile);
    figure();
    subplot(2,2,1)
    imagesc(data,[0 255]);
    colormap(gray);
    title('transmitted')
    subplot(2,2,2)
    imagesc(rx_data,[0 255]);
    colormap(gray);
    title('received - noisy')
    subplot(2,2,3)
    imagesc(rx_data_clean,[0 255]);
    colormap(gray);
    title('received - ideal')
    subplot(2,2,4)
    imagesc(sign(abs(rx_data - rx_data_clean)));
    colormap(gray);
    title('error locations')
end
% starting point for the lab
if 1
    %now lets get some stats on the number of errors
    output_wavfile0 = AWGNchannel(wavfile,0);
    output_wavfile3 = AWGNchannel(wavfile,3);
    output_wavfile6 = AWGNchannel(wavfile,6);
    output_wavfile9 = AWGNchannel(wavfile,9);
    output_wavfile12 = AWGNchannel(wavfile,12);
    rx_data0 = receiveOFDM(output_wavfile0);
    rx_data3 = receiveOFDM(output_wavfile3);
    rx_data6 = receiveOFDM(output_wavfile6);
    rx_data9 = receiveOFDM(output_wavfile9);
    rx_data12 = receiveOFDM(output_wavfile12);
    figure();
    [BitsInError0,BER0] = biterror(rx_data0,rx_data_clean);
    [BitsInError3,BER3] = biterror(rx_data3,rx_data_clean);
    [BitsInError6,BER6] = biterror(rx_data6,rx_data_clean);
    [BitsInError9,BER9] = biterror(rx_data9,rx_data_clean);
    [BitsInError12,BER12] = biterror(rx_data12,rx_data_clean);

    theSNR = [0, 3, 6, 9]
    link_spectral_efficiency = 1; %units (bit/s)/Hz
    EboverN0 = 10.^(theSNR/10)./link_spectral_efficiency
    BERa = 0.5*erfc(sqrt(EboverN0))

    BERm = [BER0,BER3,BER6,BER9];

    figure();
    semilogy(theSNR,BERa,theSNR,BERm)
    legend('theoretical','sim')
end %if

%Then do the parity check thing
% which can ARQ anything it spots.

data = imread('greytee.png');

data_parity = data .* 0;

for ii = 1:size(data,1)
    for jj = 1:size(data,2)

        Nbits = 0;
        %sacrifice the least significant bit to be a parity bit
        for Nbit = 2:8  %count number of 1s in the rest of the byte
           Nbits = Nbits + bitget(data(ii,jj),Nbit);   
        end %for Nbit
        %choose even parity
        data_parity(ii,jj) = bitset(data(ii,jj),1,mod(Nbits,2)); 
    end
end

wavfile='testparity.wav';
transmitOFDM(data_parity,wavfile);

%clean ideal channel (never happens in real life!!)
rx_data_clean = receiveOFDM(wavfile);

SNR = 5; %dB
dirty_wavfile = AWGNchannel(wavfile,5);
rx_data_dirty = receiveOFDM(dirty_wavfile);

clean_parity_error = rx_data_clean .* 0;
rx_data_clean_stripped_parity = rx_data_clean .* 0;

for ii = 1:size(rx_data_clean,1)
    for jj = 1:size(rx_data_clean,2)

        Nbits = 0;
        %sacrifice the least significant bit to be a parity bit
        for Nbit = 2:8  %count number of 1s in the rest of the byte
           Nbits = Nbits + bitget(rx_data_clean(ii,jj),Nbit);   
        end %for Nbit
        clean_parity_error(ii,jj) =mod(Nbits,2); %0 if correctly rx'd even parity
        rx_data_clean_stripped_parity(ii,jj) = ...
            bitand(rx_data_clean(ii,jj),254);
    end
end

dirty_parity_error = rx_data_dirty .* 0;
rx_data_dirty_stripped_parity = rx_data_dirty .* 0;

for ii = 1:size(rx_data_dirty,1)
    for jj = 1:size(rx_data_dirty,2)

        Nbits = 0;
        %sacrifice the least significant bit to be a parity bit
        for Nbit = 2:8  %count number of 1s in the rest of the byte
           Nbits = Nbits + bitget(rx_data_dirty(ii,jj),Nbit);   
        end %for Nbit
        dirty_parity_error(ii,jj) = mod(Nbits,2); %0 if correctly rx'd even parity
        rx_data_dirty_stripped_parity(ii,jj) = ...
                bitand(rx_data_dirty(ii,jj),254); 
        if dirty_parity_error(ii,jj) > 0
            %assume we ARQ until we get error-free transmission and
            % that we don't encounter any multibit errors that could 
            % trick us into thinking we had good data (life isn't that
            % easy, but to model imperfect retransmission would add 
            % extra complexity to this model. Try it if you like....
            rx_data_dirty_ARQ(ii,jj) = ...
                rx_data_clean_stripped_parity(ii,jj);
        else
            rx_data_dirty_ARQ(ii,jj) = ...
                rx_data_dirty_stripped_parity(ii,jj);
        end
                      
    end
end

figure()
subplot(2,2,1)
imagesc(data_parity,[0 255])
colormap('gray')
title('transmitted 7-bit image')
subplot(2,2,2)
imagesc(rx_data_clean_stripped_parity,[0 255])
colormap('gray')
title('received over ideal channel')
subplot(2,2,3)
imagesc(rx_data_dirty_stripped_parity,[0 255])
colormap('gray')
title('received over 5dB SNR AWGN channel')
subplot(2,2,4)
imagesc(rx_data_dirty_ARQ,[0 255])
title('ideal ARQ over 5dB SNR AWGN channel')
colormap('gray')