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
% reception
function rx_data = receiveOFDM(received_signal_file)

    % construct filename for datafile output with extension .mat
    dotspot = max(strfind(received_signal_file,'.'));
    if length(dotspot) == 0
        dotspot = length(received_signal_file);
    end
    rx_datafile = strcat(received_signal_file(1:dotspot-1),'_rx_data.mat');

    % modulation specific information - these values are interdependent so 
    % don't bother changing them unless you REALLY want to have some fun.
    % also, note that these values are possibly hardcoded in some places
    % so take to debug with this in mind, if you want to change the modulation.
    nIFFT = 2048;
    nCyclic = nIFFT/4;
    nTotal = (nIFFT+nCyclic)*320; %not used elsewhere but change for arb length
    pilot_amplitude = 2;
    npilots = 16;
    pilotdistance = 256 / npilots;
    offset = nIFFT*3; % to have some delay before it starts
    bytes_per_symbol = 257;  
    %%% end of channel info %%%

    signal = wavread(received_signal_file);  

    %now we work out how many symbols were sent 
    %(you can't do this in real life, you just wait til it finishes...)
    nSymbols = (length(signal)-offset)./(nIFFT+nCyclic)./2; %note downsampling
    rx_data = zeros(nSymbols,bytes_per_symbol);
    % also note in real life we have no a-priori information on the start time
    % so here we are still getting an easy life on the synchronisation problem.

    nSymbol = nIFFT+nCyclic;  %note that we have two similarly named variables

    % now we downsample by factor 2
    s = 1;
    rxindex = 1+offset;
    for a = 1:(length(signal)/2-offset/2)
      realpart = s * signal(rxindex);
      rxindex = rxindex + 1;
      imagpart = s * signal(rxindex);
      rxindex = rxindex + 1;
      rx_complex_signal(a) = realpart + i * imagpart;
      s = s * -1;
      end

    % loop for the number of symbols we will receive
    for y = 1:nSymbols

      rx_symbol = ...
        rx_complex_signal((y-1)*nSymbol+nCyclic+1:(y-1)*nSymbol+nCyclic+nIFFT);

      % negative time
      rx_symbol_dft(nIFFT/2+1:nIFFT) = rx_symbol(1:nIFFT/2);
      % positive time
      rx_symbol_dft(1:nIFFT/2) = rx_symbol(nIFFT/2+1:nIFFT);

      isymbol = fft(rx_symbol_dft);

      rand ('seed',1); 
      %if you change machines between transmitting and receiving, you will 
      %not be able to recover the original bitstream, because you get a 
      %different sequence of random numbers on different machines.
      %you can of course change machines if you do like the real world and
      %agree a PRBS in advance. How? Consider a simple description that
      % allows  a long PRBS to be described compactly, e.g. shift register
      % and XOR feedback - just like EE1X Christmas lights (except you will
      % need more stages). 
      % Historical trivia - the
      % creators of the original Elite game on the BBC micro had 
      % insufficient
      % memory to store their vast universe of galaxies, so they used 
      % Fibonacci sequences to predictably choose galaxy details from
      % custom-made 'dictionaries';  

      % we start at frequency index 512
      k = 1024+512-npilots/2;

      % we loop through one line in the image
      for x = 1:bytes_per_symbol

        % decode one byte from 4 bytes in the FFT
        % we first create an array which contains the bits in separate rows
        bitstream = zeros(8,1);
        % loop through four bytes in the fft 
        for cnum = 1:4

          % need to skip the pilots
          if (mod(k,pilotdistance)==(npilots/2))
            k = k + 1;
              if (k >= nIFFT)
          k = 1;
              end
           end

          % first bit is in the real part of the coefficient
          bitstream(cnum*2-1) = heaviside(real(isymbol(k)));
          % second bit is in the imag part of the coefficient
          bitstream(cnum*2) = heaviside(imag(isymbol(k)));
          % get the next FFT coefficient
          k = k + 1;
          if (k >= nIFFT)
        k = 1;
          end
        end

        % now let's assemble the bits into into a proper byte by
        % using bit-wise or
        greyvalue = 0;

        % let's loop through the bits
        for bit = 7:-1:0
          mask = bitshift(1,bit);
          if (bitstream(8-bit) == 1)
        greyvalue = bitor(mask,greyvalue);
          end %if
        end %for

        % de-scramble the byte
        r = floor( rand(1) * 255 );
        greyvalue = bitxor(greyvalue,r);

        % store it in the image, trimming the parity bit
        rx_data(y,x) = greyvalue;

        end % for x, bytes_per_symbol

    end % for y, nSymbols

    save(rx_datafile,'rx_data');

end % receiveOFDM
