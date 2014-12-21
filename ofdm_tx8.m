%	OFDM transmission and reception with energy dispersal
%
%	Copyright (C) 2014 Bernd Porr <mail@berndporr.me.uk>

%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU Lesser General Public License as published by
%	the Free Software Foundation; either version 2 of the License, or
%	(at your option) any later version.

%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.

%	You should have received a copy of the GNU Lesser General Public License
%	along with this program; if not, write to the Free Software
%	Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

%       using more than one pilots to eliminate ambiguities.

% load the cheesy image
a = imread('greytee.png');

% the number of lines in the image
ymax = size(a)(1);
xmax = size(a)(2);

% the plan is that we generate for every line of the image
% one symbol, thus, one IFFT

% our inverse FFT has 2048 frequencies
nIFFT = 2048;

% distance between pilots
pilot_distance = 16;

% amplitudes of the pilot carrier at the beginning
pilot_amplitude = 2;

% to have some delay before it starts
offset = nIFFT*3;

% total number of complex samples in the timedomain
nTotal = nIFFT*ymax+offset;

% first k index used
k_start = 1024+512-1024/pilot_distance/2;

% some dummy bytes before we start transmission
signal = zeros(offset,1);

% this is our complex data stream
tx_complex_signal = [];

% loop for the y coordinate, line by line
for y = 1:ymax

  % get a line from the image
  row = a(y,:);

  % create an empty spectrum with all complex frequency values set to zero
  spectrum = zeros(nIFFT,1);

  % we start with a negative frequency and then 
  % work ourselves up to positive ones
  % we have 2048 frequency samples and 1024 frequencies we use for the
  % image
  k = k_start;

  % set the random number generator to a known start value
  % will generate always the same sequence from this start value
  % We xor its value with the grey values from the image to
  % generate a pseudo random sequence which is called "engery dispersal".
  rand ('seed',1);

  % pilot signal, make it stronger than the payload data
  % so that it's easy to recognise
  % we can use it to fine tune the symbol start in the receiver
  % However, only one pilot won't be enough in a real receiver
  % because of its periodic nature. We need to scatter them
  % over the spectrum.

  % counter for the pilots
  pilot_counter = pilot_distance/2;

  % we loop through one line in the image
  for x = 1:xmax

    % lets convert one grey value from the image into a bitstream
    bitstream = zeros(8,1);
    % get the grey value
    greyvalue = row(x);
    % generate the random number
    r = floor( rand(1) * 255 );
    % xor the grey value with the random number
    greyvalue = bitxor(greyvalue,r);

    % create the bitstream
    for bit = 7:-1:0
      m = bitshift(1,bit);
      testbit = bitand(m,greyvalue);
      if ( testbit > 0 )
	bitstream(8-bit) = 1;
      end %if
    end %for

    % now we have 8 bits which we distibute over four frequency samples
    % with 4-QAM / QPSK coding
    for cnum = 1:4

      pilot_counter = pilot_counter - 1;
      if (pilot_counter == 0)
	spectrum(k) = pilot_amplitude;
	k = k + 1;
	if (k >= nIFFT)
	  % avoiding DC
	  k = 2;
        end	
	pilot_counter = pilot_distance;
      end

      spectrum(k) = bitstream(cnum*2-1)*2-1 + i*bitstream(cnum*2)*2-i;
      % increase the frequency index
      k = k + 1;
      % wrap to positive frequencies once we have reached the last index
      if (k >= nIFFT)
	k = 1;
      end
    end

  end % for x

  % create one symbol by transforming our frequency samples into
  % complex timedomain samples
  complex_symbol = ifft(spectrum);

  % add the complex symbol to our complex data stream
  tx_complex_signal = [ tx_complex_signal ; complex_symbol ];

end % for y


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reception


% our image
rx_image = zeros(ymax,xmax);

% number of samples of the symbol
nSymbol = nIFFT;

rx_complex_signal = tx_complex_signal;

% loop for the y coordinate
for y = 1:ymax

  % lets retrieve one complex symbol
  rx_symbol = rx_complex_signal((y-1)*nSymbol+1:(y-1)*nSymbol+nIFFT);
  
  % perform a FFT to get the frequency samples which code our signal as QPSK pairs
  isymbol = fft(rx_symbol);

  % set the random number generator to the same value as in the transmitter so that
  % we have exactly the same sequence
  rand ('seed',1);

  % we start at frequency index 512
  k = k_start;

  % counter for the pilots
  pilot_counter = pilot_distance/2;

  % we loop through one line in the image
  for x = 1:xmax

    % decode one byte from 4 bytes in the FFT
    % we first create an array which contains the bits in separate rows
    bitstream = zeros(8,1);
    % loop through four bytes in the fft 
    for cnum = 1:4

      % test for pilots and ignore
      pilot_counter = pilot_counter - 1;
      if (pilot_counter == 0)
	k = k + 1;
	if (k >= nIFFT)
          % hopping over the DC
	  k = 2;
        end	
	pilot_counter = pilot_distance;
      end

      % first bit is in the real part of the coefficient
      bitstream(cnum*2-1) = heaviside(real(isymbol(k)));
      % second bit is in the imag part of the coefficient
      bitstream(cnum*2) = heaviside(imag(isymbol(k)));
      % get the next FFT coefficient
      k = k + 1;
      % we wrap to positive frequencies
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

    % store it in the image
    rx_image(y,x) = greyvalue;

    % now we have 8 bits which we distibute over 4 Fourier coefficients
    % with the QPSK coding

  end % for x

end % for y

imagesc(rx_image);
colormap(gray);
