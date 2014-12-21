%
% every 128 Bytes will generate a symbols, which is means every line
% will be transimitted by 2 symbol.
%
%
%
%
%% Transmitter
% clear all;
% close all;
% clc;

% imagestream = imread('greytee.png');
% colu = 256;
% row = 100;
% rowDoThisTime =100;
% fs=8000;
% 
% figure(1);
% subplot(211);
% imshow(imagestream);
% title('Original Picture');
% 
% 
% % convert the pic matrix to bit stream, row is the same with 
% % the row in pic, data conver to bin
% imagebitstream = zeros(rowDoThisTime*2,colu/2*8);
% 
% for i=1:rowDoThisTime
%     for j=1:colu/2
%         imagebitstream(2*i-1,(8*(j-1)+1):8*j) = dec2bin(imagestream(i,j),8) - '0';
%     end
%     for j=1:colu/2
%         imagebitstream(2*i,(8*(j-1)+1):8*j) = dec2bin(imagestream(i,j+colu/2),8) - '0';
%     end
% end
% 
% % encode the bitstream
% % besides the extra 4 bits at ending (encoding from adding '0 0'), 4 more 
% % '0' bits adding to make it a Byte
% imagebitEncode = zeros(rowDoThisTime*2,colu/2*8*2+8);
% 
% for i=1:rowDoThisTime*2
%    imagebitEncode(i,1:colu/2*8*2+4) = encode(imagebitstream(i,:)); 
% end
% 
% % convert back to Byte
% imageByteEncode = zeros(rowDoThisTime*2,colu+1);
% for i=1:rowDoThisTime*2
%     for j=1:colu+1
%         imageByteEncode(i,j) = bi2de(imagebitEncode(i,1+8*(j-1):8*j),'left-msb');
%     end
% end
% imageByteEncode = cast(imageByteEncode,'uint8');
% 
% % generate sound wav
% transmitOFDM(imageByteEncode,'OFDMtest100_8000_R1.wav');

%% Receiver Side

% output_wavfile = AWGNchannel('OFDMtest2.wav',0.2);
close all;
% receive sound wav
nIFFT = 2048;
nCyclic = nIFFT/4;
npilots = 16;
bytes_per_symbol = 257;
pilotdistance = 256 / npilots;
[signal,fs] = wavread('OFDM100_8000_R1.wav');
signal = signal/max(abs(signal)); % normalize the signal



%% start bit detect with cross correlation (time-saving)

% cyclic prefix parameters
chunk1 = zeros(nCyclic*2,1);
chunk2 = zeros(nCyclic*2,1);
s = 0;
c1 = zeros(1,length(signal)-(nIFFT*2+nCyclic*2-1));

% detect all symbol start in for loop
for s=1:(length(signal)-(nIFFT*2+nCyclic*2-1))
    chunk1 = signal(s:(s+nCyclic*2-1));
    chunk2 = signal((s+nIFFT*2):(s+nIFFT*2+nCyclic*2-1));
    c1(s) = (chunk1' * chunk2)^2;  
end

%% start bit detect with cross correlation 
sbuffer = [];
prefixTH = 200;
getflag = 0;
getavoid = 0;
temp = 0;
temp1 = 0;
temp2 = 0;
cvalue = 0;

for s=1:(length(signal)-(nIFFT*2+nCyclic*2-1))

    cvalue = c1(s);
    
    if getavoid == 0
        if ((cvalue >= prefixTH) && (getflag == 0))
            temp2 = s;
            getavoid = 1000;
            getflag = 1;
        end

        if ((cvalue <= prefixTH) && (getflag == 1))
            [temp1,temp] = max(c1(temp2:s));
            temp = temp + temp2 - 1;
            sbuffer = cat(2,sbuffer,temp-80,temp+80);
            getavoid = 4000;
            getflag = 0;        
        end
    else
        getavoid = getavoid -  1;
    end
end

h=figure();
plot(linspace(0,1000*length(c1)/fs,length(c1)),c1);
axis([0 1000*length(c1)/fs min(c1) max(c1)]);
title('Cross-Correlation Result');
xlabel('time/ms');
ylabel('amplitude');
saveas(h,'DigitalComm_AS_CCResult','eps');

%% start bit detect with pilot tone

% detect the precise start point
% parameters for this procedure
fftbuffer = zeros(1,nIFFT*2);
fftdownsample = zeros(1,nIFFT);
fftdownsample_ts = zeros(1,nIFFT);
fftResult = zeros(1,nIFFT);
imagsum = 0;
imagsumbuffer = [];
counter = 0;
temp = 0;
Startbit = [];

for s = 1:length(sbuffer)/2
    while (sbuffer(2*s-1) + counter) <= sbuffer(2*s) 
        % % % load the original chunk
        fftbuffer = signal((sbuffer(2*s-1) + nCyclic * 2 + counter):(sbuffer(2*s-1) + nCyclic * 2 + counter + 2 * nIFFT - 1)); 
        
        % % % down sample it
        temp = 1; 
        for a = 1:(nIFFT)
          realpart = temp * fftbuffer(2*a-1);
          imagpart = temp * fftbuffer(2*a);
          fftdownsample(a) = realpart + 1i * imagpart;
          temp = temp * -1;
        end
        
        % % % fft
        fftResult = fft(fftdownsample); % do the fft to the chunk
        
        % % % get the summation of imaginary part
        imagsum = 0;
        
        k = 1024+512-npilots/2;
        
        for pilotnum = 1:69
            imagsum = imagsum + abs(imag(fftResult(k))); 
            k = k + pilotdistance;
               
            if (k >= nIFFT)
                k = 8;
            end
        end       
        imagsumbuffer = cat(2,imagsumbuffer,imagsum); % save the value to the buffer
        counter = counter + 1; 
    end
    [temp1, temp] = min(imagsumbuffer); % get the minimum number's index number
    Startbit = cat(2,Startbit,(sbuffer(2*s-1) + temp -1 + nCyclic*2)); % get the index for start bit (not including CP)
    
    if (s<=4)
        h = figure;
        subplot(211);
        plot(linspace(0,1000*length(imagsumbuffer)/fs,length(imagsumbuffer)),imagsumbuffer);
        axis([0 1000*length(imagsumbuffer)/fs min(imagsumbuffer) max(imagsumbuffer)]);
        title('Start bit detection with pilot tone (summation of pilot tone imaginary part)');
        xlabel('time/ms');
        ylabel('amplitude');
        subplot(212);
        plot(linspace(0,1000*length(imagsumbuffer)/fs,length(imagsumbuffer)),imagsumbuffer);
        axis([temp/8-4 temp/8+4 240 1100]);
        title('Start bit detection with pilot tone with detail');
        xlabel('time/ms');
        ylabel('amplitude');
        saveas(h,'DigitalComm_AS_StartPilot','eps');
    end
    
    imagsumbuffer = [];
    counter = 0;
end

%% angle shift for every symbol

% parameters
nSymbol = size(Startbit,2); % get the number of symbol
fftdownsample = zeros(1,nIFFT);
pilotshiftResult = zeros(nSymbol,nIFFT);
diff = 0;

for i=1:nSymbol
    % % % load symbol
    fftbuffer = signal(Startbit(i):(Startbit(i) + 2 * nIFFT -1)); 
    
    % % % down sample it
    temp = 1; 
    for a = 1:(nIFFT)
      realpart = temp * fftbuffer(2*a-1);
      imagpart = temp * fftbuffer(2*a);
      fftdownsample(a) = realpart + 1i * imagpart;
      temp = temp * -1;
    end
    
    % % % time shift
    % negative time
    fftdownsample_ts(nIFFT/2+1:nIFFT) = fftdownsample(1:nIFFT/2);
    % positive time
    fftdownsample_ts(1:nIFFT/2) = fftdownsample(nIFFT/2+1:nIFFT);
    
    % % % fft
    fftResult = fft(fftdownsample_ts); % do the fft to the symbol
    pilotshiftResult(i,:) = fftResult; % initial pilotshiftResult first
    
    % % % shift angle according to the piolt tone shifting
    k = 1024+512-npilots/2;
    for pilotnum = 1:(68 - 1)
        if (k + pilotdistance) > nIFFT
            k = 8;
        end
        
        % for pilot tone
        absValue = abs(pilotshiftResult(i,k));
        angleValue = 0;
        pilotshiftResult(i,k) = absValue * exp(1i * angleValue); 
        
        % for others
        diff = (angle(fftResult(k+pilotdistance)) - angle(fftResult(k)))/pilotdistance;
        
        for j = 1:(pilotdistance-1)   
            absValue = abs(pilotshiftResult(i,k+j));
            angleValue = angle(fftResult(k+j)) + diff;
            pilotshiftResult(i,k+j) = absValue * exp(1i * angleValue);            
        end
        k = k + pilotdistance;
    end
end

%% get the transmit data

rx_data = zeros(nSymbol,bytes_per_symbol);

for i = 1:nSymbol
    
    rand ('seed',1);
    % start with setting point
    k = 1024+512-npilots/2;
    
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
          bitstream(cnum*2-1) = heaviside(real(pilotshiftResult(i,k)));
          % second bit is in the imag part of the coefficient
          bitstream(cnum*2) = heaviside(imag(pilotshiftResult(i,k)));
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
        rx_data(i,x) = greyvalue;

    end % for x, bytes_per_symbol
end

%% decoding
% convert to bit
rowReceiveThisTime = size(rx_data,1)/2;
RXimagebit = zeros(2*rowReceiveThisTime,colu * 8 + 8);
buffer = zeros(1,8);

for i=1:rowReceiveThisTime*2
    for j=1:size(rx_data,2)
        buffer(:) = dec2bin(rx_data(i,j),8) - '0';
        RXimagebit(i,(8*(j-1)+1):8*j) = buffer(:);
    end    
end

% decode signal
RXimagebitDecode = zeros(rowReceiveThisTime*2,colu/2 * 8);
RXimageByteDecode = zeros(rowReceiveThisTime*2,colu/2);
Bytebuf = zeros(1,8);

for i=1:rowReceiveThisTime*2
    RXimagebitDecode(i,:) = decode(RXimagebit(i,1:(size(RXimagebit,2)-4)));
end

RXimagebitDecode = cast(RXimagebitDecode,'uint8');
    
for i=1:rowReceiveThisTime*2 
    for j=1:colu/2
        RXimageByteDecode(i,j) = bi2de(RXimagebitDecode(i,1+8*(j-1):8*j),'left-msb');
    end
end

RXimage = zeros(rowReceiveThisTime,colu);

for i=1:rowReceiveThisTime
    RXimage(i,1:128) = RXimageByteDecode(2*i-1,:);
    RXimage(i,129:256) = RXimageByteDecode(2*i,:);
end

RXimage = cast(RXimage,'uint8');

figure();
subplot(212);
imshow(RXimage);
h=figure(1);
saveas(h,'DigitalComm_AS_ImgRecover','eps');

