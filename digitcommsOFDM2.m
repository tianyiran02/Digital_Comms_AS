%
% every 128 Bytes will generate a symbols, which is means every line
% will be transimitted by 2 symbol.
%
%
%
%
%% Transmitter
clear all;
close all;
clc;

imagestream = imread('greytee.png');
colu = 256;
row = 100;
rowDoThisTime = 50;

figure(1);
subplot(211);
imshow(imagestream);
title('Original Picture');


% convert the pic matrix to bit stream, row is the same with 
% the row in pic, data conver to bin
imagebitstream = zeros(rowDoThisTime*2,colu/2*8);

for i=1:rowDoThisTime
    for j=1:colu/2
        imagebitstream(2*i-1,(8*(j-1)+1):8*j) = dec2bin(imagestream(i,j),8) - '0';
    end
    for j=1:colu/2
        imagebitstream(2*i,(8*(j-1)+1):8*j) = dec2bin(imagestream(i,j+colu/2),8) - '0';
    end
end

% encode the bitstream
% besides the extra 4 bits at ending (encoding from adding '0 0'), 4 more 
% '0' bits adding to make it a Byte
imagebitEncode = zeros(rowDoThisTime*2,colu/2*8*2+8);

for i=1:rowDoThisTime*2
   imagebitEncode(i,1:colu/2*8*2+4) = encode(imagebitstream(i,:)); 
end

% convert back to Byte
imageByteEncode = zeros(rowDoThisTime*2,colu+1);
for i=1:rowDoThisTime*2
    for j=1:colu+1
        imageByteEncode(i,j) = bi2de(imagebitEncode(i,1+8*(j-1):8*j),'left-msb');
    end
end
imageByteEncode = cast(imageByteEncode,'uint8');

% generate sound wav
transmitOFDM(imageByteEncode,'OFDMtest.wav');
%% Receiver Side

 dirty_signal = AWGNchannel('OFDMtest.wav',10);

% receive sound wav
rx_data = receiveOFDM(dirty_signal);

% convert to bit
RXimagebit = zeros(2*rowDoThisTime,colu * 8 + 8);
rowReceiveThisTime = size(rx_data,1)/2;
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
subplot(212)
imshow(RXimage)
