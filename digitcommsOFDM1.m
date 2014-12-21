clear all;
close all;
clc;

imagestream = imread('greytee.png');
colu = 256;
row = 100;
rowDoThisTime = 4;

flowest = 20000;
fstep = 100;

fmax = flowest + 100 * (2048 - 1);
fs = 2 * (fmax + flowest);

fnlowest = flowest/fs;
fnstep = fstep/fs;

figure(1);
subplot(211);
imshow(imagestream);
title('Original Picture');

% add random to the pic
for i=1:row
    for j=1:colu
        imagestream(i,j) = bitxor(imagestream(i,j) , floor(rand(1) * 255));
    end
end

% convert the pic matrix to bit stream, row is the same with 
% the row in pic, data conver to bin
imagebitstream = zeros(rowDoThisTime,colu * 8);
buffer = zeros(colu,8);

for i=1:rowDoThisTime
    buffer = dec2bin(imagestream(i,:)) - '0';
    for j=1:colu
        imagebitstream(i,(8*(j-1)+1):8*j) = buffer(j,:);
    end
end

% add cyclic bit to the bit stream
cyclic_bit = zeros(1,uint32(colu/4)*8);
imbit_stream = zeros(rowDoThisTime,length(imagebitstream(1,:)) + uint32(colu/4)*8);

for i=1:rowDoThisTime
    cyclic_bit = imagebitstream(i,(colu*8 - uint32(colu/4)*8 + 1):colu*8);
    imbit_stream(i,:) = [cyclic_bit imagebitstream(i,:)];
end

% generate Fourier coefficient
% the size of it is 1coe/2bits + 1tone/16bits
% txcoe = zeros(rowDoThisTime,length(imbit_stream)/2 + uint32(length(imbit_stream)/16));
txcoe = zeros(rowDoThisTime,2048);
k = 1;
for i=1:rowDoThisTime
   for j=1:length(txcoe(1,:))
       if mod(j,8) == 0
           txcoe(i,j) = 1;
       elseif (2*k) <= length(imbit_stream(1,:))
           txcoe(i,j) = (-1)^(imbit_stream(i,(2*k-1)) + 1) + 1i * (-1)^(imbit_stream(i,2*k) + 1);
           k = k + 1;
       end
   end
   k = 1;
end

% apply IDFT
txBuf = zeros(rowDoThisTime,length(txcoe(1,:)));
for i=1:rowDoThisTime
   txBuf(i,:) = ifft(txcoe(i,:)); 
end