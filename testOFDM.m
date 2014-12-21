close all
clear all
tx_data = imread('greytee.png');
wavfile='testOFDM.wav';
transmitOFDM(tx_data,wavfile);
%channel model goes here
%SNR = 5;% dB. This is pretty noisy.
%wavfile = AWGNchannel(wavfile,SNR);
rx_data = receiveOFDM(wavfile);
figure();
subplot(2,2,1)
imagesc(tx_data,[0 255]);
colormap('gray')
title('transmitted image')
subplot(2,2,2)
imagesc(rx_data,[0 255]);
colormap('gray')
title('received image')
subplot(2,2,3)
error_loc = abs(sign(double(tx_data)-rx_data));
imagesc(error_loc,[0 1]);
colormap('gray')
title('error locations shown in white')

