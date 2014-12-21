function hd = hamming(a,b)
% hd = hamming(a,b)
% hd is the Hamming distance between a and b, where a and b are 
% 1D vectors of 1s and 0s. This is NOT a bitlevel function.
% example usage
% a = [0,0];
% b = [0,1];
% hd = hamming(a,b);
% 
% (hd should be equal to 1 in this case)

% Tim Drysdale 19 Nov 2014

hd = sum(xor(a,b));