function data_out = encode(data_in)
% new encoder!
% Tim Drysdale 19 Nov 2014
n0 = 2;
k0 = 1;
L =  3;
G = [ 1, 0, 1;
      1, 1, 1;];

data_in = [data_in, 0, 0];  
  
mem = [0,0,0];

data_out = []; 

sum = 0;

for ii=1:length(data_in)
    % shift 
    mem(3) = mem(2);
    mem(2) = mem(1);
    mem(1) = data_in(ii);
    %output

    sum = 0;
    for j = 1:3
        sum = sum + G(1,j) * mem(j);
    end
    U1 = mod(sum,2);
    
    sum = 0;
    for j = 1:3
        sum = sum + G(2,j) * mem(j);
    end
    U2 = mod(sum,2);    
    
    data_out = cat(1,data_out,U1,U2);
end

    
  