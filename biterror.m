function [BitsInError,BER] = biterror(rx_data,rx_data_clean)
% BitsInError,BER = biterror(rx_data,rx_data_clean)

% T Drysdale Nov 2014
    BitsInError = zeros(1,8);
    [X,Y] = find(sign(abs(rx_data - rx_data_clean))>0);
    for ii = 1:length(X)
        bitdiff = bitxor(rx_data(X(ii),Y(ii)),rx_data_clean(X(ii),Y(ii)));
        Nbits = 0;
        for Nbit = 1:8
           Nbits = Nbits + bitget(bitdiff,Nbit);   
        end %for Nbit
        BitsInError(Nbits)=BitsInError(Nbits)+1;
    end
    BER = sum(BitsInError)./ (size(rx_data,1)*size(rx_data,2)*8);
end