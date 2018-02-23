%% template for psk2char.m
function [str,bits] = psk2char(y,M)
% psk2char.m
% Function to decode M-PSK symbols to a string of ASCII text.
% (See also char2psk.m)
%
% Inputs:
% M      size of PSK constellation (2 or 4; or 'bpsk or 'qpsk')
% y      list of PSK IQ points (complex numbers)
% Outputs:
% strbin string variable with ASCII text
% bits   list of 0/1 bits (as logical variables)
%
% Digital Communication Laboratory
% Autumn 2014

%% error checks
if(nargin ~= 2)
    error('Error: psk2char.m requires two input arguments')
end
if (isnumeric(y) ~= 1 || isempty(y))
    error('Error: the input y must be numeric.')
end

%% symbol slicing
% first, symbol string to string of integer labels: 0,1,...,log2(M)
% labels are integers; use double data type
  switch M
     case {2,'bpsk','BPSK'}
        M = 2;
        labels = (real(y) >= 0);
     case {4,'qpsk','QPSK'}
        M = 4;
        for i=1:length(y)
        if real(y(i)) > 0 && imag(y(i))< 0
            labels(i)= 2;
        end
        if real(y(i)) < 0 && imag(y(i)) < 0
            labels(i)= 0;
        end
        if real(y(i)) > 0 && imag(y(i)) > 0
            labels(i)= 3;
        end
        if real(y(i)) < 0 && imag(y(i))> 0
            labels(i)= 1;
        end
        end
        
  end
% second, labels into string of bits
  N = length(labels)*log2(M)/8;%number of characters
  tmp = dec2bin(labels,log2(M));%label to binary
  strbin = reshape(tmp.',8,N).';%array, labels to 8 bits/character
% third, convert from char to logical:
  bits = logical(bin2dec(reshape(tmp.',1,N*8).').');
% fourth, convert binary to decimal ASCII to character
  str=char(bin2dec(strbin)).';%
% end of function
