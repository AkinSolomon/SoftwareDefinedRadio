function [] = eyediagram(y_up,L,N,IQ)
% eyediagram.m
% Function to plot an eye diagram
%
% y_up    up-sampled matched filter outputs
% L       downsampling factor (integer)
% N       number of symbols
% IQ      string, 'complex', if both I and Q are desired 
%
% Note that first sample of y_up is assumed to occur at a symbol time.

% Digital Communication Laboratory
% Autumn 2014

%% error checks
if(nargin < 2)
    error('Error: eyediagram.m requires two input arguments')
end
if(nargin == 3)
    IQ='real';
end
if(strcmp(IQ,'complex'))
    IQ=1;%plot both I and Q
else
    IQ=0;
end
%N-2 segments for N symbols

%% extract symbol period length segments
% start of first full interval
start = 1 + ceil(L/2);
Y_up = reshape(y_up(start:start+(N-2)*L-1),L,N-2); % extract segments

%% plots
if(IQ == 0)
    % plot eye diagram, I channel only
    figure;
    plot(-floor(L/2)/L+(0:L-1)/L,real(Y_up));% superimpose them
    title('Eye Diagram','fontsize',13);
    xlabel('relative symbol index','fontsize',13);
else
    % plot eye diagram, I & Q channels
    figure;
    subplot(211);plot(-floor(L/2)/L+(0:L-1)/L,real(Y_up));
    title('Eye Diagram','fontsize',13);
    subplot(212);plot(-floor(L/2)/L+(0:L-1)/L,imag(Y_up));
    xlabel('relative symbol index','fontsize',13);
end


