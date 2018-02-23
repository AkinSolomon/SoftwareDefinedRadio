function [ rx_Trim ] = packetdetect( yup , L,  pkt_len)
%PACKETDETECT  energy detector to approximately locate a data packet;
%              trims a received upsampled signal
%
% [rx_Trim] = packetdetect(yup,L,pkt_len)
%
% yup       fractionally sampled matched filter output 
% L         upsampling factor
% pkt_len   packet length (including both header and payload)

% Digital Communications Laboratory
% Autumn 2014

Mask=ones(1,pkt_len*L);
Energy_Corr=conv(abs(yup).^2,Mask);
Energy_Corr=Energy_Corr(length(Mask):end);
[~,Max_Energy_Index]=max(Energy_Corr);
End_Packet_Index=min(length(yup),Max_Energy_Index+(pkt_len+3)*L);
rx_Trim=yup(Max_Energy_Index-3*L-round(L*rand(1,1)) : ...
    End_Packet_Index);
%for educational use:
%randomized offset to test frame timing;
%presumes at least 4L noise samples precede packet
end