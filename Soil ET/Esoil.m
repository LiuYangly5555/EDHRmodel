function obj = Esoil(SWC, SWCr, SWCfc, Va,  RH, T)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
Beta=min(1, (SWC-SWCr)./(SWCfc-SWCr));
obj=0.4.*Va.*Beta.*(1-RH./100).*6.112.*exp(17.62.*T./(243.12+T))/461.5./T*1000/10000*600; %单位cm/10min

% es = 6.112 * exp((17.67 * T)/(T + 243.5));  saturation vapor pressure in mb;
% e = es * (RH/100.0);  vapor pressure in mb;
% p=e/Rv/T; Rv is the specific gas constant for water vapor, 461.5 J kg−1 K−1

