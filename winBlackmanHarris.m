function w = winBlackmanHarris(dB,M)
% Ventana Blackman Harris de M términos
% --------------------------------
% dB = 62 dB, 70dB, 74dB, 92 dB
% M = Número de muestras
switch dB
  case 62
    a0=0.44859;
    a1=0.49364;
    a2=0.05677;
    a3=0;
  case 70
    a0=0.42323;
    a1=0.49755;
    a2=0.07922;
    a3=0;
  case 74
    a0=0.402217;
    a1=0.49703;
    a2=0.09892;
    a3=0.00188;
  case 92
    a0=0.35875;
    a1=0.48829;
    a2=0.14128;
    a3=0.01168;
  otherwise
    a0=0;
    a1=0;
    a2=0;
    a3=0;
end
n = 0 : M-1;    % Dominio temporal
w=(a0-a1*cos(2*n*pi/M)+a2*cos(4*n*pi/M)-a3*cos(6*n*pi/M));
w = w';