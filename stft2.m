function [y,yfft] = stft2(x,N_Win,N_fft,Hop,Win)
% Función de transformada STFT
% -------------------------------------------------------------------------
% y = señal de salida en dominio del tiempo
% yfft = señal de salida en dominio de frecuencia
% x = señal de entrada en dominio temporal
% N_Win = orden de la ventana a utilizar
% N_fft = orden de la transformada de fourier
% Hop = Factor de solapamiento en %
% Win = Tipo de Ventana a Utilizar (Rect,Hann,Hamm,BH62,BH70,BH74,BH90)
% -------------------------------------------------------------------------
% ***************** Cálculo de parámetros de trabajo **********************
H = round((Hop/100)*N_Win);                       % Solapamiento en Vector
cont_fin = length(x) - N_Win;                     % Ventaneo Final de STFT
% ****************** Inicializar Variables de Salida **********************
x = [zeros(1,N_Win/2-H-1),x];      % Agregar ceros al inicio ajustar vector
y = zeros(1,length(x));            % Vector de salida vacio
yfft = zeros(1,N_fft);            % Vector de salida vacio
% ******************* Escoger Ventana de Análisis *************************
if strcmp(Win,'Rect')
  w = window(@rectwin,N_Win);                     % Ventana Rectangular
elseif strcmp(Win,'Hann')
  w = window(@hann,N_Win);                        % Ventana Hanning
elseif strcmp(Win,'Hamm')
  w = window(@hamming,N_Win);                     % Ventana Hamming
elseif strcmp(Win,'BH62')
  w = winBlackmanHarris(62,N_Win);                % Blackman-Harris 62dB
elseif strcmp(Win,'BH70')
  w = winBlackmanHarris(70,N_Win);                % Blackman-Harris 70dB
elseif strcmp(Win,'BH74')
  w = winBlackmanHarris(74,N_Win);                % Blackman-Harris 70dB
elseif strcmp(Win,'BH90')
  w = window(@blackmanharris,N_Win);              % Blackman-Harris 90dB
end
w = w';
cont = 1;                                         % Iterador de STFT
while cont < cont_fin
  segmentoSenal = x(cont : cont+N_Win-1);   % Segmento de señal a ventanear
  segmentoFFT = segmentoSenal.*w;          % Ventanear el segmento de señal
  X = fft(segmentoFFT,N_fft);
  yfft = yfft+X;
  segmentoFFT = real(ifft(X));
  y(cont : cont+N_Win-1) = y(cont : cont+N_Win-1) + segmentoFFT(1:N_Win);
  cont = cont + H;
end
end