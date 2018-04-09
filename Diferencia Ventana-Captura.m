% Diferencia entre técnicas de ventana para captura
clear; close all; clc;
N = 65;                          % Orden de las ventanas
w1 = window(@rectwin,N);         % Ventana Rectangular
w2 = window(@hann,N);            % Ventana Hanning
w3 = window(@hamming,N);         % Ventana Hamming
w4 = winBlackmanHarris(62,N);    % Blackman-Harris 62dB
w5 = winBlackmanHarris(70,N);    % Blackman-Harris 70dB
w6 = winBlackmanHarris(74,N);    % Blackman-Harris 74dB
w7 = window(@blackmanharris,N);  % Blackman-Harris 90dB
% Ventanas en dominio de la frecuencia
nfft = 513;                      % Orden de la FFT
w1f = fft(w1,nfft); w2f = fft(w2,nfft);
w3f = fft(w3,nfft); w4f = fft(w4,nfft);
w5f = fft(w5,nfft); w6f = fft(w6,nfft);
w7f = fft(w7,nfft); 
n = 0 : N-1;                     % Dominio temporal discreto
fmax = N-1;                      % Limite para graficar en bins
f = 0 : fmax-1;                  % Dominio frecuencial en bins
% Cálculo de Magnitudes en dB de cada Ventana
magW1 = 10*log10(abs(w1f(1:fmax)));
magW2 = 10*log10(abs(w2f(1:fmax)));
magW3 = 10*log10(abs(w3f(1:fmax)));
magW4 = 10*log10(abs(w4f(1:fmax)));
magW5 = 10*log10(abs(w5f(1:fmax)));
magW6 = 10*log10(abs(w6f(1:fmax)));
magW7 = 10*log10(abs(w7f(1:fmax)));
% Gráficas en dominio temporal
figure(); subplot(331); plot(n,w1); grid; axis([0 N-1 0 1]);
xlabel('[n]'); ylabel('y[n]'); title('Rectangular');
subplot(332); plot(n,w2); grid; axis([0 N-1 0 1]);
xlabel('[n]'); ylabel('y[n]'); title('Hanning');
subplot(333); plot(n,w3); grid; axis([0 N-1 0 1]);
xlabel('[n]'); ylabel('y[n]'); title('Hamming');
subplot(334); plot(n,w4); grid; axis([0 N-1 0 1]);
xlabel('[n]'); ylabel('y[n]'); title('Blk-Harris 62dB');
subplot(335); plot(n,w5); grid; axis([0 N-1 0 1]);
xlabel('[n]'); ylabel('y[n]'); title('Blk-Harris 70dB');
subplot(336); plot(n,w6); grid; axis([0 N-1 0 1]);
xlabel('[n]'); ylabel('y[n]'); title('Blk-Harris 74dB');
subplot(337); plot(n,w7); grid; axis([0 N-1 0 1]);
xlabel('[n]'); ylabel('y[n]'); title('Blk-Harris 90dB');
saveas(gcf,'Ej11_DominioT.png')
% Gráficas en dominio frecuencial
figure(); subplot(331); plot(f,magW1); grid; axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Rectangular');
subplot(332); plot(f,magW2); grid; axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Hanning');
subplot(333); plot(f,magW3); grid; axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Hamming');
subplot(334); plot(f,magW4); grid; axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Blk-Harris 62dB');
subplot(335); plot(f,magW5); grid; axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Blk-Harris 70dB');
subplot(336); plot(f,magW6); grid; axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Blk-Harris 74dB');
subplot(337); plot(f,magW7); grid; axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Blk-Harris 90dB');
saveas(gcf,'Ej11_DominioF.png')
% Gráfica comparativa en frecuencia
figure(); plot(f,magW1); grid; hold on;
plot(f,magW2); plot(f,magW3); plot(f,magW4); plot(f,magW5);
plot(f,magW6); plot(f,magW7); axis([0 fmax-1 -inf inf]);
xlabel('Bins [k]'); ylabel('dB'); title('Comparativa Ventanas FFT');
legend('Rect','Hann','Hamm','B-H 62dB','B-H 70dB','B-H 74dB','B-H 90dB');
saveas(gcf,'Ej11_Comparativa.png')