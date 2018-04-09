% Ejercicio de Ejemplo de Uso de la Función STFT
clear; close all; clc;
f = 11000;
t = 0 : 1/f : 1;
f = 0 : 5/512 : 5-5/512;
x = sin(2*pi*5*t);
M = 64;
N = 512;
hop = 50;
metodo = 'BH90';
[y,yfft] = stft2(x,M,N,hop,metodo);
subplot(311); plot(t,x); xlabel('tiempo [n]'); ylabel('y[n]');
title('Señal Original');
subplot(312); plot(f,abs(yfft)); xlabel('frecuencia [Hz]'); ylabel('H[w]');
title('STFT');
subplot(313); plot(t,y); xlabel('tiempo [n]'); ylabel('y[n]');
title(strcat('Señal Sintentizada Usando Ventana: ',metodo));
saveas(gcf,'STFT_42.png');