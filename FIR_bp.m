function [audio_out,Hz] = FIR_bp(audio_in,n_audio_out,fs1,fp1,fp2,fs2,As)
% Diseño y aplicación de filtro Pasa-Banda FIR a partir de un audio
% -----------------------------------------
% audio_out = Señal de salida aplicada el filtro
% Hz = Función de transferencia en dominio z
% audio_in = Nombre del audio de entrada que se leerá
% n_audio_out = Nombre del archivo de audio de salida que se escribirá
% fp1 = frecuencia de corte de banda de paso 1 [Hz]
% fs1 = frecuencia de banda de rechazo 1 [Hz]
% fp2 = frecuencia de corte de banda de paso 2 [Hz]
% fs2 = frecuencia de banda de rechazo 2 [Hz]
% As = atenuación  en la banda de rechazo en dB
  [x,fms] = audioread(audio_in);         % Leer archivo entrada
  Ts = 1/fms;                            % Obtención de tiempo de muestreo
  OmegaP1=2*pi*fp1; OmegaS1=2*pi*fs1;    % Frecuencia en rad/seg HP
  OmegaP2=2*pi*fp2; OmegaS2=2*pi*fs2;    % Frecuencia en rad/seg LP
  wp1 = OmegaP1/fms; ws1 = OmegaS1/fms;  % Frecuencia digital en rad
  wp2 = OmegaP2/fms; ws2 = OmegaS2/fms;  % Frecuencia digital en rad
  delta_wHP = wp1 - ws1;                 % banda de transición HP[rad]
  delta_wLP = ws2 - wp2;                 % banda de transición LP [rad]
  delta_w = min(delta_wLP,delta_wHP);    % banda de transición Global
  wc1 = (ws1+wp1)/2;                     % Frecuencia de corte ideal
  wc2 = (ws2+wp2)/2;                     % Frecuencia de corte ideal
  if(As >= 0 && As <= 21)                % Aplicar ventana Rectangular
    fprintf('Ventana Aplicada: Rectangular');
    M = ceil(1.8*pi/delta_w);            % Calcular orden de ventana
    win = (boxcar(M))';
  elseif(As > 21 && As <= 25)            % Aplicar ventana Barlett
    fprintf('Ventana Aplicada: Bartlett');
    M = ceil(6.1*pi/delta_w);            % Calcular orden de ventana
    win = (bartlett(M))';
  elseif(As > 25 && As <= 44)            % Aplicar ventana Hanning
    fprintf('Ventana Aplicada: Hanning');
    M = ceil(6.2*pi/delta_w);            % Calcular orden de ventana
    win = (hann(M))';
  elseif(As > 44 && As <= 53)            % Aplicar ventana Hamming
    fprintf('Ventana Aplicada: Hamming');
    M = ceil(6.6*pi/delta_w);            % Calcular orden de ventana
    win = (hamming(M))';
  elseif(As > 53 && As <= 74)            % Aplicar ventana Blackman
    fprintf('Ventana Aplicada: Blackman');
    M = ceil(11*pi/delta_w);             % Calcular orden de ventana
    win = (blackman(M))';
  end
  hd_HP = ideal_hp(wc1,M);                 % Generar filtro pasa altos
  hd_LP = ideal_lp(wc2,M);                 % Generar filtro pasa bajos
  h_HP = hd_HP.*win;                       % Producto punto wind y filtro1
  h_LP = hd_LP.*win;                       % Producto punto wind y filtro2
  Hz1 = filt(h_HP,1,Ts);                   % F. Transferencia Digital 1
  Hz2 = filt(h_LP,1,Ts);                   % F. Transferencia Digital 2
  Hz = Hz1*Hz2                             % F. Transferencia Digital BP
  audio_out = filter(cell2mat(Hz.num),1,x);% Aplicación de filtro a entrada
  audiowrite(n_audio_out,audio_out,fms);   % Escritura de archivo de salida
  [db,mag,fase,grd,w] = freqz_m(cell2mat(Hz.num),1);
  % ********************* Impresión de Datos ******************************
  fprintf('\n* Datos de %s utilizando Filtro FIR Pasa-Banda',audio_in);
  fprintf('\n* Frecuencia de muestreo fs = %d [Hz]',fms);
  fprintf('\n* tiempo de muestreo Ts = %.3e [seg]',Ts);
  fprintf('\n* wp1 = %.3f [rad] y ws1 = %.3f [rad]',wp1,ws1);
  fprintf('\n* wp2 = %.3f [rad] y ws2 = %.3f [rad]',wp2,ws2);
  fprintf('\n* OmegaP1=%.3f [rad/s] y OmegaS1=%.3f [rad/s]',OmegaP1,OmegaS1);
  fprintf('\n* OmegaP2=%.3f [rad/s] y OmegaS2=%.3f [rad/s]',OmegaP2,OmegaS2);
  fprintf('\n* f_pass1 = %.0f [Hz] y f_stop1 = %.0f [Hz]',fp1,fs1);
  fprintf('\n* f_pass2 = %.0f [Hz] y f_stop2 = %.0f [Hz]',fp2,fs2);
  % *************** Diagrama Polos y Ceros en z ***************************
  figure(); zplane(cell2mat(Hz.num),1);   % Gráfico de Polos y ceros en 'z'
  title('Diagrama Polos y Ceros H(z)'); grid on; 
  saveas(gcf,'polos_zeros_Hz.png');
  % *************** Bode de Filtro Discreto *******************************
  figure(); subplot(221); plot(w/pi,mag); grid on; axis([0 1 -0.1 1.1]);
  title('Magnitud en abs'); xlabel('w/pi'); ylabel('|H(z)|');
  subplot(222); plot(w/pi,fase*180/pi); grid on;
  title('Fase en grados'); xlabel('w/pi'); ylabel('Fase H(z)');
  subplot(223); plot(w/pi,db); grid on; axis([0 1 -100 1]);
  title('Magnitud en dB'); xlabel('w/pi'); ylabel('|H(z)| [dB]');
  subplot(224); plot(w/pi,grd); grid on;
  title('Retardo de grupo'); xlabel('w/pi'); ylabel('retardo de grupo');
  saveas(gcf,'Bode_digital.png');
  % *************** Bode de Filtro Continuo *******************************
  figure(); subplot(211); plot((w*fms)/(2*pi),mag); grid on;
  axis([0 15000 -0.1 1.1]); title('Magnitud en valor absoluto');
  xlabel('Hz'); ylabel('|H(z)|');
  subplot(212); plot((w*fms)/(2*pi),fase*180/pi); grid on;
  axis([0 10000 -190 190]); title('Fase en grados');
  xlabel('Hz'); ylabel('Fase H(z)');
  saveas(gcf,'Bode_Hz.png');
end