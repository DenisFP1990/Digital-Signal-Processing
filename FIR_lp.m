function [audio_out,Hz] = FIR_lp(audio_in,n_audio_out,wp,ws,As)
% Diseño y aplicación de filtro Pasa-Bajos FIR a partir de un audio
% -----------------------------------------
% audio_out = Señal de salida aplicada el filtro
% Hz = Función de transferencia en dominio z
% audio_in = Nombre del audio de entrada que se leerá
% n_audio_out = Nombre del archivo de audio de salida que se escribirá
% wp = frecuencia de corte de banda de paso digital [rad]
% ws = frecuencia de banda de rechazo digital [rad]
% As = atenuación  en la banda de rechazo en dB
  [x,fms] = audioread(audio_in);         % Leer archivo entrada
  Ts = 1/fms;                            % Obtención de tiempo de muestreo
  OmegaP = wp*fms; OmegaS = ws*fms;      % Frecuencias analógicas [rad/s]
  fp=OmegaP/(2*pi); fs=OmegaS/(2*pi);    % Frecuencias analógicas [rad/s]
  delta_w = ws - wp;                     % Ancho banda de transición [rad]
  wc = (ws+wp)/2;                        % Frecuencia de corte ideal
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
  hd = ideal_lp(wc,M);                     % Generar filtro pasa bajos
  h = hd.*win;                             % Producto punto wind y filtro
  audio_out = filter(h,1,x);               % Aplicación de filtro a entrada
  audiowrite(n_audio_out,audio_out,fms);   % Escritura de archivo de salida
  [db,mag,fase,grd,w] = freqz_m(h,1);
  % ********************* Impresión de Datos ******************************
  fprintf('\n* Datos de %s utilizando Filtro FIR Pasa-Bajo',audio_in);
  fprintf('\n* Frecuencia de muestreo fs = %d [Hz]',fms);
  fprintf('\n* tiempo de muestreo Ts = %.3e [seg]',Ts);
  fprintf('\n* OmegaP = %.3f [rad/seg] y OmegaS = %.3f [rad/seg]',OmegaP,OmegaS);
  fprintf('\n* f_pass = %.0f [Hz] y f_stop1 = %.0f [Hz]',fp,fs);
  Hz = filt(h,1,Ts)                       % F. Transferencia Digital
  % *************** Diagrama Polos y Ceros en z ***************************
  figure(); zplane(h,1);                  % Gráfico de Polos y ceros en 'z'
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
  axis([0 10000 -0.1 1.1]); title('Magnitud en valor absoluto');
  xlabel('Hz'); ylabel('|H(z)|');
  subplot(212); plot((w*fms)/(2*pi),fase*180/pi); grid on;
  axis([0 10000 -190 190]); title('Fase en grados');
  xlabel('Hz'); ylabel('Fase H(z)');
  saveas(gcf,'Bode_Hz.png');
end