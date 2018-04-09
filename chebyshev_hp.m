function [audio_out,Hs,Hz] = chebyshev_hp(audio_in,n_audio_out,fp,fs,Rp,As)
% Diseño y aplicación de filtro Pasa-Altos Chebyshev-1 a partir de un audio
% -----------------------------------------
% audio_out = chebyshev1_lp(audio_in,n_audio_out,wp,ws,Rp,As,metodo);
% audio_out = Señal de salida aplicada el filtro
% Hs = Función de transferencia en dominio s
% Hz = Función de transferencia en dominio z
% audio_in = Nombre del audio de entrada que se leerá
% n_audio_out = Nombre del archivo de audio de salida que se escribirá
% fp = frecuencia de corte de banda de paso [Hz]
% fs = frecuencia de banda de rechazo [Hz]
% Rp = Rizado de banda pasante en dB
% As = atenuación  en la banda de rechazo en dB
  [x,fms] = audioread(audio_in);         % Leer archivo entrada
  Ts = 1/fms;                            % Obtención de tiempo de muestreo
  % Transformar parámetros a frecuencias angulares
  OmegaP = 2*pi*fp; OmegaS = 2*pi*fs;
  % Cálculo de frecuencias digitales
  wp = 2*atan(OmegaP*Ts/2); ws = 2*atan(OmegaS*Ts/2);
  % Parámetros y cálculo del prototipo pasa-bajos
  wplp = 0.2*pi;
  alpha = -(cos((wplp+wp)/2))/(cos((wplp-wp)/2));
  wslp = angle(-(exp(-1j*ws)+alpha)/(1+alpha*exp(-1j*ws)));
  OmegalP = (2/Ts)*tan(wplp/2);
  OmegalS = (2/Ts)*tan(wslp/2);
  [ns,ds] = afd_chb1(OmegalP,OmegalS,Rp,As);  % Cálculo de filtro analógico
  [nlpz,dlpz] = bilinear(ns,ds,fms);          % Transformación bilinear
  % Transformación a pasa-altos usando Z-Mapping
  Nz = -[alpha,1]; Dz = [1,alpha];
  [nhpz,dhpz] = zmapping(nlpz,dlpz,Nz,Dz);
  % Impresión de datos
  fprintf('\n* Datos de archivo %s utilizando método Bilineal',audio_in);
  fprintf('\n* Frecuencia de muestreo fs = %d [Hz]',fms);
  fprintf('\n* tiempo de muestreo Ts = %.2e [seg]',Ts);
  fprintf('\n* OmegaP = %.2f [rad/seg] y OmegaS = %.2f [rad/seg]',OmegaP,OmegaS);
  fprintf('\n* f_pass = %.2f [Hz] y f_stop = %.2f [Hz]',fp,fs);
  fprintf('\n* wp = %.3f [rad] y ws = %.3f [rad]',wp,ws);
  Hz = filt(nhpz,dhpz,Ts)                % Función de transferencia en 'z'
  Hs = d2c(Hz,'tustin')                  % Función de transferencia en 's'
    audio_out = filter(nhpz,dhpz,x);       % Aplicación de filtro a entrada
  audiowrite(n_audio_out,audio_out,fms)   % Escritura de archivo de salida
  % Extracción de parámetros para gráficos de diagramas de magnitud y fase
  [db,mag,fase,grd,w] = freqz_m(nhpz,dhpz);
  fmax = 15e3; wmax = 2*pi*fmax;
  [db2,mag2,fase2,w2] = freqs_m(cell2mat(Hs.num),cell2mat(Hs.den),wmax);
  % *************** Gráficas pedidas para Trabajo *************************
  % *************** Diagrama Polos y Ceros en s ***************************
  figure(); pzmap(Hs);                    % Gráfico de Polos y ceros en 's'
  title('Diagrama Polos y Ceros H(s)');
  saveas(gcf,'polos_zeros_Hs.png');
  % *************** Diagrama Polos y Ceros en z ***************************
  figure(); zplane(nhpz,dhpz);                % Gráfico de Polos y ceros en 'z'
  title('Diagrama Polos y Ceros H(z)'); grid on; 
  saveas(gcf,'polos_zeros_Hz.png');
  % *************** Bode de Filtro Discreto *******************************
  figure(); subplot(221); plot(w/pi,mag); grid on;
  title('Magnitud en abs'); xlabel('w/pi'); ylabel('|H(z)|');
  subplot(222); plot(w/pi,fase*180/pi); grid on;
  title('Fase en grados'); xlabel('w/pi'); ylabel('Fase H(z)');
  subplot(223); plot(w/pi,db); grid on;
  title('Magnitud en dB'); xlabel('w/pi'); ylabel('|H(z)| [dB]');
  subplot(224); plot(w/pi,grd); grid on;
  title('Retardo de grupo'); xlabel('w/pi'); ylabel('retardo de grupo');
  saveas(gcf,'Bode_digital.png');
  % *************** Bode de Filtro Continuo *******************************
  figure(); subplot(211); plot(w2/(2*pi),mag2); grid on;
  title('Magnitud en valor absoluto'); xlabel('Hz'); ylabel('|H(z)|');
  subplot(212); plot(w2/(2*pi),fase2*180/pi); grid on;
  title('Fase en grados'); xlabel('Hz'); ylabel('Fase H(z)');
  saveas(gcf,'Bode_Hz.png');
end