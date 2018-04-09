function [audio_out,Hs,HZbp] = chebyshev_bp(audio_in,n_audio_out,fs1,fp1,fp2,fs2,Rp,As)
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
  OmegaP1 = 2*pi*fp1; OmegaS1 = 2*pi*fs1;
  OmegaP2 = 2*pi*fp2; OmegaS2 = 2*pi*fs2;
  % Cálculo de frecuencias digitales
  wp1 = 2*atan(OmegaP1*Ts/2); ws1 = 2*atan(OmegaS1*Ts/2);
  wp2 = 2*atan(OmegaP2*Ts/2); ws2 = 2*atan(OmegaS2*Ts/2);
  % **************** 1) Cálculo del filtro pasa-bajos *********************
  [ns1,ds1] = afd_chb1(OmegaP2,OmegaS2,Rp,As);% Cálculo de filtro analógico
  [nlpz1,dlpz1] = bilinear(ns1,ds1,fms);      % Transformación bilinear
  % **************** 2) Cálculo de Filtro pasa-altos **********************
  wplp = 0.2*pi;
  alpha = -(cos((wplp+wp1)/2))/(cos((wplp-wp1)/2));
  wslp = angle(-(exp(-1j*ws1)+alpha)/(1+alpha*exp(-1j*ws1)));
  OmegalP = (2/Ts)*tan(wplp/2);
  OmegalS = (2/Ts)*tan(wslp/2);
  [ns2,ds2] = afd_chb1(OmegalP,OmegalS,Rp,As);% Cálculo de filtro analógico
  [nlpz,dlpz] = bilinear(ns2,ds2,fms);        % Transformación bilinear
  % Transformación a pasa-altos usando Z-Mapping
  Nz = -[alpha,1]; Dz = [1,alpha];
  [nhpz2,dhpz2] = zmapping(nlpz,dlpz,Nz,Dz);
  % *************** 3) Cálculo de Filtro pasa-banda ***********************
  Hzlp = filt(nlpz1,dlpz1,Ts); % F. transf Filtro pasa-bajo
  Hzhp = filt(nhpz2,dhpz2,Ts); % F. transf Filtro pasa-bajo
  HZbp = Hzlp*Hzhp;            % F. Transf. Filtro pasa-banda
  nbpz = cell2mat(HZbp.num);
  dbpz = cell2mat(HZbp.den);
  % Impresión de datos
  fprintf('\n* Datos de archivo %s utilizando método Bilineal',audio_in);
  fprintf('\n* Frecuencia de muestreo fs = %d [Hz]',fms);
  fprintf('\n* tiempo de muestreo Ts = %.2e [seg]',Ts);
  fprintf('\n* OmegaP1 = %.2f [rad/seg] y OmegaS1 = %.2f [rad/seg]',OmegaP1,OmegaS1);
  fprintf('\n* OmegaP2 = %.2f [rad/seg] y OmegaS2 = %.2f [rad/seg]',OmegaP2,OmegaS2);
  fprintf('\n* f_pass1 = %.2f [Hz] y f_stop1 = %.2f [Hz]',fp1,fs1);
  fprintf('\n* f_pass2 = %.2f [Hz] y f_stop2 = %.2f [Hz]',fp2,fs2);
  fprintf('\n* wp1 = %.3f [rad] y ws1 = %.3f [rad]',wp1,ws1);
  fprintf('\n* wp2 = %.3f [rad] y ws2 = %.3f [rad]',wp2,ws2);
  HZbp
  Hs = d2c(HZbp,'tustin')                 % Función de transferencia en 's'
  audio_out = filter(nbpz,dbpz,x);        % Aplicación de filtro a entrada
  audiowrite(n_audio_out,audio_out,fms)   % Escritura de archivo de salida
  % Extracción de parámetros para gráficos de diagramas de magnitud y fase
  [db,mag,fase,grd,w] = freqz_m(nbpz,dbpz);
  fmax = 15e3; wmax = 2*pi*fmax;
  [db2,mag2,fase2,w2] = freqs_m(cell2mat(Hs.num),cell2mat(Hs.den),wmax);
  % *************** Gráficas pedidas para Trabajo *************************
  % *************** Diagrama Polos y Ceros en s ***************************
  figure(); pzmap(Hs);                    % Gráfico de Polos y ceros en 's'
  title('Diagrama Polos y Ceros H(s)');
  saveas(gcf,'polos_zeros_Hs.png');
  % *************** Diagrama Polos y Ceros en z ***************************
  figure(); zplane(nbpz,dbpz);            % Gráfico de Polos y ceros en 'z'
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