function [audio_out,Hs,Hz] = chebyshev1_lp(audio_in,n_audio_out,wp,ws,Rp,As,metodo)
% Diseño y aplicación de filtro Pasa-Bajos Chebyshev-1 a partir de un audio
% -----------------------------------------
% audio_out = chebyshev1_lp(audio_in,n_audio_out,wp,ws,Rp,As,metodo);
% audio_out = Señal de salida aplicada el filtro
% Hs = Función de transferencia en dominio s
% Hz = Función de transferencia en dominio z
% audio_in = Nombre del audio de entrada que se leerá
% n_audio_out = Nombre del archivo de audio de salida que se escribirá
% wp = frecuencia de corte de banda de paso digital []
% ws = frecuencia de banda de rechazo digital []
% Rp = Rizado de banda pasante en dB
% As = atenuación  en la banda de rechazo en dB
% metodo = metodo de transformacion: 'inv_impulso' 'bilinear'
  [x,fs] = audioread(audio_in);         % Leer archivo entrada
  fprintf('\n* Datos de archivo %s utilizando método %s',audio_in,metodo);
  Ts = 1/fs;                            % Obtención de tiempo de muestreo
  % Transformar parámetros a continuo dependiendo de método utilizado
  if strcmp(metodo,'inv_impulso')
    OmegaP = wp*fs;
    OmegaS = ws*fs;
  elseif strcmp(metodo,'bilinear')
    OmegaP = (2/Ts)*tan(wp/2);
    OmegaS = (2/Ts)*tan(ws/2);
  end
  fpass = OmegaP/(2*pi); fstop = OmegaS/(2*pi);
  fprintf('\n* Frecuencia de muestreo fs = %d [Hz]',fs);
  fprintf('\n* tiempo de muestreo Ts = %.2e [seg]',Ts);
  fprintf('\n* OmegaP = %.2f [rad/seg] y OmegaS = %.2f [rad/seg]',OmegaP,OmegaS);
  fprintf('\n* f_pass = %.2f [Hz] y f_stop = %.2f [Hz]',fpass,fstop);
  [ns,ds]=afd_chb1(OmegaP,OmegaS,Rp,As);  % Cálculo de filtro analógico
  % Transformar parámetros a continuo dependiendo de método utilizado
  if strcmp(metodo,'inv_impulso')
    [nz,dz] = imp_invr(ns,ds,Ts);        % Transformación invariancia imp.
    nz = nz'; dz = dz';
  elseif strcmp(metodo,'bilinear')
    [nz,dz] = bilinear(ns,ds,fs);        % Transformación bilinear
  end
  Hs = tf(ns,ds)                         % Función de transferencia en 's'
  Hz = filt(nz,dz,Ts)                    % Función de transferencia en 'z'
  audio_out = filter(nz,dz,x);           % Aplicación de filtro a entrada
  audiowrite(n_audio_out,audio_out,fs)   % Escritura de archivo de salida
  % Extracción de parámetros para gráficos de diagramas de magnitud y fase
  [db,mag,fase,grd,w] = freqz_m(nz,dz);
  fmax = 10e3; wmax = 2*pi*fmax;
  [db2,mag2,fase2,w2] = freqs_m(ns,ds,wmax);
  % *************** Gráficas pedidas para Trabajo *************************
  % *************** Diagrama Polos y Ceros en s ***************************
  figure(); pzmap(Hs);                    % Gráfico de Polos y ceros en 's'
  title('Diagrama Polos y Ceros H(s)');
  saveas(gcf,'polos_zeros_Hs.png');
  % *************** Diagrama Polos y Ceros en z ***************************
  figure(); zplane(nz,dz);                % Gráfico de Polos y ceros en 'z'
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