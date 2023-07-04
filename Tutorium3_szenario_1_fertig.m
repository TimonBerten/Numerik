%% HAUPTFUNKTION
% Finite-Volumen-Code zur Loesung der 1D Flachwassergleichungen
% Q_t + F_x(Q) = 0
% mit Q = [h, h*u] und F(Q) = [h*u, h*u^2 + 1/2*g*h^2]

%% HAUPTFUNKTION
% =========================================================================
function Flachwasser_FV()

close all

% Anfangswert-Szenario:
% 1 = Wellenproblem (Gauß'sche Glockenkurve)
% 2 = Dammbruchproblem (Stufenprofil)
szenario = 1;

% Modellparameter
g   = 9.81;% Erdbeschleunigung
h_0 = 1;% Wasserhoehe Ruhelage
u_0 = 0;% Horizontalgeschwindigkeit Ruhelage

% Raum-Diskretisierung
x_links  = 0;% Rechengebiet linker Rand
x_rechts = 1;% Rechengebiet rechter Rand
N        = 100;% Anzahl Gitterzellen

% Zeit-Diskretisierung
CFL   = 0.5;% CFL-Zahl
t_0   = 0;% Startzeit
t_end = 0.15;% Endzeit

% Variablen deklarieren
Q   = zeros(2,N); % Loesungsvektor (wird in jedem Zeitschritt ueberschrieben):
% erste Zeile = h, zweite Zeile = h*u
dx  = (x_rechts-x_links)/N; % Raum-Schrittweite
x   = ((dx/2):(dx):(x_rechts-dx/2))% Koordinaten Zellmittelpunkte

% Anfangswerte
if szenario == 1  % Wellenproblem
  % Amplitude der Anfangsstoerung
  dh = 0.1;
  % Gauss'sche Glockenkurve
  for i = 1:N
    Q(1,i) = h_0 + dh * exp(-(50*(x(i)-0.5))^2);
    Q(2,i) = Q(1,i)*u_0;
  end
elseif szenario == 2  % Dammbruchproblem
  % Amplitude der Anfangsstoerung
  % ...
else  % Szenario nicht implementiert
  error('Szenario nicht implementiert - waehle 1 (Wellenproblem) oder 2 (Dammbruchproblem)!')
end
  
% Initialisierung Laufvariablen
t = t_0;
spotted = false;  % Flag, ob bereits Welle detektiert wurde

% Visualisierung Anfangswerte
figure(1)
Qprim = cons2prim(Q);
visualize(Qprim,x,t,szenario,1)

% Zeitschleife ----------------------------------------------------------
while t<t_end

  % Berechnung Zeitschritt
  lambda_max = get_lambdamax(Q,g,N);
  dt = CFL*dx/lambda_max;

  % Anpassung Zeitschritt, falls Endzeitpunkt ueberschritten
  if t+dt > t_end
    dt = t_end - t;
  end 

  % Loesung zum naechsten Zeitschritt (ueberschreibe Loesungsvariable)
  Q = Euler_explizit(Q,dt,dx,lambda_max,N,g,szenario);
  
  % Zeitinkrement
  t = t + dt;

  % Anwendung Wave Spotter
  Qprim = cons2prim(Q);
  spotted = wave_spotter(Qprim,t,x,szenario,spotted);

  % Live-Visualisierung
  visualize(Qprim,x,t,szenario,1)
  drawnow   % aktualisiere Plot zur Laufzeit
  pause(0.01)

end % -------------------------------------------------------------------

end


%% UNTERFUNKTIONEN
% =========================================================================

% Zeitintegration 1.Ordnung mittels explizitem Euler-Verfahren ------------
function [Q_neu] = Euler_explizit(Q,dt,dx,lambda_max,N,g,szenario)
  % Raumoperator auswerten
  % ...
  R = FV(Q,dx,lambda_max,N,g,szenario);
  

  % explizites Euler-Verfahren anwenden
  % ...
  Q_neu = Q+dt*R;
end % ---------------------------------------------------------------------


% FV-Verfahren ------------------------------------------------------------
function [dQ] = FV(Q,dx,lambda_max,N,g,szenario)
  
  dQ = zeros(2,N);    % Loesungsinkrement
  
  % Ghost States an linkem und rechtem Rand des Rechengebiets
  if szenario == 1       % Gauss-Puls: Ausflussraender
    Qghost_L = Q(:,1);
    Qghost_R = Q(:,N);
  elseif szenario == 2   % Dammbruch: reflektierende Waende
    Qghost_L = 0;
    Qghost_R = 0;
  end
  
  % FV-Raumoperator zellweise zusammensetzen
  for i=1:N
    % numerische Fluesse am linken (G_L = G_{i-1/2}) und rechten (G_R = G_{i+1/2}) Zellrand
    if i==1       % linker Rand des Rechengebiets
      % ...
      G_L = LaxFriedrichs(Qghost_L,Q(:,i),lambda_max,g);
      G_R = LaxFriedrichs(Q(:,i),Q(:,i+1),lambda_max,g)
    elseif i==N   % rechter Rand des Rechengebiets
      % ...
      G_R = LaxFriedrichs(Q(:,N),Qghost_R,lambda_max,g);
      G_L = LaxFriedrichs(Q(:,i-1),Q(:,i),lambda_max,g);
    else          % innere Zellen
      G_L = LaxFriedrichs(Q(:,i-1),Q(:,i),lambda_max,g);
      G_R = LaxFriedrichs(Q(:,i),Q(:,i+1),lambda_max,g)
    end
    % FV-Raumoperator fuer aktuelle Zelle
    dQ(:,i) = - 1/dx * (G_R - G_L);
  end
  
end % ---------------------------------------------------------------------


% Numerische Flussfunktion: Lax-Friedrichs --------------------------------
function [G] = LaxFriedrichs(Q_L,Q_R,lambda_max,g)
  
  % Auswertung physikalischer Fluss
  F_L = flux(Q_L,g);
  F_R = flux(Q_R,g);
  
  % globaler Lax-Friedrichs-Fluss
  G = 0.5 * (F_L + F_R) - 0.5 * lambda_max * (Q_R - Q_L);
end % ---------------------------------------------------------------------


% Physikalische Flussfunktion der Flachwasser-Gleichungen -----------------
function [F] = flux(Q,g) 
  F = zeros(size(Q));
  
  % h*u
  F(1,:) = Q(2,:);
  % h*u^2 + 1/2*g*h^2
  F(2,:) = (Q(2,:)^2/Q(1,:)+0.5*g*Q(1,:)^2);
end % ---------------------------------------------------------------------


% Umrechnung konservative in primitive Variablen --------------------------
function [Qprim] = cons2prim(Q)
  Qprim = zeros(size(Q));
  
  % h
  Qprim(1,:) = Q(1,:);
  % u = h*u / h
  Qprim(2,:) = Q(2,:)./Q(1,:);
end % ---------------------------------------------------------------------


% Berechnung maximaler Eigenwert ------------------------------------------
function [lambda_max] = get_lambdamax(Q,g,N)
  a = zeros(1,N);
  Qprim = cons2prim(Q);
  
  for i=1:N

    u = Qprim(2,i);
    h = Qprim(1,i);
    lambda1 = u - sqrt(g*h);
    lambda2 = u + sqrt(g*h);
    
    a(i) = max(abs(lambda1),abs(lambda2))
  end
  
  % maximaler Eigenwert global
  lambda_max = max(a);
end % ---------------------------------------------------------------------



%% HILFSFUNKTIONEN
% =========================================================================

% Wave Spotter ------------------------------------------------------------
function [spotted_out] = wave_spotter(Qprim,t,x,szenario,spotted_in)
  
  % Definition Spotter
  x_spotter   = 0.8;  % Spotter-Position
  h_threshold = 1.01; % Schwellwert Hoehe
  
  % naehester Knotenpunkt zu Spotter-Position
  [~,i_spotter] = min(abs(x_spotter-x));
  
  % detektiere, wenn Hoehe h(x_spotter,t) Schwellwert zum ersten Mal ueberschreitet
  if not(spotted_in) && (Qprim(1,i_spotter) > h_threshold)
    
    % Ausgabe Zeitpunkt
    spotted_out = true;
    fprintf('Welle detektiert zum Zeitpunkt t = %8.6f\n',t)
    
    % Visualisierung in separater Figure
    visualize(Qprim,x,t,szenario,2)
    subplot(2,1,1),hold on
    plot(x_spotter*ones(1,2),ylim,'k--')   % Spotter-Position
    legend('Loesung','Wave Spotter','location','northwest')

  else
    
    % Welle bereits detektiert oder Schwellwert nicht erreicht
    spotted_out = spotted_in;
  end
end % ---------------------------------------------------------------------


% Visualisierung ----------------------------------------------------------
function visualize(Qprim,x,t,szenario,figure_index)
  figure(figure_index)

  % Hoehe
  subplot(2,1,1)
  plot(x,Qprim(1,:),'-ob')
  xlabel('x')
  ylabel('Hoehe h')
  title(sprintf('t = %5.3f',t))
  xlim([x(1),x(end)])
  if szenario == 1  % Wellenproblem
    ylim([0.9,1.1])
  else              % Dammbruchproblem
    ylim([0.9,1.4])
  end

  % Geschwindigkeit
  subplot(2,1,2)
  plot(x,Qprim(2,:),'-or')
  xlabel('x')
  ylabel('Geschwindigkeit u')
  xlim([x(1),x(end)])
  if szenario == 1  % Wellenproblem
    ylim([-0.1,0.1])
  else              % Dammbruchproblem
    ylim([-0.1,0.5])
  end

end % ---------------------------------------------------------------------