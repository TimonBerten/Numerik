% Bestimme Naeherungsloesung fuer zweidimensionale, stationaere Waermeleitungsgleichung mittels FD-Verfahren
% d^2 T(x,y) / dx^2 + d^2 T(x,y) / dy^2 = 0
% mit -a/2 <= x <= a/2, 0 <= y <= b

%% HAUPTFUNKTION
% =========================================================================
function Waermeleitung()

% Wechseln zwischen Randbedingungen:
% Dirichlet (aufgabe=1) oder Neumann/Robin (aufgabe=2)
aufgabe = 2;

% Loesung des resultierenden LGS:
% Backslash-Operator (jacobi=false) oder Jacobi-Verfahren (jacobi=true)
jacobi = true;

% Rechengebiet
a = 0.001;
b = 0.005;

% Dirichlet-Randbedingung
T_A = 303.15;
T_B = 333.15;
% Neumann-Randbedingung
q_B    = 25000;
lambda = 240;
% Robin-Randbedingung
alpha = 80;
T_inf = 293.15;

% Diskretisierung
M = 11;
N = 51;

% Rechengitter
dx = a/(M-1);
dy = b/(N-1);


% Systemmatrix und Vektor der rechten Seite
A = zeros(M*N, M*N);
r = zeros(M*N,1);

% -------------------------------------------------------------------------
for j=1:N
  for i=1:M
    % Differenzenstern: Index-Mapping (i,j) -> k
    ij   =  M*(j-1)+i;% (i,j)
    im1j =  ij-1;% (i-1,j)
    ip1j =  ij+1;% (i+1,j)
    ijm1 =  ij-M;% (i,j-1)
    ijp1 =  ij+M;% (i,j+1)
    
    if (1<i)&&(i<M) && (1<j)&&(j<N)  % innere Punkte
      % innere Punkte
      A(ij,ij)   = -2*(1/dx^2+1/dy^2);
      A(ij,im1j) = 1/dx^2;
      A(ij,ip1j) = 1/dx^2;
      A(ij,ijm1) = 1/dy^2;
      A(ij,ijp1) = 1/dy^2;
    
    else  % Randpunkte
      
      if aufgabe==1 % .....................................................
        if (j==1)       % unterer Rand: Dirichlet

            A(ij,ij) = 1;
            r(ij)    = T_B;
        elseif (j==N)   % oberer Rand: Dirichlet
            
            A(ij,ij) = 1;
            r(ij)    = T_A;
        elseif (i==1)   % linker Rand: Dirichlet
           
            A(ij,ij) = 1;
            r(ij)    = T_A;
        elseif (i==M)   % rechter Rand: Dirichlet
           
            A(ij,ij) = 1;
            r(ij)    = T_A;
        end
        
      elseif aufgabe==2 % .................................................
        if (j==1)       % unterer Rand: Neumann
             
            A(ij,ij) = -1/dy;
            A(ij,ijp1) = 1/dy;
            r(ij) = -q_B/lambda;

        elseif (j==N)   % oberer Rand: Robin  

            A(ij,ij) = 1/dy+alpha/lambda;
            A(ij,ijm1) = -1/dy;
            r(ij) = alpha/lambda *T_inf;

        elseif (i==1)   % linker Rand: Robin

           A(ij,ij) = -(1+alpha*dx/lambda)/dx;
           A(ij,ip1j) = 1/dx;
           r(ij) = -alpha/lambda *T_inf;

        elseif (i==M)   % rechter Rand: Robin
    
           A(ij,ij) = (1+alpha*dx/lambda)/dx;
           A(ij,im1j) = -1/dx;
           r(ij) = alpha/lambda *T_inf;

        end
        
      end % ...............................................................
    
    end % innere Punkte / Raender
  end % j
end % i
% -------------------------------------------------------------------------

% Loesung des resultierenden LGS
if jacobi
  T = Jacobi(A,r);
else
  T = A\r;
end

% Umkehrung des Index-Mapping k -> (i,j)
x     = linspace(-a/2,a/2,M);
y     = linspace(0,b,N);
[X,Y] = meshgrid(x,y);
T_mat = reshape(T,[M,N])';

% numerische Loesung plotten
figure(),hold on
surf(X,Y,T_mat)
xlabel('x [m]')
ylabel('y [m]')
zlabel('T [K]')

% Ausgabe gefragte Temperaturen

T_Kern=T(6)

end


%% UNTERFUNKTIONEN
% =========================================================================

% Jacobi-Verfahren zur Loesung von A x = r
% -------------------------------------------------------------------------
function [x_out] = Jacobi(A,rhs)
  % Ueberpruefung Input-Argumente
  if size(A,1)~=size(A,2)
    error('Matrix M muss symmetrisch sein')
  else
    MN  = size(A,1);
    xp  = zeros(size(rhs)); % alte Iterierte x^(p)
    xpp = xp;               % neue Iterierte x^(p+1)
  end
  
  % Parameter definieren
  eps_max   = 10^(-4);% Toleranz Residuum
  iter_max  = 500;% maximale Anzahl Iterationen
  iter      = 0;% Iterationsindex

  % Iterationsschleife bis Abbruchbedingung erfuellt
  while 1
    % Bestimmung von neuer Iterierten .....................................
    for l=1:MN
      sum1=0;
      for k=1:(l-1)
          sum1=sum1+A(l,k)*xp(k);
      end
      sum2=0;
      for k=(l+1):MN
          sum2=sum2+A(l,k)*xp(k);
      end

      xpp(l) = (rhs(l) - sum1 - sum2)/A(l,l);
     
    end % .................................................................
    % Residuum
    eps = max(abs(xpp-xp));
    % Inkrementierung
    xp = xpp;
    iter = iter+1
    % Abbruchkriterium
    %if (iter>=iter_max) || (eps<=eps_max)
    %  break;
    %end

    if (iter>=iter_max)
      break;
    end
  end
  
  % Ausgabe
  x_out = xpp;
  if iter==iter_max
    fprintf('Residuum = %.4E - nach %d Iterationen abgebrochen\n', eps, iter);
  else
    fprintf('Residuum < %1.0E in %d Iterationen erreicht\n', eps, iter);
  end
end


% exakte Loesung Aufgabe 1
% -------------------------------------------------------------------------
function [out] = T_exakt(x,y,a,b,T_A,T_B)
  if size(x)~=size(y)
    error('x und y muessen selbe Groesse haben!')
  else
    Theta = zeros(size(x));
    x_tilde = x / (0.5*a);
    y_tilde = y / (0.5*a);
    b_tilde = b / (0.5*a);
  end
  % Entdimensionalisierung
  Theta_B = (T_B-T_A)/T_A;
  N_inf   = 20;
  % Reihenloesung
  for n = 1:N_inf
    lambda_n = (2*n-1)/2*pi;
    C_n      = Theta_B/lambda_n*2*sin(lambda_n);
    Theta    = Theta + C_n * cos(lambda_n*x_tilde) .* (cosh(lambda_n*y_tilde) - 1/tanh(lambda_n*b_tilde)*sinh(lambda_n*y_tilde));
  end
  % Redimensionalisierung
  out = T_A*(Theta+1);
end