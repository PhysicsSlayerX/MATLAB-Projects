

hold(20)

function ho1d(N)
    
    clf     % clear figure
    f2f     % = figure(gcf)
    
    K = 2*N+1;  % number of total basis functions
    
    dx = 0.4;   % spacing between basis functions
    alpha=1;    % 1/spread of basis fns
    
    xvals = -4:0.01:4;        % the x-axis
    xmin = -N*dx;             % minimum x value, i.e., center of basis function
    
    S = zeros(K,K);           % the overlap matrix of basis functions
    H = zeros(K,K);           % the Hamiltonian matrix
    
    for A=1:K                          % loop over all basis functions
        xA = xmin + (A-1)*dx;          % the center of basis function A
        for B=1:K                      % loop over all other basis functions
            xB = xmin + (B-1)*dx;      % the center of basis function B
            
            % Equation for SAB is given in Szabo and Ostlund, p47
            % Equation for HAB is the typical bra-ket integral for the harm osc
            S(A,B) = sqrt(0.5*pi/alpha)*exp(-0.5*alpha*(xA-xB)^2);
            H(A,B) = 0.5*S(A,B)*(alpha - (alpha^2)*(xA-xB)^2 ...
                + 0.25*(1/alpha + (xA+xB)^2) );
            
        end
    end
    
    SinvH = inv(S) * H;           % S^(-1) H
    [c D] = eig(SinvH);           % Let MATLAB compute (S^(-1) H) c = E c
    [c D] = sorteig(c,D)          
    
    diag(D)
    
    % Normalize wavefunctions
    % The Gaussians we began with aren't normalized, so the wavefunctions
    % shouldn't be either. Wavefunctions must be normalized to be used  
    % to create probability distributions (Pauling and Wilson, p104).
    for i=1:K
       norm=0;
       
       for A=1:K
           for B=1:K
               norm = norm+c(A,i)*c(B,i)*S(A,B);
           end
       end
       norm = sqrt(norm);
       c(:,i) = c(:,i)/norm;
        
    end
    
    
    %plot wavefunctions
    szx = size(xvals);    % Size (shape, really) of x values array.
    nel = szx(2);         % The actual number of x elements
    psi0 = zeros(1,nel);  % Initialize the wavefunction array
    c0 = c(:,1);          % wavefunction, ground state
    
    psi1 = zeros(1,nel);
    c1 = c(:,2);          % wavefunction, 1st excited state
    
    psi2 = zeros(1,nel);
    c2 = c(:,3);          % wavefunction, 2nd excited state
    
    psi9 = zeros(1,nel);
    c9 = c(:,10);         % wavefunction, 9th excited state
    
    
    for A=1:K
        xA = xmin + (A-1)*dx;
        psi0 = psi0 + c0(A)*exp(-alpha*(xvals-xA) .^ 2);
        psi1 = psi1 + c1(A)*exp(-alpha*(xvals-xA) .^ 2);
        psi2 = psi2 + c2(A)*exp(-alpha*(xvals-xA) .^ 2);
        psi9 = psi9 + c9(A)*exp(-alpha*(xvals-xA) .^ 2);
    end
    
    hold on
    
    subplot(2,2,1)
    hold on
    plot(xvals,-psi0)
    % check0 = correct ground state
    check0 = (pi^(-0.25))*exp(-0.5* xvals .^ 2);
    plot(xvals,check0,'r');
    
    subplot(2,2,2)
    hold on
    plot(xvals,-psi1)
    hold off % checkl = correct 1st excited state 
    check1 = (pi^(-0.25))*exp(-0.5* xvals .^ 2)*sqrt(2) .* xvals;
    plot(xvals,check1,'r')
    
    subplot(2,2,3)
    plot(xvals,psi2)
    
    hold on
    % check2 = correct 2nd excited state
    check2 = (pi^(-0.25))*exp(-0.5* xvals .^ 2) .* (2*xvals .^2 - 1)/sqrt(2);
    plot(xvals,check2,'r')
    
    
    subplot(2,2,4)
    hold on
    plot(xvals,psi9)
    % check9 = correct 9th excited state. Note, usage of Hermite polynomials.
    check9 = (pi^(-0.25))*exp(-0.5* xvals .^ 2) .* ...
        (512*xvals .^9 -9216*xvals .^7 +48384*xvals .^5 ...
        -80640*xvals .^3 + 30240*xvals)/(2304*sqrt(35));
    plot(xvals,check9,'r')
    
    
    energies = diag(D);
end