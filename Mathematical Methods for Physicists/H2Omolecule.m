  K=7; % number of contracted Gaussians, 1 per orbital type (e.g., 1s, 2s)
    L=2; % number of primitive Gaussians per contracted Gaussian
    
    Nelec=10; % number of electrons
    Nnuc=3;   % number of nuclei
    
    nucchg = [1 1 8]; % nuclear charge
    
    spreads = zeros(L,K);  % spreads for each primitive Gaussian
    d = zeros(L,K);        % "contraction coefficient", p155
    centers = zeros(L,K,3);% centers of each primitive Gaussian
    
    
    % zeta represents the size of the electron distribution on an atom
    zetaH = 1.24;
    zetaH2 = zetaH^2;
    zetaO = 7.66;
    zetaO2 = zetaO^2;
    
    
    % alphas were calculated by linear fitting to Gaussian functions to
    % Slater-type functions, p157 eq. 3.221
    
    %    1s orbital on hydrogen #1:
    d(1,1) = 0.164964;
    d(2,1) = 0.381381;
    spreads(1,1) = 0.151623*zetaH2;
    spreads(2,1) = 0.851819*zetaH2;
    
    %    1s orbital on hydrogen #2:
    d(1,2) = 0.164964;
    d(2,2) = 0.381381;
    spreads(1,2) = 0.151623*zetaH2;
    spreads(2,2) = 0.851819*zetaH2;
    
    %    1s orbital on oxygen:
    d(1,3) = 0.164964;
    d(2,3) = 0.381381;
    spreads(1,3) = 0.151623*zetaO2;
    spreads(2,3) = 0.851819*zetaO2;
    
    %    2s orbital on oxygen:
    d(1,4) = 0.168105;
    d(2,4) = 0.0241442;
    spreads(1,4) = 0.493363;
    spreads(2,4) = 1.945230;
    
    %    2px orbital on oxygen:
    d(1,5) =  1;
    d(2,5) = -1;
    spreads(1,5)=0.9;
    spreads(2,5)=0.9;
    
    %    2py orbital on oxygen:
    d(1,6) =  1;
    d(2,6) = -1;
    spreads(1,6)=0.9;
    spreads(2,6)=0.9;
    
    %    2pz orbital on oxygen:
    d(1,7) =  1;
    d(2,7) = -1;
    spreads(1,7)=0.9;
    spreads(2,7)=0.9;
    
    
    % Save the centers of the nuclei:
    dist = 1.809 % the current distance being computed
    theta = 104.52;
    Rnuc = zeros(Nnuc,3);
    Rnuc(1,1) =  dist*cosd(90-theta/2);
    Rnuc(1,2) =  dist*sind(90-theta/2);
    Rnuc(2,1) = -dist*cosd(90-theta/2);
    Rnuc(2,2) =  dist*sind(90-theta/2);
    
    % Center the 1s basis element for hydrogen #1:
    centers(1,1,1) = Rnuc(1,1);
    centers(2,1,1) = Rnuc(1,1);
    centers(1,1,2) = Rnuc(1,2);
    centers(2,1,2) = Rnuc(1,2);
    
    % Center the 1s basis element for hydrogen #2:
    centers(1,2,1) = Rnuc(2,1);
    centers(2,2,1) = Rnuc(2,1);
    centers(1,2,2) = Rnuc(2,2);
    centers(2,2,2) = Rnuc(2,2);
    
    % Note: 1s, 2s orbitals for oxygen remain at the origin.
    
    offset=0.1;
    % 2px, offset the 2 Gaussian functions:
    centers(1,5,1) =  offset;
    centers(2,5,1) = -offset;
    
    % 2py, offset the 2 Gaussian functions:
    centers(1,6,2) =  offset;
    centers(2,6,2) = -offset;
    
    % 2pz, offset the 2 Gaussian functions:
    centers(1,7,3) =  offset;
    centers(2,7,3) = -offset;
    
    % one-electron matrix elements (see eq 3.153, p141)
    S = zeros(K,K);
    T = zeros(K,K);
    V = zeros(K,K);
    
    
    %Calculate the parts of the Fock matrix hamiltonian:
        % overlap:       3.228
        % kinetic:       3.151,     A.11
        % potential:     3.153,     A.33, A.32
    for mu=1:K                   % loop over each basis element
        for p=1:L                % loop over each of the 3 primitives per basis element
            for nu = 1:K         % loop over each basis element
                for q=1:L        % loop over each of the 3 primitives per basis element
                    
                    RA = [centers(p,mu,1) centers(p,mu,2) centers(p,mu,3)];
                    RB = [centers(q,nu,1) centers(q,nu,2) centers(q,nu,3)];
                    alpha = spreads(p,mu);
                    beta = spreads(q,nu);
                    
                    S(mu,nu) = S(mu,nu) + d(p,mu)*d(q,nu)*...
                        overlap(RA,RB,alpha,beta);
                    
                    T(mu,nu) = T(mu,nu) + d(p,mu)*d(q,nu)*...
                        kinetic(RA,RB,alpha,beta);
                    
                    for i=1:Nnuc
                        RC = Rnuc(i,:);
                        V(mu,nu) = V(mu,nu) + d(p,mu)*d(q,nu)*...
                            nucchg(i)*elec_nuc(RA,RB,alpha,beta,RC);
                    end
                    
                end
            end
        end
    end
    
    % Two-electron matrix elements
    % Calculate eq 3.211 w/ primitive Gaussians 3.212:
    two_elec = zeros(K,K,K,K);
    
    for mu=1:K
        for nu=1:K
            for lambda=1:K
                for sigma=1:K
                    for p=1:L
                        for q=1:L
                            for s=1:L
                                for t=1:L
                                    
       RA = [centers(p,mu,1) centers(p,mu,2) centers(p,mu,3)];
       RB = [centers(q,nu,1) centers(q,nu,2) centers(q,nu,3)];
       RC = [centers(s,lambda,1) centers(s,lambda,2) centers(s,lambda,3)];
       RD = [centers(t,sigma,1) centers(t,sigma,2) centers(t,sigma,3)];
       
       alpha = spreads(p,mu);
       beta = spreads(q,nu);
       gamma = spreads(s,lambda);
       delta = spreads(t,sigma);
       
       two_elec(mu,nu,lambda,sigma) = two_elec(mu,nu,lambda,sigma) + ...
           d(p,mu)*d(q,nu)*d(s,lambda)*d(t,sigma)* ...
           elec_elec(RA,RB,RC,RD,alpha,beta,gamma,delta);
       
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    % Solve one-electron problem to get H^{core}(mu,nu)  (eq 3.149, 3.150)
    h = T+V;
    Sinv = inv(S);
    [c D] = eig(Sinv*h);
    [c D] = sorteig(c,D);
    
    size(c)
    
    nsteps=10;
    for step=1:nsteps
        
        % Normalize via "canonical orthonormalization", p144
        for i=1:K
            norm = 0;
            for mu=1:K
                for nu=1:K
                    norm = norm + c(mu,i)*c(nu,i)*S(mu,nu);
                end
            end
            c(:,i) = c(:,i)/sqrt(norm);
        end
        
        % Calculate the density matrix P(mu,nu), (eq 3.145)
        P = zeros(K,K);
        for i=1:Nelec/2
            for mu=1:K
                for nu=1:K
                    P(mu,nu) = P(mu,nu) + 2*c(mu,i)*c(nu,i);
                end
            end
        end
    
        % Calculate the G(mu,nu) part of the Fock matrix, (eq 3.154)
        F = h;          % begin by adding the 1-electron part
        for mu=1:K
            for nu=1:K
                for lambda=1:K
                    for sigma=1:K
                        F(mu,nu) = F(mu,nu) + ...
                            P(lambda,sigma)*( two_elec(mu,nu,lambda,sigma)...
                            - 0.5*two_elec(mu,lambda,nu,sigma) );
                        
                    end
                end
            end
        end
          
        energy = 0;
        % Calculate the total electronic energy:
        for mu=1:K
            for nu=1:K
                energy = energy + 0.5*P(mu,nu)*...
                    (h(mu,nu) + F(mu,nu));
            end
        end
        
        % Add the nuclear-nuclear repuslion term (eq 3.185) to get total energy: 
        % (this needs to get changed when doing the water part)
        internucEnergy = 8/dist + 8/dist + 1/(2*Rnuc(1,1));
        energy = energy + internucEnergy % compare with Table 3.13, 192
        
        [c D] = eig(Sinv*F);
        [c D] = sorteig(c,D);    
    end
    
    % Note: Eigenvalues here are approximations to ionization energies of
    %       the system.
    
    
    %plot wavefunction
    xvals = -1:0.1:3;
    yvals = -2:0.1:2;
    
    nel = numel(xvals);
    
    psi = zeros(nel,nel);
    for i=1:nel
        x = xvals(i);
        for j=1:nel
            y = yvals(j);
            for mu=1:K
                for p=1:L
                    xA = centers(p,mu,1);
                    yA = centers(p,mu,2);
                                   
                    alpha = spreads(p,mu);
                    
        psi(i,j) = psi(i,j) + d(p,mu)*c(mu,1)* ...
            exp(-alpha*(x-xA)^2-alpha*(y-yA)^2);     
                    
                   
                    
                end
            end
            
        end
    end
    
    contour(xvals,yvals,psi,15)
    axis equal
    drawnow
    
function TAB = kinetic(RA,RB,alpha,beta)
    
    apb = alpha+beta;
    
    abfac = alpha*beta/apb;
    
    dRAB = RA-RB;
    dRAB2 = dRAB*dRAB';
    
    TAB = abfac*(3 - 2*abfac*dRAB2)* ...
        ((pi/apb)^1.5)*exp(-abfac*dRAB2);
end

function SAB = overlap(RA,RB,alpha,beta)
    % RA and RB are row vectors
    
    apb = alpha+beta;
    
    abfac = alpha*beta/apb;
    
    dRAB = RA-RB;
    dRAB2 = dRAB*dRAB';
    
    SAB = ((pi/apb)^1.5)*exp(-abfac*dRAB2);
end

function VAB = elec_nuc(RA,RB,alpha,beta,RC)
    
    apb = alpha+beta;
    
    abfac = alpha*beta/apb;
    
    dRAB = RA-RB;
    dRAB2 = dRAB*dRAB';
    
    RP = (alpha*RA + beta*RB)/apb;
    dRPC = RP-RC;
    dRPC2 = dRPC*dRPC';
    
    
    VAB = -(2*pi/apb)*exp(-abfac*dRAB2)* ...
        F0(apb*dRPC2);
end

function F0 = F0(x)
    
    if (x < 1e-5)
        F0 = 1-x/3;
    else
        F0 = (sqrt(pi)/2)*erf(sqrt(x))/sqrt(x);
    end
end

function elec_elec = elec_elec(RA,RB,RC,RD,alpha,beta,gamma,delta)
    
    apb = alpha+beta;
    gpd = gamma + delta;
    
    abfac = alpha*beta/apb;
    gdfac = gamma*delta/gpd;
    
    dRAB = RA-RB;
    dRAB2 = dRAB*dRAB';
    
    dRCD = RC-RD;
    dRCD2 = dRCD*dRCD';
    
    RP = (alpha*RA + beta*RB)/apb;
    RQ = (gamma*RC + delta*RD)/gpd;
    dRPQ = RP-RQ;
    dRPQ2 = dRPQ*dRPQ';
    
    elec_elec = ((2*pi^2.5)/(apb*gpd*sqrt(apb+gpd))) * ...
        exp(-abfac*dRAB2-gdfac*dRCD2) * ...
        F0((apb*gpd/(apb+gpd))*dRPQ2);
end

function [Usort Dsort] = sorteig(U,D)

sz = size(U);
nel = sz(2);

d = diag(D);

[dsort index] = sort(d);

Dsort = zeros(nel,nel);

for i=1:nel
    
    Usort(:,i) = U(:,index(i));
    Dsort(i,i) = dsort(i);
    
end
end