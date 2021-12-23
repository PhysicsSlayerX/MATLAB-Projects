function C = convmat(A, P, Q, R)

% CONVMAT    Build Convolution Matrix
% EMPossible
%
% C = convmat(A, P);        for 1D Problems
% C = convmat(A, P, Q);     for 2D Problems
% C = convmat(A, P, Q, R);  for 3D Problems
%
% This MATLAB function construct convolution matrices from
% a real-space grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DETERMINE SIZE OF A
[Nx, Ny, Nz] = size(A);

% HANDLE NUMBER OF HARMONICS FOR ALL DIMENSIONS
if nargin == 2
    Q = 1;
    R = 1;
elseif nargin == 3
    R = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD CONVOLUTION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE INDICES OF SPATIAL HARMONICS
NH = P*Q*R;
p   = [-floor(P/2): + floor(P/2)];
q   = [-floor(Q/2): + floor(Q/2)];
r   = [-floor(R/2): + floor(R/2)];

% COMPUTE FOURIER COEFFICIENTS OF A 
A = fftshift(fftn(A))/(Nx*Ny*Nz);

% COMPUTE ARRAY INDICES OF ZERO-ORDER HARMONIC
p0 = 1 + floor(Nx/2);
q0 = 1 + floor(Ny/2);
r0 = 1 + floor(Nz/2);

% INITIALIZE CONVOLUTION MATRIX
C = zeros(NH, NH);

% CONSTRUCT CONVOLUTION MATRIX
for rrow = 1 : R
    for qrow = 1 : Q
        for prow = 1 : P
            row = (rrow - 1)*Q*P + (qrow - 1)*P  + prow;
            for rcol = 1 : R
                for qcol = 1 : Q
                    for pcol = 1 : P
                        col = (rcol - 1)*Q*P + (qcol - 1)*P  + pcol;
                        
                        pfft = p(prow) - p(pcol);
                        qfft = q(qrow) - q(qcol);
                        rfft = r(rrow) - r(rcol);
                        
                        C(row, col) = A(p0+pfft, q0+qfft, r0+rfft);

                    end
                end
            end
            
        end
    end
end



% imagesc(real(A));
% colorbar;
















