function S = star(SA, SB)
% STAR   Redheffer Star Product
%
% S = star(SA, SB)
% 
% INPUT ARGUMENTS
% ==========================================
% 
% SA            First Scattering Matrix
%   .S11        S11 Scattering Parameter
%   .S12        S12 Scattering Parameter
%   .S21        S21 Scattering Parameter
%   .S21        S22 Scattering Parameter
%
% SB            Second Scattering Matrix
%   .S11        S11 Scattering Parameter
%   .S12        S12 Scattering Parameter
%   .S21        S21 Scattering Parameter
%   .S21        S22 Scattering Parameter
%
% OUTPUT ARGUMENTS
% ==========================================
% 
% S             Combined Scattering Matrix
%   .S11        S11 Scattering Parameter
%   .S12        S12 Scattering Parameter
%   .S21        S21 Scattering Parameter
%   .S21        S22 Scattering Parameter
%

% CONSTRUCT IDENTITY MATRIX
[M,N]  = size(SA.S11);
I      = eye(M,N);

% COMPUTE COMMON TERMS
D = SA.S12/(I - SB.S11*SA.S22);
F = SB.S21/(I - SA.S22*SB.S11);

% COMPUTE COMBINED SCATTERING MATRIX
S.S11 = SA.S11 + D*SB.S11*SA.S21;
S.S12 = D*SB.S12;
S.S21 = F*SA.S21;
S.S22 = SB.S22 + F*SA.S22*SB.S12;








