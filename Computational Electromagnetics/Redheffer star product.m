function [ SAB11,SAB12,SAB21,SAB22 ] = redhefferstar( SA11,SA12,SA21,SA22,SB11,SB12,SB21,SB22)

% redhefferstar product combines two scattering matrices to form a overall
% scattering matrix. It is used in forming scattering matrices of
% dielectric stacks in transfer matrix method
% SA and SB are  scattering matrices of two different layers
% and this function outputs
% SAB which is the combined scaterring matrix of two layers

N=length(SA11);
I=eye(N);
SAB11=SA11+(SA12*((I-(SB11*SA22))^-1)*SB11*SA21);
SAB12=SA12*((I-(SB11*SA22))^-1)*SB12;
SAB21=SB21*((I-(SA22*SB11))^-1)*SA21;
SAB22=SB22+(SB21*((I-(SA22*SB11))^-1)*SA22*SB12);

end

