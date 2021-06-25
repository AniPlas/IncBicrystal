%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code computes the stress fields in an infinite bi-crystal with a planar interface. 
% The bi-crystal is submitted to a macroscopic homogeneous and remotely applied stress as well as to piecewise uniform plastic strains. 
% The bi-crystal can have different grain volume fractions.
% The two grains are assumed perfectly bonded  along the planar interface whose normal is along e2. 
% The code provides also the effective stiffness tensor, the macroscopic
% strain tensor and the effective plastic strain tensor.
% References for the model can be found in:
% T. Richeton, S. Berbenni, Eur. J. Mech. A. Solids 37 (2013) 231–247
% T. Richeton, S. Berbenni, Int. J. Sol. Struct. 51 (2014) 794-807
% T. Richeton, I. Tiba, S. Berbenni, O. Bouaziz, Phil. Mag. 95 (2015) 12-31
% I. Tiba, T. Richeton, C. Motz, H. Vehoff, S. Berbenni, Acta Mater. 83 (2015) 227-238
% T. Richeton, Crystals 7 (2017) 203 (1-14)
% T. Richeton, Scripta Mater. 169 (2019) 14-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear;

x = [1;0;0];
y = [0;1;0];
z = [0;0;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volume fractions (to be changed)
f1 = 0.5 % Volume fraction of material I 
f2 = 1-f1 % Volume fraction of material II 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macroscopic stress tensor in the global (EBSD) frame (to be changed)
SIGMA=zeros(3);
SIGMA(1,1)=500;
SIGMA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slip on the highest resolved systems in both materials (to be changed)
gammaI=0.04
gammaII=0.04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler angles in degree (to be changed)
% Material I
phi1_1 = 61;
phi_1 = 46.26;
phi2_1 = 64.71;

% Material II
phi1_2 = 179.09;
phi_2 = 18.71;
phi2_2 = 2.24;

% Degree --> Radian
phi1_1 = phi1_1*pi/180;
phi_1 = phi_1*pi/180;
phi2_1 = phi2_1*pi/180;
phi1_2 = phi1_2*pi/180;
phi_2 = phi_2*pi/180;
phi2_2 = phi2_2*pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orientation of the interface in the EBSD frame (to be changed)
psi = 161; % Angle between the trace of the interface at the surface and the X direction of the EBSD
psi = psi*pi/180;
epsilon = 0*pi/180; % Inclination of the interface (0° --> interface perpendicular to the surface)
epsilon = epsilon*pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic stiffnesses tensor of Material I in MPa (to be changed)
% Iron (FCC austenite)
C11=197500;
C12=125000;
C44=122000;

DLOCAL1=zeros(6);
for j=1:3
    DLOCAL1(j,j)=C11;
    for i=1:3
        if (i~=j)
            DLOCAL1(i,j)=C12;
        end
    end
    DLOCAL1(j+3,j+3)=C44;
end
DLOCAL1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic stiffnesses tensor of Material II in MPa (to be changed)
% Iron (BCC ferrite)
C11=232000;
C12=135000;
C44=116000;

DLOCAL2=zeros(6);
for j=1:3
    DLOCAL2(j,j)=C11;
    for i=1:3
        if (i~=j)
            DLOCAL2(i,j)=C12;
        end
    end
    DLOCAL2(j+3,j+3)=C44;
end
DLOCAL2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP SYSTEMS
nbsysFCC = 12;
load CFCSYST.DAT; 
f = [CFCSYST(:,2) CFCSYST(:,3) CFCSYST(:,4) CFCSYST(:,5) CFCSYST(:,6) CFCSYST(:,7)];
% sn: glide plane normal   % sb: slip direction
snFCC=zeros(nbsysFCC,3);
sbFCC=zeros(3,nbsysFCC);
for ss = 1 : nbsysFCC
    snFCC(ss,1) = f(ss,1);
    snFCC(ss,2) = f(ss,2);
    snFCC(ss,3) = f(ss,3);
    sbFCC(1,ss) = f(ss,4);
    sbFCC(2,ss) = f(ss,5);
    sbFCC(3,ss) = f(ss,6);
end
% Normalization
for ss = 1 : nbsysFCC 
    snFCC(ss,:) = snFCC(ss,:)/norm(snFCC(ss,:));
    sbFCC(:,ss) = sbFCC(:,ss)/norm(sbFCC(:,ss));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIP SYSTEMS
nbsysBCC = 24;
load CC24SYST.DAT; 
f = [CC24SYST(:,2) CC24SYST(:,3) CC24SYST(:,4) CC24SYST(:,5) CC24SYST(:,6) CC24SYST(:,7)];
% sn: glide plane normal   % sb: slip direction
snBCC=zeros(nbsysBCC,3);
sbBCC=zeros(3,nbsysBCC);
for ss = 1 : nbsysBCC
    snBCC(ss,1) = f(ss,1);
    snBCC(ss,2) = f(ss,2);
    snBCC(ss,3) = f(ss,3);
    sbBCC(1,ss) = f(ss,4);
    sbBCC(2,ss) = f(ss,5);
    sbBCC(3,ss) = f(ss,6);
end
% Normalization
for ss = 1 : nbsysBCC 
    snBCC(ss,:) = snBCC(ss,:)/norm(snBCC(ss,:));
    sbBCC(:,ss) = sbBCC(:,ss)/norm(sbBCC(:,ss));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate vectors of the GB frame expressed in the EBSD frame
y_GB=[cos(epsilon)*sin(psi);cos(epsilon)*cos(psi);-sin(epsilon)];
z_GB=[sin(epsilon)*sin(psi);sin(epsilon)*cos(psi);cos(epsilon)];
x_GB=cross(y_GB,z_GB);

BM=zeros(3);
BL=zeros(3);
for i = 1 : 3
    BM(i,1)=x_GB(i);
    BM(i,2)=y_GB(i);
    BM(i,3)=z_GB(i);
    BL(i,1)=x(i);
    BL(i,2)=y(i);
    BL(i,3)=z(i);    
end
Rgb = BL*BM^-1; % crossing matrix from from EBSD frame to GB frame
Tgb = Rgb^-1; % crossing matrix from GB frame to EBSD frame
    
% Crossing matrices from GB frame to crystal frames (R1 and R2)    
[R11] = mrot(z,phi1_1);
[R21] = mrot(x,phi_1);
[R31] = mrot(z,phi2_1);
R1 = R31*R21*R11*Tgb;
T1 = R1'; % crossing matrix from crystal FCC frame to GB frame
[R12] = mrot(z,phi1_2);
[R22] = mrot(x,phi_2);
[R32] = mrot(z,phi2_2);
R2 = R32*R22*R12*Tgb;
T2 = R2'; % crossing matrix from crystal BCC frame to GB frame
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passing the macroscopic stress tensor into the grain boundary frame
Sig = Rgb*SIGMA*Tgb;
 
% Macroscopic stress in vector notation in GB frame
Sig_Vect = [Sig(1,1);Sig(2,2);Sig(3,3);Sig(2,3);Sig(3,1);Sig(1,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passing the 6x6 elastic stiffnesses matrix of crystal I into the grain boundary frame
% conversion T1 (3x3) --> T1p (6x6)
T1p=zeros(6,6);

for j=1:3
    j1=1+floor(j/3);
    j2=2+floor(j/2);
    for i=1:3
        i1=1+floor(i/3);
        i2=2+floor(i/2);
% convention 12 --> 4, 31 --> 5, 23 --> 6 
%         T1p(i,j)=T1(i,j)^2;
%         T1p(i,j+3)=2*T1(i,j1)*T1(i,j2);
%         T1p(i+3,j)=T1(i1,j)*T1(i2,j);
%         T1p(i+3,j+3)=T1(i1,j1)*T1(i2,j2)+T1(i1,j2)*T1(i2,j1);
% convention 23 --> 4, 31 --> 5, 12 --> 6         
        T1p(i,j)=T1(i,j)^2;
        T1p(i,7-j)=2*T1(i,j1)*T1(i,j2);
        T1p(7-i,j)=T1(i1,j)*T1(i2,j);
        T1p(7-i,7-j)=T1(i1,j1)*T1(i2,j2)+T1(i1,j2)*T1(i2,j1);   
    end
    
end

C1=zeros(6,6);
for i=1:6
    for j=1:6
        for k=1:6
            for l=1:6
                C1(i,j)=C1(i,j)+DLOCAL1(k,l)*T1p(i,k)*T1p(j,l);
            end
        end

    end
end
% 4--6x6 elastic compliances matrix of crystal 1
S1=C1^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passing the 6x6 elastic stiffnesses matrix of crystal II into the grain boundary frame
% conversion T2 (3x3) --> T2p (6x6)
T2p=zeros(6,6);

for j=1:3
    j1=1+floor(j/3);
    j2=2+floor(j/2);
    for i=1:3
        i1=1+floor(i/3);
        i2=2+floor(i/2);
% convention 12 --> 4, 31 --> 5, 23 --> 6         
%         T2p(i,j)=T2(i,j)^2;
%         T2p(i,j+3)=2*T2(i,j1)*T2(i,j2);
%         T2p(i+3,j)=T2(i1,j)*T2(i2,j);
%         T2p(i+3,j+3)=T2(i1,j1)*T2(i2,j2)+T2(i1,j2)*T2(i2,j1);  
% convention 23 --> 4, 31 --> 5, 12 --> 6         
        T2p(i,j)=T2(i,j)^2;
        T2p(i,7-j)=2*T2(i,j1)*T2(i,j2);
        T2p(7-i,j)=T2(i1,j)*T2(i2,j);
        T2p(7-i,7-j)=T2(i1,j1)*T2(i2,j2)+T2(i1,j2)*T2(i2,j1);
    end
    
end

C2=zeros(6,6);
for i=1:6
    for j=1:6
        for k=1:6
            for l=1:6
                C2(i,j)=C2(i,j)+DLOCAL2(k,l)*T2p(i,k)*T2p(j,l);
            end
        end
    end
end
% 4--6x6 elastic compliances matrix of crystal 2
S2=C2^(-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s_tilde
ss = f2*S1 + f1*S2;

% Jump of S
saut_S = S2 - S1;

% Denominators
den=ss(1,1)*ss(3,5)^2+ss(3,3)*ss(1,5)^2+ss(5,5)*ss(1,3)^2-ss(1,1)*ss(3,3)*ss(5,5)-2*ss(1,3)*ss(1,5)*ss(3,5);

% Tensors G
G=zeros(6,6);
G(1,1)=(ss(3,5)^2-ss(3,3)*ss(5,5))/den; 
G(3,3)=(ss(1,5)^2-ss(1,1)*ss(5,5))/den;
G(5,5)=(ss(1,3)^2-ss(1,1)*ss(3,3))/den;
G(1,3)=(ss(1,3)*ss(5,5)-ss(1,5)*ss(3,5))/den; 
G(1,5)=(ss(1,5)*ss(3,3)-ss(1,3)*ss(3,5))/den;
G(3,5)=(ss(1,1)*ss(3,5)-ss(1,3)*ss(1,5))/den;
G(3,1)=G(1,3);
G(5,1)=G(1,5);
G(5,3)=G(3,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incompatibility stresses in elasticity (in GB frame)
sig_vectI = Sig_Vect + f2*G*saut_S*Sig_Vect;
sig_vectII = Sig_Vect - f1*G*saut_S*Sig_Vect;

% Stress tensors from stress vectors (in GB frame)
sigmaI = [[sig_vectI(1) sig_vectI(6) sig_vectI(5)];
          [sig_vectI(6) sig_vectI(2) sig_vectI(4)];
          [sig_vectI(5) sig_vectI(4) sig_vectI(3)]];
      
sigmaII = [[sig_vectII(1) sig_vectII(6) sig_vectII(5)];
          [sig_vectII(6) sig_vectII(2) sig_vectII(4)];
          [sig_vectII(5) sig_vectII(4) sig_vectII(3)]]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slip systems in GB frame  
snM1=zeros(nbsysFCC,3);
sbM1=zeros(3,nbsysFCC);
snM2=zeros(nbsysBCC,3);
sbM2=zeros(3,nbsysBCC);
LSF_I=zeros(nbsysFCC,1);
LSF_II=zeros(nbsysBCC,1);
for ksys = 1 : nbsysFCC        
    snM1(ksys,:) = snFCC(ksys,:)*R1;
    sbM1(:,ksys) = T1*sbFCC(:,ksys);            
end
for ksys = 1 : nbsysBCC        
    snM2(ksys,:) = snBCC(ksys,:)*R2;
    sbM2(:,ksys) = T2*sbBCC(:,ksys);            
end

% Resolved shear stresses
for ksys = 1 : nbsysFCC
    [SF_I] = cission(sigmaI,sbM1(:,ksys),snM1(ksys,:));
    LSF_I(ksys) = SF_I;
end
for ksys = 1 : nbsysBCC
    [SF_II] = cission(sigmaII,sbM2(:,ksys),snM2(ksys,:));
    LSF_II(ksys) = SF_II;
end
[tri_I,pos_I] = sort(abs(LSF_I));
[tri_II,pos_II] = sort(abs(LSF_II));
RSS_I = LSF_I(pos_I(12));
RSS_II = LSF_II(pos_II(12));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective elastic moduli tensor in GB frame
% <S>
S_aver = f1*S1 + f2*S2;

% Seff in GB frame
S_eff = S_aver - f1*f2*saut_S*G*saut_S;

% Ceff in GB frame
C_eff = S_eff^-1;

% Passing the 6x6 effective elastic stiffnesses matrix into the global frame
% conversion T (3x3) --> Tp (6x6)
Tp=zeros(6,6);

for j=1:3
    j1=1+floor(j/3);
    j2=2+floor(j/2);
    for i=1:3
        i1=1+floor(i/3);
        i2=2+floor(i/2);
% convention 12 --> 4, 31 --> 5, 23 --> 6         
%         Tp(i,j)=Tgb(i,j)^2;
%         Tp(i,j+3)=2*Tgb(i,j1)*Tgb(i,j2);
%         Tp(i+3,j)=Tgb(i1,j)*Tgb(i2,j);
%         Tp(i+3,j+3)=Tgb(i1,j1)*Tgb(i2,j2)+Tgb(i1,j2)*Tgb(i2,j1);  
% convention 23 --> 4, 31 --> 5, 12 --> 6         
        Tp(i,j)=Tgb(i,j)^2;
        Tp(i,7-j)=2*Tgb(i,j1)*Tgb(i,j2);
        Tp(7-i,j)=Tgb(i1,j)*Tgb(i2,j);
        Tp(7-i,7-j)=Tgb(i1,j1)*Tgb(i2,j2)+Tgb(i1,j2)*Tgb(i2,j1);
    end
end

Ceff=zeros(6,6);
for i=1:6
    for j=1:6
        for k=1:6
            for l=1:6
                Ceff(i,j)=Ceff(i,j)+C_eff(k,l)*Tp(i,k)*Tp(j,l);
            end
        end
    end
end
Ceff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Incompatibility stresses in plasticity (in GB frame)
if (RSS_I<0)
    gammaI=-gammaI;
end
if (RSS_II<0)
    gammaII=-gammaII;
end

% Plastic strain tensors
epsilonp_I=zeros(3);
epsilonp_II=zeros(3);
for i=1:3
    for j=1:3
        epsilonp_I(i,j) = gammaI*0.5*(sbM1(i,pos_I(12))*snM1(pos_I(12),j) + sbM1(j,pos_I(12))*snM1(pos_I(12),i));
        epsilonp_II(i,j) = gammaII*0.5*(sbM2(i,pos_I(12))*snM2(pos_I(12),j) + sbM2(j,pos_I(12))*snM2(pos_I(12),i));
    end
end

% Plastic strains in vector notation
epsilon_p1=[epsilonp_I(1,1);epsilonp_I(2,2);epsilonp_I(3,3);epsilonp_I(2,3);epsilonp_I(3,1);epsilonp_I(1,2)];
epsilon_p2=[epsilonp_II(1,1);epsilonp_II(2,2);epsilonp_II(3,3);epsilonp_II(2,3);epsilonp_II(3,1);epsilonp_II(1,2)];

% Jumps of plastic strain
saut_ep = epsilon_p2 - epsilon_p1;

% Stress vectors (in GB frame)
sig_vectI = Sig_Vect + f2*G*(saut_S*Sig_Vect + saut_ep);
sig_vectII = Sig_Vect - f1*G*(saut_S*Sig_Vect + saut_ep);

% Stress tensors from stress vectors (in GB frame)
sigmaI = [[sig_vectI(1) sig_vectI(6) sig_vectI(5)];
          [sig_vectI(6) sig_vectI(2) sig_vectI(4)];
          [sig_vectI(5) sig_vectI(4) sig_vectI(3)]];
      
sigmaII = [[sig_vectII(1) sig_vectII(6) sig_vectII(5)];
          [sig_vectII(6) sig_vectII(2) sig_vectII(4)];
          [sig_vectII(5) sig_vectII(4) sig_vectII(3)]]; 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macroscopic strain tensor in vector notation in GB frame
strain_vect = f1*(S1*sig_vectI+epsilon_p1) + f2*(S2*sig_vectII+epsilon_p2);

% Effective plastic strain tensor in vector notation in GB frame
ep_vect = strain_vect - S_eff*Sig_Vect;

% Macroscopic strain tensor from strain vector in GB frame
strain = [[strain_vect(1) 0.5*strain_vect(6) 0.5*strain_vect(5)];
          [0.5*strain_vect(6) strain_vect(2) 0.5*strain_vect(4)];
          [0.5*strain_vect(5) 0.5*strain_vect(4) strain_vect(3)]];

% Effective plastic strain tensor in the GB frame
ep = [[ep_vect(1) 0.5*ep_vect(6) 0.5*ep_vect(5)];
          [0.5*ep_vect(6) ep_vect(2) 0.5*ep_vect(4)];
          [0.5*ep_vect(5) 0.5*ep_vect(4) ep_vect(3)]];      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In global frame
% Passing the stress tensors into the global frame
SIGMAI = Tgb*sigmaI*Rgb
SIGMAII = Tgb*sigmaII*Rgb

% Passing the macroscopic strain tensor into the global frame
STRAIN = Tgb*strain*Rgb

% Macroscopic stress in vector notation in the global frame
SIGMA_Vect = [SIGMA(1,1);SIGMA(2,2);SIGMA(3,3);SIGMA(2,3);SIGMA(3,1);SIGMA(1,2)];

% Macroscopic strain in vector notation in the global frame
STRAIN_Vect = [STRAIN(1,1);STRAIN(2,2);STRAIN(3,3);2*STRAIN(2,3);2*STRAIN(3,1);2*STRAIN(1,2)];

% Effective plastic strain tensor in vector notation in the global frame
Ep_Vect = STRAIN_Vect - (Ceff^-1)*SIGMA_Vect;

% Effective plastic strain tensor in the global frame 
Ep = [[Ep_Vect(1) 0.5*Ep_Vect(6) 0.5*Ep_Vect(5)];
      [0.5*Ep_Vect(6) Ep_Vect(2) 0.5*Ep_Vect(4)];
      [0.5*Ep_Vect(5) 0.5*Ep_Vect(4) Ep_Vect(3)]]
  

  







