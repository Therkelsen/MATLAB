clc
clear
close all
%% Parametre
fa = 1000; % Afskæringsfrekvens [Hz]
fs = 8000; % Samplingsfrekvens [Hz]
T = 1/fs;  % Sampleinterval [s]
%% 1. Normerede filter
s = tf('s');
H1 = 0.49417/(s+0.49417);
H2 = 0.99421/(s^2+0.49417*s+0.99421);
%% 2. Poler for normerede filter
sigma1_nom = pole(H1);
pole2 = pole(H2);
sigma2_nom = real(pole2(1));
omega2_nom = imag(pole2(1));
%% 3. Poler for denormerede filter
sigma1 = sigma1_nom*2*pi*fa;
sigma2 = sigma2_nom*2*pi*fa;
omega2 = omega2_nom*2*pi*fa;
H1_den = -sigma1/(s-sigma1);
H2_den = (sigma2^2+omega2^2)/((s-(sigma2+1i*omega2))*(s-(sigma2-1i*omega2)));
%% 4. Bestem den digitale overføringsfunktions koefficienter
z = tf('z',T);
% Koefficienter for G1
b11 = -exp(sigma1*T);
a01 = 1+b11;
G1 = a01/(z+b11);
G1t = z*G1;
% Koefficienter for G2
b12 = -2*exp(sigma2*T)*cos(omega2*T);
b22 = exp(2*sigma2*T);
a02 = 1+b12+b22;
G2 = a02/(z^2+b12*z+b22);
G2t = z^2*G2;

%% Bode Plot 
figure
bode(H1_den,G1,G1t)
figure
bode(H2_den,G2,G2t)
figure
bode(H1_den*H2_den,G1*G2,G1t*G2t)
%% Verifikation
c2d(H1_den,T,'matched')
c2d(H2_den,T,'matched')
c2d(H1_den*H2_den,T,'matched')
c2d(H1_den,T,'impulse')
c2d(H2_den,T,'impulse')