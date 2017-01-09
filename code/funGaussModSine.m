function F = funGaussModSine(x0,xdata)

% The matching pursuit approach based on the modulated Gaussian pulse for
% efficient guided-wave damage inspection Jin-Chul Hong, Kyung Ho Sun and
% Yoon Young Kim1 Published 4 May 2005 

A      = x0(1); % amplitude
fc     = x0(2); % driving frequency
u      = x0(3); % u
sigma  = x0(4); % k
phi    = x0(5); % phi

F = exp(-1/2*((xdata-u)./(sigma)).^2).*cos(2*pi*fc*(xdata-u)+ phi);
z = hilbert(real(F));
r = (real(z).^2+imag(z).^2).^0.5;
F = A/max(r)*real(F);