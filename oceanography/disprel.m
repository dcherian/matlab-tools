% Makes mega dispersion relation chart

% Limits & parameters
k = -100:100;
N2 = 1e-5;
f0 = 1e-4;
beta = 2e-13;
g = 9.81;
H = 5000;
figure

% surface gravity waves

% Internal waves
%wigw1 = sqrt(f0^2 + N2*(k.^2)./(k.^2 + (pi/H)^2));

% Rossby waves

for n=1:5
	wr(n,:) = -beta.*k./(k.^2 + (n*pi/(N2*H^2/f0^2)));
end


% Equatorial waves

% Continental shelf waves

% Edge waves

% Kelvin Waves
wkl = sqrt(g*H)*k;

plot(k,wkl);
hold on
plot(k,wr);