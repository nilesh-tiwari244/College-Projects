function [h]=func_4_channel_coeff_march28trial_twoantenna_in_room(length,theta,gamma,R,psi,fov)

h = zeros(length); 
for i=1:1:length
    for j=1:1:length
if abs(theta(i,j)) < fov/2
    h(i,j)=((gamma+1)*(1/(R(i,j)^2))*(0.5*10e-4)*(1/pi)* cos(theta(i,j))*(cos(psi(i,j)))^(gamma));
else
    h(i,j)=0;
end
    end
end