clc;
close all;

Aov = pi*(80/180);
ax=20;
ay=20;
S1x=cos(pi*(ax/180));
S2y=cos(pi*(ay/180));
S3z=1-sqrt(S1x^2 + S2y^2);
func_cone_plotting(Aov,S1x,S2y,S3z)
