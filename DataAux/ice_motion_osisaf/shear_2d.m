function [curlz, shearz] = shear_2d(u,v)
%%%%% adapted from Matlab curv function

u = double(u);
v = double(v);

[~, py] = gradient(u);
[qx, ~] = gradient(v);


curlz = 0.5*(qx-py); 
shearz= 0.5*(qx+py); 