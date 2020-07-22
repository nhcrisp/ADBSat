folderpath = 'C:\Users\mbgm9nc4\Documents\MATLAB\ADBSat_v2.1\inou\results\SOAR_Fin20_Counter_Default_sentman_16315720190228547\aedb';
file_name = 'SOAR_Fin20_Counter_Default';
modname = 'SOAR_Fin20_Counter_Default';
param = 'cp';

aoa = 10*pi/180;
aos = 10*pi/180;

% Load model mesh
eval(['load ', modname,'.mat']);
eval(['x = ', modname, '.XData;']);
eval(['y = ', modname, '.YData;']);
eval(['z = ', modname, '.ZData;']);

% Load results for indicated aoa, aos

%eval(['load ', folderpath, '/',file_name,'_a',mat2str(aoa*180/pi),'s',mat2str(aos*180/pi),'.mat ',param, ' Lref']);
%%
L_bw = dcmbody2wind(aoa, aos);
%[cos(aoa)*cos(aos) -cos(aoa)*sin(aos) -sin(aoa);
 %       sin(aos)          cos(aos)                0   ;
  %      sin(aoa)*cos(aos) -sin(aoa)*sin(aos)  cos(aoa)];
    
L_gb = [1 0 0; 0 -1 0; 0 0 -1];    

L_gw = L_gb*L_bw;

L_fb = [-1 0 0; 0 1 0; 0 0 -1];

axlenght = 0.7*1;

x1 = [0;0;0];
y1 = x1;
z1 = x1;

g = [1;1;1];

b = L_gb' * g;
f = L_fb * b;
w = L_gw' * g;

figure 
% Wind Axes
quiver3(x1,y1,z1,w,angle2dcm(0,pi/2,0,'ZYX')*w,angle2dcm(pi/2,0,0,'ZYX')*w,axlenght)
hold on
% Body Axes
quiver3(x1,y1,z1,[b(1);0;0],[0;b(2);0],[0;0;b(3)],axlenght)
hold on
% Geometric
quiver3(x1,y1,z1,[g(1);0;0],[0;g(2);0],[0;0;g(3)],axlenght)
% Flight
quiver3(x1,y1,z1,[f(1);0;0],[0;f(2);0],[0;0;f(3)],axlenght)

legend('Wind', 'Body', 'Geometric', 'Flight', 'Location','NorthWest')

h = patch(x, y, z, cp);
axis equal
axis off
