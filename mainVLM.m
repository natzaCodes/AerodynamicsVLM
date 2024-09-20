%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following script executes a VLM to study the aerodynamic properties
% of a finite swept wing.
% Authors: Albert Canovas Cots & Natalia Zalewska
% Date: 10.10.2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
%% Figure's positions
% ----------------------------------------------------------------------- %
scrsz = get(0,'ScreenSize');
scW=scrsz(3); % Screen Width (px)
scH=scrsz(4); % Screen Height (px)
% figure('Name','NameFigure','Position',[x0 y0 width height])
%figure('Name','Wing','Position',[1 scH/4 scW/2 scH/2.5])
%figure('Name','Distributions','Position',[scW/2 scH/4 scW/2 scH/2.5])

% Plots configuration
font=15;  % font size
lw=2;     % linewidth
szp=100;  % scatter size plot


%% Mesh Configuration
% ----------------------------------------------------------------------- %
mesh.npx = 4; % Number of panels in the streamwise direction
mesh.npy = 32; % Number of panels in the SEMI-spanwise direction
mesh.npt = mesh.npx*mesh.npy;

%% Flight Conditions
% ----------------------------------------------------------------------- %

flow.alpha = deg2rad(4); % Angle of Attack (rad)
flow.tas  = 1;       % Wind Speed (m/s)
flow.rho   = 1.225;    % Air Density (kg/m3)

%% Wing Geometry
% ----------------------------------------------------------------------- %
% Model 14 Naca TR 2445
% geo.sweep = deg2rad(0);            % Leading edge sweep angle (rad)
% geo.phi = deg2rad(0);  % Dihedral angle (rad)
% geo.b2   = 0.762;            % SEMI-span (m)
% geo.cr  = 0.254;            % Root chord (m)
% geo.ct  = 0.254;  % Tip chord (m)
% geo.sw   = (geo.ct + geo.cr)*geo.b2;            % Wing surface (m)
% geo.ar  = ((2*geo.b2)^2)/geo.sw;            % Wing Aspect Ratio (-)

% B787
geo.sweep = deg2rad(35);            % Leading edge sweep angle (rad)
geo.phi = deg2rad(11);  % Dihedral angle (rad)
geo.b2   = 60.12/2;            % SEMI-span (m)
geo.cr  = 11.9;            % Root chord (m)
geo.ct  = 2.37;  % Tip chord (m)
geo.sw   = (geo.ct + geo.cr)*geo.b2;            % Wing surface (m)
geo.ar  = ((2*geo.b2)^2)/geo.sw;            % Wing Aspect Ratio (-)

% nx = [1,2,4,8,16];
% ny = [1,2,4,8,16,32,64,128];
% for i = 1:size(nx,2)
% 	for j = 1:size(ny,2)
% 		mesh.npx = nx(i); % Number of panels in the streamwise direction
% 		mesh.npy = ny(j); % Number of panels in the SEMI-spanwise direction
% 		mesh.npt = mesh.npx*mesh.npy;
% 		[coef,mesh] = VLM(geo,flow,mesh);
% 		CD(i,j) = coef.CD;
% 		CL(i,j) = coef.CL;
% 		mesh.npt
% 	end
% end
%  
% [NX,NY] = meshgrid(nx,ny);
% figure
% surf(NX,NY,CD')
% xlabel("Number of panels chordwise")
% ylabel("Number of panels spanwise")
% zlabel("Induced drag coefficient")
% 
% figure
% surf(NX,NY,CL')
% xlabel("Number of panels chordwise")
% ylabel("Number of panels spanwise")
% zlabel("Lift coefficient")

[coef,mesh] = VLM(geo,flow,mesh);

alpha = deg2rad(linspace(0,10,10)); %AoA
Clvec = zeros(1,10);
a=1;

for i=1:length(alpha)
    flow.alpha = alpha(i);
    [coef2, mesh2] = VLM(geo,flow,mesh);
    Clvec(1,a) = coef2.CL;
    a=a+1;
    
end

flow.alpha = deg2rad(4);
sweep = deg2rad(linspace(0,50,20));
Clsweep = zeros(1,20);
a=1;
for i=1:length(sweep)
    geo.sweep = sweep(i);
    [coef3, mesh3] = VLM(geo,flow,mesh);
    Clsweep(1,a) = coef3.CL;
    a=a+1;
    
end

geo.sweep = deg2rad(35);
AR = linspace(1,geo.ar,20);
Clar = zeros(1,20);
a=1;
for i=1:length(AR)
    geo.ar = AR(i);
    [coef4, mesh4] = VLM(geo,flow,mesh);
    Clar(1,a) = coef4.CL;
    a=a+1;
    
end

geo.b2   = 60.12/2;        
geo.cr  = 11.9;      
geo.ct  = 2.37;
geo.sw   = (geo.ct + geo.cr)*geo.b2;
geo.ar  = ((2*geo.b2)^2)/geo.sw;            % Wing Aspect Ratio (-)
Phi = deg2rad(linspace(1,20,20));
Clphi = zeros(1,20);
a=1;
for i=1:length(Phi)
    geo.phi = Phi(i);
    [coef5, mesh5] = VLM(geo,flow,mesh);
    Clphi(1,a) = coef5.CL;
    a=a+1;
    
end

%% Plots and figures
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)
set(groot,'defaulttextFontSize',12)
set(groot,'defaultLegendFontSize',12)

% Mesh representation
figure
scatter3(reshape(mesh.coor.p1(:,:,2),[],1),reshape(mesh.coor.p1(:,:,1),[],1),reshape(mesh.coor.p1(:,:,3),[],1),'b','filled')
hold on;
scatter3(reshape(mesh.coor.cp(:,:,2),[],1),reshape(mesh.coor.cp(:,:,1),[],1),reshape(mesh.coor.cp(:,:,3),[],1),'r','*')
scatter3(reshape(mesh.coor.p2(:,:,2),[],1),reshape(mesh.coor.p2(:,:,1),[],1),reshape(mesh.coor.p2(:,:,3),[],1),'b','filled')

surf(mesh.coor.corners(:,:,2),mesh.coor.corners(:,:,1),mesh.coor.corners(:,:,3),'FaceColor','none')
set(gca, 'YDir','reverse')
set(gca,'DataAspectRatio',[1 1 1])
xlabel('y (m)')
ylabel('x (m)')
zlabel('z (m)')
title('Representation of the mesh')
legend('Bound vortex end points','Control points')

% Force vectors
figure
surf(mesh.coor.corners(:,:,2),mesh.coor.corners(:,:,1),mesh.coor.corners(:,:,3),'FaceColor','none')
hold on;
quiver3(mesh.coor.cpbb(:,:,2),mesh.coor.cpbb(:,:,1),mesh.coor.cpbb(:,:,3),coef.F(:,:,2),coef.F(:,:,1),coef.F(:,:,3),10)
set(gca, 'YDir','reverse')
set(gca,'DataAspectRatio',[1 1 1])
xlabel('y (m)')
ylabel('x (m)')
zlabel('z (m)')
title('Force acting on each panel')

% Pressure distribution
figure
surf(mesh.coor.corners(:,:,2),mesh.coor.corners(:,:,1),mesh.coor.corners(:,:,3),coef.dcp,'EdgeColor','none')
caxis([min(min(coef.dcp)) max(max(coef.dcp))])
colormap(jet)
colorbar;
set(gca, 'YDir','reverse')
set(gca,'DataAspectRatio',[1 1 1])
xlabel('y (m)')
ylabel('x (m)')
zlabel('z (m)')
title('Pressure coefficient difference distribution')

% Lift and drag coefficient distribution along the span
figure
plot(mesh.coor.cp(1,:,2),coef.Cly,'linewidth',1.5)
hold on;
grid on;
xlabel(' y (m)')
ylabel(' $C_l$')
title('Lift distribution')
figure
plot(mesh.coor.cp(1,:,2),coef.Cdy,'linewidth',1.5)
hold on;
grid on;
xlabel(' y (m)')
ylabel(' $C_{Di}$')
title('Induced drag distribution')


% Lift coefficient evolution for AoA
figure
plot(rad2deg(alpha),Clvec(1,:),'linewidth',1.5)
hold on;
grid on;
xlabel('$\alpha$ (deg)')
ylabel(' $C_{L}$')
title('Lift coefficient over angle of attack')

% Lift coefficient evolution for Sweep angle
figure
plot(rad2deg(sweep),Clsweep(1,:),'linewidth',1.5)
hold on;
grid on;
xlabel('sweep angle(deg)')
ylabel(' $C_{L}$')
title('Lift coefficient over sweep angle')

% Lift coefficient evolution for AR
figure
plot(AR,Clar(1,:),'linewidth',1.5)
hold on;
grid on;
xlabel('AR')
ylabel(' $C_{L}$')
title('Lift coefficient over aspect ratio changes')

% Lift coefficient evolution with the dihedral angle
figure
plot(rad2deg(Phi),Clphi(1,:),'linewidth',1.5)
hold on;
grid on;
xlabel('$\phi$ (deg)')
ylabel(' $C_{L}$')
title('Lift coefficient evolution with the dihedral angle')

