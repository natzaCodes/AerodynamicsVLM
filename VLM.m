function [coef,mesh] = VLM(geo,flow,mesh)

	%% Panels´ Geometrical Properties
	% ----------------------------------------------------------------------- %
	mesh.coor.p1 = zeros(mesh.npx,mesh.npy,3);
	mesh.coor.p2 = zeros(mesh.npx,mesh.npy,3);
	mesh.coor.cp = zeros(mesh.npx,mesh.npy,3);
	mesh.coor.cpbb = zeros(mesh.npx,mesh.npy,3);
    mesh.coor.corners = zeros(mesh.npx+1,mesh.npy+1,3);
	dy = geo.b2 / mesh.npy;
	tglam = tan(geo.sweep);
	tgphi = tan(geo.phi);
	
	%% Matrices of circulations and coefficients
	A = zeros(mesh.npt,mesh.npt);
	RHS = zeros(mesh.npt,1);
	

	tic
	%% Geometry discretization
	chc = zeros(mesh.npy,1); %vector of chords spanwise, needed for local Cl
	for j = 1:mesh.npy
		y0 = dy * (j-1);% leftmost position y of the panel
		
		chl = Interpolate(y0,0,geo.b2,geo.cr,geo.ct);
		chr = Interpolate(y0+dy,0,geo.b2,geo.cr,geo.ct);
		chc(j) = Interpolate(y0+dy/2,0,geo.b2,geo.cr,geo.ct);
		dxl = chl / mesh.npx; %delta x of panel in the left
		dxr = chr / mesh.npx; %delta x of panel in the right
		dxc = chc(j) / mesh.npx; %delta x of panel in the center
		for i = 1:mesh.npx
			% Point 1 at c/4
			mesh.coor.p1(i,j,1) = y0 * tglam + dxl*(i-0.75);
			mesh.coor.p1(i,j,2) = y0;
			mesh.coor.p1(i,j,3) = y0 * tgphi;
		
			%Point 2 at c/4
			mesh.coor.p2(i,j,1) = (y0+dy) * tglam + dxr*(i-0.75);
			mesh.coor.p2(i,j,2) = y0+dy;
			mesh.coor.p2(i,j,3) = (y0+dy) * tgphi;
		
			%Control point at 3c/4
			mesh.coor.cp(i,j,1) = (y0+dy/2) * tglam + dxc*(i-0.25);
			mesh.coor.cp(i,j,2) = y0+dy/2;
			mesh.coor.cp(i,j,3) = (y0+dy/2) * tgphi;
			
			%Control point at c/4
			mesh.coor.cpbb(i,j,1) = (y0+dy/2) * tglam + dxc*(i-0.75);
			mesh.coor.cpbb(i,j,2) = y0+dy/2;
			mesh.coor.cpbb(i,j,3) = (y0+dy/2) * tgphi;
            
            % Corners for plotting the panel boundaries
            mesh.coor.corners(i,j,1) = y0 * tglam + dxl*(i-1);
			mesh.coor.corners(i,j,2) = y0;
			mesh.coor.corners(i,j,3) = y0 * tgphi;
            if(i == mesh.npx)
                mesh.coor.corners(i+1,j,1) = y0 * tglam + dxl*(i);
                mesh.coor.corners(i+1,j,2) = y0;
                mesh.coor.corners(i+1,j,3) = y0 * tgphi;
            end
            if(j == mesh.npy)
                mesh.coor.corners(i,j+1,1) = (y0+dy) * tglam + dxr*(i-1);
                mesh.coor.corners(i,j+1,2) = (y0+dy);
                mesh.coor.corners(i,j+1,3) = (y0+dy) * tgphi;
            end
            if(i == mesh.npx && j == mesh.npy)
                mesh.coor.corners(i+1,j+1,1) = (y0+dy) * tglam + dxr*(i);
                mesh.coor.corners(i+1,j+1,2) = (y0+dy);
                mesh.coor.corners(i+1,j+1,3) = (y0+dy) * tgphi;
            end
		end
	end
	% Time enlapsed in the geometric loop
%	T_geom_loop=toc % (s)  
	
	cp=[0,0,0];%Control point
	cps=[0,0,0];%Control point at the symmetry
	p1=[0,0,0];%Point 1 of the horseshoe
	p2=[0,0,0];%Point 2 of the horseshoe
	v_in = [0,0,0];%Induced velocity
	v_ins = [0,0,0];%Induced velocity by symmetry
	v_int = [0,0,0];%total induced velocity
	k = [0,-sin(geo.phi),cos(geo.phi)];
	countij=0;
	countmn=0;
	for j = 1:mesh.npy
		for i = 1:mesh.npx
			countij = countij + 1;
			countnm = 0;
			cp(1) = mesh.coor.cp(i,j,1); %cpx
			cp(2) = mesh.coor.cp(i,j,2); % cpy
			cp(3) = mesh.coor.cp(i,j,3);
			% Symmetric point 
			cps = cp;
			cps(2) = -cp(2); %Symmetric point has negative y
			for m = 1:mesh.npy
				for n = 1:mesh.npx
					countnm = countnm + 1;
					p1(1) = mesh.coor.p1(n,m,1);
					p1(2) = mesh.coor.p1(n,m,2);
					p1(3) = mesh.coor.p1(n,m,3);
					
					p2(1) = mesh.coor.p2(n,m,1);
					p2(2) = mesh.coor.p2(n,m,2);
					p2(3) = mesh.coor.p2(n,m,3);
					v_in = VORHS(cp,p1,p2,1);
					
					% Account for symmetry
					v_ins = VORHS(cps,p1,p2,1);
					
					v_int(1) = v_in(1) + v_ins(1);
					v_int(2) = v_in(2) - v_ins(2);
					v_int(3) = v_in(3) + v_ins(3);
					
					A(countij,countnm) = dot(v_int,k);	% TODO revise indices order
				end
			end
			
		end
	end
%	T_coef_loop=toc
	%% Right hand side of the equation
	RHS = RHS - flow.tas * sin(flow.alpha) * cos(geo.phi);
	
	%% Solving the system of equations
	GAMMA = linsolve(A,RHS);
%	T_solve_loop=toc
    
    %% Induced velocities
    u = zeros(mesh.npx,mesh.npy,3);% Total seen velocity by each panel
	vin = zeros(mesh.npx,mesh.npy,3);% Induced velocity in each panel
	alphain = zeros(mesh.npx,mesh.npy);% Induced angle of attack at each panel
	countnm=0;
    for j = 1:mesh.npy
		for i = 1:mesh.npx
			countnm = 0;
			cp(1) = mesh.coor.cpbb(i,j,1);
			cp(2) = mesh.coor.cpbb(i,j,2);
			cp(3) = mesh.coor.cpbb(i,j,3);
			% Symmetric point
			cps = cp;
			cps(2) = -cp(2);
			for m = 1:mesh.npy
				for n = 1:mesh.npx
					countnm = countnm + 1;
					p1(1) = mesh.coor.p1(n,m,1);
					p1(2) = mesh.coor.p1(n,m,2);
					p1(3) = mesh.coor.p1(n,m,3);
					
					p2(1) = mesh.coor.p2(n,m,1);
					p2(2) = mesh.coor.p2(n,m,2);
					p2(3) = mesh.coor.p2(n,m,3);
					v_in = VOTRAIL(cp,p1,p2,GAMMA(countnm));
					
					% Account for symmetry
					v_ins = VOTRAIL(cps,p1,p2,GAMMA(countnm));
					
					vin(i,j,1) = vin(i,j,1) + v_in(1) + v_ins(1);
                    vin(i,j,2) = vin(i,j,2) + v_in(2) - v_ins(2);
                    vin(i,j,3) = vin(i,j,3) + v_in(3) + v_ins(3);
				end
			end
			u(i,j,1) = vin(i,j,1) + flow.tas * cos(flow.alpha);
			u(i,j,3) = vin(i,j,3) + flow.tas * sin(flow.alpha);
			alphain(i,j) = vin(i,j,3) / flow.tas;
		end
	end
%	T_indv_loop=toc
	
	%% Total forces and lift spanwise distribution
	Ly = zeros(mesh.npy,1);
	coef.Cly = zeros(mesh.npy,1);
	Dy = zeros(mesh.npy,1);
	coef.Cdy = zeros(mesh.npy,1);
	coef.dcp = zeros(mesh.npx,mesh.npy);% pressure coefficient difference
	coef.F = zeros(mesh.npx,mesh.npy,3);% pressure coefficient difference
	countij = 0;
	dL = 0;
	L = 0;
	dD = 0;
	D = 0;
    dyvec = [0,dy,0];
	dyvec2 = [0,0,0];
    
    
	for j = 1:mesh.npy
		for i = 1:mesh.npx
			countij = countij + 1;
			dL = GAMMA(countij) * flow.rho * dy * flow.tas;
			dD = - GAMMA(countij) * flow.rho * dy * vin(i,j,3);
			up(1) = u(i,j,1);
			up(2) = u(i,j,2);
			up(3) = u(i,j,3);
			dyvec2 = [mesh.coor.p2(i,j,:) - mesh.coor.p1(i,j,:)];
			dyvec(1) = dyvec2(1,1,1);
			dyvec(2) = dyvec2(1,1,2);
			dyvec(3) = dyvec2(1,1,3);
			dF = flow.rho * GAMMA(countij) * cross(up,dyvec);
			coef.F(i,j,:) = dF;
			dL = dF(3)*cos(flow.alpha)-dF(1)*sin(flow.alpha);
			dD = dF(3)*sin(flow.alpha)+dF(1)*cos(flow.alpha);
			coef.dcp(i,j) = dL / (dy * (mesh.coor.cp(i,j,1)-mesh.coor.cpbb(i,j,1)) * 2 * 0.5 * flow.rho * flow.tas * flow.tas);
			Ly(j) = Ly(j) + dL/dy;			
			L = L + 2 * dL;
			Dy(j) = Dy(j) + dD/dy;			
			D = D + 2 * dD;
			
		end
		coef.Cly(j) = 2 * Ly(j) / (flow.rho * flow.tas * flow.tas * chc(j));
		coef.Cdy(j) = 2 * Dy(j) / (flow.rho * flow.tas * flow.tas * chc(j));
	end
	
	coef.CL = 2 * L / (flow.rho * flow.tas * flow.tas * geo.sw);
	coef.CD = 2 * D / (flow.rho * flow.tas * flow.tas * geo.sw);
%	T_coef_loop=toc
end