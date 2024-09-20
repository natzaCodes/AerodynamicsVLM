function [vin] = VOTRAIL(cp,p1,p2,gamma)
%VOTRAIL calculates the velocity induced at point cp by trailing vortices of HS vortex defined by
% points p1 and p2
%   @param cp: Control point
%	@param p1: point 1
%	@param p2: point 2
%	@return v: induced velocity
	
	p1inf = p1;
	p1inf(1) = p1inf(1) + 1e10;
	
	p2inf = p2;
	p2inf(1) = p2inf(1) + 1e10;
	vin = VORTXL(cp,p1inf,p1,gamma) + VORTXL(cp,p2,p2inf,gamma);	
end

