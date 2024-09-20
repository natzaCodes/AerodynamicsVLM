function vin = VORTXL(cp,p1,p2,gamma)
%VORTXL calculates the velocity induced at point cp by line vortex defined by
% points p1 and p2
%   @param cp: Control point
%	@param p1: point 1
%	@param p2: point 2
%	@return v: induced velocity
	r0 = p2-p1;
	r1 = cp-p1;
	r2 = cp-p2;   
	cros = cross(r1,r2);
	crosnorm = norm(cros);
	vin = gamma / (4 * pi) * cros / (crosnorm^2) * (dot(r0,r1)/norm(r1)-dot(r0,r2)/norm(r2));
end

