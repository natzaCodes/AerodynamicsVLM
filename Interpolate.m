function xvec = Interpolate(x,x0,x1,y0,y1)
if (x0<x1)
    if (x<x0 || x>x1)
        fprintf('\n Warning! Interpolation point (%d) is outside the given interval [%d,%d] !!!!\n',x,x0,x1);
    end
end
if (x1<x0)
    if (x>x0 || x<x1)
        fprintf('\n Warning! Interpolation point (%d) is outside the given interval [%d,%d] !!!!\n',x,x0,x1);
    end
end
xvec = y0 + (x-x0) * (y1 - y0) / (x1-x0);
end