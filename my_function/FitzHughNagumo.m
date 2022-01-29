function dXdt = FitzHughNagumo(X, varargin)
    x    = X(1);
    w    = X(2);
    
    if length(varargin)==1    
        par  = varargin{1};
    else
        par  = varargin;
    end

    Iext = par{1};
    tau  = par{2};
    a    = par{3};
    b    = par{4};
    
    dxdt = x - (x.^3) - w +  Iext;
    dwdt = 1/tau .* (x - a - b .* w);

    dXdt = [dxdt, dwdt];
end