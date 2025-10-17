function [n_eff, status] = te_slab_solver(t_m, lambda_m, n_core, n_clad, mode)

%Default
n_eff  = NaN;
status = "unset";

%Helpers
a  = t_m/2;                 %half thickness
k0 = (2*pi)/lambda_m;       %free-space wavenumber
V  = k0*a*sqrt(n_core^2 - n_clad^2);   %normalized frequency

%Check TE1 existence
mode = upper(string(mode));
if mode == "TE1" && V < (pi/2)
    status = "TE1 below cutoff";
    return;
end

%Dispersion function definition
if mode == "TE0"
    f = @(u) u.*tan(u) - sqrt(max(V^2 - u.^2, 0));
else % "TE1"
    f = @(u) -u.*cot(u) - sqrt(max(V^2 - u.^2, 0));
end

%ALL OF THIS IS CHATGPT FROM HERE TO NUMERIC ROOT (line 27-> line 43)
delta = 1e-3;           %safety delta
if mode == "TE0"
    if V <= (pi/2)
        uL = delta; uR = max(V - delta, 2*delta);
    else
        uL = delta; uR = (pi/2) - delta;
    end
else % TE1
    uL = (pi/2) + delta;
    uR = min(V - delta, pi - delta);
    if uR <= uL
        status = "no convergence (narrow TE1 bracket)";
        return;
    end
end
u0 = 0.5*(uL + uR);     %single starting point (midpoint)

%Find u roots
try
    u_star = fzero(f, u0);
catch
    status = "no convergence";
    return;
end
if ~(isfinite(u_star) && u_star > 0 && u_star < V)
    status = "no convergence";
    return;
end

%Find effective refraction
w   = sqrt(max(V^2 - u_star^2, 0));
kx  = u_star / a;                      %core transverse k
beta = sqrt(max((k0*n_core)^2 - kx^2, 0));  %propagation constant
n_eff = beta / k0;                     %effective index

status = "ok";
end
