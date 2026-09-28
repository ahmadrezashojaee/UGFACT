function [poro, perm] = generateReservoirPropertiesFFTMA(NX, NY, NZ, ...
                Lx, Ly, Lz, corrX, corrY, corrZ, ...
                poro_median, poro_std, poro_min, ...
                perm_median, perm_std, rho, seed)
% FFT-MA generator for porosity & permeability (mD) with correct correlation scaling

if ~isempty(seed)
    rng(seed);
end

% --- Physical spacing ---
dx = Lx/NX; 
dy = Ly/NY; 
dz = Lz/NZ;

% --- Coordinates in meters ---
[x,y,z] = ndgrid((0:NX-1)*dx, (0:NY-1)*dy, (0:NZ-1)*dz);

% --- Periodic wrap distances ---
hx = min(x, Lx-x);
hy = min(y, Ly-y);
hz = min(z, Lz-z);

% --- Exponential covariance ---
C = exp(-sqrt((hx./corrX).^2 + (hy./corrY).^2 + (hz./corrZ).^2));

% --- FFT of covariance ---
CF = abs(fftn(C));
CF = sqrt(CF / max(CF(:))); % normalize

% --- Random Gaussian field ---
Z = randn(NX,NY,NZ);
Zf = fftn(Z);

% Multiply in Fourier space
field = real(ifftn(Zf .* CF));

% Normalize
field = (field - mean(field(:))) / std(field(:));

% --- Porosity ---
poro = poro_median + poro_std*field;
poro = max(poro, poro_min);
poro = min(poro, 1.0);

% --- Permeability (lognormal with scatter) ---
noiseField = randn(NX,NY,NZ);
Z_perm = rho*field + sqrt(1-rho^2)*noiseField;
Z_perm = (Z_perm - mean(Z_perm(:))) / std(Z_perm(:));

logk_median = log10(perm_median);
logk_std    = log10(perm_median + perm_std) - logk_median;
logk = logk_median + logk_std*Z_perm;
perm = 10.^logk; % mD

% Flatten
poro = poro(:);
perm = perm(:);

end
