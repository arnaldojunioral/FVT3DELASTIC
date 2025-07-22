% A THREE_DIMENSIONAL FINITE-VOLUME THEORY CODE FOR STRESS ANALYSIS IN
% CONTINUUM ELASTIC STRUCTURES

function FVT3DELASTIC(nx, ny, nz)
%% ____________________________________________________________USER-DEFINED
[L, H, B] = deal(500, 100, 100);                                           % Cantelever beam dimensions (mm)
[E, nu]  = deal(150e3, 0.3);                                               % material's Young modulus (N/mm2: MPa) and Poisson's ratio
P  = -2e3;                                                                 % applied load (N)
[pb, amp] = deal('flexure', 1e3);                                          % pb: 'flexure', 'torsion', 'torsion-flexure'; amp: deformation amplification

%% _____________________________________________________BOUNDARY CONDITIONS
[l, h, b] = deal(L/nx, H/ny, B/nz);                                        % subvolume dimensions
[nodes, dof, ndof, iK, jK] = DOFassembly(nx, ny, nz);                      % degrees of freedom

% Essential boundary conditions
supp = unique(dof(1:nx:end, 1:3));                                         % fixed DOFs at x = 0
fdof = setdiff(dof(:), supp(:));                                           % free DOFs

% Default: full material
rho = ones(nx, ny, nz);

switch pb
    case 'flexure'
        F = sparse(dof(nx:nx:end, 5), 1, P / (ny * nz), ndof, 1);          % Load in y-direction
        cp = [1, 6];                                                       % sigma_11 and sigma_12

    case 'torsion'
        if H ~= B
            error('Limited to square section: no warping!');
        end

        J = (B * H / 12) * (B^2 + H^2);
        val = P / J * b * h;

        [jpos, kpos] = meshgrid(0:ny-1, 0:nz-1);
        y = jpos(:) * h - H/2 + h/2;
        z = kpos(:) * b - B/2 + b/2;

        load = reshape([-val * y, val * z]', [], 1);
        F = sparse(dof(nx:nx:end, 5:6)', 1, load, ndof, 1);                % Load in y and z directions
        cp = [5, 6];                                                       % sigma_13 and sigma_12

    case 'torsion-flexure'
        % Convert wall thickness to number of subvolumes
        t_y = max(1, ceil((B/2) / h));
        t_z = max(1, ceil((B/2) / b));

        % Hollow out the central region to form a C-shaped cross-section
        rho(:, t_y+1:end-t_y, t_z+1:end) = 1e-6;  % inner region set to minimum density

        % Identify subvolumes with material at the face x = L
        mask = squeeze(rho(nx, :, :) > 1e-5);                              % material presence at last x-slice
        [j, k] = find(mask);
        sv = sub2ind([nx, ny, nz], nx * ones(size(j)), j, k);              % linear indices of solid subvolumes

        F = sparse(ndof, 1);
        F(dof(sv, 5)) = P / numel(sv);                                     % Load in y-direction
        cp = [1, 6];                                                       % sigma_11 and sigma_12

    otherwise
        error('Unsupported problem type.');
end

%%

% FINITE-VOLUME THEORY ANALYSIS
[KL, ab, Ab, C] = LocalStiffMatrix(E, nu, l, h, b);                        % local stiffness matrix
sK = KL(:) * rho(:)';                                                      % stiffness interpolation
K = sparse(iK, jK, sK); K = (K + K')/2;                                    % assemblage of the stiffness matrix
U = Tikhonov(fdof, K, F);                                                  % compute global displacements employing Tikhonov strategy
[u, str] = NodalResponse(nx, ny, nz, l, h, b, ab, Ab, C, nodes, dof, U, rho);   % compute nodal displacements and stress fields

% PLOT DEFORMED STRUCTURE AND STRESS FIELDS
Plot3D(nx, ny, nz, l, h, b, nodes, u, amp, str, cp, rho);

%% _____________________________________________ORDENING DEGREES OF FREEDOM
function [nodes, dof, ndof, iK, jK] = DOFassembly(nx, ny, nz)

% Indices for subvolumes:
[i, j, k] = ndgrid(1:nx, 1:ny, 1:nz);

% Compute the base index for node positioning
base = (j - 1) * ((nx + 1) * (nz + 1)) + (k - 1) * (nx + 1);

% Compute the node indices efficiently in a vectorized manner
nodes = [i(:)   + base(:), ...                                             % Node 1
    i(:)+1 + base(:), ...                                                  % Node 2
    i(:)+1 + (j(:) * ((nx+1) * (nz+1))) + (k(:) - 1) * (nx+1), ...         % Node 3
    i(:)   + (j(:) * ((nx+1) * (nz+1))) + (k(:) - 1) * (nx+1), ...         % Node 4
    i(:)   + base(:) + (nx+1), ...                                         % Node 5
    i(:)+1 + base(:) + (nx+1), ...                                         % Node 6
    i(:)+1 + (j(:) * ((nx+1) * (nz+1))) + k(:) * (nx+1), ...               % Node 7
    i(:)   + (j(:) * ((nx+1) * (nz+1))) + k(:) * (nx+1)];                  % Node 8

% Number of horizontal faces:
Nhf = nx * nz * (ny + 1);

% Number of vertical faces in the x-direction:
Nvfx = (nx + 1) * nz * ny;

% Subvolume faces:
fc1 = Nhf + i(:) + (k(:) - 1) * (nx + 1) + (j(:) - 1) * nz * (nx + 1);     % Left lateral face
fc2 = Nhf + i(:) + 1 + (k(:) - 1) * (nx + 1) + (j(:) - 1) * nz * (nx + 1); % Right lateral face
fc3 = i(:) + (k(:) - 1) * nx + (j(:) - 1) * nx * nz;                       % Bottom face
fc4 = i(:) + (k(:) - 1) * nx + j(:) * nx * nz;                             % Top face
fc5 = Nhf + Nvfx + i(:) + (k(:) - 1) * nx + (j(:) - 1) * nx * (nz + 1);    % Back face
fc6 = Nhf + Nvfx + i(:) + k(:) * nx + (j(:) - 1) * nx * (nz + 1);          % Front face

% Final combination of faces:
faces = [fc1, fc2, fc3, fc4, fc5, fc6];

% Degrees of fereedom
dof = zeros(nx * ny * nz, 18);
dof(:, 1:3:16) = 3 * faces - 2;
dof(:, 2:3:17) = 3 * faces - 1;
dof(:, 3:3:18) = 3 * faces;

ndof = max(dof(:));                                                        % total number of degrees of freedom
iK = kron(dof, ones(18, 1))';                                              % iK - line mapping indice to the global stiffness location
jK = kron(dof, ones(1, 18))';                                              % jK - column mapping indice to the global stiffness location

%% ______________________________________________________KINEMATIC RELATION
function E = KinematicRelation(x, y, z)

% strain = E * Wi - Kinematic Relation
E = [1, 0, 0, 3*x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*y, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*z;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*z, 0, 1, 0, 0, 3*y, 0;
    0, 0, 1, 0, 0, 3*z, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3*x, 0, 0;
    0, 1, 0, 0, 3*y, 0, 1, 0, 0, 3*x, 0, 0, 0, 0, 0, 0, 0, 0];

%% ____________LOCAL STIFFNESS MATRIX FOR ZEROTH ORDER FINITE VOLUME THEORY
function [KL, ab, Ab, C0] = LocalStiffMatrix(E, nu, l, h, b)

% Constitutive matrix
C0 = constitutive_matrix(E, nu);

% Normal vectors
N = zeros(18,36);
N1 = [-1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -1; 0, 0, 0, 0, -1, 0];
N3 = [0, 0, 0, 0, 0, -1; 0, -1, 0, 0, 0, 0; 0, 0, 0, -1, 0, 0];
N5 = [0, 0, 0, 0, -1, 0; 0, 0, 0, -1, 0, 0; 0, 0, -1, 0, 0, 0];

N(1:3, 1:6)     = N1;
N(4:6, 7:12)    = -N1;
N(7:9, 13:18)   = N3;
N(10:12, 19:24) = -N3;
N(13:15, 25:30) = N5;
N(16:18, 31:36) = -N5;

% Kinematic matrices:
a = repmat(eye(3), 6, 1);
A = sparse(18, 18);

A(1,1)  = -l/2;    A(1,4)  = l^2/4;
A(2,7)  = A(1,1);  A(2,10) = A(1,4);
A(3,13) = A(1,1);  A(3,16) = A(1,4);
A(4,1)  = -A(1,1); A(4,4)  = A(1,4);
A(5,7)  = -A(1,1); A(5,10) = A(1,4);
A(6,13) = -A(1,1); A(6,16) = A(1,4);

A(7,2)   = -h/2;    A(7,5)   = h^2/4;
A(8,8)   = A(7,2);  A(8,11)  = A(7,5);
A(9,14)  = A(7,2);  A(9,17)  = A(7,5);
A(10,2)  = -A(7,2); A(10,5)  = A(7,5);
A(11,8)  = -A(7,2); A(11,11) = A(7,5);
A(12,14) = -A(7,2); A(12,17) = A(7,5);

A(13,3)  = -b/2;     A(13,6)  = b^2/4;
A(14,9)  = A(13,3);  A(14,12) = A(13,6);
A(15,15) = A(13,3);  A(15,18) = A(13,6);
A(16,3)  = -A(13,3); A(16,6)  = A(13,6);
A(17,9)  = -A(13,3); A(17,12) = A(13,6);
A(18,15) = -A(13,3); A(18,18) = A(13,6);

% Strain/Displacement operator per face
E1 = KinematicRelation(-l/2, 0, 0);
E2 = KinematicRelation(l/2, 0, 0);
E3 = KinematicRelation(0, -h/2, 0);
E4 = KinematicRelation(0, h/2, 0);
E5 = KinematicRelation(0, 0, -b/2);
E6 = KinematicRelation(0, 0, b/2);
E = [E1; E2; E3; E4; E5; E6];

% Constitutive matrix: 36 by 36 elements
C = zeros(36, 36);

% Use a loop to place C0 in the diagonal blocks
idx = 1:6:36;  % Indices for each block
for i = idx
    C(i:i+5, i:i+5) = C0;  % Assign C0 to each block
end

% Operator that relates the local surface-averaged tractions x Wi
B = N * C * E;

% Eficcienty matrix inverse
identI = eye(18, 18);
invA = A \ identI;

% Summation components of matrix (B * S) per face
sumB = (B(1:3, :) + B(4:6, :)) * b * h + ...
    (B(7:9, :) + B(10:12, :)) * l * b + ...
    (B(13:15, :) + B(16:18, :)) * l * h;

% Matrices a_barra and A_barra
ab = (sumB * invA * a) \ (sumB * invA);
Ab = invA * (identI - a * ab);

% Pseudo local stiffness matrix
KL = B * Ab;

% Create an 18x18 matrix of zeros for the subvolume faces area
S = zeros(18, 18);

% Fill the first 6x6 block with b * h on the diagonal
S(1:6, 1:6) = b * h * eye(6);

% Fill the next 6x6 block with l * b on the diagonal
S(7:12, 7:12) = l * b * eye(6);

% Fill the last 6x6 block with l * h on the diagonal
S(13:18, 13:18) = l * h * eye(6);

% Modified local stiffness matrix
KL = S * KL;

%% __________________________________________TIKHONOV REGULARIZATION SCHEME
function Unew = Tikhonov(fdof, K, F)

% Define tolerance
tol = 1e-6;

% Number of degrees of freedom
ndof = length(K);

% Initial regularization parameter
lambda0 = 10.^(-12);
lambda = lambda0 .* trace(K) / ndof;
max_iter = length(lambda);

% Initialize variables
U = zeros(ndof, 1); Unew = U;
error = 1; iter = 1;

while error > tol && iter <= max_iter

    % Tikhonov regularization
    Knew = K + lambda(iter) * speye(ndof);

    % Solve the system of equations
    Unew(fdof) = Knew(fdof, fdof) \ F(fdof);

    % Compute the relative error
    error = norm(K(fdof, fdof) * Unew(fdof) - F(fdof)) / norm(F(fdof));

    % Update the displacement vector
    U(fdof) = Unew(fdof);
end
fprintf('\n');

%% ____________________________COMPUTE NODAL DISPLACEMENTS AND STRESS FIELD
function [u, str] = NodalResponse(nx, ny, nz, l, h, b, ab, Ab, C, nodes, dof, U, rho)

% Number of nodes
nnodes = numel(unique(nodes(:)));

% Initialize displacement and stress fields
u = zeros(nnodes, 3);
str = zeros(nnodes, 6);

% Local coordinates for each node within an subvolume
x = [-1, 1, 1, -1, -1, 1, 1, -1] * 0.5 * l;
y = [-1, -1, 1, 1, -1, -1, 1, 1] * 0.5 * h;
z = [-1, -1, -1, -1, 1, 1, 1, 1] * 0.5 * b;

% Precompute constants for quadratic terms of the displacement field
cx = 0.5 * (3 * (x.^2) - l^2 / 4);
cy = 0.5 * (3 * (y.^2) - h^2 / 4);
cz = 0.5 * (3 * (z.^2) - b^2 / 4);

% Loop over subvolumes
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            % Subvolume index
            q = i + (j - 1) * nx + (k - 1) * nx * ny;

            % Coefficients of the displacement field
            W0 = ab * U(dof(q, :));                                        % Zeroth-order coefficients
            Wi = Ab * U(dof(q, :));                                        % Higher-order coefficients

            % Compute displacements
            u1 = W0(1) + x(:) * Wi(1) + y(:) * Wi(2) + z(:) * Wi(3) + ...
                cx(:) * Wi(4) + cy(:) * Wi(5) + cz(:) * Wi(6);

            u2 = W0(2) + x(:) * Wi(7) + y(:) * Wi(8) + z(:) * Wi(9) + ...
                cx(:) * Wi(10) + cy(:) * Wi(11) + cz(:) * Wi(12);

            u3 = W0(3) + x(:) * Wi(13) + y(:) * Wi(14) + z(:) * Wi(15) + ...
                cx(:) * Wi(16) + cy(:) * Wi(17) + cz(:) * Wi(18);

            u(nodes(q, :), 1) = u(nodes(q, :), 1) + u1(:);
            u(nodes(q, :), 2) = u(nodes(q, :), 2) + u2(:);
            u(nodes(q, :), 3) = u(nodes(q, :), 3) + u3(:);

            % Compute stress field
            for m = 1:8
                E = KinematicRelation(x(m), y(m), z(m));
                val = C * rho(q) * E * Wi;
                str(nodes(q, m), :) = str(nodes(q, m), :) + val';
            end
        end
    end
end

% Count occurrences of each node index
freq = accumarray(nodes(:), 1);

% Normalize each component using element-wise division
u = u ./ freq;
str = str ./ freq;

%% ________________________________PLOT DEFORMED STRUCTURE AND STRESS FIELD
function Plot3D(nx, ny, nz, l, h, b, nodes, u, amp, str, comp, rho)


[Vertices, Faces] = MapIndices(nx, ny, nz, l, h, b, rho);
[Deformed, stress] = StressAndDisplacement(nx, ny, nz, l, h, b, nodes, u, amp, str, rho);

% Swap columns [2 3] function
sw = @(V) V(:, [1 3 2]);

% Plot deformed structure
figure;
patch('Faces', Faces, 'Vertices', sw(Deformed), 'FaceColor', 'w', 'LineWidth', 1.1);

% Configuration of the axes:
scrsz = get(0,'ScreenSize');
set(gcf,'Outerposition',[1 1 0.4 * scrsz(3) 0.6 * scrsz(4)]);

set(gcf, 'Name', 'Deformed Structure - FVT', 'NumberTitle', 'off', 'Color', 'w');
axis equal; axis tight; axis off;
view([30, 20]);

% Plot stress field
labels = {'\sigma_{11}', '\sigma_{22}', '\sigma_{33}', '\sigma_{23}', ...
    '\sigma_{13}', '\sigma_{12}'};

for i = 1:numel(comp)

    figure;
    patch('Faces', Faces, 'Vertices', sw(Vertices), 'FaceVertexCData',...
        stress(:, comp(i)), 'FaceColor', 'interp', 'EdgeColor', 'none');

    % Configure figure appearance
    scrsz = get(0,'ScreenSize');
    set(gcf,'Outerposition',[1 1 0.4 * scrsz(3) 0.6 * scrsz(4)]);
    set(gcf, 'Name', 'Stress Field - FVT', 'NumberTitle', 'off', 'Color', 'w');
    colormap jet;
    cbar = colorbar;
    set(cbar, 'FontSize', 20, 'Location', 'eastoutside');
    title(cbar, labels{comp(i)}, 'FontSize', 18);
    axis equal; axis tight; axis off;
    view([-30, 30]);
end

% Mapping face's and vertice's indices
function [Vertices, Faces] = MapIndices(nx, ny, nz, l, h, b, rho)

% Total number of subvolumes
nsv = nx * ny * nz;

% Local face definitions (relative to the 8 nodes of an subvolume)
local_faces = [1,5,8,4; 2,6,7,3; 1,2,6,5; 4,3,7,8; 1,4,3,2; 5,8,7,6];

% Preallocate global matrices for efficiency
Vertices = zeros(nsv * 8, 3);
Faces    = zeros(nsv * 6, 4);

% Initialize vertex and face indices
idV = 0;
idF = 0;

% Iterate through the 3D grid to compute deformed geometry
for k = 1:nz
    for j = 1:ny
        for i = 1:nx

            % Compute subvolume vertices
            coord = GenerateSubvolumeCoordinates(i, j, k, l, h, b);

            if (rho(i, j, k) == 1)

                % Store results in global arrays:
                Vertices(idV + (1:8), :)  = coord;
                Faces(idF + (1:6), :) = local_faces + idV;

                % Update indices
                idV = idV + 8;
                idF = idF + 6;
            end
        end
    end
end

% Removing subvolumes with densities lower than the predefined threshold
Vertices = Vertices(1:idV, :);
Faces = Faces(1:idF, :);

function [Deformed, stress] = StressAndDisplacement(nx, ny, nz, l, h, b, nodes, u, amp, str, rho)

% Total number of subvolumes
nsv = nx * ny * nz;

% Preallocate global matrices for efficiency
Deformed = zeros(nsv * 8, 3);
stress   = zeros(nsv * 8, 6);

% Initialize vertex and face indices
idV = 0;
idF = 0;

% Iterate through the 3D grid to compute deformed geometry
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            % Compute global index of the current subvolume
            q = i + (j - 1) * nx + (k - 1) * nx * ny;

            % Global node indices
            idnode = nodes(q, :);

            % Compute subvolume vertices
            coord = GenerateSubvolumeCoordinates(i, j, k, l, h, b);

            % Apply displacement
            local_deformed = [coord(:,1) + amp * u(idnode, 1), ...
                coord(:,2) + amp * u(idnode, 2), ...
                coord(:,3) + amp * u(idnode, 3)];

            if (rho(i, j, k) == 1)

                % Store results in global arrays: deformed and stress
                Deformed(idV + (1:8), :)  = local_deformed;
                stress(idV + (1:8), :) = str(idnode, :);

                % Update indices
                idV = idV + 8;
                idF = idF + 6;
            end
        end
    end
end

% Removing subvolumes with densities lower than the predefined threshold
stress = stress(1:idV, :);
Deformed = Deformed(1:idV, :);

% Constitutive matrix for 3D linear elasticity (isotropic)
function C = constitutive_matrix(E, nu)
C = E / ((1+nu)*(1-2*nu)) * ...
    [1-nu,  nu,    nu,    0,       0,       0;
    nu,   1-nu,  nu,    0,       0,       0;
    nu,   nu,    1-nu,  0,       0,       0;
    0,    0,     0,     (1-2*nu)/2, 0,    0;
    0,    0,     0,     0,  (1-2*nu)/2,   0;
    0,    0,     0,     0,       0, (1-2*nu)/2];

function coord = GenerateSubvolumeCoordinates(i, j, k, l, h, b)
% Local vertex positions in reference cube
vx = [0 1 1 0 0 1 1 0];
vy = [0 0 1 1 0 0 1 1];
vz = [0 0 0 0 1 1 1 1];

% Compute coordinates in physical space
coord = [(i - 1 + vx') * l, ...
    (j - 1 + vy') * h, ...
    (k - 1 + vz') * b];

%% _____________________________________________________________________END