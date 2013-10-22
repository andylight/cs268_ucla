function hw2b()
%HW2B
%   Recall: 
%   spirit1983.png := 1
%   spirit1706.png := 2
%   spirit1433.png := 3

%halp();
DO_MASKS = false;
SHOW_CORRS = true;

focal = 14.67 * 1e-3;       % 14.67 mm
mx = sqrt(12) * 1e-6;       % sqrt(12) micrometers per pixel
%mx = 11.87 * 1e-6;
my = mx;
u0 = 843 / 2;
v0 = 843 / 2;
%u0 = 1024 / 2;
%v0 = 1024 / 2;
K = [[focal/mx 0 u0];
     [0 focal/my v0];
     [0 0 1]];
Kinv = inv(K);

[imgpaths, imgs, ptspaths, pts] = load_data();

newpts = cell([1, length(pts)]);
for i=1:length(pts)
    A = pts{i};
    Acalib = zeros(size(A));
    for idx=1:size(A, 1)
        curpt = [A(idx,:) 1];
        pt_calibrated = Kinv*curpt';
        Acalib(idx, :) = pt_calibrated(1:2)';
    end
    newpts{i} = Acalib;
end
pts = newpts;
    
% G_21 := 3x3 affine matrix mapping spirit1983.png -> spirit1706.png
% G_32 := 3x3 affine matrix mapping spirit1706.png -> spirit1433.png
% PtMats := {pts_12, pts_23}, where each pts_ij is a cell array:
%               {pts_i, pts_j}    pts_i := nx2 matrix of points
%           such that each point k in pts_i,pts_j correspond.
[G_21, G_32, PtMats] = hw2b_getTs('K', K);
if SHOW_CORRS
    [~,~,PtMats_uncalib] = hw2b_getTs();
    figure;
    idx1983 = indexof(imgpaths, 'spirit1983');
    idx1706 = indexof(imgpaths, 'spirit1706');
    idx1433 = indexof(imgpaths, 'spirit1433');
    showMatchedFeatures(imgs{idx1983}, imgs{idx1706}, PtMats_uncalib{1}{1}, PtMats_uncalib{1}{2}, 'montage');
    title(sprintf('Image %s to: %s', imgpaths{idx1983}, imgpaths{idx1706}));
    legend('imgs1', 'imgs2');
    figure;
    showMatchedFeatures(imgs{idx1706}, imgs{idx1433}, PtMats_uncalib{2}{1}, PtMats_uncalib{2}{2}, 'montage');
    title(sprintf('Image %s to: %s', imgpaths{idx1706}, imgpaths{idx1433}));
    legend('imgs2', 'imgs3');
end

if DO_MASKS
    disp('==== Doing MaSKS homography2Motion ====');
    p1 = pts2maskspts(PtMats{1}{1}); % 3xN, spirit1983
    p2 = pts2maskspts(PtMats{1}{2}); %      spirit1706
    [Sols2, H2] = homography2Motion(p1,p2);
    
    p3 = pts2maskspts(PtMats{2}{1}); %      spirit1706
    p4 = pts2maskspts(PtMats{2}{2}); %      spirit1433
    [Sols3, H3] = homography2Motion(p3,p4);
    disp('Done computing homographies.');
    sols2 = filtersols(Sols2);
    sols3 = filtersols(Sols3);
else
    p1 = pts2maskspts(PtMats{1}{1}); % 3xN, spirit1983
    p2 = pts2maskspts(PtMats{1}{2}); %      spirit1706   
    p3 = pts2maskspts(PtMats{2}{1}); %      spirit1706
    p4 = pts2maskspts(PtMats{2}{2}); %      spirit1433
    G_21 = compute_H_special(p1, p2);
    G_32 = compute_H_special(p3,p4);
    
    H2_norm = fix_sign(normalize_H(G_21), PtMats{1}{1}(1,:), PtMats{1}{2}(1,:));
    H3_norm = fix_sign(normalize_H(G_32), PtMats{2}{1}(1,:), PtMats{2}{2}(1,:));
    sols2 = decompose_H(H2_norm);
    sols3 = decompose_H(H3_norm);
end
disp('==== spirit1983 -> spirit1706 ====');
for i=1:length(sols2)
    sol = sols2{i};
    R = sol{1};     N = sol{2};     T = sol{3};
    H = (R + T*N');
    H = H / H(3,3);
    pt1 = H*[PtMats{1}{1}(1,:) 1]';
    pt1 = pt1 / pt1(3);
    pt1 = pt1(1:2)';
    pt2 = PtMats{1}{2}(1,:);
    err = norm(pt1-pt2, 2);
    disp(sprintf('err: %f', err));
    T_norm = T / norm(T,2);
    %[~, Rvert] = align_vector([0 0 -1], N', true);
    %Tworld = Rvert * T_norm;
    %motion_2 = Tworld / (Tworld(3)/-277);
    %myZ = R*[0 0 -1]';
    %[myZ_align, ~] = align_vector(myZ', N');
    %[Tworld, ~] = align_vector(T_norm', myZ_align');
    %Tworld = Rvert * T_norm;
    %motion_2 = Tworld / (Tworld(3)/-277);
    motion_2 = T_norm / (T_norm(3)/-277);
    v_2 = norm(motion_2 / 3.75, 2);
    disp(sprintf('    motion: %f %f %f', motion_2));
    disp(sprintf('    velocity: %f', v_2));
    % Plot plane
    %[xx,yy] = ndgrid(-10:10, -10:10);
    %z = (-N(1)*xx - N(2)*yy) / N(3);
    %figure;
    %surf(xx,yy,z); hold on;
    %draw_frame_scaled([R T; 0 0 0 1], 4); hold on;
end
disp('==== spirit1706 -> spirit1433 ====');
for i=1:length(sols3)
    sol = sols3{i};
    R = sol{1};     N = sol{2};     T = sol{3};
    H = (R + T*N');
    pt1 = H*[PtMats{2}{1}(1,:) 1]';
    pt1 = pt1 / pt1(3);
    pt1 = pt1(1:2)';
    pt2 = PtMats{2}{2}(1,:);
    err = norm(pt1-pt2, 2);
    disp(sprintf('err: %f', err));    
    T_norm = T / norm(T,2);
    motion_2 = T_norm / (T_norm(3)/-273);
    v_2 = norm(motion_2 / 3.75, 2);
    disp(sprintf('    motion: %f %f %f', motion_2));
    disp(sprintf('    velocity: %f', v_2));
    % Plot plane
    %[xx,yy] = ndgrid(-10:10, -10:10);
    %z = (-N(1)*xx - N(2)*yy) / N(3);
    %figure;
    %surf(xx,yy,z); hold on;
    %draw_frame_scaled([R T; 0 0 0 1], 4); hold on;
end

T_2 = sols2{1}{3};
T_3 = sols3{1}{3};

%%%%%%%% secs      meters
data = {{-25        1983},
        {-21.25     1706},
        {-17.5      1433}};

C_1 = [0 0 1983];
C_2 = C_1 + (T_2 / (T_2(3)/-277))';
C_3 = C_2 + (T_3 / (T_3(3)/-273))';

disp(sprintf('C_2 is: %s', str_vec(C_2)));
disp('    Last entry should be: 1706 meters');
disp(sprintf('C_3 is: %s', str_vec(C_3)));
disp('    Last entry should be: 1433 meters');

v_12 = (C_2 - C_1) / (data{2}{1} - data{1}{1});
v_23 = (C_3 - C_2) / (data{3}{1} - data{2}{1});

v_1 = norm(v_12, 2);
v_2 = norm(v_23, 2);

disp(sprintf('v_1 := %.4f m/s v_2 := %.4f m/s', v_1, v_2));
disp(sprintf('    v_avg: %.4f m/s', (v_1+v_2)/2));

end

function H_norm = normalize_H(H)
% Note: The sign of H_norm is ambiguous. It is up to the caller to
% disambiguate the sign, by enforcing that for two corresonding points on
% two different images:
%       (x_i)' * H_norm * (y_i) > 0
[U,S,V] = svd(H);
sigma_2 = S(2,2);   % second-biggest singular value
H_norm = H / sigma_2;
end

function H_out = fix_sign(H_norm, x_i, y_i)
%FIX_SIGN Return +H_norm or -H_norm depending on the constraint:
%   (x_i) * H_norm * (y_i)' > 0
x_i = [x_i 1];
y_i = [y_i 1];
if x_i * H_norm * y_i' > 0
    H_out = H_norm;
else
    H_out = -H_norm;
end
end

function sols = decompose_H(H)
%DECOMPOSE_H H must be normalized and sign-adjusted, i.e. normalize_H and
%fix_sign MUST be called prior to calling decompose_H.
[U, S, V] = svd(H'*H);
detU = det(U);
if (detU <= (-1 + 1e-4)) && (detU >= (-1 - 1e-4))
    % If det(U) == -1, replace both U,V with -U,-V
    U = -U;
    V = -V;
end
sig1 = S(1,1);  sig2 = S(2,2);  sig3 = S(3,3);
u1 = (sqrt(1 - sig3^2)*U(:,1) + sqrt(sig1^2 - 1)*U(:,3)) / sqrt(sig1^2 - sig3^2);
u2 = (sqrt(1 - sig3^2)*U(:,1) - sqrt(sig1^2 - 1)*U(:,3)) / sqrt(sig1^2 - sig3^2);
U1 = [U(:,2) u1 cross(U(:,2), u1)];
U2 = [U(:,2) u2 cross(U(:,2), u2)];
W1 = [H*U(:,2) H*u1 cross(H*U(:,2), H*u1)];
W2 = [H*U(:,2) H*u2 cross(H*U(:,2), H*u2)];
% Compute 4 possible solutions
R1 = W1*U1';     N1 = cross(U(:,2), u1);     T1 = (H - R1)*N1;
R2 = W2*U2';     N2 = cross(U(:,2), u2);     T2 = (H - R2)*N2;
R3 = R1;         N3 = -N1;                   T3 = -T1;
R4 = R2;         N4 = -N2;                   T4 = -T2;
possible_sols = {{R1 N1 T1} {R2 N2 T2} {R3 N3 T3} {R4 N4 T4}};
sols = {};      solidxs = {};
% Physical Constraint: points must lie in front of camera
%       n3 > 0
for i=1:length(possible_sols)
    cursol = possible_sols{i};
    R = cursol{1};  N = cursol{2};  T = cursol{3};
    if N(3) > 0
        sols{length(sols)+1} = {R, N, T};
        solidxs{length(solidxs)+1} = i;
    end
end
if length(sols) > 1
    %disp(sprintf('Found %d feasible solutions', length(sols)));
end
end

function halp()
%HALP Why are there two feasible sols?
R = [[cos(pi/10) 0 sin(pi/10)];
     [0 1 0];
     [-sin(pi/10) 0 cos(pi/10)]];
T = [2 0 0];
N = [1 0 2];
d = 5; lambda = 4;

HL = lambda * (R + (1/d)*T'*N);
disp('HL is:');
HL

HL_norm = normalize_H(HL);

sols = decompose_H(HL_norm);

R1 = sols{1}{1};    T1 = sols{1}{2};    N1 = sols{1}{3};
R2 = sols{2}{1};    T2 = sols{2}{2};    N2 = sols{2}{3};
end

function p = pts2maskspts(pts)
pts = [pts ones([size(pts,1),1])];
pts = pts';
p = pts;
end

function sols = filtersols(Sols)
%FILTERSOLS: Filters out physically impossible g=(R,T),N. For MaSKS.
%   Sols := 3x5x4 matrix, where:
%       Sols(:, 1:3, i) := Rotation matrix of i-th solution
%       Sols(:, 4, i)   := Translation of i-th solution
%       Sols(:, 5, i)   := normal vector N of i-th solution
sols = {};
solidxs = {};
for i=1:size(Sols, 3)
    R = Sols(:, 1:3, i);
    T = Sols(:, 4, i);
    N = Sols(:, 5, i);
    if N(3) > 0
        sols{length(sols)+1} = {R, N, T};
        solidxs{length(solidxs)+1} = i;
    end
end
end
