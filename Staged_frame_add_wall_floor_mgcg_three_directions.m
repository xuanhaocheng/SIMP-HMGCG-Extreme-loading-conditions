%function top3D125(nelx,nely,nelz,volfrac,penal,rmin,ft,ftBC,eta,beta,move,maxit)
clear;
nl=8;
nelx=3*2^(nl-1);
nely=3*2^(nl-1);
nelz=4*2^(nl-1);  %2.5
% volfrac=0.0;
volfrac0=0.12;
%penal=3;
penal=5;
rmin=sqrt(3)*2;
ft=3;
ftBC='N';
eta=0.4;
move=0.10;
maxit=30;
beta=2;

% ---------------------------- PRE. 1) MATERIAL AND CONTINUATION PARAMETERS
E0 = 1;                                                                    % Young modulus of solid
Emin = 1e-9;                                                               % Young modulus of "void"
nu = 0.3;                                                                  % Poisson ratio
penalCnt = { 1, 1, 25, 0.25 };                                             % continuation scheme on penal
betaCnt  = { 1, 1, 25,    2 };                                             % continuation scheme on beta
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end                      % filter BC selector
% ----------------------------------------- PRE. 2) DISCRETIZATION FEATURES
nEl = nelx * nely * nelz;                                                  % number of elements          #3D#
nodeNrs = int32( reshape( 1 : ( 1 + nelx ) * ( 1 + nely ) * ( 1 + nelz ), ...
    1 + nely, 1 + nelz, 1 + nelx ) );                                      % nodes numbering             #3D#
cVec = reshape( 3 * nodeNrs( 1 : nely, 1 : nelz, 1 : nelx ) + 1, nEl, 1 ); %                             #3D#
% cMat = cVec+int32( [0,1,2,3*(nely+1)*(nelz+1)+[0,1,2,-3,-2,-1],-3,-2,-1,3*(nely+...
%    1)+[0,1,2],3*(nely+1)*(nelz+2)+[0,1,2,-3,-2,-1],3*(nely+1)+[-3,-2,-1]]);% connectivity matrix         #3D#
nDof = ( 1 + nely ) * ( 1 + nelz ) * ( 1 + nelx ) * 3;                     % total number of DOFs        #3D#
Ke = 1/(1+nu)/(2*nu-1)/144 *( [ -32;-6;-6;8;6;6;10;6;3;-4;-6;-3;-4;-3;-6;10;...
    3;6;8;3;3;4;-3;-3; -32;-6;-6;-4;-3;6;10;3;6;8;6;-3;-4;-6;-3;4;-3;3;8;3;...
    3;10;6;-32;-6;-3;-4;-3;-3;4;-3;-6;-4;6;6;8;6;3;10;3;3;8;3;6;10;-32;6;6;...
    -4;6;3;10;-6;-3;10;-3;-6;-4;3;6;4;3;3;8;-3;-3;-32;-6;-6;8;6;-6;10;3;3;4;...
    -3;3;-4;-6;-3;10;6;-3;8;3;-32;3;-6;-4;3;-3;4;-6;3;10;-6;6;8;-3;6;10;-3;...
    3;8;-32;-6;6;8;6;-6;8;3;-3;4;-3;3;-4;-3;6;10;3;-6;-32;6;-6;-4;3;3;8;-3;...
    3;10;-6;-3;-4;6;-3;4;3;-32;6;3;-4;-3;-3;8;-3;-6;10;-6;-6;8;-6;-3;10;-32;...
    6;-6;4;3;-3;8;-3;3;10;-3;6;-4;3;-6;-32;6;-3;10;-6;-3;8;-3;3;4;3;3;-4;6;...
    -32;3;-6;10;3;-3;8;6;-3;10;6;-6;8;-32;-6;6;8;6;-6;10;6;-3;-4;-6;3;-32;6;...
    -6;-4;3;6;10;-3;6;8;-6;-32;6;3;-4;3;3;4;3;6;-4;-32;6;-6;-4;6;-3;10;-6;3;...
    -32;6;-6;8;-6;-6;10;-3;-32;-3;6;-4;-3;3;4;-32;-6;-6;8;6;6;-32;-6;-6;-4;...
    -3;-32;-6;-3;-4;-32;6;6;-32;-6;-32]+nu*[ 48;0;0;0;-24;-24;-12;0;-12;0;...
    24;0;0;0;24;-12;-12;0;-12;0;0;-12;12;12;48;0;24;0;0;0;-12;-12;-24;0;-24;...
    0;0;24;12;-12;12;0;-12;0;-12;-12;0;48;24;0;0;12;12;-12;0;24;0;-24;-24;0;...
    0;-12;-12;0;0;-12;-12;0;-12;48;0;0;0;-24;0;-12;0;12;-12;12;0;0;0;-24;...
    -12;-12;-12;-12;0;0;48;0;24;0;-24;0;-12;-12;-12;-12;12;0;0;24;12;-12;0;...
    0;-12;0;48;0;24;0;-12;12;-12;0;-12;-12;24;-24;0;12;0;-12;0;0;-12;48;0;0;...
    0;-24;24;-12;0;0;-12;12;-12;0;0;-24;-12;-12;0;48;0;24;0;0;0;-12;0;-12;...
    -12;0;0;0;-24;12;-12;-12;48;-24;0;0;0;0;-12;12;0;-12;24;24;0;0;12;-12;...
    48;0;0;-12;-12;12;-12;0;0;-12;12;0;0;0;24;48;0;12;-12;0;0;-12;0;-12;-12;...
    -12;0;0;-24;48;-12;0;-12;0;0;-12;0;12;-12;-24;24;0;48;0;0;0;-24;24;-12;...
    0;12;0;24;0;48;0;24;0;0;0;-12;12;-24;0;24;48;-24;0;0;-12;-12;-12;0;-24;...
    0;48;0;0;0;-24;0;-12;0;-12;48;0;24;0;24;0;-12;12;48;0;-24;0;12;-12;-12;...
    48;0;0;0;-24;-24;48;0;24;0;0;48;24;0;0;48;0;0;48;0;48 ] );             % elemental stiffness matrix  #3D#
Ke0( tril( ones( 24 ) ) == 1 ) = Ke';
Ke0 = reshape( Ke0, 24, 24 );
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );                                    % recover full matrix
% ----------------------------- PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS
[ pasS, pasV ] = deal( [], [] );                                            % passive solid and void elements
elmNrs = int32( reshape( 1 : nelx*nely*nelz ,nely, nelz, nelx));            % elements numbering             #3D#

width_column = round(0.04444* nelx);                                       %Calculate the position of holes and columns
x_column1 = round(0.35436*nelx);
y_column1 = round(0.4556*nely);
x_column2 = round(0.68898*nelx);
y_column2 = round(0.53763*nely);
x_stair = round(0.39716 * nelx);
y_stair = round(0.67162 * nelx);
wall_width = round(300 / 9000 * nelx);
z_beam_height = round(800 / 12000 * nelz);
z_board_height = round(300 / 12000 * nelz);
% x_stair_length =round(0.29315*nelx);
zz = round(1 / 3* nelz);
pasV=elmNrs(wall_width + 1:nely-wall_width,[1:zz-z_beam_height,zz+1 : 2*zz-z_beam_height, 2*zz + 1: nelz - z_beam_height ],wall_width + 1:nelx-wall_width);

pasSSx=elmNrs(y_column1 + 1:y_column1+width_column,1:nelz,x_column1 + 1:x_column1+width_column);
pasSSy=elmNrs(y_column2 + 1: y_column2 + width_column,1:nelz,x_column2 + 1: x_column2 + width_column);
pasSS = union(pasSSx,pasSSy);
pasSS1x = elmNrs(y_column1 + 2:y_column1+width_column-1,1:nelz,x_column1 + 2:x_column1+width_column-1);
pasSS1y = elmNrs(y_column2 + 2: y_column2 + width_column-1,1:nelz,x_column2 + 2: x_column2 + width_column-1);
pasSS1 = union(pasSS1x,pasSS1y);
intersection = intersect(pasSS, pasSS1);
pasSS = setdiff(pasSS, intersection);

intersection = intersect(pasSS, pasV);
pasV = setdiff(pasV, intersection);
pasS = double(reshape(pasSS, [], 1));
pasV = double(reshape(pasV, [], 1));

fixed1=3*nodeNrs( 1:nely+1, 1, 1:nelx+1)-2;
fixed2=3*nodeNrs( 1:nely+1, 1, 1:nelx+1)-1;
fixed3=3*nodeNrs( 1:nely+1, 1, 1:nelx+1);
fixed = double(union(fixed3,union(fixed1,fixed2)));

free = setdiff( 1 : nDof, fixed );                                         % set of free DOFs
act = setdiff( ( 1 : nEl )', union( pasS, pasV ) );                        % set of active d.v.
% [y_indices, z_indices, x_indices] = ind2sub([nely, nelz, nelx],pasS) ;
% 
% % Display the 3D coordinates of each element in pasV1
% pasV1_index = [x_indices, y_indices, z_indices];
% 
% % Use the scatter3 function to draw the corresponding units
% scatter3(pasV1_index(:, 1), pasV1_index(:, 2), pasV1_index(:, 3), 'filled');
% xlabel('i');
% ylabel('j');
% zlabel('k');
% title('Visualization of pasV1\_index');
% grid on;
% --------------------------------------- PRE. 4) DEFINE IMPLICIT FUNCTIONS
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));                                   % projection
deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .* ...
    sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );                % projection eta-derivative 
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative
cnt = @(v,vCnt,l) v+(l>=vCnt{1}).*(v<vCnt{2}).*(mod(l,vCnt{3})==0).*vCnt{4};
% -------------------------------------------------- PRE. 5) PREPARE FILTER
[dy,dz,dx]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,...
    -ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1 );
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 + dz.^2 ) );                        % conv. kernel                #3D#
Hs = imfilter( ones( nely, nelz, nelx ), h, bcF );                         % matrix of weights (filter)  #3D#
dHs = Hs;
% ------------------------ PRE. 6) ALLOCATE AND INITIALIZE OTHER PARAMETERS
[ x, dsK, dV ] = deal( zeros( nEl, 1 ) );                                  % initialize vectors
x( act ) = volfrac0;
x( pasS ) = 1.0;                                                            % set x = 0.05 on pasS set
x( pasV ) = 1e-6;
F0=zeros(nDof,1);
F = zeros(nDof,1);
volfrac_wall = 0.5;
volfrac_floor = 0.3;
% Obtain the force on each element node based on the density, quality, and overall force situation of the elements
aerocrete_density = 680.0;                                                  %  680kg/m3
concrete_density = 2400.0;                                              %  C30 2400 kg/m3
steel_density = 7850.0;                                                 %  7850.0 kg/m3

one_element_concrete = concrete_density * power((9.0 / nelx ),3);
one_element_aerocrete = aerocrete_density * power((9.0 / nelx ),3);
one_element_steel = steel_density * power((9.0 / nelx ),3);
force_wall = 121 * one_element_aerocrete * volfrac_wall;
force_floor = 16 * one_element_concrete * volfrac_floor; 
one_element_g_force_wall = force_wall / one_element_steel * 10.0;
one_element_g_force_floor = force_floor / one_element_steel * 10.0;
dead_load_live_load = 400 * 1.2+200 * 1.4;                              %  760 kg / m2
concrete_dead_load_live_load = dead_load_live_load * power((9.0 / nelx ),2);
force_g_dead_load_live_load = concrete_dead_load_live_load / one_element_steel * 10.0; % Convert the force into the load that each unit bears on the surface of the board

actual_add_one_element_g_force_wall = one_element_g_force_wall - force_g_dead_load_live_load;
actual_add_g_force_DL_floor = force_g_dead_load_live_load + one_element_g_force_floor;
x_stair_length =round(0.29315*nelx);

x_wall1_1 = round(x_column2 + width_column+1);
x_wall1_2 = round(nelx - width_column);
y_wall1_1 = round(y_column2 + width_column - wall_width+1);
y_wall1_2 = round(y_column2 + width_column);

x_wall2_1 = round(wall_width +1);
x_wall2_2 = round(x_column1);
y_wall2_1 = round(y_column1 + 1);
y_wall2_2 = round(y_column1 + wall_width);

x_wall3_1 = round(wall_width +1);
x_wall3_2 = x_stair;
y_wall3_1 = y_stair+1;
y_wall3_2 = y_stair+wall_width;

x_wall4_1 = round(0.21528 * nelx+ 1);
x_wall4_2 = round(x_wall4_1 + wall_width -1);
y_wall4_1 = round(y_column1 + wall_width + 1);
y_wall4_2 = round(y_stair);

x_wall5_1 = round(x_stair +1 - wall_width);
x_wall5_2 = round(x_stair);
y_wall5_1 = round(y_stair + 1 + wall_width);
y_wall5_2 = round(nely - wall_width);

x_wall6_1 = round(x_column2 +1);
x_wall6_2 = round(x_column2 + wall_width);
y_wall6_1 = round(wall_width + 1);
y_wall6_2 = round(0.23349*nely);

x_wall7_1 = round( x_column2+1);
x_wall7_2 = round(x_column2 + wall_width);
y_wall7_1 = round(y_column2 + 1 + width_column);
y_wall7_2 = round(nely -  wall_width);

act_wall1 = nodeNrs(y_wall1_1:y_wall1_2+1,[zz+1, zz *2+1,nelz+1],x_wall1_1:x_wall1_2+1);
act_wall2 = nodeNrs(y_wall2_1:y_wall2_2+1,[zz+1, zz *2+1,nelz+1],x_wall2_1:x_wall2_2+1);
act_wall3 = nodeNrs(y_wall3_1:y_wall3_2+1,[zz+1, zz *2+1,nelz+1],x_wall3_1:x_wall3_2+1);
act_wall4 = nodeNrs(y_wall4_1:y_wall4_2+1,[zz+1, zz *2+1,nelz+1],x_wall4_1:x_wall4_2+1);
act_wall5 = nodeNrs(y_wall5_1:y_wall5_2+1,[zz+1, zz *2+1,nelz+1],x_wall5_1:x_wall5_2+1);
act_wall6 = nodeNrs(y_wall6_1:y_wall6_2+1,[zz+1, zz *2+1,nelz+1],x_wall6_1:x_wall6_2+1);
act_wall7 = nodeNrs(y_wall7_1:y_wall7_2+1,[zz+1, zz *2+1,nelz+1],x_wall7_1:x_wall7_2+1);
act_wall = union(act_wall7,union(act_wall6,union(act_wall5,union(act_wall4,union(act_wall3,union(act_wall1,act_wall2))))));
act_wall = double(reshape(act_wall,[],1));

lcDof11 = 3 * nodeNrs( 1 +wall_width: nely + 1 -wall_width , nelz+1, 1+wall_width:nelx + 1 -wall_width); % -2 is horizontal, and not subtracting 2 is vertical.
lcDof12 = 3 * nodeNrs( 1 +wall_width: nely + 1-wall_width, zz *2+1, 1+wall_width:nelx + 1 -wall_width);
lcDof13 = 3 * nodeNrs( 1 +wall_width: nely + 1-wall_width, zz+1, 1+wall_width:nelx + 1 -wall_width);
lcDof14 = 3 * act_wall;

% [y_indices, z_indices, x_indices] = ind2sub([nely+1, nelz+1, nelx+1],act_wall) ;
% 
% % Display the 3D coordinates of each element in pasV1
% pasV1_index = [x_indices, y_indices, z_indices];
% 
% % Use the scatter3 function to draw the corresponding units
% scatter3(pasV1_index(:, 1), pasV1_index(:, 2), pasV1_index(:, 3), 'filled');
% xlabel('i');
% ylabel('j');
% zlabel('k');
% title('Visualization of pasV1\_index');
% grid on;

% intersection = intersect(act_wall, pasV);
% act_wall = setdiff(act_wall,intersection);

% lcDof = nodeNrs()
direction3 = 3;
g3 = 10.0;
direction2 = 2;
g2 = 0.2 * 10.0 / sqrt(2.0);
direction1 = 1;
g1 = 0.2 * 10.0 / sqrt(2.0);
%Assembly_F0_ Before loop
Assembly_F0_1(F0,pasS, cVec,g3, nelx, nely, nelz,direction3);
Assembly_F0_1(F0,pasS, cVec,g2, nelx, nely, nelz,direction2);
Assembly_F0_1(F0,pasS, cVec,g1, nelx, nely, nelz,direction1);
F0(lcDof11) = actual_add_g_force_DL_floor;
F0(lcDof12) = actual_add_g_force_DL_floor;
F0(lcDof13) = actual_add_g_force_DL_floor;
F0(lcDof14) = actual_add_one_element_g_force_wall;

volfrac=(numel(act)*volfrac0+numel(pasS)*1.0+numel(pasV)*1e-6)/nEl;
dV( act, 1 ) = 1/nEl/volfrac; 
[ xPhys, xOld, ch, loop, U ] = deal( x, 1, 1, 0, zeros( nDof, 1 ) );       % old x, x change, it. counter, U
% ================================================= START OPTIMIZATION LOOP
nu=0.3; nswp=16; printLev=1; hx=1; hy=1; hz=1; cgtol=1e-5; cgmax=400;
while ch > 1e-8 && loop < maxit
  tic
  loop = loop + 1;                                                         % update iter. counter
  % ----------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD (AND ETA IF PROJECT.)
  xTilde = imfilter( reshape( x, nely, nelz, nelx ), h, bcF ) ./ Hs;       % filtered field              #3D#
  xPhys( act ) = xTilde( act );                                            % reshape to column vector
  if ft > 1                                                                % compute optimal eta* with Newton
      f = ( mean( prj( xPhys, eta, beta ) ) - volfrac )  * (ft == 3);      % function (volume)
      while abs( f ) > 1e-6           % Newton process for finding opt. eta
          eta = eta - f / mean( deta( xPhys, eta, beta ) );
          f = mean( prj( xPhys, eta, beta ) ) - volfrac;
      end
      dHs = Hs ./ reshape( dprj( xPhys, eta, beta ), nely, nelz, nelx );   % sensitivity modification    #3D#
      xPhys = prj( xPhys, eta, beta );                                     % projected (physical) field
  end
  ch = norm( xPhys - xOld ) ./ nEl;
  xOld = xPhys;
  % -------------------------- RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS
  dsK( act ) = -penal * ( E0 - Emin ) * xPhys( act ) .^ ( penal - 1 ); 
  F(:) = F0(:);
  Assembly_F_1(F,act, xPhys,cVec,g3, nelx, nely, nelz,direction3);
  Assembly_F_1(F,act, xPhys,cVec,g2, nelx, nely, nelz,direction2);
  Assembly_F_1(F,act, xPhys,cVec,g1, nelx, nely, nelz,direction1);
  F(fixed)=0.0;
  CX=Emin+xPhys.^penal*(E0-Emin);
  nfd=length(fixed);
  [cgres,cgiters]=mgcg9(Ke0,F,U,CX,fixed',nfd,nelx,nely,nelz,nl,cgtol,cgmax,nswp, printLev, nu, hx, hy, hz);
  if cgiters<0, break; end
  fprintf('MGCG relres: %4.2e iters: %4i \n',cgres,cgiters);
  % ------------------------------------------ RL. 3) COMPUTE SENSITIVITIES
  % dc = dsK .* sum( ( U( cMat ) * Ke0 ) .* U( cMat ), 2 );                  % derivative of compliance
  dc = sum_UKU7(Ke0,dsK, U, cVec, nelx, nely, nelz);                       % derivative of compliance
  dc = imfilter( reshape( dc, nely, nelz, nelx ) ./ dHs, h, bcF );         % filter objective sens.      #3D#
  dV0 = imfilter( reshape( dV, nely, nelz, nelx ) ./ dHs, h, bcF );        % filter compliance sens.     #3D#
  % ----------------- RL. 4) UPDATE DESIGN VARIABLES AND APPLY CONTINUATION
  xT = x( act );
  [ xU, xL ] = deal( xT + move, xT - move );                               % current upper and lower bound
  ocP = xT .* sqrt( - dc( act ) ./ dV0( act ) );                           % constant part in resizing rule
  l = [ 0, mean( ocP ) / volfrac ];                                        % initial estimate for LM
  while ( l( 2 ) - l( 1 ) ) / ( l( 2 ) + l( 1 ) ) > 1e-4                   % OC resizing rule
      lmid = 0.5 * ( l( 1 ) + l( 2 ) );
      x( act ) = max( max( min( min( ocP / lmid, xU ), 1 ), xL ), 0 );
      if mean( x ) > volfrac, l( 1 ) = lmid; else, l( 2 ) = lmid; end
  end
  [penal,beta] = deal(cnt(penal,penalCnt,loop), cnt(beta,betaCnt,loop));   % apply conitnuation on parameters
  % -------------------------- RL. 5) PRINT CURRENT RESULTS AND PLOT DESIGN
  fprintf( 'It.:%5i C:%6.5e V:%7.3f ch.:%0.2e penal:%7.2f beta:%7.1f eta:%7.2f lm:%0.2e \n', ...
      loop, F'*U, mean(xPhys(:)), ch, penal, beta, eta, lmid );
  toc
end
% Specify the file name to save (file path can be customized)
filename = 'xphys_frame_data.mat';
% Use the save function to save variables to a file
save(filename, 'xPhys','-v7.3');
isovals = shiftdim( reshape( xPhys, nely, nelz, nelx ), 2 );
isovals = smooth3( isovals, 'box', 1 );
patch(isosurface(isovals, .5),'FaceColor',[192, 192, 192]/255,'EdgeColor','none'); %'b'
patch(isocaps(isovals, .5),'FaceColor','r','EdgeColor','none');
drawnow; view( [ 145, 25 ] ); axis equal tight off; camlight;
savefig(strcat('QNC_V1_staged_frame',num2str(loop),'.fig'));
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by F. Ferrari, O. Sigmund                   %
% Dept. of Solid Mechanics-Technical University of Denmark,2800 Lyngby (DK)%
% Please send your comments to: feferr@mek.dtu.dk                          %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper Ferrari, F. Sigmund, O. - A new generation 99 %
% line Matlab code for compliance Topology Optimization and its extension  %
% to 3D, SMO, 2020                                                         %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
