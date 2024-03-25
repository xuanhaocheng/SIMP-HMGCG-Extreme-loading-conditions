%function top3D125(nelx,nely,nelz,volfrac,penal,rmin,ft,ftBC,eta,beta,move,maxit)
clear;
nl=7; %7
nelx=25*2^(nl-1);
nely=4*2^(nl-1);
nelz=3*2^(nl-1);  %2.5
% volfrac=0.0;
volfrac0=0.04; %0.04
%penal=3;
penal=5;
rmin=sqrt(3)*2;
ft=3;
ftBC='N';
eta=0.4;
move=0.1;
maxit=50;
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
%-[ sI, sII ] = deal( [ ] );
%-for j = 1 : 24
%-    sI = cat( 2, sI, j : 24 );
%-    sII = cat( 2, sII, repmat( j, 1, 24 - j + 1 ) );
%-end
%-[ iK , jK ] = deal( cMat( :,  sI )', cMat( :, sII )' );
%Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' ); clear iK jK              % reduced assembly indexing
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
fixed1=3*nodeNrs( [1,nely+1], 1, [1,nelx+1])-2;
fixed2=3*nodeNrs( [1,nely+1], 1, [1,nelx+1])-1;
fixed3=3*nodeNrs( [1,nely+1], 1, [1,nelx+1]);
fixed = double(union(fixed3,union(fixed1,fixed2)));

z1 = round(nelz*0.1278 + 2);
z2 = round(nelz - nelz*0.1278);
y1 = round(nely*0.0548 + 1);
y2 = round(nely - nely*0.0548);
zz = round(nelz*0.1278 + 1);
pasS1 = elmNrs(y1:y2,zz,:);
pasS_double = double(pasS1);
total_rows = size(pasS_double, 1) * size(pasS_double, 2)* size(pasS_double, 3);
pasS = reshape(pasS_double, total_rows, 1);

pasV1=elmNrs(y1:y2,z1:z2,:);
pasV_double = double(pasV1);
total_rows1 = size(pasV_double, 1) *size(pasV_double, 2)* size(pasV_double, 3);
pasV = reshape(pasV_double, total_rows1, 1);

free = setdiff( 1 : nDof, fixed );                                         % set of free DOFs
act = setdiff( ( 1 : nEl )', union( pasS, pasV ) );                        % set of active d.v.
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
x( pasS ) = 0.001;                                                            % set x = 0.05 on pasS set
x( pasV ) = 1e-6;
F0=zeros(nDof,1);
F = zeros(nDof,1);

% Obtain the force on each element node based on the density, quality, and overall force situation of the elements
% aerocrete_density = 680.0;                                                  %  680kg/m3
concrete_density = 2400.0;                                              %  C30 2400 kg/m3
steel_density = 7850.0;                                                 %  7850.0 kg/m3
line_load1 = 1050;  % 1050 kg/m
line_load2 = 36000;  % 36000 kg/m
roll_length = round(0.6/ 14.6 * nely);
floor_y = nely - 2 * round(nely*0.0548);

one_element_concrete = concrete_density * power((91.25 / nelx ),3);
one_element_steel = steel_density * power((91.25 / nelx ),3);
dead_load = 10.0 * (1.2 * 0.25 / (91.25 / nelx) * one_element_concrete) /  one_element_steel;
live_load1_line = 10.0 * (1.4 * line_load1 * 91.25/ (2 * (roll_length + 1)* (nelx + 1))) /  one_element_steel;
live_load2_line = 10.0 * (1.4 * line_load2 / (2 * (floor_y + 1))) /  one_element_steel;
direction3 = 3;
g4 = 30.0;
g3 = 10.0;
line1_1_y1 = round(1.9 / 14.6 * nely + 1);    % The position of force in space
line1_1_y2 = line1_1_y1 + 1+ roll_length;
line1_2_y1 = round(3.7 / 14.6 * nely + 1);
line1_2_y2 = line1_2_y1 + 1+ roll_length;
line1_3_y1 = round(6.1 / 14.6 * nely + 1);
line1_3_y2 = line1_3_y1 + 1+ roll_length;
line1_4_y1 = round(7.9 / 14.6 * nely + 1);
line1_4_y2 = line1_4_y1 + 1+ roll_length;
line1_5_y1 = round(10.3 / 14.6 * nely + 1);
line1_5_y2 = line1_5_y1 + 1+ roll_length;
line1_6_y1 = round(12.1 / 14.6 * nely + 1);
line1_6_y2 = line1_6_y1 + 1+ roll_length;

middle_line = round(nelx / 2);

lcdof_line1 = 3 * nodeNrs([line1_1_y1:line1_1_y2,line1_2_y1:line1_2_y2,line1_3_y1:line1_3_y2,...
    line1_4_y1:line1_4_y2,line1_5_y1: line1_5_y2,line1_6_y1:line1_6_y2],zz+1 , 1:nelx+1);
lcdof_line2 = 3 * nodeNrs(y1: y2+1,zz+1 , middle_line:middle_line+1);

% lcdof_line1 = nodeNrs([line1_1_y1:line1_1_y2,line1_2_y1:line1_2_y2,line1_3_y1:line1_3_y2,...
%     line1_4_y1:line1_4_y2,line1_5_y1: line1_5_y2,line1_6_y1:line1_6_y2],zz+1 , 1:nelx+1);
% lcdof_line2 = nodeNrs(y1: y2+1,zz+1 , middle_line:middle_line+1);
% lcdof_line = union(lcdof_line1, lcdof_line2);

% [y_indices, z_indices, x_indices] = ind2sub([nely+1, nelz+1, nelx+1],reshape(lcdof_line,[],1)) ;
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

% direction2 = 2;
% g2 = 0.2 * 10.0 / sqrt(2.0);
% direction1 = 1;
% g1 = 0.2 * 10.0 / sqrt(2.0);
F0(lcdof_line1) = live_load1_line;
F0(lcdof_line2) = live_load2_line;
Assembly_F0_1(F0,pasS, cVec,dead_load, nelx, nely, nelz,direction3);
volfrac=(numel(act)*volfrac0+numel(pasS)*0.001+numel(pasV)*1e-6)/nEl;
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
  Assembly_F_1(F,act, xPhys,cVec,g3, nelx, nely, nelz,direction3); %Apply gravity load based on element density in each loop
  F(fixed)=0.0;
  CX=Emin+xPhys.^penal*(E0-Emin);
  nfd=length(fixed);
  [cgres,cgiters]=mgcg9(Ke0,F,U,CX,fixed',nfd,nelx,nely,nelz,nl,cgtol,cgmax,nswp, printLev, nu, hx, hy, hz);
  if cgiters<0, break; end
  fprintf('MGCG relres: %4.2e iters: %4i \n',cgres,cgiters);
  
  % ------------------------------------------ RL. 3) COMPUTE SENSITIVITIES
  % dc = dsK .* sum( ( U( cMat ) * Ke0 ) .* U( cMat ), 2 );                  % derivative of compliance
  dc = sum_UKU7(Ke0,dsK, U, cVec, nelx, nely, nelz);                        % derivative of compliance
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
isovals = shiftdim( reshape( xPhys, nely, nelz, nelx ), 2 );
isovals = smooth3( isovals, 'box', 1 );
patch(isosurface(isovals, .5),'FaceColor',[192, 192, 192]/255,'EdgeColor','none'); %'b'
patch(isocaps(isovals, .5),'FaceColor','r','EdgeColor','none');
drawnow; view( [ 145, 25 ] ); axis equal tight off; camlight;
savefig(strcat('QNC_V1_bridge_suxiang',num2str(loop),'.fig'));
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
