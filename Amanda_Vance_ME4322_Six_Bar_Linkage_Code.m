% six bar linkage 
% Static Equilibrium 

clc;
clear;

% define the joints 
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43 32 0];
G = [45 17 0];

% Define length of bars 
AB = norm(B - A);
BC = norm(C - B);
CD = norm(D - C);
DE = norm(E - D);
EF = norm(F - E);
FG = norm(G - F);

% Weight of Each Link 
WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];

%Center of mass of each link 
S1 = (A + B) / 2;
S2 = (B + C) / 2;
S3 = (C + D) / 2;
S4 = (D + E) / 2;
S5 = (F + G) / 2;

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin 


ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [0 0 Tin];

% Applied Force 
AppliedForce = [50 0 0];

% Static Equilibrium Conditions for Link AB 

% Sum of forces = 0
%Fa + Fb + WeightofAB = 0
eqn1 = ForceA + ForceB + WAB == 0;

% Sum of moments = 0 with respect to CoM of Link AB 
% S1A x FA + S1B x FB + InputTorque = 0

eqn2 = cross(A-S1,ForceA) + cross(B-S1, ForceB) + InputTorque == 0;

% Equations for Link BEC
%Sum of Forces = 0
% -FB + FC +FE +WBEC = 0

eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;

% Sum of Moments = 0
% Sum of Moments = 0 with respect to CoM of LInk BEC 

eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE) == 0;

% Equations for Link CD 
% Sum of Forces = 0 for Link CD
% -fc +fd +WCD = 0
eqn5 = -ForceC + ForceD + WCD == 0;

eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;

% Equations for Link EF 
eqn7 = -ForceE +ForceF + WEF == 0;

eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;

eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

% Sum of Moments = 0 with respect to CoM of Link FG

eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) == 0;

% Solving the 10 equations 

eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

% Extracting the solution for each force and torque
Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);
Force_Cx = double(StaticSolution.FCx);
Force_Cy = double(StaticSolution.FCy);
Force_Dx = double(StaticSolution.FDx);
Force_Dy = double(StaticSolution.FDy);
Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);
Force_Fx = double(StaticSolution.FFx);
Force_Fy = double(StaticSolution.FFy);
Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);
Input_Torque = double(StaticSolution.Tin);

disp('ForceA');
disp([Force_Ax, Force_Ay]);
disp('Force_B');
disp([Force_Bx, Force_By]);
disp('Force_C');
disp([Force_Cx, Force_Cy]);
disp('Force_D');
disp([Force_Dx, Force_Dy]);
disp('Force_E');
disp([Force_Ex, Force_Ey]);
disp('Force_F');
disp([Force_Fx, Force_Fy]);
disp('Force_G');
disp([Force_Gx, Force_Gy]);
disp('Input Torque');
disp(Input_Torque);




% Angular Velocity Calculations 
% Loop ABCDA

syms wBEC wCD
omega_AB = [0 0 1];
omega_BEC = [0 0 wBEC];
omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB,B-A) + cross(omega_BEC,C-B) + cross(omega_CD,D-C) == 0;

loop1Solution = solve(eqn11,[wBEC wCD]); 

% Extract angular velocities from the loop solution 
angularVelocity_BEC = double(loop1Solution.wBEC)
angularVelocity_CD = double(loop1Solution.wCD)

omegaBEC = [0 0 angularVelocity_BEC]
omegaCD = [0 0 angularVelocity_CD]

% Loop DCEFGD

syms wEF wFG
omega_EF = [0 0 wEF];
omega_FG = [0 0 wFG];
eqn12 = cross(omega_CD, C-D) + cross(omega_BEC,E-C) + cross(omega_EF, F-E) + cross(omega_FG, G-F) == 0;


loop2Solution = solve(eqn12, [wEF wFG]);
loop2Solution_sub = subs(loop2Solution, [wBEC, wCD], [angularVelocity_BEC, angularVelocity_CD]);
angularVelocity_EF = double(loop2Solution_sub.wEF)
angularVelocity_FG = double(loop2Solution_sub.wFG)


% Angular Acceleration 
% Loop 1 ABCDA
% Angular Acceleration for Loop 1 ABCDA
syms aBEC aCD
alpha_AB = [ 0 0 0];
alpha_BEC = [0 0 aBEC];
alpha_CD = [0 0 aCD];
% Angular Acceleration for Loop 1 ABCDA
a_B_A = cross(alpha_AB, B-A) + cross(omega_AB,cross(omega_AB,B-A));
a_C_B = cross(alpha_BEC, C-B) + cross(omega_BEC,cross(omega_BEC,C-B));
a_D_C = cross(alpha_CD, D-C) + cross(omega_CD,cross(omega_CD,D-C));

eqn13 = a_B_A + a_C_B + a_D_C == 0;

loop1AccSolution = solve(eqn13,[aBEC aCD]);

loop1AccSolution_sub = subs(loop1AccSolution, [wBEC, wCD], [angularVelocity_BEC, angularVelocity_CD]);
alphaBEC = double(loop1AccSolution_sub.aBEC)
alphaCD  = double(loop1AccSolution_sub.aCD)

alphaBEC_vector = [0 0 alphaBEC];
alphaCD_vector = [0 0 alphaCD];

% Loop 2 DCEFGD

syms aEF aFG
% Angular Acceleration for Loop 2 DCEFGD
alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

% a_C_D + a_E_C + a_F_E + a_G_F = 0 
a_C_D = cross(alphaCD_vector,C-D) + cross(omegaCD,cross(omegaCD,C-D));
a_E_C = cross(alphaBEC_vector,E-C) + cross(omegaBEC,cross(omegaBEC,E-C));

angVel_EF = [0 0 angularVelocity_EF];
angVel_FG = [0 0 angularVelocity_FG];

a_F_E = cross(alpha_EF,F-E) + cross(angVel_EF, cross(angVel_EF, F-E));

a_G_F = cross(alpha_FG, G-F) + cross(angVel_FG,cross(angVel_FG,G-F));

eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;

loop2AccSolution = solve(eqn14,[aEF aFG]);

alphaEF = double(loop2AccSolution.aEF)
alphaFG = double(loop2AccSolution.aFG)
 

% S1 = center of mass of link AB
S1 = (A + B) / 2;

% S2 = center of mass of rigid triangular link BEC
S2 = (B + C + E) / 3;

% S3 = center of mass of link CD
S3 = (C + D) / 2;

% S4 = center of mass of link EF
S4 = (E + F) / 2;

% S5 = center of mass of link FG
S5 = (F + G) / 2;

% VELOCITIES AT ALL JOINTS

% Fixed joints
vA = [0 0 0];
vD = [0 0 0];
vG = [0 0 0];

% Velocity of Joint B


vB_A = cross(omega_AB,B-A);

vB = vA + vB_A;

% Velocity of Joint C

% Calculate C from B using rigid body BEC
vC_B = vB + cross(omegaBEC,C-B);

% Calculate C from D using rigid body CD
vC_D = vD + cross(omegaCD,C-D);

% Use B calculation as the velocity of C
vC = vC_B;

% Velocity of Joint E

% Calculate E from B using rigid body BEC
vE_B = vB + cross(omegaBEC,E-B);

% Calculate E from C as a check
vE_C = vC + cross(omegaBEC,E-C);

% Use B calculation as the velocity of E
vE = vE_B;

% Velocity of Joint F

% Calculate F from E using link EF
vF_E = vE + cross(angVel_EF,F-E);

% Calculate F from G using link FG
vF_G = vG + cross(angVel_FG,F-G);

% Use E calculation as the velocity of F
vF = vF_E;


% Display Joint Velocities

disp('Velocity of A = ');
disp(vA);

disp('Velocity of B = ');
disp(vB);

disp('Velocity of C = ');
disp(vC);

disp('Velocity of D = ');
disp(vD);

disp('Velocity of E = ');
disp(vE);

disp('Velocity of F = ');
disp(vF);

disp('Velocity of G = ');
disp(vG);

% Center of Mass S1 - Link AB

vS1 = vA + cross(omega_AB,S1-A);

% Center of Mass S2 - Link BEC


vS2 = vB + cross(omegaBEC,S2-B);


% Center of Mass S3 - Link CD

vS3 = vD + cross(omegaCD,S3-D);

% Center of Mass S4 - Link EF


% Calculate from E
vS4_E = vE + cross(angVel_EF,S4-E);

% Calculate from F as a check
vS4_F = vF + cross(angVel_EF,S4-F);

% Use E calculation
vS4 = vS4_E;

% Center of Mass S5 - Link FG
vS5 = vG + cross(angVel_FG,S5-G);


% Display COM Velocities

disp('Velocity of S1 = ');
disp(vS1);

disp('Velocity of S2 = ');
disp(vS2);

disp('Velocity of S3 = ');
disp(vS3);

disp('Velocity of S4 = ');
disp(vS4);

disp('Velocity of S5 = ');
disp(vS5);



% ANGULAR ACCELERATION VECTORS


% Your code already calculated:
%
% alphaBEC
% alphaCD
% alphaEF
% alphaFG

alpha_AB = [0 0 0];

alphaBEC_vector = [0 0 alphaBEC];

alphaCD_vector = [0 0 alphaCD];

alphaEF_vector = [0 0 alphaEF];

alphaFG_vector = [0 0 alphaFG];


% ANGULAR VELOCITY VECTORS


omega_AB = [0 0 1];

omegaBEC = [0 0 angularVelocity_BEC];

omegaCD = [0 0 angularVelocity_CD];

angVel_EF = [0 0 angularVelocity_EF];

angVel_FG = [0 0 angularVelocity_FG];


% ACCELERATIONS AT ALL JOINTS

% Fixed joints
aA = [0 0 0];
aD = [0 0 0];
aG = [0 0 0];


% Acceleration of Joint B

aB_A = cross(alpha_AB,B-A) + cross(omega_AB,cross(omega_AB,B-A));

aB = aA + aB_A;


% Acceleration of Joint C

% Calculate C from B using BEC
aC_B = cross(alphaBEC_vector,C-B) + cross(omegaBEC,cross(omegaBEC,C-B));

aC_from_B = aB + aC_B;


% Calculate C from D using CD
aC_D = cross(alphaCD_vector,C-D) ...
     + cross(omegaCD,cross(omegaCD,C-D));

aC_from_D = aD + aC_D;

% Use B calculation
aC = aC_from_B;

% Acceleration of Joint E

% Calculate E from B using BEC
aE_B = cross(alphaBEC_vector,E-B) + cross(omegaBEC,cross(omegaBEC,E-B));

aE_from_B = aB + aE_B;


% Calculate E from C as a check
aE_C = cross(alphaBEC_vector,E-C) + cross(omegaBEC,cross(omegaBEC,E-C));

aE_from_C = aC + aE_C;

% Use B calculation
aE = aE_from_B;

% Acceleration of Joint F


% Calculate F from E using EF
aF_E = cross(alphaEF_vector,F-E) + cross(angVel_EF,cross(angVel_EF,F-E));

aF_from_E = aE + aF_E;


% Calculate F from G using FG
aF_G = cross(alphaFG_vector,F-G) + cross(angVel_FG,cross(angVel_FG,F-G));

aF_from_G = aG + aF_G;

% Use E calculation
aF = aF_from_E;


% Display Joint Accelerations


disp('Acceleration of A = ');
disp(aA);

disp('Acceleration of B = ');
disp(aB);

disp('Acceleration of C = ');
disp(aC);

disp('Acceleration of D = ');
disp(aD);

disp('Acceleration of E = ');
disp(aE);

disp('Acceleration of F = ');
disp(aF);

disp('Acceleration of G = ');
disp(aG);


% ACCELERATION OF CENTER OF MASS OF EACH LINK

% Acceleration of S1 - Link AB


aS1 = aA + cross(alpha_AB,S1-A) + cross(omega_AB,cross(omega_AB,S1-A));


% Acceleration of S2 - Link BEC

aS2 = aB + cross(alphaBEC_vector,S2-B) + cross(omegaBEC,cross(omegaBEC,S2-B));


% Acceleration of S3 - Link CD


aS3 = aD + cross(alphaCD_vector,S3-D) + cross(omegaCD,cross(omegaCD,S3-D));


% Acceleration of S4 - Link EF

% Calculate from E
aS4_E = aE + cross(alphaEF_vector,S4-E) + cross(angVel_EF,cross(angVel_EF,S4-E));

% Calculate from F as a check
aS4_F = aF + cross(alphaEF_vector,S4-F) + cross(angVel_EF,cross(angVel_EF,S4-F));

% Use E calculation
aS4 = aS4_E;

% Acceleration of S5 - Link FG


aS5 = aG + cross(alphaFG_vector,S5-G) + cross(angVel_FG,cross(angVel_FG,S5-G));


% Display COM Accelerations


disp('Acceleration of S1 = ');
disp(aS1);

disp('Acceleration of S2 = ');
disp(aS2);

disp('Acceleration of S3 = ');
disp(aS3);

disp('Acceleration of S4 = ');
disp(aS4);

disp('Acceleration of S5 = ');
disp(aS5);


% MAGNITUDES OF VELOCITIES


disp('|vB| = ');
disp(norm(vB));

disp('|vC| = ');
disp(norm(vC));

disp('|vE| = ');
disp(norm(vE));

disp('|vF| = ');
disp(norm(vF));

disp('|vS1| = ');
disp(norm(vS1));

disp('|vS2| = ');
disp(norm(vS2));

disp('|vS3| = ');
disp(norm(vS3));

disp('|vS4| = ');
disp(norm(vS4));

disp('|vS5| = ');
disp(norm(vS5));


% MAGNITUDES OF ACCELERATIONS

disp('|aB| = ');
disp(norm(aB));

disp('|aC| = ');
disp(norm(aC));

disp('|aE| = ');
disp(norm(aE));

disp('|aF| = ');
disp(norm(aF));

disp('|aS1| = ');
disp(norm(aS1));

disp('|aS2| = ');
disp(norm(aS2));

disp('|aS3| = ');
disp(norm(aS3));

disp('|aS4| = ');
disp(norm(aS4));

disp('|aS5| = ');
disp(norm(aS5));

