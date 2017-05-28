% fprintf('Pair (1):\n');
% load('Popt_108-0p1230.mat');
% printP(P);
% P2 = repeatP(P,2);
% printP(P2);
% P3 = repeatP(P,3);
% printP(P3);
% P4 = repeatP(P,4);
% printP(P4);
% 
% fprintf('Pair (2):\n');
% load('popu50de3sp1.mat');
% P = P.opt;
% printP(P);
% P2 = repeatP(P,2);
% printP(P2);

% fprintf('Pair (3):\n');
% load('5DE1SPinit.mat');
% P = popu{3}.P;
printP(P);
P2 = repeatP(P,2);
printP(P2);
P3 = repeatP(P,3);
printP(P3);
P4 = repeatP(P,4);
printP(P4);