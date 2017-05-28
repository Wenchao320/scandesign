%% initialize
% Ranges are: 
% rng.de.tr = [17.5 Inf];
% rng.sp.tr = [11.8 Inf];
% rng.de.aex = [1 60] * (pi/180);                
% rng.sp.aex = [1 40] * (pi/180); 
load('Popt_263p7-0p2654.mat');
Praw = P;
figure;
hold on;
tic;

% %% change de.aex(1)
% P = Praw;
% subplot(2,4,1);
% x = linspace(2, 58, 100);
% f1 = ones(1, length(x));
% for i = 1:length(x)
%     P.de.aex(1) = x(i) * (pi/180);
%     try
%         f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
%     catch ME
%         f1(i) = NaN;
%     end
% end
% plot(x, f1);
% title('change de.aex(1)');
% axis([min(x), max(x), min(f1), 0.005]);
% 
% %% change de.aex(2)
% P = Praw;
% subplot(2,4,2);
% x = linspace(2, 58, 100);
% f1 = ones(1, length(x));
% for i = 1:length(x)
%     P.de.aex(2) = x(i) * (pi/180);
%     try
%         f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
%     catch ME
%         f1(i) = NaN;
%     end
% end
% plot(x, f1);
% title('change de.aex(2)');
% axis([min(x), max(x), min(f1), 0.005]);
% 
% %% change de.aex(3)
% P = Praw;
% subplot(2,4,3);
% x = linspace(2, 58, 100);
% f1 = ones(1, length(x));
% for i = 1:length(x)
%     P.de.aex(3) = x(i) * (pi/180);
%     try
%         f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
%     catch ME
%         f1(i) = NaN;
%     end
% end
% plot(x, f1);
% title('change de.aex(3)');
% axis([min(x), max(x), min(f1), 0.005]);

%% change sp.aex(1)
P = Praw;
subplot(1,2,1);
x = linspace(2, 38, 100);
f1 = ones(1, length(x));
for i = 1:length(x)
    P.sp.aex(1) = x(i) * (pi/180);
    try
        f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
    catch ME
        f1(i) = NaN;
    end
end
plot(x, f1);
title('change sp.aex(1)');
% axis([min(x), max(x), 0, 0.0003]);

% %% change de.tr(1)
% P = Praw;
% subplot(2,4,5);
% x = linspace(17.5, 50, 100);
% f1 = ones(1, length(x));
% for i = 1:length(x)
%     P.de.tr(1) = x(i);
%     try
%         f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
%     catch ME
%         f1(i) = NaN;
%     end
% end
% plot(x, f1);
% title('change de.tr(1)');
% axis([min(x), max(x), min(f1), 0.005]);
% 
% %% change de.tr(2)
% P = Praw;
% subplot(2,4,6);
% x = linspace(17.5, 50, 100);
% f1 = ones(1, length(x));
% for i = 1:length(x)
%     P.de.tr(2) = x(i);
%     try
%         f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
%     catch ME
%         f1(i) = NaN;
%     end
% end
% plot(x, f1);
% title('change de.tr(2)');
% axis([min(x), max(x), min(f1), 0.005]);
% 
% %% change de.tr(3)
% P = Praw;
% subplot(2,4,7);
% x = linspace(17.5, 50, 100);
% f1 = ones(1, length(x));
% for i = 1:length(x)
%     P.de.tr(3) = x(i);
%     try
%         f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
%     catch ME
%         f1(i) = NaN;
%     end
% end
% plot(x, f1);
% title('change de.tr(3)');
% axis([min(x), max(x), min(f1), 0.005]);

%% change sp.tr(1)
P = Praw;
subplot(1,2,2);
x = linspace(11.8, 50, 100);
f1 = ones(1, length(x));
for i = 1:length(x)
    P.sp.tr(1) = x(i);
    try
        f1(i) = dess_spgr_2comp_cost(P, subArg.cost{:});
    catch ME
        f1(i) = NaN;
    end
end
plot(x, f1);
title('change sp.tr(1)');
% axis([min(x), max(x), 0, 0.0003]);

%% done
t = toc;
fprintf('Program runs %0.2f minutes.\n', t / 60);