function Pout = repeatP(P, times)

Pout.de.aex = [];
Pout.sp.aex = [];
Pout.de.tr = [];
Pout.sp.tr = [];
for i = 1:times
    Pout.de.aex = [Pout.de.aex; P.de.aex];
    Pout.sp.aex = [Pout.sp.aex; P.sp.aex];
    Pout.de.tr = [Pout.de.tr; P.de.tr];
    Pout.sp.tr = [Pout.sp.tr; P.sp.tr];
end