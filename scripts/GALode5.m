function dydt = GALode5( t, y, param )
% 2017.09.07 GAL model v5.0
% Major update - add the 12th species HXT into the equations, similar
% formula of G2
y = max(y,0);  % to solve the negative value problem

% define variables
G1 = y(1);
G2 = y(2);
G3 = y(3);
G4 = y(4);
G80 = y(5);
G3s = y(6);
Mig1tot = y(7);
C83 = y(8);
C84 = y(9);
glu = y(10);
gal = y(11);
HXT = y(12);

% unwrap parameters
a1 = param.a1;
a2 = param.a2;
a3 = param.a3;
a4 = param.a4;
a80 = param.a80;
aR = param.aR;
a0HXT = param.a0HXT;
aHXT = param.aHXT;
ag1 = param.ag1;
ag2 = param.ag2;
ag3 = param.ag3;
ag4 = param.ag4;
ag80 = param.ag80;
d = param.d;
dsugar = param.dsugar;
kf3 = param.kf3;
kr3 = param.kr3;
kf83 = param.kf83;
kr83 = param.kr83;
kf84 = param.kf84;
kr84 = param.kr84;
KG1 = param.KG1;
KG2 = param.KG2;
KG3 = param.KG3;
KG80 = param.KG80;
KR1 = param.KR1;
KR3 = param.KR3;
KR4 = param.KR4;
KRs = param.KRs;    % the affinity of glucose binding to free Mig1
KHXT = param.KHXT;  % the hill constant of G4 inhibiting HXT synthesis
kG2 = param.kG2;        % G2 transportation rate
rcat = param.rcat;      % kHXT/kG2              the ratio of the max uptake rate
rG2 = param.rG2;        % KGgal/KGglu           the ratio of G2 binding affinity
rHXT = param.rHXT;       % KHXTgal/KHXTglu      the ratio of HXT binding affinity
KGglu = param.KGglu;    % the binding affinity between G2 and glucose
KHXTglu = param.KHXTglu;    % the binding affinity between HXT and glucose
n1 = param.n1;
n2 = param.n2;
n3 = param.n3;
n80 = param.n80;
nR1 = param.nR1;
nR3 = param.nR3;
nR4 = param.nR4;
nRs = param.nRs;    % the hill coefficient of glucose binding to free Mig1
nHXT = param.nHXT;  % the hill coefficient of G4 inhibiting HXT synthesis
exglu = param.exglu;
exgal = param.exgal;
kHXT = rcat * kG2;
KGgal = rG2 * KGglu;
KHXTgal = rHXT * KHXTglu;
Mig1s = glu^nRs/(KRs^nRs + glu^nRs) * Mig1tot;

G1_synthesis_term = G4^n1/(G4^n1 + KG1^n1) * KR1^nR1/(KR1^nR1 + Mig1s^nR1);
G2_synthesis_term = G4^n2/(G4^n2 + KG2^n2);
G3_synthesis_term = G4^n3/(G4^n3 + KG3^n3) * KR3^nR3/(KR3^nR3 + Mig1s^nR3);
G4_synthesis_term = KR4^nR4/(KR4^nR4 + Mig1s^nR4);
G80_synthesis_term = G4^n80/(G4^n80 + KG80^n80);
HXT_synthesis_term = KHXT^nHXT/(KHXT^nHXT + G4^nHXT);

G3star_f = kf3 * gal * G3;
G3star_r = kr3 * G3s;

C83_f = kf83 * G3s * G80;
C83_r = kr83 * C83;

C84_f = kf84 * G4 * G80;
C84_r = kr84 * C84;

dG1 = a1 + ag1 * G1_synthesis_term - d * G1;
dG2 = a2 + ag2 * G2_synthesis_term - d * G2;
dG3 = a3 + ag3 * G3_synthesis_term - G3star_f + G3star_r - d * G3;
dG4 = a4 + ag4 * G4_synthesis_term - C84_f + C84_r - d * G4;
dG80 = a80 + ag80 * G80_synthesis_term - C83_f + C83_r - C84_f + C84_r - d * G80;
dG3s = G3star_f - G3star_r - C83_f + C83_r - d * G3s;
dMig1tot = aR - d * Mig1tot;
dC83 = C83_f - C83_r - d * C83;
dC84 = C84_f - C84_r - d * C84;
dglu = kG2 * G2 * exglu/(1/rG2 * exgal + exglu + KGglu) + kHXT * HXT * exglu/(1/rHXT * exgal + exglu + KHXTglu) - dsugar * glu;
dgal = kG2 * G2 * exgal/(rG2 * exglu + exgal + KGgal) + kHXT * HXT * exgal/(rHXT * exglu + exgal + KHXTgal) - G3star_f + G3star_r  - dsugar * gal;
dHXT = a0HXT + aHXT * HXT_synthesis_term - d * HXT;

dydt = ...
    [ dG1 ...
    ; dG2 ...
    ; dG3 ...
    ; dG4 ...
    ; dG80 ...
    ; dG3s ...
    ; dMig1tot ...
    ; dC83 ...
    ; dC84 ...
    ; dglu ...
    ; dgal ...
    ; dHXT];

end