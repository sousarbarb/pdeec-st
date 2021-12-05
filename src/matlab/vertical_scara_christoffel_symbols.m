%% INITIALISATION
scara_inertia_matrix
q = [q1 q2 q3];

%% CHRISTOFFEL SYMBOLS
% k=1
c111 = ChristoffelSymbols(D,q,1,1,1)
c121 = ChristoffelSymbols(D,q,1,2,1)
c131 = ChristoffelSymbols(D,q,1,3,1)
c211 = ChristoffelSymbols(D,q,2,1,1)
c221 = ChristoffelSymbols(D,q,2,2,1)
c231 = ChristoffelSymbols(D,q,2,3,1)
c311 = ChristoffelSymbols(D,q,3,1,1)
c321 = ChristoffelSymbols(D,q,3,2,1)
c331 = ChristoffelSymbols(D,q,3,3,1)

% k=2
c112 = ChristoffelSymbols(D,q,1,1,2)
c122 = ChristoffelSymbols(D,q,1,2,2)
c132 = ChristoffelSymbols(D,q,1,3,2)
c212 = ChristoffelSymbols(D,q,2,1,2)
c222 = ChristoffelSymbols(D,q,2,2,2)
c232 = ChristoffelSymbols(D,q,2,3,2)
c312 = ChristoffelSymbols(D,q,3,1,2)
c322 = ChristoffelSymbols(D,q,3,2,2)
c332 = ChristoffelSymbols(D,q,3,3,2)

% k=3
c113 = ChristoffelSymbols(D,q,1,1,3)
c123 = ChristoffelSymbols(D,q,1,2,3)
c133 = ChristoffelSymbols(D,q,1,3,3)
c213 = ChristoffelSymbols(D,q,2,1,3)
c223 = ChristoffelSymbols(D,q,2,2,3)
c233 = ChristoffelSymbols(D,q,2,3,3)
c313 = ChristoffelSymbols(D,q,3,1,3)
c323 = ChristoffelSymbols(D,q,3,2,3)
c333 = ChristoffelSymbols(D,q,3,3,3)