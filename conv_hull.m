%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Constructing a convex hull %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1. The first P represents the stabilizer states on the Bloch sphere.
%     Later P wil be used as a variable to store all extreme points.
% 2. (a,b,c) represents the point in the Bloch corresponding to some state
%     whose convex polytope is to be constructed using stabilizer states 
%     and given state's clifford equivalent states.
%   ** You can change (a,b,c) according to your state
%   ** Using this you can construct two convex sets one for the initial
%   state and the other for the final state. If convex polytope of state 1
%   lies inside the convex polytope of state 2, then we can convert state 2
%   to state 1.


P = [0 0 1; 0 1 0; 1 0 0; 0 -1 0; -1 0 0; 0 0 -1];
[k,av]= convhull(P);
trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','blue','Edgecolor','blue')


hold on

%{
a = 1/(sqrt(12));
b = sqrt(11/24);
P = [0 0 1; 0 1 0; 1 0 0; 0 -1 0; -1 0 0; 0 0 -1;
     a b b;  b a b;  b b a;
     a b -b; b a -b; b b -a;
     a -b b; b -a b; b -b a;
     -a b b; -b a b; -b b a;
     a -b -b; b -a -b; b -b -a;
     -a b -b; -b a -b; -b b -a;
     -a -b b; -b -a b; -b -b a;
     -a -b -b; -b -a -b; -b -b -a];
[k,av]= convhull(P)
trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','none', 'EdgeColor', 'red')
%}
%hold on

%
a = 0.344944711374856;%1/(sqrt(3));
b = 0.854087851690644;%1/(sqrt(3));%2/sqrt(14);
c = 0.389290492677355;%1/(sqrt(3));%3/sqrt(14);
P = [0 0 1; 0 1 0; 1 0 0; 0 -1 0; -1 0 0; 0 0 -1;
     a b c;  c a b;  b c a;
     a c -b; b a -c; c b -a;
     a -c b; b -a c; c -b a;
     -a c b; -b a c; -c b a;
     a -b -c; c -a -b; b -c -a;
     -a b -c; -c a -b; -b c -a;
     -a -b c; -c -a b; -b -c a;
     -a -c -b; -b -a -c; -c -b -a];
[k,av]= convhull(P);
trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','none', 'EdgeColor', 'black')
%}
%hold on

%
a = 0.108866544180631;%0.032133301801681;%2/3;
b = 0.564657410189317;%0.426628556971787;%2/3;
c = 0.556156528934490;%0.903855920648580;%1/3;
P = [0 0 1; 0 1 0; 1 0 0; 0 -1 0; -1 0 0; 0 0 -1;
     a b c;  c a b;  b c a;
     a c -b; b a -c; c b -a;
     a -c b; b -a c; c -b a;
     -a c b; -b a c; -c b a;
     a -b -c; c -a -b; b -c -a;
     -a b -c; -c a -b; -b c -a;
     -a -b c; -c -a b; -b -c a;
     -a -c -b; -b -a -c; -c -b -a];
[k,av]= convhull(P);
trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','green', 'EdgeColor', 'red')

hold off




