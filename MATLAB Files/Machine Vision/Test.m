% Simpel algebra
2+2

2*2

% Variabel manipulation
x=1:5

y=x*2

% Funktioner
f=@(x) x.^3+x;

f(2)

plot(f(1:5))
xlabel('x')
ylabel('y')

% Matricer
A=[1 2 3; 4 5 6; 7 8 9]

A(1,2)

 % Største værdi i hver kolonne
max(A)

 % Største værdi i hele matrixen
max(max(A))

 % Mindste værdi i hver kolonne
min(A)
 
 % Mindste værdi i hele matrixen
min(min(A))

 % Find det element som har værdien x
find(A==3)



