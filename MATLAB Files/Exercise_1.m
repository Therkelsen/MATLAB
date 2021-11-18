clear all
clc

3+4
240*124

2123/124*234
3/(5+6)

x=25;
y=30;

z=x*y

disp(['Dette er værdien af z: ',num2str(z)])

c=input('Skriv din mors vægt: ');
disp(['Din mor vejer ',num2str(c),' tons!'])

%% Matricer

A=[1 2 3]

B=[2 -3 4;3 4 -5]

C=2*B

B(2,3)

D=B.*C


G = zeros(4,5)

G = ones(4,5)

G(2,2)=8

B(end+1,:)=A

G = rand(4,5)

P = rand(4,5,3,9)

v = 1:10;
x = zeros(1,10);
y = flip(v);
z = rand(1,10);

R = [v; x; y; z]

x = -10*(1:4);

R(:,end+1)=x

R(:,11)

R(:,4)

R(:,8)=10*R(:,4)

R(:,2)=ones(4,1)

h=2:6

max(h)

max(h,pi)

max(B)

B

transpose(B)

B'

H = [2.3 -3.5 4.9; 3.3 4.4 -5.2; 3.3 4.2 -5.7]

floor(H)

H = floor(H)

size(P)

%% Lineære systemer

syms x
eqn1 = 2*x + 2 == 10;

xeqn1 = solve(eqn1)

eqn2 = 4*x^2+3 == 20;

xeqn2 = solve(eqn2)

syms x y z
eqn3 = 2*x + 3*y + z == 5;
eqn4 = -x + y - 4*z == 4;
eqn5 = x + 2*y + 3*z == -9;

[A,B] = equationsToMatrix([eqn3, eqn4, eqn5],[x, y, z])

xeqn345 = linsolve(A,B)

xeqn345 = A\B

%% Plots

x = -10:10

y = x.^3+2*x.^2-8

plot(x,y)

axis([-10 7 -20 20])

x=-10:0.2:10

y = x.^3+2*x.^2-8

plot(x,-y)
hold on
plot(x,y)
hold off
%%
figure(1)
plot(x,y)

figure (2)
plot(x,y)
hold on
plot(x,-y)
hold off

f=@(x) (x.^2+x+3)
xlabel('x')
ylabel('y')
plot(x,y, 'r-.')
hold on
plot(x,f(x))
hold off
print('figurnavn','-djpeg')

r = 2;
xc = 4;
yc = 3;

figure(1)
plot(x,y)
theta = linspace(0,2*pi);
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;
figure(2)
plot(x,y)
axis equal

%% Funktioner og plots i flere variable

x = linspace(-10,10,200)

y = linspace(-2,2,100)

[X,Y] = meshgrid(x,y)

g = cos(X).*exp(-(Y).^2)

g(80,15)

surf(g)

h = pcolor(X,Y,g);
set(h,'edgecolor','none')

