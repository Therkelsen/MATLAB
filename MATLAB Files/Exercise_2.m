clc
clear all

save('data')

%%
clc
clear all

x = input('Angiv en x værdi: ');; % Det er denne vi vil tjekke.
minVal = 2;
maxVal = 6;

if (x >= minVal) && (x <= maxVal)
disp('værdien er inden for intervallet.')
elseif (x > maxVal)
disp('Værdien er større end maximum værdien.')
else
disp('Værdien er mindre end minimum værdien.')
end

%% Løkker
clc
clear all

z=[]

for n = 1:6
z(n) = (n+3)*4;
end

z

h=[];

for n = 4:20
h(n) = n * (n - 1) + 24 - n;
end

h

h=[];

for n = 20:-1:4
h(n) = n * (n - 1) + 24 - n;
end

h

for i = 1:10000
    i
    if i>500
        break
    end
end

n = 1;

z=[]

while n < 10
n = n + 1
z(n)=(n+3)*4;
end

z

k=[]

while n < 20
n = n + 1
k(n)=n*(n-1)-5;
end

k

%% Numeriske funktioner

f = @(x) (x.^3+2*x.^2-8);
[x0,fopt] = fzero(f,-3);

[x0,fopt]

plot(f(-5:5))

%% Definer egne funktioner
clc
clear all

z=[];

for n = 1:6
z(n) = (n+3)*4;
end

z

%% Ekstraopgave
clc
clear all

n=10;
k=0;
for i = 1:n
    for j = 1:n
        k=k+i+j;
    end
end
k

n=100;
k=100;
for i = 1:n
    for j = 1:n
        k=(i*j);
    end
end
k

n=321;
k=64;
for i = 1:n
    for j = 1:n
        k=(i+j);
    end
end
k

for i = 2:100
    for j = 2:100
        if(~mod(i,j))
            break;
        end
    end
    if(j > (i/j))
        fprintf('%d er et primtal\n', i);
    end
end



%% Custom funktioner
function [ z ] = minfunktion( n )
% minfunktion Skriv noget forklarende tekst her
for m = 1:n
z(m) = (m+3)*4;
end;
end

function [ z ] = minfunktion2( n )
% minfunktion Skriv noget forklarende tekst her
n=1
while n<10
n = n+1;
z(n) = (n+3)*4;
end;
end

