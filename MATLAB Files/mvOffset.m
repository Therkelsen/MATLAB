x = [-50 0 50 100 150 200 250 300 350 400 450 500 550 600 650 700]
y = [-50 0 50 100 150 200 250 300 350 400 450 500 550 600 650 700]
z = transpose(testafvision)
nan = zeros(1,48)
nan(:,:) = missing
z = [z nan]
z = reshape(z, [16,16])
figure(1)
stem3(x, y, z)
title('Offset')
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Offset [mm]')
grid on
xv = linspace(min(x), max(x), 20);
yv = linspace(min(y), max(y), 20);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);
figure(2)
surf(X, Y, Z);
grid on
set(gca, 'ZLim',[-2 10])
shading interp
title('Offset')
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Offset [mm]')