function [] = visualize(K)

% draw a unit sphere
n = 25;
r = ones(n, n); % radius is 5
[th, phi] = meshgrid(linspace(0, pi, n), linspace(0, 2*pi, n));
[x,y,z]= sph2cart(th, phi, r);
surface(x,y,z,'FaceColor', 'none');
hold on;

for k = 1:K
   data_file = strcat('comp',num2str(k),'.dat');
   M = load(data_file);
   x = M(:,1);
   y = M(:,2);
   z = M(:,3);
   colors = rand(1,3);
   plot3(x,y,z,'.','Color',colors);
end    

legend('unit sphere','comp1','comp2');
end