%[theta phi] = meshgrid(0:0.1:359.9,0:0.1:179.9);
%M = load('heat_map.data');
M = load('bins_res1_ideal');
mesh(M);
%image([0 180], [0 360], M, 'CDataMapping', 'scaled');