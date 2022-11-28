clear; clc; close all;

% Description: Plots the velocity profiles for the fundamental solution of
% stoke's flow

% Load the data
xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');
xv = load('v_x_mesh.txt');
yv = load('v_y_mesh.txt');
xp = load('p_x_mesh.txt');
yp = load('p_y_mesh.txt');

% Domain size
Lx = max(xu,[],'all');
Ly = max(yv,[],'all');

uFile = dir(strcat('u_0','*'));
vFile = dir(strcat('v_0','*'));
PFile = dir(strcat('P_0','*'));
FFile = dir(strcat('F_0','*'));

nFiles = length(uFile);
u = load(uFile(nFiles).name);
v = load(vFile(nFiles).name);


% Store the values at cell centers
Nx = size(u,2)-1;
Ny = size(u,1)-2;
uc = zeros(Ny,Nx);
vc = zeros(Ny,Nx);

umag = sqrt(uc.^2+vc.^2);

% Total number of files (Same for velocity components, cilia and particles)
nFiles = length(uFile);

% Visualize cilia motion over velocity field
colormap(jet)
vid = VideoWriter('ibm_c.avi','Uncompressed AVI');
open(vid);
figure(1)
fig = gcf;
fig.Position = [1 1 1920 961];
for iFile = 1:nFiles
    %subplot(2,2,1)
    clc
    fprintf('Processing file # %d of %d\n',iFile,nFiles)
    hold on
    
    %% Plot velocity field
    u = load(uFile(iFile).name);
    v = load(vFile(iFile).name);
    % Convert velocity to cell center values
    [uc,vc] = cellcenter(uc,vc,u,v,Nx,Ny);
    umag = sqrt(uc.^2+vc.^2);
    
    contourf(xp,yp,umag,50,'edgecolor','none')
    
     %% Plot mesh
     mesh(xp,yp,0*xp,'FaceAlpha','0.0','EdgeColor','w','LineStyle','-','EdgeAlpha','0.25')
% view(90,0)
    colorbar
    axis equal
    title(uFile(iFile).name)

    %% Plot FEM mesh
    P = load(PFile(iFile).name)';
    M = load(FFile(iFile).name)';
    for i = 1:2:size(P,2)
        triplot(M,P(:,i),P(:,i+1),'k-','linewidth',5)
    end
%     pause(0.5)
    writeVideo(vid,getframe(gca));
    if iFile ~= nFiles
        clf
    end
end
close(vid)

% Function to convert staggered velocity to cell-centered velocity
function [uc,vc] = cellcenter(uc,vc,u,v,Nx,Ny)
    for j = 2:Nx+1
        for i = 1:Ny
            vc(i,j-1) = v(i,j); 
        end
    end

    for j = 1:Nx
        for i = 2:Ny+1
            uc(i-1,j) = u(i,j); 
        end
    end
end
