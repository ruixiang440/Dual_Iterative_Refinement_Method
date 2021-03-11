%%
clc;
close all;
clear;
addpath(genpath('.'));
% files = [1,4,13,16];
files = [5];
spec_dim = 1000; % max spectrum dimension
spec_dim_cut = 960; % max spectrum dimension of landmarks
low = 0.1; % local distortion lower bound
gap = 5; % local distortion gaps
iter_number = 10; % total iteration number
th = 0.5-(0.5-low)/gap:-(0.5-low)/gap:low; % local distortion
% e.g here, there will be totally 10 iterations with local distortion at
% each iteration equals to [0.42, 0.34, 0.26, 0.18, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
option1.nb_iter_max = 30; % iteration number needed for fast marching
option2.nb_iter_max = 120;
%%
res = zeros(length(files),1000);
ctime = zeros(length(files),1);
tic;
[surf2.pt,surf2.trg] = ReadOFF('mesh000.off');
S2 = MESH.MESH_IO.read_shape('mesh000.off');
surf2.Phi = MESH.compute_LaplacianBasis(S2, spec_dim);
surf2.n = length(surf2.pt);
num2 = surf2.n;
opts.shot_num_bins = 10; % number of bins for shot
opts.shot_radius = 5; % percentage of the diameter used for shot
Xdesc = calc_shot(surf2.pt', surf2.trg', 1:num2, opts.shot_num_bins, opts.shot_radius*surf2.Phi.sqrt_area/100, 3)';
vertex2 = surf2.pt';
faces2 = surf2.trg';
pre_time = toc;
s = randsample(num2,1000);
option3.nb_iter_max = inf;
D3 = perform_fast_marching_mesh(vertex2, faces2, s, option3);


for iter = 1:length(files)
    tic
    if files(iter)<10
        [surf1.pt,surf1.trg] = ReadOFF(['mesh00' num2str(files(iter)) '.off']);
        S1 = MESH.MESH_IO.read_shape(['mesh00' num2str(files(iter)) '.off']);
    else
        [surf1.pt,surf1.trg] = ReadOFF(['mesh0' num2str(files(iter)) '.off']);
        S1 = MESH.MESH_IO.read_shape(['mesh0' num2str(files(iter)) '.off']);
    end
    surf1.Phi = MESH.compute_LaplacianBasis(S1, spec_dim);
    true_P = 1:num2;
    num = length(surf1.pt);
    surf1.n = num;
    MA = full(diag(surf1.Phi.A));
    xdesc = calc_shot(surf1.pt', surf1.trg', 1:num, opts.shot_num_bins, opts.shot_radius*surf1.Phi.sqrt_area/100, 3)';
    pertF = knnsearch(Xdesc, xdesc,'NSMethod','kdtree');
    ini = pertF;
    vertex1 = surf1.pt';
    faces1 = surf1.trg';
    cnt = 0;
    e = zeros(num,1);
    D1 = perform_fast_marching_mesh(vertex1, faces1, 1:num, option1);
    D2 = perform_fast_marching_mesh(vertex2, faces2, 1:num, option2);
    R_max = max(max(D2));
    landmarks = [];
    for kk = 1:iter_number
        DD1 = cell(num,1);
        DD2 = cell(num,1);
        idx = cell(num,1);
        idx2 = cell(num,1);
        idx1 = cell(num,1);
        ee = 0;
        good = 1:num;
        D1T = D1(:,good);
        D2T = D2(:,pertF(good));
        for i = 1:length(good)
            idx{i} = find(D1T(:,i) ~= 0);
            DD1{i} = D1T(idx{i},i);
            DD2{i} = D2T(idx{i},i);
            DD2{i}(DD2{i}==0) = R_max;
            r = max(DD1{i});
            ee(i) = sum(((abs(DD1{i}-DD2{i}))/r).*MA(idx{i}))/sum(MA(idx{i}));
        end
        e(good)=ee;
        if kk<=length(th)
            landmarks = find(e<th(kk));
        else
            landmarks = find(e<th(end));
        end
        sub_landmarks = landmarks;
        good = setdiff(1:num,sub_landmarks);
        goodF = setdiff(1:num,pertF(sub_landmarks));
        H = surf1.Phi.evecs(sub_landmarks,:)'*surf2.Phi.evecs(pertF(sub_landmarks),:);
        [U,D,V] = svd(H);
        dim = findK(diag(D));
        if dim>spec_dim_cut
            cnt = cnt + 1;
        end
        if cnt>3
            break;
        end
        specB = surf2.Phi.evecs(:,1:dim);
        specA = surf1.Phi.evecs(:,1:dim);
        H = specA(sub_landmarks,:)'*specB(pertF(sub_landmarks),:);
        [U,D,V] = svd(H);
        C = U*V';
        specBB = specB(goodF,:);
        specAA = specA(good,:);
        p = knnsearch(specBB*C',specAA,'NSMethod','kdtree');
        pertF(good) = goodF(p);
    end
    ctime(iter) = toc + pre_time;
    D3T = D3(pertF,:);
    geoerror = zeros(1000,1);
    for i=1:length(s)
        a=true_P(s(i));
        dis1 = D3T(a,i);
        geoerror(i)=dis1/surf1.Phi.sqrt_area;
    end
    hold on
    res(iter,:)=full(geoerror);
    disp(['mesh00' num2str(files(iter)) '.off mapping' ' is finished,' 'time is ' num2str(ctime(iter)) 'seconds.'])
    cdfplot(geoerror)
    title('Geo error')
    ylabel('%Correspondence')
    xlim([0,0.25])
end