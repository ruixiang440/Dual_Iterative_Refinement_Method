function S = compute_LaplacianBasis(S,numEigs)
%Compute the eigen vectors of the cotangent Laplacian matrix
if nargin < 2, numEigs = 200; end
fprintf('Computing %d Eigen-functions...',numEigs); tic;

% A: mixed voronoi area weights
[W, A] = cotLaplacian(get_mesh_vtx_pos(S), S.surface.TRIV);
% A: one-ring neighbor area weights
% A = vertexAreas(get_mesh_vtx_pos(S), S.surface.TRIV);
W =((W+W')/2);
A = sparse(1:length(A), 1:length(A), A);
W(W==Inf) = 0;
W(W==-Inf) = 0;
try
    [evecs, evals] = eigs(W, A, numEigs, 1e-6);
catch
    % In case of trouble make the laplacian definite
    [evecs, evals] = eigs(W - 1e-8*speye(S.nv), A, numEigs, 'sm');
end
evals = diag(evals);
if ~isreal(evecs)
    evecs(1:2:end) = real(evecs(1:2:end));
    evecs(2:2:end) = imag(evecs(1:2:end));
end
[evals, order] = sort(abs(evals),'ascend');
S.evals = evals;
evecs = evecs(:,order);
S.evecs = evecs;
S.A = A;
S.W = W;
S.area = diag(S.A);
S.sqrt_area = sqrt(sum(S.area));
t =toc; fprintf('done:%.4fs\n',t);
end