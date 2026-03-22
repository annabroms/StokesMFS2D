function [u,lambda_x,lambda_y] = getRecompletionFlowStokes1BEnhanced(q,F,T,geom)
%GETRECOMPLETIONFLOWSTOKES1BENHANCED Completion flow for 1B mobility solves.
%
% Syntax:
%   [u,lambda_x,lambda_y] = getRecompletionFlowStokes1BEnhanced(q,F,T,geom)
%
% Inputs:
%   q    - Complex particle centers.
%   F    - Px2 matrix of net forces.
%   T    - Px1 vector of torques.
%   geom - Geometry struct from prepareStokes1BEnhanced.
%
% Outputs:
%   u        - Stacked completion flow on the collocation grid.
%   lambda_x  - Stacked x-strengths of the completion sources.
%   lambda_y  - Stacked y-strengths of the completion sources.

q = q(:);
T = T(:);
P = numel(q);
if size(F,1) ~= P || size(F,2) ~= 2 || numel(T) ~= P
    error('getRecompletionFlowStokes1BEnhanced:badLoads', ...
        'F must be Px2 and T must have one entry per particle.');
end

lambda_x = zeros(geom.total_source_count,1);
lambda_y = zeros(geom.total_source_count,1);

for k = 1:P
    Ksrc = geom.Ksrc{k};
    Fk = [F(k,:)'; T(k)];
    lambda0k = Ksrc*((Ksrc'*Ksrc)\Fk);

    ns = geom.source_count(k);
    sx = geom.source_x_ind{k};
    lambda_x(sx) = lambda0k(1:ns);
    lambda_y(sx) = lambda0k(ns+1:end);
end

u = -getVelocityField(geom.rvec_in,geom.rout,lambda_x,lambda_y,[],[],[],[],[],[],geom.use_fmm);

end
