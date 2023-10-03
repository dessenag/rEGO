function [MAC_MATRIX] = compute_mac(U1,U2)
%Computes MAC
%   Computes MAC value between two sets of mode shapes
% mac

MAC_MATRIX=zeros(size(U1,2),size(U2,2));

for i=1:size(U1,2)
    for j=1:size(U2,2)
        MAC_MATRIX(i,j)=(dot((U1(:,i)),U2(:,j))^2)/dot(dot(U1(:,i),U1(:,i)),dot(U2(:,j),U2(:,j)));
    end
end


