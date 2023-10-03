function yy = frame_u_opti(x,mat,type)

mat2 = frame_u(x,mat(2,:));
freq = mat(1,:); ffreq = mat2(1,:);
fi = mat(3:end,:); ffi = mat2(3:end,:);
if strcmp(type,'rmse_freq')
    yy = RMSE(ffreq,freq);
elseif strcmp(type,'rmse')
    yy = RMSE(ffreq,freq)+RMSE(ffi,fi);
elseif strcmp(type,'mu_mac')
    mac = diag(compute_mac(fi,ffi));
    yy=mean(1-mac);
elseif strcmp(type,'tmac')
    mac = diag(compute_mac(fi,ffi));
    yy=1 - prod(mac);
elseif strcmp(type,'mtmac')
    mac = diag(compute_mac(fi,ffi))';
    yt = mac.*((1+abs((ffreq-freq)./(ffreq+freq))).^(-1));
    yy=1-prod(yt);
else
    disp('error')
end
end

function y=RMSE(x,x0)
y=sqrt(mean((x - x0).^2,'all'));
end
function [MAC_MATRIX] = compute_mac(U1,U2)
% Compute MAC between to ModeShape matrices U1 and U2

MAC_MATRIX=zeros(size(U1,2),size(U2,2));

for i=1:size(U1,2)
    for j=1:size(U2,2)
        MAC_MATRIX(i,j)=(dot((U1(:,i)),U2(:,j))^2)/dot(dot(U1(:,i),U1(:,i)),dot(U2(:,j),U2(:,j)));
    end
end
end