function value = frame_u(x,zi)

m_1 = 6.279+.143;
m_2 = 6.279+.286;
m_4 = 6.532+.143;
med = [0.965673613318987,1.14306334469746,1.25282420940468,...
    1.14692873288260,1.46307302627545,1.70301279337578,1.67849727191707];
k = med(5:7).*4*68.167e3.*x(5:7);
% k = 4*100e3.*ones(1,3);
z = [0 zi];
M = diag(med(1:4).*[m_1,m_2,m_2,m_4].*x(1:4));
K =[[1+k(1) -k(1) 0 0];
    [-k(1) k(1)+k(2) -k(2) 0];
    [0 -k(2) k(3)+k(2) -k(3)]
    [0 0 -k(3) k(3)]];

[fi,la] = eig(M \ K);

% mode shapes
for ii = 1:4
    fi(:,ii) = fi(:,ii) / sqrt( fi(:,ii).' * M * fi(:,ii) );
end
fi = fi./max(abs(fi));
%natural frequencies
wn = sqrt(diag(la));
wd = wn' .* sqrt(1 - z.^2);
[freq,I] = sort(wd/(2*pi));
mat = [real(freq);z;fi(:,I)];
value = mat(:,2:end);
end
