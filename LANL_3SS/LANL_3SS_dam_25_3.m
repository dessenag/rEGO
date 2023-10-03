clear all
close all
%% correction factors - from update
med = [0.965673613318987,1.14306334469746,1.25282420940468,...
    1.14692873288260,1.46307302627545,1.70301279337578,1.67849727191707];
%% Init Mass matrix
m_1 = 6.279+.143;
m_2 = 6.279+.286;
m_4 = 6.532+.143;
M = diag(med(1:4).*[m_1,m_2,m_2,m_4]);
%% Init Stiffness matrix
k = med(5:7).*4*68.167e3;
k(3) = k(3)*.75;
K =[[1+k(1) -k(1) 0 0];
    [-k(1) k(1)+k(2) -k(2) 0];
    [0 -k(2) k(3)+k(2) -k(3)]
    [0 0 -k(3) k(3)]];

z = [0 0.06 0.02 0.008]; %damping

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
damaged = [real(freq);z;fi(:,I)];

save LANL_3SS_dam_25_3 damaged  