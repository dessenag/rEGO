function value = frame_u(x,zi)
%% Please cite the works under "References" when using this program.
% This code can be found at https://doi.org/10.5281/zenodo.8406030
% G. Dessena, rEGO – A tutorial on the refined Efficient Global 
% Optimisation, Zenodo, Oct. 10, 2023. doi: 10.5281/zenodo.8406030
%% Disclaimer
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.

%% Credits
% Implementation of the refined Efficient Global Optimisation:
% Gabriele Dessena
% gdessena@ing.uc3m.es
% Universidad Carlos III de Madrid
% 2023/10/02 v0.1 - pre release
% 2023/10/10 v1.0b - release
%% References
%% Please cite the works under "References" when using this program
% [1] G. Dessena, D. I. Ignatyev, J. F. Whidborne, and L. Zanotti 
% Fragonara, ‘A global–local meta-modelling technique for model updating’, 
% Computer Methods in Applied Mechanics and Engineering, vol. 418. Elsevier 
% BV, p. 116511, Jan. 2024. doi: 10.1016/j.cma.2023.116511.
% [2] G. Dessena, rEGO – A tutorial on the refined Efficient Global 
% Optimisation, Zenodo, Oct. 10, 2023. doi: 10.5281/zenodo.8406030
% [3] E. Figueiredo, G. Park, J. Figueiras, C. Farrar, and K. Worden, 
% ‘Structural health monitoring algorithm comparisons using standard data 
% sets’, Office of Scientific and Technical Information (OSTI), Mar. 2009. 
% doi: 10.2172/961604.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
