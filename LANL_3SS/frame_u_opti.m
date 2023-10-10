function yy = frame_u_opti(x,mat,type)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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