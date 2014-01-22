# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:46:21 2013

@author: Will
"""

def Field_Realistic(x,y,z,
    B_xbias,B_ybias,B_zbias, fin_horz_params, fin_vert_params, fin_norm_params):
## ----------Defs
mu=(4*pi)*10^-7 #N/A^2
#Beta=y.^2+z.^2
gridsize = size(x)

## Finite Horizontal wires

#Unpack parameters
LL = fin_horz_params(1,:)
LR = fin_horz_params(2,:)
Y0_fin_horz = fin_horz_params(3,:)
Z0_fin_horz = fin_horz_params(4,:)
I_Guide = fin_horz_params(5,:)

n_fin_horz = size(fin_horz_params)n_fin_horz = n_fin_horz(2)

#preallocate resources
const_G = zeros(size(I_Guide))
B_G = zeros([gridsize,n_fin_horz])
B_Gy = B_G
B_Gz_horz = B_G

for i=1:n_fin_horz
    Beta = (z-Z0_fin_horz(i)).^2 + (y-Y0_fin_horz(i)).^2
    const_G(i)=mu.*I_Guide(i)./(4*pi)
    B_G(:,:,i)=const_G(i).*((x-LL(i))./(Beta.*sqrt(Beta+...
        (x-LL(i)).^2))-(x-LR(i))./(Beta.*sqrt(Beta+(x-LR(i)).^2)))
    B_Gy(:,:,i)=B_G(:,:,i).*(Z0_fin_horz(i)-z)
    B_Gz_horz(:,:,i)=B_G(:,:,i).*(y-Y0_fin_horz(i))
end

## Finite Vertical Wires

#unpack parameters
YD = fin_vert_params(1,:)
YU = fin_vert_params(2,:)
X0_fin_vert = fin_vert_params(3,:)
Z0_fin_vert = fin_vert_params(4,:)
I_Guide = fin_vert_params(5,:)

n_fin_vert = size(fin_vert_params)n_fin_vert = n_fin_vert(2)

#preallocate resources
const_G = zeros(size(I_Guide))
B_G = zeros([gridsize,n_fin_vert])
B_Gx = B_G
B_Gz_vert = B_G

for i=1:n_fin_vert
    Beta = (z-Z0_fin_vert(i)).^2 + (x-X0_fin_vert(i)).^2
    const_G(i)=mu.*I_Guide(i)./(4*pi)
    B_G(:,:,i)=const_G(i).*((y-YD(i))./(Beta.*sqrt(Beta+...
        (y-YD(i)).^2))-(y-YU(i))./(Beta.*sqrt(Beta+(y-YU(i)).^2)))
    B_Gx(:,:,i)=B_G(:,:,i).*(z-Z0_fin_vert(i))
    B_Gz_vert(:,:,i)=B_G(:,:,i).*(X0_fin_vert(i)-x)
end

## Finite Normal Wires

#unpack parameters
ZD = fin_norm_params(1,:)
ZU = fin_norm_params(2,:)
X0 = fin_norm_params(3,:)
Y0 = fin_norm_params(4,:)
I = fin_norm_params(5,:)

n_fin_norm = size(fin_norm_params)n_fin_norm = n_fin_norm(2)

#preallocate resources
const_G = zeros(size(I))
B_G = zeros([gridsize,n_fin_norm])
B_Gx_norm = B_G
B_Gy_norm = B_G

for i=1:n_fin_norm
    Beta = (x-X0(i)).^2 + (y-Y0(i)).^2
    const_G(i)=mu.*I(i)./(4*pi)
    B_G(:,:,i)=const_G(i).*((z-ZD(i))./(Beta.*sqrt(Beta+...
        (z-ZD(i)).^2))-(z-ZU(i))./(Beta.*sqrt(Beta+(z-ZU(i)).^2)))
    B_Gx_norm(:,:,i)=B_G(:,:,i).*(Y0(i)-y)
    B_Gy_norm(:,:,i)=B_G(:,:,i).*(x-X0(i))
end

B_Gy_tot=sum(B_Gy,3)+sum(B_Gy_norm,3)
B_Gz_tot=sum(B_Gz_horz,3)+sum(B_Gz_vert,3)
B_Gx_tot=sum(B_Gx,3)+sum(B_Gx_norm,3)
B_tot=sqrt((B_Gx_tot + B_xbias).^2+(B_Gz_tot...
    +B_zbias).^2+(B_Gy_tot+ B_ybias).^2)