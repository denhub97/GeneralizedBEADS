% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
%-------------------------------------------------------------------------%
% Define the spherical harmonics
% Angle grid

theta=linspace(0,pi,res_theta);
phi=linspace(0,2*pi,res_phi);
[phi,theta]=meshgrid(phi,theta);


%The complex Spherical Harmomics
% rank 0
Y00_ = 1/2*sqrt(1/pi)*ones(res_theta,res_phi);

% rank 1
Y1n1_ = sqrt(3/8/pi).*sin(theta).*exp(-1i*phi);         
Y10_ = sqrt(3/4/pi).*cos(theta);
Y11_ = -sqrt(3/8/pi).*sin(theta).*exp(1i*phi);

% rank 2
Y2n2_ = sqrt(15/(32*pi)).*sin(theta).^2.*exp(-2i*phi);
Y2n1_ = sqrt(15/(8*pi)).*sin(theta).*cos(theta).*exp(-1i*phi);
Y20_ = sqrt(5/(16*pi)).*(3*cos(theta).^2-1);
Y21_ = -sqrt(15/(8*pi)).*sin(theta).*cos(theta).*exp(1i*phi);
Y22_ = sqrt(15/(32*pi)).*sin(theta).^2.*exp(2i*phi);

% rank 3
Y3n3_ = sqrt(35/(64*pi))*sin(theta).^3.*exp(-3i*phi);
Y3n2_ = sqrt(105/(32*pi))*sin(theta).^2.*cos(theta).*exp(-2i*phi);
Y3n1_ = sqrt(21/(64*pi))*sin(theta).*(5*cos(theta).^2-1).*exp(-1i*phi);
Y30_ = sqrt(7/(16*pi))*(5*cos(theta).^3-3*cos(theta));
Y31_ = -sqrt(21/(64*pi))*sin(theta).*(5*cos(theta).^2-1).*exp(1i*phi);
Y32_ = sqrt(105/(32*pi))*sin(theta).^2.*cos(theta).*exp(2i*phi);
Y33_ = -sqrt(35/(64*pi))*sin(theta).^3.*exp(3i*phi);

%-------------------------------------------------------------------------%

%The real Spherical Harmomics
% rank 0
Y00 = Y00_;

% rank 1
Y1n1 = 1i * sqrt(1/2)*(Y1n1_ + Y11_);
Y10 = Y10_;
Y11 = sqrt(1/2) * (Y1n1_ - Y11_);

% rank 2
Y2n2 = 1i * sqrt(1/2)*(Y2n2_ - Y22_);
Y2n1 = 1i * sqrt(1/2)*(Y2n1_ + Y21_);
Y20 = Y20_;
Y21 = sqrt(1/2) * (Y2n1_ - Y21_);
Y22 = sqrt(1/2)*(Y2n2_ + Y22_);

% rank 3
Y3n3 = 1i * sqrt(1/2)*(Y3n3_ + Y33_);
Y3n2 = 1i * sqrt(1/2)*(Y3n2_ - Y32_);
Y3n1 = 1i * sqrt(1/2)*(Y3n1_ + Y31_);
Y30 = Y30_;
Y31 = sqrt(1/2) * (Y3n1_ - Y31_);
Y32 = sqrt(1/2)*(Y3n2_ + Y32_);
Y33 = sqrt(1/2)*(Y3n3_ - Y33_);

%-------------------------------------------------------------------------%