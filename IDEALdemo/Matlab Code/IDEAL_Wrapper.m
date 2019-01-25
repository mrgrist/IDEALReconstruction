function IDEAL_Wrapper(TestFlag)

%%

%This is a wrapper function to run multi-slice, multi-phase cartesian IDEAL
%reconstruction.

%% James Grist and Christian Mariager 25/05/2018

%%
%Test flag = 0,1. If 0, no tests will run.

%%

if nargin < 1
    TestFlag = 0;
end

%% Load in data function

if TestFlag == 1

[RE,IM,NP.NB,NT,HDR] = Test_Load_IDEAL;

end

%% Standard reconstruction functions (EPI, phasing  etc)


%% IDEAL reconstruction section

nx = size();
ny = size();
nz = size();
nmetabolites = size();
ntime = size();


%Pre-allocate memory for deconvolved images, and B0 map
Images = zeros(nx,ny,nz,nmetabolites,ntime);

B0Map = zeros(nx,ny,nz,ntime);

Tolerance = 40; %(Hz) If initial B0 map is greater than this, will not iterate.

%Run this over every voxel, slice, and phase

for x = 1:nx
    
    for y = 1:ny
        
        for z = 1:nz
            
            for numtime = 1:ntime
                
                [Rho_Real_New,Rho_Imag_New,DeltaPhi,InitialDeconv, cnj, dnj]...
                 = IDEAL_FirstIteration...
                 (fidne,te,ne,interTE,delay2,frequencies);
                
                %Put a check in here to discard voxels with too big an error
                if abs(DeltaPhi) < Tolerance
                    
                    [B0_Final, Complex_Rho_Final] =...
                    IDEAL_LoopingIterations...
                    (Rho_Real_New, Rho_Imag_New, DeltaPhi, InitialDeconv,...
                    cnj, dnj);
                    
                    B0Map = B0_Final;
                    
                    Images(nx,ny,nz,:,ntime) = ...
                        Complex_Rho_Final;
                    
                else
                    
                    B0Map(nx,ny,nz,ntime) = 0;
                    
                    Images(nx,ny,nz,:,ntime) = 0;
                    
                end
                
            end
            
        end
        
    end
    
end

%% DICOM saving




end
