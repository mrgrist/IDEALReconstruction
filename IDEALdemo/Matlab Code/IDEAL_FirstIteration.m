function [Rho_Real_New,Rho_Imag_New,DeltaPhi] = IDEAL_FirstIteration...
    (fidne,te,ne,interTE,delay2,frequencies)

%First iteration of IDEAL reconstruction, for one pizel


%%
%Inputs:
%Echo data (time,echos)
%te = echo time ()
%ne = number of echoes 
%interTE = Echo spacing ()
%delay2 = 
%frequencies = Reconstruction frequencies (Hz)

%%
%Outputs:


%% James Grist and Christian Mariager 25/05/2018

%%

% IDEAL recon

counter = 1;

%Pre-allocate memory
cnj = zeros(size(frequencies,2));

dnj = zeros(1,ne);

t = linspace(te,te+(ne-1).*(interTE+delay2),ne)'; % timepoints for the ne echoes
    
for k = mtssfpindex(1):mtssfpindex(2) % frame loop

    % estimate A matrix
    for x = 1:size(frequencies,2)
       
        for y = 1:ne
           
            cnj(x,y) = cos(2.*pi.*frequencies(x).*t(y)); % calculate coefficients for A
            
            dnj(y,x) = sin(2.*pi.*frequencies(x).*t(y));
        
        end
        
    end
    
    sizeA = 1;
    
    for x = 1:size(frequencies,2)
        
        A(:,sizeA:sizeA+1) = [cnj(:,x) ,  -dnj(:,x);  dnj(:,x) , cnj(:,x)];
        
        sizeA = size(A,2) + 1;
    
    end
    
    pinv_A = pinv(A'*A)*A'; % Moore-Penrose pseudoinverse of A
    
    InitialDeconv = (squeeze([squeeze(real(fidne(:,counter)));squeeze(imag(fidne(:,counter)))])); % sort echo in real and imag components to match pinv_A
    
    rho_estimate = pinv_A*InitialDeconv; % rho for time evolution of one pixel, structure: [real_1;imag_1;real_2;imag_2; ...]
    
    newrhoreal(1:size(frequencies,2),1) = rho_estimate(1:2:end);  % real components of initial rho estimate
    newrhoimag(1:size(frequencies,2),1) = rho_estimate(2:2:end);  % imag components of initial rho estimate
    
    % calculate g matrix for and create B matrix
    for x = 1:size(frequencies,2)
        
        sumrealsum(:,x) = (-newrhoreal(x).*dnj(:,x) - newrhoimag(x).*cnj(:,x)); % real sum over M species for g coefficients
        
        sumimagsum(:,x) = (newrhoreal(x).*cnj(:,x) - newrhoimag(x).*dnj(:,x));  % imag sum over M species for g coefficients
    
    end
    
    sumreal = 2.*pi*t.*sum(sumrealsum,2);                                   % real components of g matrix
    
    sumimag = 2.*pi*t.*sum(sumimagsum,2);                                   % imag components of g matrix
    
    gmatrix = [sumreal; sumimag];
    
    B = [gmatrix, A];
    
    % calculate Moore-Penrose pseudoinverse of B using the initial image data and new g coefficients (contained in B) and find y for each pixel
    y = pinv(B'*B)*B'*InitialDeconv; % first iteration of y
    
    DeltaPhi = y(1);   % first iteration of field map offset

    SumDeltaRhoImag = y(3:2:end,1); % sum all iterations of delta rho
    
    SumDeltaRhoReal = y(2:2:end,1); % sum all iterations of delta rho
    
    DeltaRhoImag = y(3:2:end,1);   % first iteration imag components of initial rho estimate using B
    
    DeltaRhoReal = y(2:2:end,1);   % first iteration real components of initial rho estimate using B

    Rho_Real_New = DeltaRhoReal + newrhoreal; % update rho before loop
    
    Rho_Imag_New = DeltaRhoImag + newrhoimag;
    
    % Update frequencies using the first field map
    delf = frequencies+DeltaPhi;
    
    % Find coefficients again using the first field map
    for x = 1:size(frequencies,2)
        
        for y = 1:ne
            
            cnj(y,x) = cos(2.*pi.*delf(x).*t(y)); % updated coefficients
            
            dnj(y,x) = sin(2.*pi.*delf(x).*t(y));
        
        end
        
    end


end