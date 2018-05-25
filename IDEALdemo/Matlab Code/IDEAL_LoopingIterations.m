function [B0_Final, Complex_Rho_Final] =...
    IDEAL_LoopingIterations...
    (Rho_Real_New, Rho_Imag_New, DeltaPhi, InitialDeconv, cnj, dnj)

%This function performs looping IDEAL reconstruction

%Inputs;

%Outputs:

%% James Grist and Christian Mariager 25/05/18
Iter = 1;

%%
while abs(DeltaPhi) > 1 && Iter < 500 %Keep the field map to less than 2Hz. --> change this to run longer!
    
    sizeA = 1;
    
    for x = 1:size(frequencies,2)
        
        A(:,sizeA:sizeA+1) = [cnj(:,x) ,  -dnj(:,x);  dnj(:,x) , cnj(:,x)];    % calculate A again in order to find next iteration of B
        
        sizeA = size(A,2) + 1;
        
    end
    
    newrhoreal = Rho_Real_New; % use updated rho to determine g matrix
   
    newrhoimag = Rho_Imag_New;
    
    % calculate g matrix for and create B matrix
    for x = 1:size(frequencies,2)
        
        sumrealsum(:,x) = (-newrhoreal(x).*dnj(:,x) - newrhoimag(x).*cnj(:,x)); % real sum over M species for g coefficients
        
        sumimagsum(:,x) = (newrhoreal(x).*cnj(:,x) - newrhoimag(x).*dnj(:,x));  % imag sum over M species for g coefficients
    
    end
    
    sumreal = 2.*pi*t.*sum(sumrealsum,2);                                   % real components of g matrix
    
    sumimag = 2.*pi*t.*sum(sumimagsum,2);                                   % imag components of g matrix
    
    gmatrix = [sumreal; sumimag];
    
    B = [gmatrix, A];
    
    % calculate Moore-Penrose pseudoinverse of B using the latest image data and new g coefficients (contained in B) and find y for each pixel
    y = pinv(B'*B)*B'*InitialDeconv; % first iteration of y
    
    DeltaPhi = y(1);   % latest iteration of field map offset
    
    DeltaPhiNew = y(1); % latest iteration of field map offset (is this used?)
    
    DeltaRhoImag = y(3:2:end,1);   % latest iteration imag components of initial rho estimate using B
    
    DeltaRhoReal = y(2:2:end,1);   % latest iteration real components of initial rho estimate using B
    
    SumDeltaRhoImag = SumDeltaRhoImag+ y(3:2:end,1); % sum all iterations of delta rho
    
    SumDeltaRhoReal = SumDeltaRhoReal+ y(2:2:end,1); % sum all iterations of delta rho
    
    B0_Final(l1,l2,counter) = B0_Final(l1,l2,counter) + DeltaPhi;    % update field map with the latest iteration field map offset
    
    Rho_Real_New = DeltaRhoReal + newrhoreal; % update rho inside loop
    Rho_Imag_New = DeltaRhoImag + newrhoimag;
    
    % Update frequencies using the latest field map
    delf = delf+DeltaPhi;
    
    % Find coefficients again using the latest field map
    for x = 1:size(frequencies,2)
        for y = 1:ne
            cnj(y,x) = cos(2.*pi.*delf(x).*t(y)); % updated coefficients for next loop round
            dnj(y,x) = sin(2.*pi.*delf(x).*t(y));
        end
    end
    
    Iter = Iter + 1;
end % iteartive while loop


% smooth final B0 map
B0_Final(l1,l2,counter) = B0_Final(l1,l2,counter);  % apply a low pass smoothing filter here!

% now calculate the signal from final field map and do final IDEAL recon
ddrealnew = 2.*pi.*B0_Final(l1,l2,counter).*t.*sum(sumrealsum,2); % calculate first term of eq. B3 (real)
ddimagnew = 2.*pi.*B0_Final(l1,l2,counter).*t.*sum(sumimagsum,2); % calculate first term of eq. B4 (imag)

% calculate second and third terms of eq. B3 and B4
for x = 1:size(frequencies,2)
    B3_1(:,x) = SumDeltaRhoReal(x).*cnj(:,x); % calculate second term of eq. B3
    B3_2(:,x) = SumDeltaRhoImag(x).*dnj(:,x); % calculate third term of eq. B3
    B4_1(:,x) = SumDeltaRhoReal(x).*dnj(:,x) + SumDeltaRhoImag(x).*cnj(:,x);  % calculate second term of eq. B4
end

B3_12sum = sum(B3_1,2) - sum(B3_2,2); % sum term two and three of eq. B3
B4_1sum = sum(B4_1,2); % sum term two of eq. B4

% calculate eq. B3 and B4
ddrealnew = ddrealnew + B3_12sum; % updated real signal (eq. B3)
ddimagnew = ddimagnew + B4_1sum; % updated imag signal (eq. B4)

InitialDeconv = [ddrealnew; ddimagnew];

% pinv_A = pinv(A'*A)*A'; % Moore-Penrose pseudoinverse of A (if it should eb updated and used in the final IDEAL calculation)
rho_estimate = pinv_A*InitialDeconv; % is this correct or should I use pinv_A where the frequencies are updated pr. the field map (no: that is even more ugly)?

Complex_Rho_Final(l1,l2,:,counter) = complex(rho_estimate(1:2:end),rho_estimate(2:2:end)); % final species estimats (k-space)


end