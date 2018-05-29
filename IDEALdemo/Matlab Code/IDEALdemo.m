% IDEALdemo.m
% Demo of the IDEAL reconstruction, using an example dataset acquired in
% the mouse liver at 9.4T using hyperpolairzed pyruvate, as well as a
% phantom dataset containing three phantoms (bicarbonate, acetate and urea).
% Place the folder containing the example data and .m files in your matlab
% path.

% Change directory to script folder
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
cd ..
% Clear workspace
clear all; close all; clc;

%%% import and use phantom data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DataPath = 'Phantom_data/s_2017012401/'; %Add check for mac or windows here
% PrintAllFigs = 1; % 1=static data/one frame, 2=dynamic data/time series
% HannFilter = 0; % use a Hanning filter on the k-space data before processing? 1=use filter, 0=do not use filter
% experimentPath = DataPath;
% C13mtssfp = 'mtssfp_';
% delay1 =   0.000000;    % [s] - correction for IDEAl using odd/even echoes
% delay2 =   0.000000;    % [s] - correction for IDEAl using all echoes
% mtssfpindex = [1 1];    % Specify first and last mtssfp data set to use for mtssfp IDEAL (index is number of files in folder)
% manual_slopeval=-0.3;  % set a manual value for phase correction
% offsetval = 0.0.*pi;    % [radians] - Linear phase correction
% interTE = 0.000608;     % [s] - delta TE
% frequencies= [-1.8404    0.0346    0.2788].*1e3; % acquired from spectral data
% C13names = {'Bicarbonate','Urea','Acetate'};     % names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% import and use in-vivo data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataPath = 'Mouse_data/s_2018022201_13C/';
frequencies=fliplr([-0.889   -0.503   -0.229    0.342].*1e3);% 0.5 0.8].*1e3;
mtssfpindex = [12 20];    % Specify first and last mtssfp data set to use for mtssfp IDEAL (index is number of files in folder)
PrintAllFigs = 2; % 1=static data / one frame, 2=dynamic data/time series
HannFilter =1;
experimentPath = DataPath;
C13mtssfp = 'mtssfp_';
delay1 =   0.000000; % [s] - correction for IDEAl using odd/even echoes
delay2 =   0.000000; % [s] - correction for IDEAl using all echoes
manual_slopeval=-0.3;
offsetval = 0.*pi;      % [radians] - Linear phase correction
interTE = 0.000708;     % [s] - delta TE
C13names = {'Lactate','Pyruvate-hy','Alanine','Pyruvate'};%,'Urea','Bicarbonate'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort image data and flip every other k-space line
mtssfpFiles = dir(strcat([experimentPath '/' C13mtssfp '*.fid']));
counter = 1;
for k=mtssfpindex(1):mtssfpindex(2)        
    mtssfpFilename = mtssfpFiles(k).name(1:end-4);
    nv = queryPP([experimentPath mtssfpFilename,'.fid'],'nv');
    ne = queryPP([experimentPath mtssfpFilename,'.fid'],'ne');
    te = queryPP([experimentPath mtssfpFilename,'.fid'],'te');
    sw = queryPP([experimentPath mtssfpFilename,'.fid'],'sw');
    pad = queryPP([experimentPath mtssfpFilename,'.fid'],'pad');
    tof = queryPP([experimentPath mtssfpFilename,'.fid'],'tof');
    [RE,IM,NP,NB,NT,HDR] = load_fid([experimentPath mtssfpFilename]);

    mtssfpfid = complex(RE,IM);      
    fid = mtssfpfid(:,1:nv*ne);   % extract each fid from raw data
    mtssfpfid(:,1:nv*ne) = [];    % remove extracted fid from raw data
    for i=1:ne                    % loop for flipping every other k-space line
        fid_org(:,:,i,counter) = fid(:,i:ne:end);
        if mod(i,2) == 1
           fidne(:,:,i,counter) = fid(:,i:ne:end);              % don't flip
        else
           fidne(1:end,:,i,counter) = flip(fid(:,i:ne:end),1);  % flip
        end
        uncorrected_mtssfp(:,:,i,counter) = fftshift(fft2(ifftshift(fid_org(:,:,i,counter)))); % unflipped echos
        flipped_mtssfp(:,:,i,counter) = fftshift(fft2(ifftshift(fidne(:,:,i,counter))));       % flipped echos
    end
    
       Profiles = abs(squeeze(sum(fidne(:,:,:,counter),2)));       % find a value to use for phase correction (make this better)
    ZipperBefore(:,:,counter) = Profiles;
    for kk=1:size(Profiles,2)
    [Pmax(kk),Ploc(kk)] = max(Profiles(:,kk));
    end
    difference = abs(mean(Ploc(1:2:end)) - mean(Ploc(2:2:end)));
    alpha = (difference./2).*2.*pi ./ (NP);
    slopeval  = -alpha./2.*pi
    save_slopeval(counter) = slopeval;
    counter = counter +1;
end


% Apply auto or manual phase correction?
if manual_slopeval==0
    slopeval = mean(save_slopeval);
else
    slopeval = manual_slopeval;
end

slopeval
% Apply phase correction
counter = 1;
for k=mtssfpindex(1):mtssfpindex(2)
    
     % use a hanning filter?
     if HannFilter ==1;
        hann = repmat(hanning(NP),1,nv);
     else
        hann=1;
     end
     
     % apply hanning filter
     for kk=1:ne
         fidne(:,:,kk,counter)=fidne(:,:,kk,counter).*hann;
     end
     
     % phase correction
     fidne_x_ky = fftshift(fft(fidne(:,:,:,counter)),1);
     offset = ones(size(fidne_x_ky,1),1).*offsetval;
     slope = ones(size(fidne_x_ky,1),1).*slopeval;
     phasecorr = zeros(size(fidne_x_ky,1),size(fidne_x_ky,2));  
   
     for d = 1:size(fidne_x_ky,1)
         arg = ([floor(-size(fidne_x_ky,1)/2+1):ceil(size(fidne_x_ky,1)/2)]...
                 *slope(d)+offset(d));        
         phasecorr = complex(cos(arg),sin(arg));
     end
     
     phasecorr = repmat(phasecorr,nv,1);
     for f = 1:ne
         if mod(f,2) == 1
            fidne_x_ky(:,:,f) = fidne_x_ky(:,:,f);
         elseif mod(f,2) == 0    
            fidne_x_ky(:,:,f) = fidne_x_ky(:,:,f) .* phasecorr';
         end
     end
       
     fidneCorr = ifft(fftshift(fidne_x_ky,1));

          % visually check phase correction         
     ProfilesAfter = abs(squeeze(sum(fidneCorr,2)));
     ZipperAfter(:,:,counter) = ProfilesAfter;
    
     fig2=figure(2);
     subplot(1,2,1)
     plot(ZipperBefore(:,1:2:end,counter),'r')
     hold on
     plot(ZipperBefore(:,2:2:end,counter),'b')
     vline([NP./2],'--k')
     xlim([1 NP])
     xlabel('kx samples')
     ylabel('Intensity')
     title('Uncorrected echo profiles')
     subplot(1,2,2)
     plot(ProfilesAfter(:,1:2:end),'r')
     hold on
     plot(ProfilesAfter(:,2:2:end),'b')
     xlim([1 NP])
     vline([NP./2],'--k')
     ylabel('Intensity')
     xlabel('kx samples')
     title('Corrected echo profiles')
     %pause(2)
     %clf(fig2)
     
     % use corrected fid's
     fidne(:,:,:,counter)=fidneCorr;
     pause(0.5)
     counter = counter +1;
end


% IDEAL recon
matrix = 64; % recon matrix
allmap_img_display=[];
allmap_img_display_initial=[];
counter = 1;
tic
%frequencies = [10000 11000 12000 13000 14000 15000];
for k=mtssfpindex(1):mtssfpindex(2) % frame loop
    clear A
    
    InitialB0Map = zeros(size(fidne));                % initial field map
    t = linspace(te,te+(ne-1).*(interTE+delay2),ne)'; % timepoints for the ne echoes

    % estimate A matrix
    for x = 1:size(frequencies,2)
        for y = 1:ne   
            cnj(y,x) = cos(2.*pi.*frequencies(x).*t(y)); % calculate coefficients for A
            dnj(y,x) = sin(2.*pi.*frequencies(x).*t(y));    
        end   
    end
    
    sizeA = 1;
    for x = 1:size(frequencies,2)    
        A(:,sizeA:sizeA+1) = [cnj(:,x) ,  -dnj(:,x);  dnj(:,x) , cnj(:,x)];    
        sizeA = size(A,2) + 1;
    end 
    
    pinv_A = pinv(A'*A)*A'; % Moore-Penrose pseudoinverse of A
    
    % pixel'vise reconstruction
    for l1=1:NP   % for each pixel do..
        for l2=1:nv   % for each pixel do..
          
          InitialDeconv = (squeeze([squeeze(real(fidne(l1,l2,:,counter)));squeeze(imag(fidne(l1,l2,:,counter)))])); % sort echo in real and imag components to match pinv_A
          
          rho_estimate = pinv_A*InitialDeconv; % rho for time evolution of one pixel, structure: [real_1;imag_1;real_2;imag_2; ...]
          
          complex_rho_initial(l1,l2,:,counter) = complex(rho_estimate(1:2:end),rho_estimate(2:2:end));  % initial IDEAL maps after one iteration and Phi_0=0
          
          
          newrhoreal(1:size(frequencies,2),1) = rho_estimate(1:2:end);  % real components of initial rho estimate
          newrhoimag(1:size(frequencies,2),1) = rho_estimate(2:2:end);  % imag components of initial rho estimate
          
%           for j = 1:ne
%               for x = 1:size(frequencies,2)
%                   sumreal(j,x) = 2*pi*t(j).*(-newrhoreal(x,1)*dnj(j,x) - newrhoimag(x,1)*cnj(j,x)); % real components of g matrix
%                   sumimag(j,x) = 2*pi*t(j).*(newrhoreal(x,1)*cnj(j,x)- newrhoimag(x,1)*dnj(j,x));   % imag components of g matrix
%               end
%           end
        
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
                
         DeltaPhiNew = y(1); % first iteration of field map offset (is this used?)
         
         SumDeltaRhoImag = y(3:2:end,1); % sum all iterations of delta rho
         
         SumDeltaRhoReal = y(2:2:end,1); % sum all iterations of delta rho
                
         DeltaRhoImag = y(3:2:end,1);   % first iteration imag components of initial rho estimate using B
                
         DeltaRhoReal = y(2:2:end,1);   % first iteration real components of initial rho estimate using B
         
         B0_Final(l1,l2,counter) = DeltaPhi;    % update field map with the first iteration field map offset
 
% calculate signal
%          ddrealnew = 2.*pi.*DeltaPhi.*t.*sum(sumrealsum,2); % calculate first term of eq. B3 (real)
%          ddimagnew = 2.*pi.*DeltaPhi.*t.*sum(sumimagsum,2); % calculate first term of eq. B4 (imag)
%          
%          % calculate second and third terms of eq. B3 and B4
%           for x = 1:size(frequencies,2)
%               B3_1(:,x) = DeltaRhoReal(x).*cnj(:,x); % calculate second term of eq. B3 
%               B3_2(:,x) = DeltaRhoImag(x).*dnj(:,x); % calculate third term of eq. B3 
%               B4_1(:,x) = DeltaRhoReal(x).*dnj(:,x) + DeltaRhoImag(x).*cnj(:,x);  % calculate second term of eq. B4 
%           end 
%           
%           B3_12sum = sum(B3_1,2) - sum(B3_2,2); % sum term two and three of eq. B3
%           B4_1sum = sum(B4_1,2); % sum term two of eq. B4
%           
%           % calculate eq. B3 and B4
%           ddrealnew = ddrealnew + B3_12sum; % updated real signal (eq. B3)
%           ddimagnew = ddimagnew + B4_1sum; % updated imag signal (eq. B4)


          Rho_Real_New = DeltaRhoReal + newrhoreal; % update rho before loop
          Rho_Imag_New = DeltaRhoImag + newrhoimag;
          
%                 
%                 for l = 1:size(dnj,1)
%                     
%                     ddrealnew(:,l) = Rho_Real_New.*cnj(l,:)' - Rho_Imag_New.*dnj(l,:)';
%                     
%                     ddimagnew(:,l) = Rho_Real_New.*dnj(l,:)' + Rho_Imag_New.*cnj(l,:)'; % check if it is correct to use old dn and cn here?
%                     
%                     
%                 end
             
                
                % Update frequencies using the first field map
                delf = frequencies+DeltaPhi;
                
                % Find coefficients again using the first field map
                for x = 1:size(frequencies,2)
                    for y = 1:ne   
                     cnj(y,x) = cos(2.*pi.*delf(x).*t(y)); % updated coefficients
                     dnj(y,x) = sin(2.*pi.*delf(x).*t(y));    
                    end   
                end
                
%                 for n = 1:ne
%                     sreal(n,1) = sum((ddrealnew(:,n).*cnj(n,:)' - ddimagnew(:,n).*dnj(n,:)')); % check if it is correct to use new dn and cn here?
%                     simag(n,1) = 1i.*(sum((ddrealnew(:,n).*dnj(n,:)' + ddimagnew(:,n).*cnj(n,:)')));
%                 end

                
                Iter = 1;
                
                %%%% START %%%%
                % iteartive part of code
        while abs(DeltaPhi) > 1 && Iter < 500 %Keep the field map to less than 2Hz. --> change this to run longer!
              clear A B gmatrix
                    
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
                    
        newrhoreal(1:size(frequencies,2),1) = rho_estimate(1:2:end);
        newrhoimag(1:size(frequencies,2),1) = rho_estimate(2:2:end);
        complex_rho_final(l1,l2,:,counter) = complex(rho_estimate(1:2:end),rho_estimate(2:2:end)); % final species estimats (k-space)
                
        end % pixel loop
    end % pixel loop
    
    %%%%% calculate images
    allmap_montage(:,:,:,counter)=fliplr(reshape(complex_rho_final(:,:,:,counter),[NP nv.*length(frequencies)]));       % print IDEAL maps
    allmap_img(:,:,:,counter) = fftshift(fft(fftshift(fft(complex_rho_final(:,:,:,counter),matrix,1),1),matrix,2),2);
    allmap_img_montage(:,:,:,counter)=(reshape(allmap_img(:,:,:,counter),[matrix matrix.*length(frequencies)]));
    allmap_img_display = [allmap_img_display;squeeze(allmap_img_montage(:,:,:,counter))];
    
    allmap_montage_initial(:,:,:,counter)=fliplr(reshape(complex_rho_initial(:,:,:,counter),[NP nv.*length(frequencies)]));       % print IDEAL maps
    allmap_img_initial(:,:,:,counter) = fftshift(fft(fftshift(fft(complex_rho_initial(:,:,:,counter),matrix,1),1),matrix,2),2);
    allmap_img_montage_initial(:,:,:,counter)=(reshape(allmap_img_initial(:,:,:,counter),[matrix matrix.*length(frequencies)]));
    allmap_img_display_initial = [allmap_img_display_initial;squeeze(allmap_img_montage_initial(:,:,:,counter))];

    counter=counter+1; % repeat for next timeframe
end % frame loop
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PrintAllFigs ==1
figure(6)

% make images of phantom data
subplot(2,4,[1 2 3])
imshow(abs(allmap_img_display),[])
C13names = fliplr(C13names);
for i=1:size(frequencies,2)
    if mod(i,2)==0
        text(matrix*(i-1),-5,C13names(i))
    else
        text(matrix*(i-1),-5,C13names(i))
    end
end
colormap(gray)
colorbar
xlabel('IDEAL maps (after iterative code)')

subplot(2,4,[5 6 7])
imshow(abs([fftshift(fft2(ifftshift(complex_rho_initial(:,:,1)))),...
    fftshift(fft2(ifftshift(complex_rho_initial(:,:,2)))),...
    fftshift(fft2(ifftshift(complex_rho_initial(:,:,3))))]),[])
colormap(gray)
colorbar
warning off
xlabel('initial IDEAL maps (no iterative component)')

subplot(2,4,[4])
imshow(B0_Final,[])
colormap(gray)
colorbar
truesize
warning off
xlabel('Final k-space field map [Hz]')

subplot(2,4,[8])
imshow(sum(abs(flipped_mtssfp),3),[])
colormap(gray)
colorbar
truesize
warning off
xlabel('Sum of all echoes')

set(gcf,'Position',get(0,'Screensize')); % large figure

elseif PrintAllFigs ==2
    
% make images of in-vivo data
figure(6)
subplot(1,3,[1])
imshow(abs(allmap_img_display),[0 10000])
C13names = fliplr(C13names);
for i=1:size(frequencies,2)
    if mod(i,2)==0
        text(matrix*(i-1),-12,C13names(i))
    else
        text(matrix*(i-1),-32,C13names(i))
    end
end
colormap(jet)
xlabel('IDEAL maps (after iterative code)') 

subplot(1,3,[2])
imshow(abs(allmap_img_display_initial),[0 10000])
C13names = C13names;
for i=1:size(frequencies,2)
    if mod(i,2)==0
        text(matrix*(i-1),-12,C13names(i))
    else
        text(matrix*(i-1),-32,C13names(i))
    end
end
colormap(jet)
xlabel('initial IDEAL maps (no iterative component)') 

subplot(1,3,[3])
fieldmap = [];
for k=1:size(B0_Final,3)
    fieldmap = [fieldmap;B0_Final(:,:,k)];
end
imshow(fieldmap,[])
colormap(jet)
colorbar
xlabel('Final field maps [Hz]') 

% subplot(1,4,[4])
% plot(freq,intensity)
% hold on;
% plot(freq(loc),intensity(loc),'*r')
% ylim([0 1.1])
% title('Sum of spectra (IDEAL frequencies)')
% xlabel('frequency [Hz]')

set(gcf,'Position',get(0,'Screensize')); % large figure
end

