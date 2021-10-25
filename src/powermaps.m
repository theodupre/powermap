classdef powermaps 
    
    properties 
        filepath;
        fs;
        start;
        stop;
            
        s;
        transform = 'afSTFT';
        hopSize = 128;
        winLen = 256;
        win;
        nfft;
        
        
        fmax = 150;
        
        S;
        f;
        K;
        t;
        tCov; 
        tMap;
        N;
        
        SHOrder;
        tDesignOrder = 15;
        
        aziSampling;
        elevSampling;
        Ynm;
        
        
        
        mapAvg = 0.8;
        mapAvgFunc;
        covAvg = 0.5;
        covAvgFunc;
        
        pwdMap;
        mvdrMap;
        cropacMap;

    end
       
    methods
        
        
        function this = powermaps(varargin)
            this = this.Set(varargin{:});
        end
        
        function this = Set(this,varargin)
            st_input = struct(varargin{:});
            
            if isfield(st_input, 'filepath')
                this.filepath = st_input.filepath;
                this.fs = audioinfo(this.filepath).SampleRate;
            else
                error("No file path provided as argument")
            end
            
            if isfield(st_input, 'start')
                this.start = st_input.start;
            else
                this.start = 0;
            end
            
            if isfield(st_input, 'stop')
                this.stop = st_input.stop;
            else
                this.stop = audioinfo(this.filepath).Duration;
            end
            
            [this.s, ~] = audioread(this.filepath, ...
                [floor(this.start*this.fs)+1, floor(this.stop*this.fs)]);
            
            
            if isfield(st_input, 'SHOrder')
                this.SHOrder = st_input.SHOrder;
                this.s = this.s(:,1:(this.SHOrder+1)^2);
            else
                this.SHOrder = sqrt(size(this.s,2)) - 1;
            end
            
            if isfield(st_input, 'transform')
                this.transform = st_input.transform;
            end
            
            if isfield(st_input, 'hopSize')
                this.hopSize = st_input.hopSize;
            end
            
            if isfield(st_input, 'winLen')
                this.winLen = st_input.winLen;
            else
                this.winLen = this.hopSize*4;
            end
            this.win = hann(this.winLen);
            this.nfft = this.winLen*2;
                        
            if isfield(st_input, 'nfft')
                this.nfft = st_input.nfft;
            end
            
            if isfield(st_input, 'fmax')
                this.fmax = st_input.fmax;
            end
            
            this.tCov = 8;
            
            [this.S, this.f, this.t] = timeFrequencyTransform(this);
            [this.K, this.N,~] = size(this.S);
            this.tMap = this.t(1:this.tCov:this.N);
            
            if isfield(st_input, 'tDesignOrder')
                this.tDesignOrder = st_input.tDesignOrder;
            end
            
            % Init sampling points and SH basis on these points
            sphSamplingCartCoor = getTdesign(this.tDesignOrder); % get Sampling from t-design algorithm
            [this.aziSampling, this.elevSampling, ~] = ...
                cart2sph(sphSamplingCartCoor(:,1), sphSamplingCartCoor(:,2), sphSamplingCartCoor(:,3)); % Convert cartesian coordinate to spherical coordinate
            this.Ynm = getSH(this.SHOrder, [this.aziSampling, pi/2 - this.elevSampling], 'real'); % get SH basis values at each sampling point, elevation coordinate is provided in inclination (pi/2 - elevation)
            
            if isfield(st_input, 'mapAvg')
                this.mapAvg = st_input.mapAvg;
            end
            this.mapAvgFunc = @(Fn, Fn_1) (1 - this.mapAvg)*Fn + this.mapAvg*Fn_1;
            
            if isfield(st_input, 'covAvg')
                this.covAvg = st_input.covAvg;
            end
            this.covAvgFunc = @(Fn, Fn_1) (1 - this.covAvg)*Fn + this.covAvg*Fn_1; 
            
            this.pwdMap = zeros(this.K, floor(this.N/this.tCov), size(this.Ynm,1), 'single');
            this.mvdrMap = zeros(this.K, floor(this.N/this.tCov), size(this.Ynm,1), 'single');
            this.cropacMap = zeros(this.K, floor(this.N/this.tCov), size(this.Ynm,1), 'single');          
        end
        function this = processMap(this)
           
            PnmCov_n = zeros(this.K,(this.SHOrder+1)^2,(this.SHOrder+1)^2, 'single');

            for n = 1:floor(this.N)/this.tCov
                tic
                for k = 1:this.K
                    % Compute covariance and smooth with previous cov
                   PnmCov_n(k,:,:) = this.covAvgFunc(covPnm_kn(this, k, (n-1)*this.tCov+1, 10), squeeze(PnmCov_n(k,:,:)));

                   % Compute PWD map
                   this.pwdMap(k,n,:) = pwdknMap(this, squeeze(PnmCov_n(k,:,:)), this.Ynm); 
                   
                   % Compute MVDR and CroPaC map
                   if trace(squeeze(PnmCov_n(k,:,:))) > 1e-8
                       % MVDR Weigths and Map
                       wMVDR_kn = mvdrWeigths(this, squeeze(PnmCov_n(k,:,:)));
                       this.mvdrMap(k,n,:) = pwdknMap(this, squeeze(PnmCov_n(k,:,:)), wMVDR_kn);
                       
                       %CroPaC Weigths and Map
                       wCroPaC_kn = croPaCWeigths(this, squeeze(PnmCov_n(k,:,:)), wMVDR_kn, k, n);
                       this.cropacMap(k,n,:) = pwdknMap(this, squeeze(PnmCov_n(k,:,:)), wCroPaC_kn);
                       
                   end
                end
                fprintf('Frame %d out of %d, processed in %f s\n', n, floor(this.N)/this.tCov, toc)
            end
        end
        
        function pwd_kn = pwdknMap(this, PnmCov_kn, w)
            % Steering matrix 
            stMat = 4*pi/length(PnmCov_kn)*w; 
            pwd_kn = real(diag(stMat*PnmCov_kn*stMat'));
        end
           
        function wMVDRk = mvdrWeigths(this, PnmkCov)
   
            wMVDRk = zeros(size(this.Ynm), 'single');
            % Apply diagonal loading
            PnmkCov = PnmkCov + 8*trace(PnmkCov)/length(PnmkCov)*eye(size(PnmkCov));
            % Option for inversion solver : Cov matrix is symmetric
            opts.SYM = true;
            for i = 1:length(this.Ynm)
                % Steering vector
                stVec = this.Ynm(i,:).';
                % Numerator of mvdr weights : Cxx^-1 * Ynm(i)
                numW = linsolve(PnmkCov, stVec, opts);
                % Denominator of mvdr weights : Ynm(i)'*Cxx^-1*Ymn(i)
                denW = stVec'*numW;
                % MVDR weights
                wMVDRk(i,:) = numW./denW;
            end
            
        end
        
        function wCroPaC_kn = croPaCWeigths(this, PnmCov_n, wMVDR_kn, k, n)
            
            wCroPaC_kn = zeros(size(wMVDR_kn), 'single');
            A = zeros(2,length(PnmCov_n), 'single');
            b = [1;0];
            lambda = 0.1;
            
            opts.SYM = true;
            % First half of the cross spectrum
            Cx_Y = this.Ynm*PnmCov_n;
            % Apply diagonal loading
            PnmCov_n = PnmCov_n + 8*trace(PnmCov_n)/length(PnmCov_n)*eye(size(PnmCov_n));
            
            for i = 1:length(this.Ynm)
                stVec = this.Ynm(i,:)';
                A(1,:) = stVec;
                A(2,:) = stVec.*diag(PnmCov_n);
                
                % w0 = (Cx^-1 * A^H) * (A * Cx^-1 * A^H)^-1 * b
                invCx_AH = linsolve(PnmCov_n, A', opts);
                A_invCx_AH = A*invCx_AH + 1e-2* trace(A*invCx_AH)/length(A*invCx_AH)*eye(size(A*invCx_AH));
                invA_invCx_AH_b = mrdivide(invCx_AH, A_invCx_AH);
                w0 = invA_invCx_AH_b*b;
                
                % Cross spectrum between static beam Y, and adaptive beam w0 (LCMV)
                Y_Cx_w0 = Cx_Y(i,:)*w0;
                g = min(abs(Y_Cx_w0), this.mvdrMap(k,n,i));
                G = sqrt(g/(this.mvdrMap(k,n,i)+eps));
                G = max(lambda, G); % spectral floor
                
                wCroPaC_kn(i,:) = wMVDR_kn(i,:)*G;
            end
        end
        
        function SknCov = covPnm_kn(this, k, n, nFrames)
            nIdx = n+(1:nFrames); % time frames to estimate the covariance
            nIdx = nIdx(nIdx<=size(this.S,2)); % make sure N doesn't contain indexes greater than length(S)
            SknCov = squeeze(this.S(k,nIdx,:))'*squeeze(this.S(k,nIdx,:));
        end
        
        function [S, f, t] = timeFrequencyTransform(this)
           
            switch this.transform
                case "STFT"
                    if contains(version, 'R2020b')
                        [S, f, t] = stft(this.s, this.fs, 'Window', this.win, ...
                            'OverlapLength', this.winLen - this.hopSize, ...
                            'FFTLength', this.nfft, 'FrequencyRange', 'onesided');
                    else
                        [S, f, t] = stft(this.s, this.fs, 'Window', this.win, ...
                            'OverlapLength', this.winLen - this.hopSize, ...
                            'FFTLength', this.nfft, 'Centered', false);
                    end
                    t = t + this.start;
                case "afSTFT"
                    f = afSTFT(this.hopSize,(this.SHOrder+1)^2,1, 'hybrid'); % init
                    f = f*this.fs/2;
                    S = afSTFT(this.s) + eps;
                    t = (0:this.hopSize:size(this.s,1)-1)'/this.fs + this.start;
                    afSTFT();
                otherwise
                    f = afSTFT(this.hopSize,(this.SHOrder+1)^2,1, 'hybrid'); % init 
            end
           
            f = f(f < this.fmax);
            S = S(1:length(f),:,:);
            
        end
              
    end
    
end