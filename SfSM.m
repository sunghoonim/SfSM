%
%  SfSM.m
%  SfSM: High Quality Structure from Small Motion for Rolling Shutter Cameras
% 
%
%  Created by Sunghoon Im on 2017. 1. 25..
%  Copyright @ 2017 Sunghoon Im. All rights reserved.
%


classdef SfSM < handle
    
    properties
       srcdir;
       type_data;
       cname;
       alpha;
       InitDepth;
    end
    
    methods
        function this = SfSM(setting)
            this.srcdir = setting.srcdir;
            this.type_data = setting.type_data;
            this.cname = setting.cname;
            this.alpha = setting.alpha;
            this.InitDepth = setting.InitDepth;
        end
        
        function SMBARS(this)
            
            mkdir(fullfile(this.srcdir, this.cname));
            
            % VIdeo Read
            VID = VideoReader(fullfile(this.srcdir, [this.cname, '.', this.type_data]));
            Ni = min(get(VID, 'NumberOfFrame'),30);
            
            Img = cell(1,Ni);
            for ni = 1:Ni
                Img{ni} = read(VID,ni);
            end
            
            % Load data
            load(fullfile(this.srcdir,'params_60D'));
            K = params.IntrinsicMatrix';
            dist_radial = zeros(1,3); dist_tangential = zeros(1,2); 
            dist_radial(1:length(params.RadialDistortion)) = params.RadialDistortion;
            dist_tangential(1:length(params.TangentialDistortion)) = params.TangentialDistortion;
            dist = [dist_radial dist_tangential];

            % Feature Tracking
            [~,v,u] = harris(rgb2gray(Img{1}),2,2,3,0);
            feats = [u,v]';
            Nf = size(feats, 2);
            Nf_prev = Nf;

            tracker = vision.PointTracker('MaxBidirectionalError', 0.1);
            initialize(tracker, feats', im2double(Img{1}));
            valid = ones(Nf, 1);
            for ni = 2:Ni
               ImgT = im2double(Img{ni});
               [feats_i, valid_i] = step(tracker, ImgT);
               feats(1+(ni-1)*2 : 2+(ni-1)*2, :) = feats_i';
               valid = valid + valid_i;
            end

            falsevalid = valid < Ni;
            feats(:, falsevalid) = [];
            Nf = size(feats,2);
            fprintf('%d number of features are remaind from %d.\n', Nf, Nf_prev);

            % 3D reconstruction
            paramters.h = F_vdist(feats, K, dist, size(Img{1}));
            paramters.K = K;
            paramters.feats = feats;
            paramters.rol = this.alpha;

            X = K\[feats(1,:); feats(2,:); ones(1,Nf)]*this.InitDepth;    
            x0 = [zeros(6*(Ni-1),1); X(:)];
             [x, ~]=HyowonLM(@Bundlecost,x0,paramters); % You can use lsqnonlin function as optimization method
%             OPTIONS = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',  'Display', 'iter','Jacobian', 'on');
%             x = lsqnonlin(@(x)Bundlecost(x,paramters),x0,[],[],OPTIONS); % Invoke optimizer
            
            Ropt = [0; 0; 0; x(1:3*(Ni-1))];
            Topt = [0; 0; 0; x(3*(Ni-1)+1:6*(Ni-1))];

            CamP = zeros(6,Ni);
            cmap=jet(Ni);
            CamP(4:6,:) = cmap';
            for ni = 1 : Ni
                thetax = Ropt(1+3*(ni-1));
                thetay = Ropt(2+3*(ni-1));
                thetaz = Ropt(3+3*(ni-1));
                ttemp = Topt(3*(ni-1)+1:3*(ni-1)+3);
                Rtemp = [1 -thetaz thetay; thetaz 1 -thetax; -thetay thetax 1];
                CamP(1:3,ni) = -Rtemp'*ttemp;
            end
            Xopt = reshape(x(6*(Ni-1)+1:end),3,[]);
            if ( nnz(Xopt(3,:)>0) < nnz(Xopt(3,:)<0) )
                Xopt = -Xopt;
                Topt = -Topt;
            end

            RTopt = [reshape(Ropt,[3,Ni]); reshape(Topt,[3,Ni])];
    
            % Feature color 
            color_feats = zeros(3,Nf);
            for ch = 1:3
                color_feats(ch,:) = interp2(im2double(Img{1}(:,:,ch)),feats(1,:),feats(2,:));
            end
            
            Xfinal = [Xopt; color_feats];

            % Save
            fn_saveply( fullfile(this.srcdir, this.cname, [this.cname, '_sparse.ply']), [Xfinal CamP] );

            % Propagation
            min_depth = min(Xfinal(3,:))/1.1;
            xopt = [feats(1:2,:); Xfinal(3,:)];
            [~,~,depth_init] = fn_depthprop(Img{1},xopt);
            
            depth_init_ = min_depth./depth_init;
            depth_init_ = min(max(depth_init_,0),1);
            imwrite(depth_init_,fullfile(this.srcdir, this.cname, [this.cname, '_depth.png']));
        end
    end
    
end


function h = F_vdist(feats, K, distort, Imgsize)

    Nf=size(feats,2);
    Ni=size(feats,1)/2;
    u_feat=feats(1:2:end,:);
    v_feat=feats(2:2:end,:);

    k1=distort(1); k2=distort(2); k3=distort(3);
    p1=distort(4); p2=distort(5);

    uv1_undist=[u_feat(:)';v_feat(:)';ones(1,Ni*Nf)];
    nxny1=K\uv1_undist;
    nx_undist=nxny1(1,:);
    ny_undist=nxny1(2,:);
    r_undist=sqrt(nx_undist.^2+ny_undist.^2);
    nx_dist=nx_undist.*(1 + k1*r_undist.^2 + k2*r_undist.^4 + k3*r_undist.^6) + (2*p1*nx_undist.*ny_undist + p2*(r_undist.^2 + 2*nx_undist.^2));
    ny_dist=ny_undist.*(1 + k1*r_undist.^2 + k2*r_undist.^4 + k3*r_undist.^6) + (p1*(r_undist.^2 + 2*ny_undist.^2) + 2*p2*nx_undist.*ny_undist);

    nxny1_dist=[nx_dist;ny_dist;ones(1,Ni*Nf)];
    uv1_dist=K*nxny1_dist;
    u_dist=reshape(uv1_dist(1,:),size(u_feat));
    v_dist=reshape(uv1_dist(2,:),size(v_feat));
    h= ((v_dist-1).*Imgsize(2)+u_dist-1)/(Imgsize(1)*Imgsize(2));
    
end


function fn_saveply(filename, X)
    out=fopen(filename,'w');
    fprintf(out,'ply\n');
    fprintf(out,'format ascii 1.0\n');
    fprintf(out,'element vertex %d\n',size(X,2));
    fprintf(out,'property float x\n');
    fprintf(out,'property float y\n');
    fprintf(out,'property float z\n');
    fprintf(out,'property uchar diffuse_red\n');
    fprintf(out,'property uchar diffuse_green\n');
    fprintf(out,'property uchar diffuse_blue\n');
    fprintf(out,'end_header\n');
    for i=1:size(X,2)
        fprintf(out,'%f %f %f %d %d %d\n',[X(1,i),X(2,i),X(3,i),min(round(X(4,i)*255),255),min(round(X(5,i)*255),255),min(round(X(6,i)*255),255)]);
    end
    fclose(out);
end
