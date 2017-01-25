%
%  Bundlecost.m
%  SfSM: High Quality Structure from Small Motion for Rolling Shutter Cameras
% 
%
%  Created by Sunghoon Im on 2017. 1. 25..
%  Copyright @ 2017 Sunghoon Im. All rights reserved.
%


function [F,J]=Bundlecost(p, params)

    fx=params.K(1,1); a=params.K(1,2); cx=params.K(1,3);
              fy=params.K(2,2); cy=params.K(2,3);
    
	scale_t=1;
    scale_X=1;
    
    n=size(params.feats,1)/2;
    ni = n-1;
    [Rvec,Tvec,X_world] = deserialize(p,n);
    rol=params.rol; 
    
    h = params.h';
    t = h*rol;
    
    nf=size(params.feats,2);
    F=zeros(n*nf*2,1);
    X_1=[X_world; ones(1,size(X_world,2))];
    X=X_world(1,:); Y=X_world(2,:); Z=X_world(3,:);
    
    r0=zeros(6*nf, 1);
    c0=zeros(6*nf, 1);
    v0=zeros(6*nf, 1);
    
    r=zeros(30*nf, ni-1);
    c=zeros(30*nf, ni-1);
    v=zeros(30*nf, ni-1);
    
    rn=zeros(18*nf, 1);
    cn=zeros(18*nf, 1);
    vn=zeros(18*nf, 1);
    
    for i=1:n
        uv_feat=params.feats(i*2-1:i*2,:);
        Rvec_i=Rvec(:,i);
        r10=Rvec_i(1); r20=Rvec_i(2); r30=Rvec_i(3);
        
        Tvec_i=Tvec(:,i);
        t10=Tvec_i(1); t20=Tvec_i(2); t30=Tvec_i(3);
        
        tt = t(:,i)';
        if( i ~= n)
            Rvec_i1=Rvec(:,i+1);
            r11=Rvec_i1(1); r21=Rvec_i1(2); r31=Rvec_i1(3);
            Tvec_i1=Tvec(:,i+1);
            t11=Tvec_i1(1); t21=Tvec_i1(2); t31=Tvec_i1(3);
            r1 = (1-tt)*r10 + tt*r11;
            r2 = (1-tt)*r20 + tt*r21;
            r3 = (1-tt)*r30 + tt*r31;
            t1 = (1-tt)*t10 + tt*t11;
            t2 = (1-tt)*t20 + tt*t21;
            t3 = (1-tt)*t30 + tt*t31;
        else
            r1 = r10*ones(1,nf);
            r2 = r20*ones(1,nf);
            r3 = r30*ones(1,nf);
            t1 = t10*ones(1,nf);
            t2 = t20*ones(1,nf);
            t3 = t30*ones(1,nf);
        end
                
        xy_i = zeros(3, nf);
        xy_i(1,:) = sum([ones(1,nf); -r3; r2; t1].*X_1);
        xy_i(2,:) = sum([r3; ones(1,nf); -r1; t2].*X_1);
        xy_i(3,:) = sum([-r2; r1; ones(1,nf); t3].*X_1);
        nxy_i = xy_i./repmat(xy_i(3,:),[3,1]);
        uv_proj=params.K*nxy_i;
        uv_proj=uv_proj(1:2,:);
        F((i-1)*nf*2+1:i*nf*2)= uv_feat(:) - uv_proj(:);
        
        %       분모 미분 +                         분자 미분
        X_i=xy_i(1,:); Y_i=xy_i(2,:); Z_i=xy_i(3,:);
        Z_i_2=Z_i.^(-2);
        
        X_i_Z_i_2=X_i.*Z_i_2;
        Y_i_Z_i_2=Y_i.*Z_i_2;
        
        dZi_dr10 = (1-tt).*Y; %
        dx_dr10 = -dZi_dr10.*X_i_Z_i_2;
        dy_dr10 = -(1-tt).*Z./Z_i -dZi_dr10.*Y_i_Z_i_2;
        
        dZi_dr20 = -(1-tt).*X; %
        dx_dr20 = (1-tt).*Z./Z_i -dZi_dr20.*X_i_Z_i_2;
        dy_dr20 = -dZi_dr20.*Y_i_Z_i_2;
        
        dx_dr30 = -(1-tt).*Y./Z_i; %
        dy_dr30 = (1-tt).*X./Z_i;
        
        dZi_dr11 = tt.*Y; %
        dx_dr11 = -dZi_dr11.*X_i_Z_i_2;
        dy_dr11 = -tt.*Z./Z_i -dZi_dr11.*Y_i_Z_i_2;
        
        dZi_dr21 = -tt.*X; %
        dx_dr21 = tt.*Z./Z_i -dZi_dr21.*X_i_Z_i_2;
        dy_dr21 = -dZi_dr21.*Y_i_Z_i_2;
        
        dx_dr31 = -tt.*Y./Z_i; %
        dy_dr31 = tt.*X./Z_i;
        
        
        dx_dt10 = (1-tt).*scale_t./Z_i; %
        dy_dt10 = zeros(1,nf);
        
        dx_dt20 = zeros(1,nf); %
        dy_dt20 = (1-tt).*scale_t./Z_i;
        
        dx_dt30 = -(1-tt).*scale_t.*X_i_Z_i_2; %
        dy_dt30 = -(1-tt).*scale_t.*Y_i_Z_i_2;
        
        
        dx_dt11 = tt.*scale_t./Z_i; %
        dy_dt11 = zeros(1,nf);
        
        dx_dt21 = zeros(1,nf); %
        dy_dt21 = tt.*scale_t./Z_i;
        
        dx_dt31 = -tt.*scale_t.*X_i_Z_i_2; %
        dy_dt31 = -tt.*scale_t.*Y_i_Z_i_2;
        
        
        dZi_dX = -scale_X*r2; %
        dx_dX = scale_X./Z_i - dZi_dX.*X_i_Z_i_2;
        dy_dX = scale_X*r3./Z_i - dZi_dX.*Y_i_Z_i_2;
        
        dZi_dY = scale_X*r1; %
        dx_dY = scale_X*(-r3)./Z_i - dZi_dY.*X_i_Z_i_2;
        dy_dY = scale_X./Z_i - dZi_dY.*Y_i_Z_i_2;
        
        dZi_dZ = scale_X; %
        dx_dZ = scale_X*(r2)./Z_i - dZi_dZ.*X_i_Z_i_2;
        dy_dZ = scale_X*(-r1)./Z_i - dZi_dZ.*Y_i_Z_i_2;

        du_dr10 = fx*dx_dr10+a*dy_dr10;
        dv_dr10 = fy*dy_dr10;
        du_dr20 = fx*dx_dr20+a*dy_dr20;
        dv_dr20 = fy*dy_dr20;
        du_dr30 = fx*dx_dr30+a*dy_dr30;
        dv_dr30 = fy*dy_dr30;
        du_dt10 = fx*dx_dt10+a*dy_dt10;
        dv_dt10 = fy*dy_dt10;
        du_dt20 = fx*dx_dt20+a*dy_dt20;
        dv_dt20 = fy*dy_dt20;
        du_dt30 = fx*dx_dt30+a*dy_dt30;
        dv_dt30 = fy*dy_dt30;
        
        
        du_dr11 = fx*dx_dr11+a*dy_dr11;
        dv_dr11 = fy*dy_dr11;
        du_dr21 = fx*dx_dr21+a*dy_dr21;
        dv_dr21 = fy*dy_dr21;
        du_dr31 = fx*dx_dr31+a*dy_dr31;
        dv_dr31 = fy*dy_dr31;
        du_dt11 = fx*dx_dt11+a*dy_dt11;
        dv_dt11 = fy*dy_dt11;
        du_dt21 = fx*dx_dt21+a*dy_dt21;
        dv_dt21 = fy*dy_dt21;
        du_dt31 = fx*dx_dt31+a*dy_dt31;
        dv_dt31 = fy*dy_dt31;
        
        
        du_dX = fx*dx_dX+a*dy_dX;
        dv_dX = fy*dy_dX;
        du_dY = fx*dx_dY+a*dy_dY;
        dv_dY = fy*dy_dY;
        du_dZ = fx*dx_dZ+a*dy_dZ;
        dv_dZ = fy*dy_dZ;
        
        
        ui=(((i-1)*nf*2+1):2:(i*nf*2))';
        vi=(((i-1)*nf*2+2):2:(i*nf*2))';
        Xi=(1:nf)';
        
        if i==1
            r0(:) = [ ui; vi;
                      ui; vi;
                      ui; vi; ];
            c0(:) = [ 6*(ni)+3*Xi-2; 6*(ni)+3*Xi-2;
                      6*(ni)+3*Xi-1; 6*(ni)+3*Xi-1;
                      6*(ni)+3*Xi; 6*(ni)+3*Xi; ];
                  
            v0(:) = [ -du_dX'; -dv_dX';
                      -du_dY'; -dv_dY';
                      -du_dZ'; -dv_dZ'; ];
                  
        elseif i == n
            rn(:) = [ ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi;
                         ui; vi; 
                         ui; vi; ];

            cn(:) = [ (i*3-2-3)*ones(nf,1); (i*3-2-3)*ones(nf,1);
                         (i*3-1-3)*ones(nf,1); (i*3-1-3)*ones(nf,1);
                         (i*3-3)*ones(nf,1); (i*3-3)*ones(nf,1);
                         3*(ni)+(i*3-2-3)*ones(nf,1); 3*(ni)+(i*3-2-3)*ones(nf,1);
                         3*(ni)+(i*3-1-3)*ones(nf,1); 3*(ni)+(i*3-1-3)*ones(nf,1);
                         3*(ni)+(i*3-3)*ones(nf,1); 3*(ni)+(i*3-3)*ones(nf,1);
                         6*(ni)+3*Xi-2; 6*(ni)+3*Xi-2;
                         6*(ni)+3*Xi-1; 6*(ni)+3*Xi-1;
                         6*(ni)+3*Xi; 6*(ni)+3*Xi; ];
                     
            vn(:) = [ -du_dr10'; -dv_dr10';
                         -du_dr20'; -dv_dr20';
                         -du_dr30'; -dv_dr30';
                         -du_dt10'; -dv_dt10';
                         -du_dt20'; -dv_dt20';
                         -du_dt30'; -dv_dt30';
                         -du_dX'; -dv_dX';
                         -du_dY'; -dv_dY';
                         -du_dZ'; -dv_dZ'; ];
                     
        else
            r(:,i-1) = [ ui; vi;
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         
                         ui; vi;
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         ui; vi; 
                         
                         ui; vi;
                         ui; vi; 
                         ui; vi; ];

            c(:,i-1) = [ (i*3-2-3)*ones(nf,1); (i*3-2-3)*ones(nf,1);
                         (i*3-1-3)*ones(nf,1); (i*3-1-3)*ones(nf,1);
                         (i*3-3)*ones(nf,1); (i*3-3)*ones(nf,1);
                         3*(ni)+(i*3-2-3)*ones(nf,1); 3*(ni)+(i*3-2-3)*ones(nf,1);
                         3*(ni)+(i*3-1-3)*ones(nf,1); 3*(ni)+(i*3-1-3)*ones(nf,1);
                         3*(ni)+(i*3-3)*ones(nf,1); 3*(ni)+(i*3-3)*ones(nf,1);
                         
                         (i*3-2-3)*ones(nf,1)+3; (i*3-2-3)*ones(nf,1)+3;
                         (i*3-1-3)*ones(nf,1)+3; (i*3-1-3)*ones(nf,1)+3;
                         (i*3-3)*ones(nf,1)+3; (i*3-3)*ones(nf,1)+3;
                         3*(ni)+(i*3-2-3)*ones(nf,1)+3; 3*(ni)+(i*3-2-3)*ones(nf,1)+3;
                         3*(ni)+(i*3-1-3)*ones(nf,1)+3; 3*(ni)+(i*3-1-3)*ones(nf,1)+3;
                         3*(ni)+(i*3-3)*ones(nf,1)+3; 3*(ni)+(i*3-3)*ones(nf,1)+3; 
                         
                         6*(ni)+3*Xi-2; 6*(ni)+3*Xi-2;
                         6*(ni)+3*Xi-1; 6*(ni)+3*Xi-1;
                         6*(ni)+3*Xi; 6*(ni)+3*Xi; ];

            v(:,i-1) = [ -du_dr10'; -dv_dr10';
                         -du_dr20'; -dv_dr20';
                         -du_dr30'; -dv_dr30';
                         -du_dt10'; -dv_dt10';
                         -du_dt20'; -dv_dt20';
                         -du_dt30'; -dv_dt30';
                         
                         -du_dr11'; -dv_dr11';
                         -du_dr21'; -dv_dr21';
                         -du_dr31'; -dv_dr31';
                         -du_dt11'; -dv_dt11';
                         -du_dt21'; -dv_dt21';
                         -du_dt31'; -dv_dt31';
                         
                         -du_dX'; -dv_dX';
                         -du_dY'; -dv_dY';
                         -du_dZ'; -dv_dZ'; ];
        end
    end
    
    J=sparse([r0(:);r(:);rn(:)],[c0(:);c(:);cn(:)],[v0(:);v(:);vn(:)],n*nf*2,length(p));
end

function [r,t,X]=deserialize(p, n)

    r=[[0;0;0],reshape(p(1:3*(n-1)),3,[])];
    t=[[0;0;0],reshape(p(3*(n-1)+1:6*(n-1)),3,[])];
    p=p(6*(n-1)+1:end);
    X=reshape(p,3,[]);
end

    