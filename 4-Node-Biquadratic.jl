#
n =[8;1;3;3]
nodes=[1 0 0;2 0.25 0;3 0.5 0;4 0 0.125;5 0.5 0.125;6 0 0.25;7 0.25 0.25;8 0.5 0.25]
elements=[1 1 3 8 6 2 5 7 4 210e9 0.3 0.025 1]
supports=[1 1 0 0;2 4 0 0;3 6 0 0]
loadings=[1 3 3.125e3 0;2 5 12.5e3 0;3 8 3.125e3 0]
nintpt=9
DoF=2
NEle=8
LDoF=NEle*n[1]
GDoF=DoF*n[1]
#

function GlobalStiffness(elements,nodes,n,NEle,nintpt,DoF,K)
  LDoF=NEle*DoF
  GDoF=n[1]*DoF
  xx = zeros(Float64,n[1])
  yy = zeros(Float64,n[1])
  for i = 1:n[1]
    xx[i]=nodes[i,2]
    yy[i]=nodes[i,3]
  end
  indice=zeros(Int64,NEle)
  coord=zeros(Float64,NEle,2)
  for e = 1:n[2]
    for i = 1:NEle
      indice[i]=elements[e,i+1]
    end
    for i = 1:NEle
      for j = 1:2
        coord[i,j]=nodes[indice[i],j+1]
      end
    end
    elementDof=zeros(Int64,LDoF)
    for i = 1:NEle
      for j = 1:DoF
        elementDof[2(i-1)+j]=indice[i]*2+(j-2)
      end
    end
    k=zeros(Float64,LDoF,LDoF)
    xi=zeros(Float64,2,nintpt)
    w′=zeros(Float64,3)
    w′=[0.555555555 0.888888888 0.55555555555]
    w=zeros(Float64,9)
    for i=1:3, j=1:3
      n=3*(i-1)+j
      w[n]=(w′[i])*(w′[j])
    end
    xi[1,1] = -0.7745966692
    xi[2,1] = xi[1,1]
    xi[1,2] = 0
    xi[2,2] = xi[1,1]
    xi[1,3] = -xi[1,1]
    xi[2,3] = xi[1,1]
    xi[1,4] = xi[1,1]
    xi[2,4] = 0
    xi[1,5] = 0
    xi[2,5] = 0
    xi[1,6] = -xi[1,1]
    xi[2,6] = 0
    xi[1,7] = xi[1,1]
    xi[2,7] = -xi[1,1]
    xi[1,8] = 0
    xi[2,8] = -xi[1,1]
    xi[1,9] = -xi[1,1]
    xi[2,9] = -xi[1,1]
    E = elements[e,10]
    ν = elements[e,11]
    t = elements[e,12]
    for intpt = 1:nintpt
      N=zeros(Float64,NEle)
      dNdxi=zeros(Float64,NEle,2)
      B=zeros(Float64,3,NEle*2)
      ξ=xi[1,intpt]
      η=xi[2,intpt]
      N[1]=((1-ξ)*(1-η))/4
      N[2]=((1+ξ)*(1-η))/4
      N[3]=((1+ξ)*(1+η))/4
      N[4]=((1-ξ)*(1+η))/4
      N[5]=(1-ξ*ξ)*(1-η)/2
      N[6]=(1+ξ)*(1-η*η)/2
      N[7]=(1-ξ*ξ)*(1+η)/2
      N[8]=(1-ξ)*(1-η*η)/2
      dNdxi[1,1] = 0.25*(1.-η)*(2.*ξ+η);
      dNdxi[1,2] = 0.25*(1.-ξ)*(ξ+2.*η);
      dNdxi[2,1] = 0.25*(1.-η)*(2.*ξ-η);
      dNdxi[2,2] = 0.25*(1.+ξ)*(2.*η-ξ);
      dNdxi[3,1] = 0.25*(1.+η)*(2.*ξ+η);
      dNdxi[3,2] = 0.25*(1.+ξ)*(2.*η+ξ);
      dNdxi[4,1] = 0.25*(1.+η)*(2.*ξ-η);
      dNdxi[4,2] = 0.25*(1.-ξ)*(2.*η-ξ);
      dNdxi[5,1] = -ξ*(1.-η);
      dNdxi[5,2] = -0.5*(1.-ξ*ξ);
      dNdxi[6,1] = 0.5*(1.-η*η);
      dNdxi[6,2] = -(1.+ξ)*η;
      dNdxi[7,1] = -ξ*(1.+η);
      dNdxi[7,2] = 0.5*(1.-ξ*ξ);
      dNdxi[8,1] = -0.5*(1.-η*η);
      dNdxi[8,2] = -(1.-ξ)*η;
      dxdxi=zeros(Float64,2,2)
      for i = 1:2, j = 1:2
        dxdxi[i,j]=0
        for l = 1:NEle
          dxdxi[i,j]=dxdxi[i,j]+coord[l,i]*dNdxi[l,j]
        end
      end
      dxidx=zeros(Float64,2,2)
      dtm=0
      dtm=dxdxi[1,1]*dxdxi[2,2]-dxdxi[1,2]*dxdxi[2,1]
      dxidx[1,1]=dxdxi[2,2]/dtm
      dxidx[2,2]=dxdxi[1,1]/dtm
      dxidx[1,2]=-dxdxi[1,2]/dtm
      dxidx[2,1]=-dxdxi[2,1]/dtm
      dNdx=zeros(Float64,NEle,2)
      for l = 1:NEle, i = 1:2
        dNdx[l,i]=0
        for j = 1:2
          dNdx[l,i]=dNdx[l,i]+dNdxi[l,j]*dxidx[j,i]
        end
      end
      if elements[e,13]==1
        D = (E/(1-ν*ν))*[1 ν 0;ν 1 0;0 0 (1-ν)/2]
      elseif elements[e,13]==2
        D = (E/((1+ν)*(1-2*ν)))*[1-ν ν 0;ν 1-ν 0;0 0 (1-2*ν)/2]
      end
      iy=0
      for i = 1:NEle
        ix=iy+1
        iy=ix+1
        B[1,ix]=dNdx[i,1]
        B[1,iy]=0
        B[2,ix]=0
        B[2,iy]=dNdx[i,2]
        B[3,ix]=dNdx[i,2]
        B[3,iy]=dNdx[i,1]
      end
      k = k + transpose(B)*D*B*w[intpt]*dtm*t
    end
    for i = 1:LDoF
      for j = 1:LDoF
        row =elementDof[i]
        column = elementDof[j]
        K[row,column]=K[row,column]+k[i,j]
      end
    end
  end
end

function fixednodes(supports,loadings,n,DoF,u,F)
  for i = 1:n[3]
    row=Int(supports[i,2]*2-1)
    for j = 1:DoF
      u[row-1+j]=supports[i,j+2]
    end
  end
  for i = 1:n[4]
    row=Int(loadings[i,2]*2-1)
    for j = 1:DoF
      F[row-1+j]=loadings[i,j+2]
    end
  end
  for i = 1:n[1]*DoF
    if(u[i]==0)
      F[i]=0
    end
  end
end

function boundaryconditions(K,u,n_no,DoF)
  for i = 1:n_no*DoF
    if(u[i]==0)
      for j = 1:n_no*DoF
        K[i,j]=0
        K[j,i]=0
        if(i==j)
          K[i,j]=1
        end
      end
    end
  end
end

function Stresses(elements,nodes,n,NEle,nintpt,DoF,K,desloc)
  LDoF=NEle*DoF
  GDoF=n[1]*DoF
  xx = zeros(Float64,n[1])
  yy = zeros(Float64,n[1])
  for i = 1:n[1]
    xx[i]=nodes[i,2]
    yy[i]=nodes[i,3]
  end
  indice=zeros(Int64,NEle)
  coord=zeros(Float64,NEle,2)
  for e = 1:n[2]
    for i = 1:NEle
      indice[i]=elements[e,i+1]
    end
    for i = 1:NEle
      for j = 1:2
        coord[i,j]=nodes[indice[i],j+1]
      end
    end
    elementDof=zeros(Int64,LDoF)
    for i = 1:NEle
      for j = 1:DoF
        elementDof[2(i-1)+j]=indice[i]*2+(j-2)
      end
    end
    k=zeros(Float64,LDoF,LDoF)
    xi=zeros(Float64,2,nintpt)
    w′=zeros(Float64,3)
    w′=[0.555555555 0.888888888 0.55555555555]
    w=zeros(Float64,9)
    for i=1:3, j=1:3
      n=3*(i-1)+j
      w[n]=(w′[i])*(w′[j])
    end
    xi[1,1] = -0.7745966692
    xi[2,1] = xi[1,1]
    xi[1,2] = 0
    xi[2,2] = xi[1,1]
    xi[1,3] = -xi[1,1]
    xi[2,3] = xi[1,1]
    xi[1,4] = xi[1,1]
    xi[2,4] = 0
    xi[1,5] = 0
    xi[2,5] = 0
    xi[1,6] = -xi[1,1]
    xi[2,6] = 0
    xi[1,7] = xi[1,1]
    xi[2,7] = -xi[1,1]
    xi[1,8] = 0
    xi[2,8] = -xi[1,1]
    xi[1,9] = -xi[1,1]
    xi[2,9] = -xi[1,1]
    E = elements[e,10]
    ν = elements[e,11]
    t = elements[e,12]
    for intpt = 1:nintpt
      N=zeros(Float64,NEle)
      dNdxi=zeros(Float64,NEle,2)
      B=zeros(Float64,3,NEle*2)
      ξ=xi[1,intpt]
      η=xi[2,intpt]
      N[1]=((1-ξ)*(1-η))/4
      N[2]=((1+ξ)*(1-η))/4
      N[3]=((1+ξ)*(1+η))/4
      N[4]=((1-ξ)*(1+η))/4
      N[5]=(1-ξ*ξ)*(1-η)/2
      N[6]=(1+ξ)*(1-η*η)/2
      N[7]=(1-ξ*ξ)*(1+η)/2
      N[8]=(1-ξ)*(1-η*η)/2
      dNdxi[1,1] = 0.25*(1.-η)*(2.*ξ+η);
      dNdxi[1,2] = 0.25*(1.-ξ)*(ξ+2.*η);
      dNdxi[2,1] = 0.25*(1.-η)*(2.*ξ-η);
      dNdxi[2,2] = 0.25*(1.+ξ)*(2.*η-ξ);
      dNdxi[3,1] = 0.25*(1.+η)*(2.*ξ+η);
      dNdxi[3,2] = 0.25*(1.+ξ)*(2.*η+ξ);
      dNdxi[4,1] = 0.25*(1.+η)*(2.*ξ-η);
      dNdxi[4,2] = 0.25*(1.-ξ)*(2.*η-ξ);
      dNdxi[5,1] = -ξ*(1.-η);
      dNdxi[5,2] = -0.5*(1.-ξ*ξ);
      dNdxi[6,1] = 0.5*(1.-η*η);
      dNdxi[6,2] = -(1.+ξ)*η;
      dNdxi[7,1] = -ξ*(1.+η);
      dNdxi[7,2] = 0.5*(1.-ξ*ξ);
      dNdxi[8,1] = -0.5*(1.-η*η);
      dNdxi[8,2] = -(1.-ξ)*η;
      dxdxi=zeros(Float64,2,2)
      for i = 1:2, j = 1:2
        dxdxi[i,j]=0
        for l = 1:NEle
          dxdxi[i,j]=dxdxi[i,j]+coord[l,i]*dNdxi[l,j]
        end
      end
      dxidx=zeros(Float64,2,2)
      dtm=0
      dtm=dxdxi[1,1]*dxdxi[2,2]-dxdxi[1,2]*dxdxi[2,1]
      dxidx[1,1]=dxdxi[2,2]/dtm
      dxidx[2,2]=dxdxi[1,1]/dtm
      dxidx[1,2]=-dxdxi[1,2]/dtm
      dxidx[2,1]=-dxdxi[2,1]/dtm
      dNdx=zeros(Float64,NEle,2)
      for l = 1:NEle, i = 1:2
        dNdx[l,i]=0
        for j = 1:2
          dNdx[l,i]=dNdx[l,i]+dNdxi[l,j]*dxidx[j,i]
        end
      end
      if elements[e,13]==1
        D = (E/(1-ν*ν))*[1 ν 0;ν 1 0;0 0 (1-ν)/2]
      elseif elements[e,13]==2
        D = (E/((1+ν)*(1-2*ν)))*[1-ν ν 0;ν 1-ν 0;0 0 (1-2*ν)/2]
      end
      iy=0
      for i = 1:NEle
        ix=iy+1
        iy=ix+1
        B[1,ix]=dNdx[i,1]
        B[1,iy]=0
        B[2,ix]=0
        B[2,iy]=dNdx[i,2]
        B[3,ix]=dNdx[i,2]
        B[3,iy]=dNdx[i,1]
      end
    end
  end
end


K=zeros(Float64,GDoF,GDoF)
GlobalStiffness(elements,nodes,n,NEle,nintpt,DoF,K)
u=ones(Float64,n[1]*DoF)
F=zeros(Float64,n[1]*DoF)
fixednodes(supports,loadings,n,DoF,u,F)
boundaryconditions(K,u,n[1],DoF)
desloc = (K\F)
σ = Stresses(elements,nodes,n,NEle,nintpt,DoF,K,desloc)
