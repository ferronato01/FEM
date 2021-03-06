#
n =[6;2;2;2]
nodes=[1 0 0;2 0.25 0;3 0.5 0;4 0 0.25;5 0.25 0.25;6 0.5 0.25]
elements=[1 1 2 5 4 210e9 0.3 0.025 1;2 2 3 6 5 210e9 0.3 0.025 1]
supports=[1 1 0 0;2 4 0 0]
loadings=[1 3 9.375e3 0;2 6 9.375e3 0]
nintpt=4
DoF=2
NEle=4
LDoF=NEle*n[1]
GDoF=DoF*n[1]
nldegree=5
increment=zeros(Float64,n[4],2+nldegree*2)
for i=1:n[4]
  increment[i,1]=loadings[i,1]
  increment[i,2]=loadings[i,2]
end
for i=1:n[4], j=1:nldegree
  k=1
  l=3
  Δqx=loadings[i,l]/nldegree
  increment[i,k*2j+1]=Δqx*j
  l+=1
  Δqy=loadings[i,l]/nldegree
  increment[i,k*2j+2]=Δqy*j
end
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
    xi=zeros(Float64,2,4)
    w=ones(Float64,4)
    xi[1,1]=-(sqrt(3))/3
    xi[2,1] = xi[1,1]
    xi[1,2] = -xi[1,1]
    xi[2,2] = xi[1,1]
    xi[1,3] = xi[1,1]
    xi[2,3] = -xi[1,1]
    xi[1,4] = -xi[1,1]
    xi[2,4] = -xi[1,1]
    E = elements[e,6]
    ν = elements[e,7]
    t = elements[e,8]
    for intpt = 1:nintpt
      N=zeros(Float64,nintpt)
      dNdxi=zeros(Float64,nintpt,2)
      B=zeros(Float64,3,8)
      ξ=xi[1,intpt]
      η=xi[2,intpt]
      N[1]=((1-ξ)*(1-η))/4
      N[2]=((1+ξ)*(1-η))/4
      N[3]=((1+ξ)*(1+η))/4
      N[4]=((1-ξ)*(1+η))/4
      dNdxi[1,1]=-(1-η)/4
      dNdxi[1,2]=-(1-ξ)/4
      dNdxi[2,1]=(1-η)/4
      dNdxi[2,2]=-(1+ξ)/4
      dNdxi[3,1]=(1+η)/4
      dNdxi[3,2]=(1+ξ)/4
      dNdxi[4,1]=-(1+η)/4
      dNdxi[4,2]=(1-ξ)/4
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
      if elements[e,9]==1
        D = (E/(1-ν*ν))*[1 ν 0;ν 1 0;0 0 (1-ν)/2]
      elseif elements[e,9]==2
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

function fixednodes(supports,n,DoF,u,F,nldegree,increment)
  δ=0
  for inc=1:nldegree
    for i = 1:n[3]
      row=Int(supports[i,2]*DoF-1)
      for j = 1:DoF
        u[row-1+j]=supports[i,j+2]
      end
    end
    for i = 1:n[4]
      row=Int(increment[i,2]*DoF-1)
      for j = 1:DoF
        F[row-1+j+δ]=increment[i,2*inc+j]
      end
    end
    δ+=n[1]*DoF
    for i = 1:n[1]*DoF
      if(u[i]==0)
        F[i]=0
      end
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

function Stress(elements,nodes,n,NEle,nintpt,DoF,desloc,σ)
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
    xi=zeros(Float64,2,4)
    w=ones(Float64,4)
    xi[1,1]=-(sqrt(3))/3
    xi[2,1] = xi[1,1]
    xi[1,2] = -xi[1,1]
    xi[2,2] = xi[1,1]
    xi[1,3] = xi[1,1]
    xi[2,3] = -xi[1,1]
    xi[1,4] = -xi[1,1]
    xi[2,4] = -xi[1,1]
    E = elements[e,6]
    ν = elements[e,7]
    t = elements[e,8]
    for intpt = 1:nintpt
      N=zeros(Float64,nintpt)
      dNdxi=zeros(Float64,nintpt,2)
      B=zeros(Float64,3,8)
      ξ=xi[1,intpt]
      η=xi[2,intpt]
      N[1]=((1-ξ)*(1-η))/4
      N[2]=((1+ξ)*(1-η))/4
      N[3]=((1+ξ)*(1+η))/4
      N[4]=((1-ξ)*(1+η))/4
      dNdxi[1,1]=-(1-η)/4
      dNdxi[1,2]=-(1-ξ)/4
      dNdxi[2,1]=(1-η)/4
      dNdxi[2,2]=-(1+ξ)/4
      dNdxi[3,1]=(1+η)/4
      dNdxi[3,2]=(1+ξ)/4
      dNdxi[4,1]=-(1+η)/4
      dNdxi[4,2]=(1-ξ)/4
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
      if elements[e,9]==1
        D = (E/(1-ν*ν))*[1 ν 0;ν 1 0;0 0 (1-ν)/2]
      elseif elements[e,9]==2
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
u=ones(Float64,n[1]*DoF*nldegree)
F=zeros(Float64,n[1]*DoF*nldegree)
fixednodes(supports,n,DoF,u,F,nldegree,increment)
boundaryconditions(K,u,n[1],DoF)
F′=zeros(Float64,n[1]*DoF)
σ=zeros(Float64,n[2]*nldegree,NEle)
δ=0
for i=1:nldegree
  for j=1:n[1]*DoF
    F′[j]=F[δ+j]
  end
  δ+=n[1]*DoF
  desloc = K\F′
#Reformular
  #Fazer o vetor desloc ser uma matriz com nldegree linhas e armazenar o deslocamento de cada K\F′ em uma linha da matriz desloc
  #Como desloc passa a ser uma matriz, na função Stress, tem que fazer um loop maior para calcular as tensões, e cada iteração percorrerá uma linha de desloc
  println(desloc)
  Stress(elements,nodes,n,NEle,nintpt,DoF,desloc,σ)
end
