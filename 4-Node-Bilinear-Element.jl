# The primary need to run this code is to have the informations about the initial mesh, and it has to be manually inserted below at the following vectors
n =[6;2;2;2]     #This is a Vector that Quantify the informations that should be inserted in the Structure and is designed to follow this sequence [Number of Nodes ; Number of Elements ; Number of Restrained Nodes ; Number of Loadings]
nodes=[1 0 0;2 0.25 0;3 0.5 0;4 0 0.25;5 0.25 0.25;6 0.5 0.25]
elements=[1 1 2 5 4 210e9 0.3 0.025 1;2 2 3 6 5 210e9 0.3 0.025 1]
supports=[1 1 0 0;2 4 0 0]
loadings=[1 3 9.375e3 0;2 6 9.375e3 0]
nintpt=4         #Number of Integration Points / As is a Linear Quadrilateral Element it's supposed to be 4 integration points per element
DoF=2            
NEle=4           
LDoF=NEle*n[1]   
GDoF=DoF*n[1]    
#
#This is a function that it's goal is to assemble the Global Stiffness Matrix (GSM)
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
    xi[1,1]=-(sqrt(3))/3            #Natural Coordinates given by Gaussian Quadrature
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
      xii=zeros(Float64,2)
      B=zeros(Float64,3,8)  
      xii[1]=xi[1,intpt]    
      xii[2]=xi[2,intpt]
      N[1]=((1-xii[1])*(1-xii[2]))/4
      N[2]=((1+xii[1])*(1-xii[2]))/4
      N[3]=((1+xii[1])*(1+xii[2]))/4
      N[4]=((1-xii[1])*(1+xii[2]))/4
      dNdxi[1,1]=-(1-xii[2])/4
      dNdxi[1,2]=-(1-xii[1])/4
      dNdxi[2,1]=(1-xii[2])/4
      dNdxi[2,2]=-(1+xii[1])/4
      dNdxi[3,1]=(1+xii[2])/4
      dNdxi[3,2]=(1+xii[1])/4
      dNdxi[4,1]=-(1+xii[2])/4
      dNdxi[4,2]=(1-xii[1])/4
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

function fixednodes(supports,loadings,n,DoF,u,F)
#This function is designed to be able to get the nodes in the Structure and write a vector that has 1 or 0 in it
#1 indicate that the node can displace/rotate in that direction, and 0 that is fixed (can't)
#If in a coordinate at Displacement Vector (u) has 1 it implie that in the same coordinate at the Force Vector (F) should be a force that is not 0
#and vice-versa
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
#This function was written to get the Displacement Vector and apply Boundary Conditions at the GSM by analyzing if that direction is free to displace
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

K=zeros(Float64,GDoF,GDoF)
GlobalStiffness(elements,nodes,n,NEle,nintpt,DoF,K)
u=ones(Float64,n[1]*DoF)
F=zeros(Float64,n[1]*DoF)
fixednodes(supports,loadings,n,DoF,u,F)
boundaryconditions(K,u,n[1],DoF)
desloc = K\F
println(desloc)