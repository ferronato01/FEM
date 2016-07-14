
#This algorithm is written to calculate a Structure from a previous inputed Mesh using 3-Nodes-Element FEM Theory and return the Displacement Vector,
#the Element Stress and Principal Stresses.
#In order to calculate this, is needed first to insert the structure, that shall be done first, then is needed
#to set the Degree of Freedom of each node, and then both Local and Global Degree of Freedom that is done automatically.
#When everything get done, the proceedure is Function 'main' calling all the functions, and if all goes right the output will be the result.
#First function to be called is one to set the Global Stiffness Matrix (GSM) by finding but not armazenating the Local Stiffness Matrix (LSM),
#when this matrix is set, we need to apply Boundary Conditions, in order to the resultant linear system get solved, the solution-vector got
#from Ax=B is the nodal displacement, and the main result, as with it everything else can be found.

#Global Variables
#Here we load the Mesh with all the data within it
n = [11;12;3;1]
nodes = [1 0 0.5;2 0.25 0.5;3 0.125 0.375;4 0 0.25;5 0.25 0.25;6 0.125 0.125;7 0 0;8 0.25 0;9 0.375 0.125;10 0.5 0.25;11 0.5 0]
elements =[1 1 3 2 210 0.3 0.025 1;2 1 4 3 210 0.3 0.025 1;3 3 5 2 210 0.3 0.025 1;4 3 4 5 210 0.3 0.025 1;5 4 6 5 210 0.3 0.025 1
            6 4 7 6 210 0.3 0.025 1;7 5 6 8 210 0.3 0.025 1;8 6 7 8 210 0.3 0.025 1;9 5 8 9 210 0.3 0.025 1;10 5 9 10 210 0.3 0.025 1
            11 8 11 9 210 0.3 0.025 1;12 9 11 10 210 0.3 0.025 1]
supports = [1 1 0 0;2 4 0 0;3 7 0 0]
loadings = [1 5 0 -12.5]
#First the Global Variables has to be setted
DoF=2              #Degree of Freedom per Node
NEle=3             #Number of Elements per Node
GDoF=DoF*n[1]      #Global Degree of Freedom
LDoF=NEle*DoF      #Local Degree of Freedom
#


function GlobalStiffness(elements, nodes,n,NEle,DoF,K)
  LDoF=NEle*DoF    #Local Degree of Freedom
  GDoF=n[1]*DoF    #Global Degree of Freedom
#Here the vectors are initialized, as will be repeated occasionally
  xx = zeros(Float64,n[1])
  yy = zeros(Float64,n[1])
#Nodes coordinates for Lenght are gotten
  for i = 1:n[1]
    xx[i]=nodes[i,2]
    yy[i]=nodes[i,3]
  end
  indice=zeros(Int64,NEle)
#Here the iterations to set the GSM starts
  for e = 1:n[2]
    for i = 1:NEle
      indice[i]=elements[e,i+1]
    end
    elementDof=zeros(Int64,LDoF)
#First each element connective has have to be found
    for i = 1:NEle
      for j = 1:DoF
        elementDof[2(i-1)+j]=indice[i]*2+(j-2)
      end
    end
#Then the Triangular-Element matrices coefficients got set, and after the constants
    A = (xx[indice[1]]*(yy[indice[2]]-yy[indice[3]]) + xx[indice[2]]*(yy[indice[3]]-yy[indice[1]]) + xx[indice[3]]*(yy[indice[1]]-yy[indice[2]]))/2
    ν = elements[e,6]
    E = elements[e,5]
    t = elements[e,7]
    p = elements[e,8]
    βi=yy[indice[2]]-yy[indice[3]]
    βj=yy[indice[3]]-yy[indice[1]]
    βm=yy[indice[1]]-yy[indice[2]]
    γi=xx[indice[3]]-xx[indice[2]]
    γj=xx[indice[1]]-xx[indice[3]]
    γm=xx[indice[2]]-xx[indice[1]]
    Β = [βi 0 βj 0 βm 0
         0 γi 0 γj 0 γm
         γi βi γj βj γm βm]/(2*A)
#Stress or Strain must be verified
    if p==1
      D = (E/(1-ν*ν))*[1 ν 0;ν 1 0;0 0 (1-ν)/2]
    elseif p==2
      D = (E/((1+ν)*(1-2*ν)))*[1-ν ν 0;ν 1-ν 0;0 0 (1-2*ν)/2]
    end
#Here we have everything we need to LSM and then GSM get set
    k = t*A*(transpose(Β))*D*Β
#Now GSM get setted
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
#1 indicate that the node can displace/rotate in that direction, and 0 that is restrained (can't)
#If in a coordinate at Displacement Vector (u) has 1 it imply that in the same coordinate at the Force Vector (F) should be a force that is not 0
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

function stress(elements, nodes,n,NEle,DoF,desloc,σ)
#This function is most like the GlobalStiffness one, but here we don't want just to find GSM (as we found before), here is wanted to get the
#Element's Stress
  LDoF=NEle*DoF
  GDoF=n[1]*DoF
  xx = zeros(Float64,n[1])
  yy = zeros(Float64,n[1])
  for e = 1:n[2]
    for i = 1:n[1]
      xx[i]=nodes[i,2]
      yy[i]=nodes[i,3]
    end
    indice=zeros(Int64,NEle)
    for i = 1:NEle
      indice[i]=elements[e,i+1]
    end
    A = (xx[indice[1]]*(yy[indice[2]]-yy[indice[3]]) + xx[indice[2]]*(yy[indice[3]]-yy[indice[1]]) + xx[indice[3]]*(yy[indice[1]]-yy[indice[2]]))/2
    ν = elements[e,6]
    E = elements[e,5]
    t = elements[e,7]
    p = elements[e,8]
    βi=yy[indice[2]]-yy[indice[3]]
    βj=yy[indice[3]]-yy[indice[1]]
    βm=yy[indice[1]]-yy[indice[2]]
    γi=xx[indice[3]]-xx[indice[2]]
    γj=xx[indice[1]]-xx[indice[3]]
    γm=xx[indice[2]]-xx[indice[1]]
    u=zeros(Float64,NEle*DoF)
#At this point is interesting to have the nodals displacement of each Element, and so is needed to call it from main and set it here
    for i = 1:NEle
      u[DoF*i-1]=desloc[indice[i]*2-1]
      u[DoF*i]=desloc[indice[i]*2]
    end
#From this on it should return to be the same as before, but the last equation, that in Finite Element Theory is a bit different
    Β = [βi 0 βj 0 βm 0
         0 γi 0 γj 0 γm
         γi βi γj βj γm βm]/(2*A)
    if p==1
      D = (E/(1-ν*ν))*[1 ν 0
                       ν 1 0
                       0 0 (1-ν)/2]
    elseif p==2
      D = (E/((1+ν)*(1-2*ν)))*[1-ν ν 0
                               ν 1-ν 0
                               0 0 (1-2*ν)/2]
    end
    σ=D*Β*u
  end
end

#Then the iterations is started by calling the functions above written
K = zeros(Float64,GDoF,GDoF)
GlobalStiffness(elements,nodes,n,NEle,DoF,K)
u=ones(Float64,n[1]*DoF)
F=zeros(Float64,n[1]*DoF)
fixednodes(supports,loadings,n,DoF,u,F)
boundaryconditions(K,u,n[1],DoF)
desloc = K\F
GlobalStiffness(elements, nodes,n,NEle,DoF,K)
cargas = K*desloc
σ=zeros(Float64,DoF*NEle)
stress(elements, nodes,n,NEle,DoF,desloc,σ)
println(desloc)

