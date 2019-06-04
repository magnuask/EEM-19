#Packages:
using CSV
using DataFrames
using JuMP
using GLPKMathProgInterface
using Gurobi

#Set working directory
cd("inputdir")

#Set true if you want to fix the investments:
fixbuilding=false

#Model:
SO = Model(with_optimizer(Gurobi.Optimizer))

#sets
setsData=CSV.read("setsData.csv")
N=1:setsData.Buildings[1,1]
H=1:setsData.Hours[1,1]
N_max=setsData.Buildings[1,1]
H_max=setsData.Hours[1,1]
#variables
buildingExp_EQ=CSV.read("BuildingExp-EQ_vnt+cnt.csv")
@variable(SO,cgrid[n in N] >= 0)
@variable(SO,cB[n in N] >= 0)
@variable(SO,cPV[n in N] >= 0)
@variable(SO,ctot >= 0)
@variable(SO,soc[n in N,h in H] >= 0)
@variable(SO,dDeltaUp[n in N,h in H] >= 0)
@variable(SO,dDeltaDown[n in N,h in H] >= 0)
@variable(SO,eexp[n in N,h in H] >= 0)
@variable(SO,eimp[n in N,h in H] >= 0)
@variable(SO,etran[h in H] >= 0)
@variable(SO,pPV[n in N,h in H] >= 0)
if fixbuilding==true
    for n in N
        JuMP.fix(cB[n],buildingExp_EQ.cB[n],force=true)
        JuMP.fix(cPV[n],buildingExp_EQ.cPV[n],force=true)
    end
end
#Parameters
#Building data
building=CSV.read("buildingData.csv")
SD=Array{Float64}(zeros(N_max))
SL=Array{Float64}(zeros(N_max))
PB=Array{Float64}(zeros(N_max))
IB=Array{Float64}(zeros(N_max))
IPV=Array{Float64}(zeros(N_max))
for n in N
    SD[n]=building.SD[n]
    SL[n]=building.SL[n]
    PB[n]=building.PB[n]
    IB[n]=building.IB[n]
    IPV[n]=building.IPV[n]
end

#building-hour data
buildingTime=CSV.read("buildingTimeData.CSV")
PVData=CSV.read("PVData.CSV")
D=zeros(N_max,H_max)
GPV=zeros(N_max,H_max)
for n in N
D[n,:]=buildingTime.Demand[buildingTime.Building.==n,:]
GPV[n,:]=PVData.EPV[:]/1000
end


#system data
systemData=CSV.read("systemData.csv")
lambda = Array{Float64}(zeros(H_max))
for h in H
lambda[h]=systemData.lambdabuy[h]
end

#DSO data
DSOData=CSV.read("DSOData.csv")
Icap=DSOData.Icap[1,1]
LDSO=DSOData.LDSO[1,1]
TAX=DSOData.TAX[1,1]
#Objective:
@expression(SO,buildingCapCost[n in N],IPV[n]*cPV[n] + IB[n]*cB[n]);

@expression(SO,buildingEnergyCost[n in N],sum((eimp[n,h]*(1+TAX)lambda[h] - eexp[n,h]*lambda[h]) for h in H));

@expression(SO,CapCostDSO,Icap*ctot);

@expression(SO,VarCostDSO,sum(etran[h]*LDSO*lambda[h] for h in H));

@expression(SO,CostDSO,CapCostDSO+VarCostDSO);

@objective(SO, Min, sum(buildingCapCost[n] + buildingEnergyCost[n] for n in N) + CostDSO);

#Constraints:
#energy balance
@constraint(SO,enBal[n in N,h in H],eimp[n,h]-eexp[n,h] == D[n,h] + dDeltaUp[n,h] - dDeltaDown[n,h] - pPV[n,h])

#flexibility balance
@constraint(SO,FlexBal[n in N,h in H],soc[n,h] == soc[n,h>1 ? h-1 : H_max]*(1-SD[n]) + dDeltaUp[n,h]*(1-SL[n]) - dDeltaDown[n,h]*(1+SL[n]))
#@constraint(SO,FlexBoundary[n=1:N],soc[n,0] == soc[n,H])

#Heat load capacity:
@constraint(SO,FlexAmount[n in N,h in H],cB[n] >= soc[n,h])

#flexibility limit up:
@constraint(SO,FlexLimUp[n in N,h in H],cB[n]*PB[n] >= dDeltaUp[n,h])

#flexibility limit down:
@constraint(SO,FlexLimDown[n in N,h in H], cB[n]*PB[n]>= dDeltaDown[n,h])

#grid capacity allocation
@constraint(SO,GridCapAll[n in N,h in H], cgrid[n] >= eimp[n,h] + eexp[n,h])

#PV generation
@constraint(SO,PVGen[n in N,h in H],cPV[n]*GPV[n,h] >= pPV[n,h])

#DSO energy transfer 1:
@constraint(SO,DSOetran1[h in H],etran[h] >= sum(eimp[n,h]-eexp[n,h] for n in N))

#DSO energy transfer 2:
@constraint(SO,DSOetran2[h in H],etran[h] >= sum(eexp[n,h]-eimp[n,h] for n in N))

#DSO total capacity:
@constraint(SO,DSOTotCap[h in H],ctot >= etran[h])

#Solve statement:
status = optimize!(SO)

#Export of building data
building.cgrid=value.(cgrid)
building.cB = value.(cB)
building.cPV = value.(cPV)
building.buildingCapCost = value.(buildingCapCost)
building.buildingEnergyCost = value.(buildingEnergyCost)
#Export of temporal building data
buildingTime.soc = zeros(N_max*H_max)
buildingTime.dDeltaUp = zeros(N_max*H_max)
buildingTime.dDeltaDown = zeros(N_max*H_max)
buildingTime.eexp = zeros(N_max*H_max)
buildingTime.eimp = zeros(N_max*H_max)
buildingTime.pPV = zeros(N_max*H_max)
for n in N
    buildingTime.soc[buildingTime.Building.==n,:] = value.(soc[n,:])
    buildingTime.dDeltaUp[buildingTime.Building.==n,:] = value.(dDeltaUp[n,:])
    buildingTime.dDeltaDown[buildingTime.Building.==n,:] = value.(dDeltaDown[n,:])
    buildingTime.eexp[buildingTime.Building.==n,:] = value.(eexp[n,:])
    buildingTime.eimp[buildingTime.Building.==n,:] = value.(eimp[n,:])
    buildingTime.pPV[buildingTime.Building.==n,:] = value.(pPV[n,:])
end

#export of system data:
DSOData.ctot = value.(ctot)
DSOData.CapCostDSO = value.(CapCostDSO)
DSOData.VarCostDSO = value.(VarCostDSO)
DSOData.CostDSO = value.(CostDSO)
DSOData.BuildingCosts   = objective_value(SO) - value.(CostDSO)
DSOData.totcost = objective_value(SO)
DSOData.vnt=0
DSOData.cnt=0
#export of temporal system data:
systemData.etran=value.(etran)

#Writing data to files:
cd("outdir")
if fixbuilding==false
    CSV.write("BuildingExp-SO.csv",building,writeheader=true)
    CSV.write("BuildingTemporalExp-SO.csv",buildingTime,writeheader=true)
    CSV.write("DSOExp-SO.csv",DSOData,writeheader=true)
    CSV.write("SystemTemporalExp-SO.csv",systemData,writeheader=true)

else
    CSV.write("BuildingExp-SO-fx.csv",building,writeheader=true)
    CSV.write("BuildingTemporalExp-SO-fx.csv",buildingTime,writeheader=true)
    CSV.write("DSOExp-SO-fx.csv",DSOData,writeheader=true)
    CSV.write("SystemTemporalExp-SO-fx.csv",systemData,writeheader=true)
end
