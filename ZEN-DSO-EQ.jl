#Packages:
using CSV
using DataFrames
using JuMP
using GLPKMathProgInterface
using Gurobi
#Set working directory
cd("C:\\Users\\magnusa\\OneDrive - SINTEF\\Publications\\EEM-19_ZEN-and-DSO\\Modeling\\DataIn")

#choose tariffs, at least one must be true
capacity=false
volumetric=true


#sets
setsData=CSV.read("setsData.csv")
N=1:setsData.Buildings[1,1]
H=1:setsData.Hours[1,1]
N_max=setsData.Buildings[1,1]
H_max=setsData.Hours[1,1]

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
building.cgrid=zeros(N_max)
building.cB=zeros(N_max)
building.cPV=zeros(N_max)
building.buildingCapCost=zeros(N_max)
building.buildingEnergyCost=zeros(N_max)
building.buildingTariffCost=zeros(N_max)

#building-hour data
buildingTime=CSV.read("buildingTimeData.CSV")
PVData=CSV.read("PVData.CSV")
D=zeros(N_max,H_max)
GPV=zeros(N_max,H_max)
for n in N
D[n,:]=buildingTime.Demand[buildingTime.Building.==n,:]
GPV[n,:]=PVData.EPV[:]/1000
end
buildingTime.soc = zeros(N_max*H_max)
buildingTime.dDeltaUp = zeros(N_max*H_max)
buildingTime.dDeltaDown = zeros(N_max*H_max)
buildingTime.eexp = zeros(N_max*H_max)
buildingTime.eimp = zeros(N_max*H_max)
buildingTime.pPV = zeros(N_max*H_max)

#system data
systemData=CSV.read("systemData.csv")
lambda = Array{Float64}(zeros(H_max))
for h in H
lambda[h]=systemData.lambdabuy[h]
end
systemData.etran=zeros(H_max)
#DSO data
DSOData=CSV.read("DSOData.csv")
Icap=DSOData.Icap[1,1]
LDSO=DSOData.LDSO[1,1]
NM=DSOData.NM[1,1]
TAX=DSOData.TAX[1,1]

counter=0
ConvergenceCriteria=0.000001
convergence=false
#initialize tariffs on first iteration

while convergence==false
    global counter+=1
    println("iteration: ",counter)
    if counter == 1
        vnt=1
        cnt=100
    end
    #And run the building problems:
    println("vnt: ",vnt," ","cnt: ",cnt)
    for n in N
        println("Consumer: ",n)
        EQ=nothing
        EQ = Model(with_optimizer(Gurobi.Optimizer))
        @variable(EQ,cgrid >= 0)
        @variable(EQ,cB >= 0)
        @variable(EQ,cPV >= 0)
        @variable(EQ,soc[h in H] >= 0)
        @variable(EQ,dDeltaUp[h in H] >= 0)
        @variable(EQ,dDeltaDown[h in H] >= 0)
        @variable(EQ,eexp[h in H] >= 0)
        @variable(EQ,eimp[h in H] >= 0)
        @variable(EQ,pPV[h in H] >= 0)

        #Objective:
        @expression(EQ,buildingCapCost,IPV[n]*cPV + IB[n]*cB);

        @expression(EQ,buildingEnergyCost,sum((eimp[h]*(1+TAX)*lambda[h] - eexp[h]*lambda[h]) for h in H));

        @expression(EQ,buildingTariffCost,sum((eimp[h]-eexp[h]*NM)*vnt for h in H) + cgrid*cnt);

        @objective(EQ, Min, buildingCapCost + buildingEnergyCost + buildingTariffCost);

        #Constraints:
        #energy balance
        @constraint(EQ,enBal[h in H],eimp[h]-eexp[h] == D[n,h] + dDeltaUp[h] - dDeltaDown[h] - pPV[h])

        #flexibility balance
        @constraint(EQ,FlexBal[h in H],soc[h] == soc[h>1 ? h-1 : H_max]*(1-SD[n]) + dDeltaUp[h]*(1-SL[n]) - dDeltaDown[h]*(1+SD[n]))
        #@constraint(EQ,FlexBoundary,soc[0] == soc[H])

        #Heat load capacity:
        @constraint(EQ,FlexAmount[h in H],cB >= soc[h])

        #flexibility limit up:
        @constraint(EQ,FlexLimUp[h in H],cB*PB[n] >= dDeltaUp[h])

        #flexibility limit down:
        @constraint(EQ,FlexLimDown[h in H], cB*PB[n]>= dDeltaDown[h])

        #grid capacity allocation
        @constraint(EQ,GridCapAll[h in H], cgrid >= eimp[h] + eexp[h])

        #PV generation
        @constraint(EQ,PVGen[h in H],cPV*GPV[n,h] >= pPV[h])

        #Solve statement:
        status = optimize!(EQ)

        #Export of building data
        global building.cgrid[n]=value.(cgrid)
        global building.cB[n] = value.(cB)
        global building.cPV[n] = value.(cPV)
        global building.buildingCapCost[n] = value.(buildingCapCost)
        global building.buildingEnergyCost[n] = value.(buildingEnergyCost)
        global building.buildingTariffCost[n] = value.(buildingTariffCost)

        #Export of temporal building data
        global buildingTime.soc[buildingTime.Building.==n,:] = value.(soc[:])
        global buildingTime.dDeltaUp[buildingTime.Building.==n,:] = value.(dDeltaUp[:])
        global buildingTime.dDeltaDown[buildingTime.Building.==n,:] = value.(dDeltaDown[:])
        global buildingTime.eexp[buildingTime.Building.==n,:] = value.(eexp[:])
        global buildingTime.eimp[buildingTime.Building.==n,:] = value.(eimp[:])
        global buildingTime.pPV[buildingTime.Building.==n,:] = value.(pPV[:])
    end

    #Import data on building decisions:
    eimpcons=zeros(N_max,H_max)
    eexpcons=zeros(N_max,H_max)
    cgridcons=zeros(N_max)
    for n in N
        eimpcons[n,:] = buildingTime.eimp[buildingTime.Building.==n,:]
        eexpcons[n,:] = buildingTime.eexp[buildingTime.Building.==n,:]
        cgridcons[n,:] = building.cgrid[building.Building.==n,:]
    end

    #Calculating need for transmission:
    etran=zeros(H_max)
    etran=broadcast(abs,sum(eimpcons,dims=1) - sum(eexpcons,dims=1))
    ctot=maximum(etran)

    #Calculating DSO costs:
    CDSOInv=Icap*ctot
    CDSOVar=sum(etran*LDSO*lambda,dims=2)
    println("DSO investment costs: ",CDSOInv," ","DSO operational costs: ",CDSOVar)
    #Calculating DSO revenues:
    RDSOCap = sum(cgridcons,dims=1)*cnt
    RDSOVar = sum(sum(eimpcons-eexpcons*NM,dims=1),dims=2)*vnt
    println("DSO capacity revenue: ",RDSOCap," ","DSO volume-based revenue: ",RDSOVar)
    #Convergence check:
    if capacity == true && volumetric == true
        if abs(RDSOCap[1,1]-CDSOInv[1,1]) < ConvergenceCriteria && abs(RDSOVar[1,1]-CDSOVar[1,1] < ConvergenceCriteria)
            global convergence=true
        else
            global convergence=false
        end
    elseif capacity == false && volumetric == true
        if abs(RDSOVar[1,1] - CDSOVar[1,1] - CDSOInv[1,1]) < ConvergenceCriteria
            global convergence = true
        else
            global convergence = false
        end
    elseif capacity == true && volumetric == false
        if abs(RDSOCap[1,1] - CDSOVar[1,1] - CDSOInv[1,1]) < ConvergenceCriteria
            global convergence = true
        else
            global convergence = false
        end
    end
    if convergence==false
        #Calculate new tariffs:
        if capacity == true && volumetric == true
            cnt=CDSOInv/sum(cgridcons,dims=1)[1,1]
            vnt=CDSOVar/sum(sum(eimpcons-eexpcons*NM,dims=1),dims=2)
            global cnt=cnt[1,1]
            global vnt = vnt[1,1]
        elseif capacity == false && volumetric == true
            global cnt=0
            vnt=(CDSOVar[1,1]+CDSOInv)/sum(sum(eimpcons-eexpcons*NM,dims=1),dims=2)[1,1]
            global vnt = vnt[1,1]
        elseif capacity == true && volumetric == false
            global vnt = 0
            cnt = (CDSOVar[1,1]+CDSOInv)/sum(cgridcons,dims=1)[1,1]
            global cnt = cnt[1,1]
        end
    elseif convergence ==true
        print("convergence after ", counter," iterations")
    end
    #Save system results:
    global systemData.etran[:] = etran[:]
    global DSOData.ctot = ctot
    global DSOData.CapCostDSO = CDSOInv[1,1]
    global DSOData.VarCostDSO = CDSOVar[1,1]
    global DSOData.CostDSO = CDSOInv[1,1] + CDSOVar[1,1]
    global DSOData.BuildingCosts = zeros(1)
    for n in N
        global DSOData.BuildingCosts[1] += building.buildingCapCost[n] + building.buildingEnergyCost[n]
    end
    global DSOData.totcost = CDSOInv[1,1] + CDSOVar[1,1] + DSOData.BuildingCosts[1]
    global DSOData.vnt=vnt
    global DSOData.cnt=cnt
end


#Writing data to files:
cd("C:\\Users\\magnusa\\OneDrive - SINTEF\\Publications\\EEM-19_ZEN-and-DSO\\Modeling\\DataExp")
CSV.write("BuildingExp-EQ.csv",building,writeheader=true)
CSV.write("BuildingTemporalExp-EQ.csv",buildingTime,writeheader=true)
CSV.write("DSOExp-EQ.csv",DSOData,writeheader=true)
CSV.write("SystemTemporalExp-EQ.csv",systemData,writeheader=true)
