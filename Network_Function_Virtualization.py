import math
import random
from multiprocessing import Pool, Array, Value
import multiprocessing
import time
from array import array
import threading
import os
import sys
import cPickle
import scipy
import numpy as np
import cvxpy as cvx
from pickle import load
from gurobipy import *
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from numpy import linalg as LA
def ismember(k,d):
    a=0
    for i in range(len(k)):
      if d==k[i]:
          a=1
    return a
def findlink(LinkInfo,i,j):
    for t in range(len(LinkInfo)):
        if LinkInfo[t]==[i,j]:
            return t+1
def eye(k):
  x=np.zeros((k,k))
  for t in range(k):
    x[t,t]=1
  return x
def findsmall(x):
    adr=len(x)*[-1]
    count=0
    for t in range(len(x)):
        if x[t] > 10**(-3):
            adr[count]=t
            count=count+1
    return adr
def integer_count(A):
    t=A.shape
    c=0
    for j in range(t[0]):
        if A[j] > .09 and A[j] < .91:
            c=c+1
def findmember(k,d):
    a=-1
    for i in range(len(k)):
      if d==k[i]:
          a=1
          return i       
fin = open('CapacityVector.txt','r')
CapacityVector1=879*[0]
CapacityVector=[]
PsumItemax=20
rho=10
for line in fin.readlines():
    CapacityVector.append( [ int (x) for x in line.split(',') ] )
for i in range(879):
    CapacityVector1[i]=CapacityVector[i][0]
LinkInfo1=[[0 for i in range(2)] for j in range(879)]
fin = open('LinkInfo.txt','r')
LinkInfo=[]
for line in fin.readlines():
    LinkInfo.append( [ int (x) for x in line.split(',') ] )
FunctionNodes=np.array([42,51,72,77,93,96,154,163,184,189,205,208])
NumOfFuncs=4
NumOfCandidates=12 # number of nodes that can provide a function
ChainLength=1 # number of functions in a function chain
NumOfFuncNodes=12
NumOfNodes=111+112
NumOfLinks=880
ProcessCap=NumOfFuncNodes*[20]
FuncNodes=FunctionNodes
NumOfFlows=4
sigma2=2*np.ones(20)
sigma2[9]=sigma2[0]*2
sigma2[10:20]=sigma2[9]
TrafficRate=NumOfFlows*[1]
NumOfR=NumOfNodes*NumOfNodes*NumOfFlows*(ChainLength+1)
NumOfX=NumOfFuncNodes+NumOfFuncNodes*NumOfFuncs+NumOfFuncNodes*NumOfFuncs*NumOfFlows
FunctionChain=[1,2,3,4]
RankStart=[[1 for i in range(112)] for j in range(2)]
for i in range(2):
    for j in range(112):
        RankStart[i][j]=(i)*112+j+1
RankStart[1][111]=0
print RankStart
Pair=[[1,112],[2,112],[3,112],[4,112]]
ChainLength=1
NodeFuncTable=[[1 for i in range(NumOfFuncs)] for j in range(NumOfFuncNodes)]

NumOfCores=2
CoreNum=[112,112]
SizeOfF=CoreNum[1]*NumOfFlows*(ChainLength+1)
Function2Node=[[0 for i in range(NumOfFuncNodes)] for j in range(NumOfFuncs)]
for i in range(NumOfFuncs):
    for j in range(NumOfFuncNodes):
        Function2Node[i][j]=i
NumOfR=NumOfLinks*NumOfFlows*(ChainLength+1)
NumOfXEach=NumOfFuncNodes*NumOfFuncs*NumOfFlows+NumOfFuncNodes*NumOfFlows+NumOfFuncNodes
NumOfX=NumOfXEach
Xtemp=np.zeros(NumOfX)
Nr=NumOfLinks*NumOfFlows
NumOfLinksFrom=np.zeros(NumOfNodes)
NumOfLinksTo=np.zeros(NumOfNodes)
NumOfLinksFromCrossCore=np.zeros(NumOfNodes)
NumOfLinksToCrossCore=np.zeros(NumOfNodes)
CapaVector=CapacityVector
CapacityMatrix=np.zeros((NumOfNodes,NumOfNodes))
for i in range(NumOfNodes):
    for j in range(NumOfNodes):
        for k in range(NumOfLinks):
            if LinkInfo[k]==[i+1,j+1]:
                CapacityMatrix[i,j]=CapaVector[k][0]

FuncNodeNum=np.zeros(NumOfCores)
for i in range(NumOfFuncNodes):
    j=FuncNodes[i]
    for Core in range(NumOfCores):
        for k in range(112):
            if RankStart[Core][k]==j+1:
                FuncNodeNum[Core] = FuncNodeNum[Core]+1

NumOfLinksFrom=[0]*112*2
NumOfLinksTo=[0]*112*2
NumOfLinksFromCrossCore=[0]*112*2
NumOfLinksToCrossCore=[0]*112*2
LinkTableFrom=[[0 for i in range(7)] for j in range(112*2)]
LinkTableTo=[[0 for i in range(7)] for j in range(112*2)]
LinkTableFromCrossCore=[[0 for i in range(7)] for j in range(112*2)]
LinkTableToCrossCore=[[0 for i in range(7)] for j in range(112*2)]
SimulationIte=0
for i in range(NumOfLinks):
  for j in range(NumOfCores):
      for k in range(112):
        if LinkInfo[i][0]==RankStart[j][k]:
          if ismember(RankStart[j],LinkInfo[i][1]) > 0:
            NumOfLinksFrom[LinkInfo[i][0]-1] = NumOfLinksFrom[LinkInfo[i][0]-1]+1
            LinkTableFrom[LinkInfo[i][0]-1][NumOfLinksFrom[LinkInfo[i][0]-1]-1]=LinkInfo[i][1]
            NumOfLinksTo[LinkInfo[i][1]-1] = NumOfLinksTo[LinkInfo[i][1]-1]+1
            LinkTableTo[LinkInfo[i][1]-1][NumOfLinksTo[LinkInfo[i][1]-1]-1]=LinkInfo[i][0]
          if ismember(RankStart[j],LinkInfo[i][1]) == 0:
            NumOfLinksFromCrossCore[LinkInfo[i][0]-1] = NumOfLinksFromCrossCore[LinkInfo[i][0]-1]+1
            LinkTableFromCrossCore[LinkInfo[i][0]-1][NumOfLinksFromCrossCore[LinkInfo[i][0]-1]-1]=LinkInfo[i][1]
            NumOfLinksToCrossCore[LinkInfo[i][1]-1] = NumOfLinksToCrossCore[LinkInfo[i][1]-1]+1
            LinkTableToCrossCore[LinkInfo[i][1]-1][NumOfLinksToCrossCore[LinkInfo[i][1]-1]-1]=LinkInfo[i][0]
count1=0
count2=0
count3=0
count4=0
ObjCoeff = np.zeros((2*len(LinkInfo)*NumOfFlows*(ChainLength+1),2))
LocM=[[0 for i in range((sum(NumOfLinksFromCrossCore)+sum(NumOfLinksToCrossCore))*(ChainLength+1))] for j in range(NumOfCores)] 
XtempM=np.zeros((NumOfCores,NumOfX))
FlowTempM = np.zeros((NumOfCores,NumOfR))
LocXM=[[0 for i in range(NumOfFlows*1)] for j in range(NumOfCores)]
pen=0
LinksFrom=[[0 for i in range(10)] for j in range(NumOfNodes)]
for i in range(NumOfNodes):
    count=0
    for j in range(NumOfLinksFrom[i]):
        count=count+1
        LinksFrom[i][count-1]=findlink(LinkInfo,i+1,LinkTableFrom[i][j])
LinksFromCrossCore=[[0 for i in range(10)] for j in range(NumOfNodes)]
for i in range(NumOfNodes):
    count=0
    for j in range(NumOfLinksFromCrossCore[i]):
        count=count+1
        LinksFromCrossCore[i][count-1]=findlink(LinkInfo,i+1,LinkTableFromCrossCore[i][j])
LinksTo=[[0 for i in range(10)] for j in range(NumOfNodes)]
for i in range(NumOfNodes):
    count=0
    for j in range(NumOfLinksTo[i]):
        count=count+1
        LinksTo[i][count-1]=findlink(LinkInfo,LinkTableTo[i][j],i+1)
LinksToCrossCore=[[0 for i in range(10)] for j in range(NumOfNodes)]
for i in range(NumOfNodes):
    count=0
    for j in range(NumOfLinksToCrossCore[i]):
        count=count+1
        LinksToCrossCore[i][count-1]=findlink(LinkInfo,LinkTableToCrossCore[i][j],i+1)
RecvLinkRate=np.zeros((NumOfR))
X=np.zeros((NumOfX))
DualXM=np.zeros((NumOfCores,NumOfX))
DualX=np.zeros((NumOfX))
DualLinkRateM=np.zeros((NumOfCores,NumOfR))
DualLinkRate=np.zeros((NumOfR))
beta=len(FuncNodes)*[1]
beta0=len(FuncNodes)*[1]
Xikf=np.zeros((20,NumOfX))
pen=100
FXM=np.zeros((NumOfCores,SizeOfF,NumOfX))
FRM=np.zeros((NumOfCores,SizeOfF,NumOfR))
FbM=np.zeros((NumOfCores,SizeOfF))
Prev=np.zeros((NumOfR))
PrevX=np.zeros((NumOfX))
count1=0
count2=0
source=[[0 for i in range(NumOfFlows)] for j in range(NumOfCores)]
dest=[[0 for i in range(NumOfFlows)] for j in range(NumOfCores)]
for Rank in range(NumOfCores):
    for i in range(NumOfFlows):
        if ismember(RankStart[Rank],Pair[i][0]) > 0:
            source[0][count1]=Pair[i][0]
            source[1][count1]=Rank
            count1=count1+1
        if ismember(RankStart[Rank],Pair[i][1]) > 0:
            dest[0][count2]=Pair[i][1]
            dest[1][count2]=Rank
            count2=count2+1
NumOfLinksM=[0]*NumOfCores
I_m_0M=np.zeros((NumOfCores,NumOfR,NumOfR))
A_m_0M=np.zeros((NumOfCores,NumOfX,NumOfX))

def LocalUpdate(I_m_0, LinkTableFrom, NumOfLinksFrom, NumOfCores, source, dest, FunctionChain, FuncNodes, NumOfFuncs, NumOfFuncNodes, beta0, beta, ObjCoeff, Rank, RecvLinkRate, DualLinkRate, FR, FX, Q, Fb, CapaVector, X, DualX, NumOfLinksFromCrossCore, NumOfLinksToCrossCore, LinkTableFromCrossCore, LinkTableToCrossCore, RankStart, SizeOfF, NumOfNodes, ChainLength, NumOfFlows, NumOfLinks, NumOfR, NumOfX):
    rho=10
    Obj=np.zeros((NumOfR))
    for i in range(NumOfR):
        for j in range(len(ObjCoeff)):
            if i==ObjCoeff[j,0]:
                Obj[i]=ObjCoeff[j,1]
    FlowTemp = cvx.Variable(NumOfR)
    Xtemp = cvx.Variable(NumOfX)
    obj = cvx.Minimize(sum(Obj*FlowTemp)+sum(DualLinkRate*np.matmul(I_m_0,RecvLinkRate)-np.matmul(DualLinkRate,I_m_0)*FlowTemp)+sum(DualX*np.matmul(A_m_0,X)-np.matmul(DualX,A_m_0)*Xtemp)+rho/2*cvx.sum_squares(np.matmul(A_m_0,X)-A_m_0*Xtemp)+rho/2*cvx.sum_squares(np.matmul(I_m_0,RecvLinkRate)-I_m_0*FlowTemp))
    constraints = [Q*FlowTemp <=CapaVector, FR*FlowTemp-FX*Xtemp==Fb, FlowTemp>=0, 1 >= Xtemp, Xtemp>=0]
    prob = cvx.Problem(obj,constraints)
    prob.solve()
    xx=np.asarray(FlowTemp.value)
    yy=np.asarray(Xtemp.value)
    return (xx,yy)
def UpdateC(sigma2, Xikf, LocC, I_m_0M, LocM, LocXM, FunctionChain, FuncNodes, NumOfFuncs, NumOfFuncNodes, beta0, beta, ObjCoeff, DualLinkRate, FR, FX, Q, Fb, CapaVector, DualX, Xtemp, FlowTemp, NumOfLinksFromCrossCore, NumOfLinksToCrossCore, LinkTableFromCrossCore, LinkTableToCrossCore, RankStart, SizeOfF, NumOfNodes, ChainLength, NumOfFlows, NumOfLinks, NumOfR, NumOfX):
    ObjC=np.zeros((NumOfR))
    for i in range(NumOfR):
        if ismember(LocC,i) > 0:
            ObjC[i]=.5
    print ObjC
    H0=np.zeros((sum(NumOfLinksFromCrossCore),NumOfR))
    Count=0
    rho=10
    capa1=np.zeros((sum(NumOfLinksFromCrossCore)))
    for Rank in range(NumOfCores):
        for i in RankStart[Rank]:
            for j in range(NumOfLinksFromCrossCore[i-1]):
                NodeTo=LinkTableFromCrossCore[i-1][j];
                for k in range(NumOfFlows):
                    for s in range(ChainLength+1):
                        H0[Count,(LinksFromCrossCore[i-1][j]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=1
                capa1[Count]=CapacityMatrix[i-1][NodeTo-1]
                Count=Count+1
    E0=np.zeros((NumOfFuncNodes,NumOfX))
    count=-1
    for Rank in range(NumOfCores):
        for i in RankStart[Rank]:
            for j in range(NumOfFuncNodes):
                if FuncNodes[j]==i:
                    count=count+1
                    for k in range(NumOfFlows):
                        if j<= NumOfFuncNodes:
                            for f in range(NumOfFuncs):
                                E0[count,(j)*NumOfFuncs*NumOfFlows+(f)*NumOfFlows+k]=TrafficRate[k]
    count=-1;
    T0=np.zeros((NumOfFlows,NumOfX))
    for k in range(NumOfFlows):
        count=count+1
        for Rank in range(NumOfCores):
            for i in RankStart[Rank]:
                for j in range(NumOfFuncNodes):
                    if FuncNodes[j]==i:
                        if j<= NumOfFuncNodes:
                            f=FunctionChain[k]
                            T0[count,(j)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k]=TrafficRate[k]
    dim=NumOfFuncNodes*NumOfFlows*NumOfFuncs+2*NumOfFlows+4*NumOfFuncNodes*NumOfFuncs
    bmat=np.zeros((dim))
    Xmat=np.zeros((dim,NumOfX))
    count=1;
    for k in range(NumOfFlows):
        f=FunctionChain[k]
        for i in range(NumOfFuncNodes):
            Xmat[count,(i-1)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k]=1
        count=count+1
    for k in range(NumOfFlows):
        f=FunctionChain[k]
        for i in range(NumOfFuncNodes):
            Xmat[count,(i-1)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k]=-1
        count=count+1
    aa=np.ones((NumOfFlows))
    bb=-np.ones((NumOfFlows))
    bmat[0:2*NumOfFlows]=np.append(aa,bb)
    for i in range(NumOfFuncNodes):
        for f in range(NumOfFuncs):
            for k in range(NumOfFlows):
                Xmat[count,(i)*NumOfFuncs*NumOfFlows+(f)*NumOfFlows+k]=1
                Xmat[count,NumOfFuncNodes*NumOfFuncs*NumOfFlows+NumOfFuncNodes+(i)*NumOfFuncs+f]=-1
                count=count+1
    for i in range(NumOfFuncNodes):
        for f in range(NumOfFuncs):
            Xmat[count,NumOfFuncNodes*NumOfFuncs*NumOfFlows+NumOfFuncNodes+(i)*NumOfFuncs+f]=1
            Xmat[count,NumOfFuncNodes*NumOfFuncs*NumOfFlows+i]=-1
            count=count+1
    for i in range(NumOfFuncNodes):
        for f in range(NumOfFuncs):
            for k in range(NumOfFlows):
                Xmat[count,(i)*NumOfFuncs*NumOfFlows+(f)*NumOfFlows+k]=TrafficRate[k]
            Xmat[count,NumOfFuncNodes*NumOfFuncs*NumOfFlows+NumOfFuncNodes+(i)*NumOfFuncs+f]=-ProcessCap[i]
            count=count+1
    for i in range(NumOfFuncNodes):
        for f in range(NumOfFuncs):
            for k in range(NumOfFlows):
                Xmat[count,(i)*NumOfFuncs*NumOfFlows+(f)*NumOfFlows+k]=TrafficRate[k]
        Xmat[count,NumOfFuncNodes*NumOfFuncs*NumOfFlows+i]=-ProcessCap[i]
        count=count+1
    CC=np.zeros((NumOfX))
    p=.4
    index=findsmall(Xikf[SimulationIte,0:NumOfFuncNodes*NumOfFuncs*NumOfFlows])
    print len(index)
    for i in range(len(index)):
        if index[i] >= 0:
            CC[i]=sigma2[SimulationIte]*p*Xikf[SimulationIte,index[i]]**(p-1)
    if SimulationIte == 0:
        x = cvx.Variable(NumOfR)
        y = cvx.Variable(NumOfX)
        obj1 = cvx.Minimize(sum(ObjC*x)+sum(np.dot(I_m_0M[0],DualLinkRateM[0])*x-np.dot(DualLinkRateM[0],I_m_0M[0])*FlowTempM[0])+sum(np.dot(A_m_0M[0],DualXM[0])*y-np.dot(DualXM[0],A_m_0M[0])*XtempM[0])+rho/2*cvx.sum_squares(A_m_0M[0]*y-np.dot(A_m_0M[0],XtempM[0]))+rho/2*cvx.sum_squares(I_m_0M[0]*x-np.dot(I_m_0M[0],FlowTempM[0]))+sum(np.dot(I_m_0M[1],DualLinkRateM[1])*x-np.dot(DualLinkRateM[1],I_m_0M[1])*FlowTempM[1])+sum(np.dot(A_m_0M[1],DualXM[1])*y-np.dot(DualXM[1],A_m_0M[1])*XtempM[1])+rho/2*cvx.sum_squares(A_m_0M[1]*y-np.dot(A_m_0M[1],XtempM[1]))+rho/2*cvx.sum_squares(I_m_0M[1]*x-np.dot(I_m_0M[1],FlowTempM[1])))
        constraints = [E0*y <= ProcessCap, 0<= x, 0<= y,H0*x <= capa1, T0*y == TrafficRate, 1 >= y]
        prob = cvx.Problem(obj1,constraints)
        prob.solve()
        xx=np.asarray(x.value)
        yy=np.asarray(y.value)
    elif SimulationIte == 1:
        x = cvx.Variable(NumOfR)
        y = cvx.Variable(NumOfX)
        obj1 = cvx.Minimize(sum(CC*y)+sum(ObjC*x)+sum(np.dot(I_m_0M[0],DualLinkRateM[0])*x-np.dot(DualLinkRateM[0],I_m_0M[0])*FlowTempM[0])+sum(np.dot(A_m_0M[0],DualXM[0])*y-np.dot(DualXM[0],A_m_0M[0])*XtempM[0])+rho/2*cvx.sum_squares(A_m_0M[0]*y-np.dot(A_m_0M[0],XtempM[0]))+rho/2*cvx.sum_squares(I_m_0M[0]*x-np.dot(I_m_0M[0],FlowTempM[0]))+sum(np.dot(I_m_0M[1],DualLinkRateM[1])*x-np.dot(DualLinkRateM[1],I_m_0M[1])*FlowTempM[1])+sum(np.dot(A_m_0M[1],DualXM[1])*y-np.dot(DualXM[1],A_m_0M[1])*XtempM[1])+rho/2*cvx.sum_squares(A_m_0M[1]*y-np.dot(A_m_0M[1],XtempM[1]))+rho/2*cvx.sum_squares(I_m_0M[1]*x-np.dot(I_m_0M[1],FlowTempM[1])))
        constraints = [Xmat*y <= bmat, E0*y <= ProcessCap, 0<= x, 0<= y,H0*x <= capa1, T0*y == TrafficRate, 1 >= y]
        prob = cvx.Problem(obj1,constraints)
        prob.solve()
        xx=np.asarray(x.value)
        yy=np.asarray(y.value)
    elif SimulationIte == 2:
        x = cvx.Variable(NumOfR)
        y = cvx.Variable(NumOfX)
        obj1 = cvx.Minimize(sum(ObjC*y)+sum(ObjC*x)+sum(np.dot(I_m_0M[0],DualLinkRateM[0])*x-np.dot(DualLinkRateM[0],I_m_0M[0])*FlowTempM[0])+sum(np.dot(A_m_0M[0],DualXM[0])*y-np.dot(DualXM[0],A_m_0M[0])*XtempM[0])+rho/2*cvx.sum_squares(A_m_0M[0]*y-np.dot(A_m_0M[0],XtempM[0]))+rho/2*cvx.sum_squares(I_m_0M[0]*x-np.dot(I_m_0M[0],FlowTempM[0]))+sum(np.dot(I_m_0M[1],DualLinkRateM[1])*x-np.dot(DualLinkRateM[1],I_m_0M[1])*FlowTempM[1])+sum(np.dot(A_m_0M[1],DualXM[1])*y-np.dot(DualXM[1],A_m_0M[1])*XtempM[1])+rho/2*cvx.sum_squares(A_m_0M[1]*y-np.dot(A_m_0M[1],XtempM[1]))+rho/2*cvx.sum_squares(I_m_0M[1]*x-np.dot(I_m_0M[1],FlowTempM[1])))
        constraints = [y[0:NumOfFuncNodes*NumOfFuncs*NumOfFlows]==Xikf[SimulationIte,0:NumOfFuncNodes*NumOfFuncs*NumOfFlows], E0*y <= ProcessCap, 0<= x, 0<= y,H0*x <= capa1, T0*y == TrafficRate, 1 >= y]
        prob = cvx.Problem(obj1,constraints)
        prob.solve()
        xx=np.asarray(x.value)
        yy=np.asarray(y.value)
    print np.shape(xx)
    print np.shape(yy)
    PP=np.zeros((NumOfR))
    VV=np.zeros((NumOfX))
    for ee in range(NumOfR):
        PP[ee]=xx[ee][0]
    for oo in range(NumOfX):
        VV[oo]=yy[oo][0]
    X=VV
    RecvLinkRate=PP
    for t in range(NumOfCores):
        DualXM[t]=DualXM[t]+rho*(np.dot(A_m_0M[t],X)-np.dot(A_m_0M[t],XtempM[t]))
        DualLinkRateM[t]=DualLinkRateM[t]+rho*(np.dot(I_m_0M[t],RecvLinkRate)-np.dot(I_m_0M[t],FlowTempM[t]))
    return (RecvLinkRate,X,DualLinkRateM,DualXM)

sum0=100000
while (sum0 >= .01):
    for Rank in range(NumOfCores):
        NumOfLinks=0
        for i in RankStart[Rank]:
            NumOfLinks=NumOfLinks+NumOfLinksFrom[i]+NumOfLinksFromCrossCore[i]+NumOfLinksToCrossCore[i]
        NumOfLinksM[Rank]=NumOfLinks
        DualLinkRate=DualLinkRateM[Rank]
        DualX=DualXM[Rank]
        del CapaVector

        CapaVector=(NumOfLinks)*[0] #####Number of links in a subnetwork
        Count=0
        Q=np.zeros((NumOfLinks,NumOfR))
        for i in RankStart[Rank]:
            for j in range(NumOfLinksFrom[i-1]):
                for k in range(NumOfFlows):
                    for s in range(ChainLength+1):
                        Q[Count,(LinksFrom[i-1][j]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=1
                CapaVector[Count]=CapacityMatrix[i-1][LinkTableFrom[i-1][j]-1]
                Count=Count+1
        for i in RankStart[Rank]:
            for j in range(NumOfLinksFromCrossCore[i-1]):
                for k in range(NumOfFlows):
                    for s in range(ChainLength+1):
                        Q[Count,(LinksFromCrossCore[i-1][j]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=1
                CapaVector[Count]=CapacityMatrix[i-1][LinkTableFromCrossCore[i-1][j]-1]
                Count=Count+1
        for i in RankStart[Rank]:
            for j in range(NumOfLinksToCrossCore[i-1]):
                for k in range(NumOfFlows):
                    for s in range(ChainLength+1):
                        Q[Count,(LinksToCrossCore[i-1][j]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=1
                CapaVector[Count]=CapacityMatrix[i-1][LinkTableToCrossCore[i-1][j]-1]
                Count=Count+1
        FX=np.zeros((SizeOfF,NumOfX))
        FR=np.zeros((SizeOfF,NumOfR))
        Fb=[0 for i in range(SizeOfF)]
        count=0
        for i in RankStart[Rank]:
            count=count+1
            Index=(count-1)*NumOfFlows*(ChainLength+1)
            for k in range(NumOfFlows):
                for s in range(ChainLength+1):
                    for l in range(NumOfLinksFrom[i-1]):
                        FR[Index+(k)*(ChainLength+1)+s,(LinksFrom[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=1
                    for l in range(NumOfLinksTo[i-1]):
                        FR[Index+(k)*(ChainLength+1)+s,(LinksTo[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=-1
                    for l in range(NumOfLinksFromCrossCore[i-1]):
                        FR[Index+(k)*(ChainLength+1)+s,(LinksFromCrossCore[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=1
                    for l in range(NumOfLinksToCrossCore[i-1]):
                        FR[Index+(k)*(ChainLength+1)+s,(LinksToCrossCore[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s]=-1
                    if s==0:
                        if Pair[k][0]==i:
                            Fb[Index+(k)*(ChainLength+1)+s]=TrafficRate[k]
                    if s==1:
                        if Pair[k][1]==i:
                            Fb[Index+(k)*(ChainLength+1)+s]=-TrafficRate[k]
                for j in range(NumOfFuncNodes):
                    if FuncNodes[j]==i:
                        if j<= NumOfFuncNodes:
                            s=0
                            f=FunctionChain[k]
                            if NodeFuncTable[j][f-1]==1:
                                FX[Index+(k)*(ChainLength+1)+s,(j)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k]=-TrafficRate[k]
                                FX[Index+(k)*(ChainLength+1)+s+1,(j)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k]=TrafficRate[k]
        FRM[Rank]=FR
        FXM[Rank]=FX
        FbM[Rank]=Fb

        Loc=(sum(NumOfLinksToCrossCore))*NumOfFlows*(ChainLength+1)*[0]
        ObjCoeff=np.zeros((len(LinkInfo)*NumOfFlows*(ChainLength+1),2))
        count=0
        count1=0
        count2=0
        for i in RankStart[Rank]:
            if i!=0:
                for k in range(NumOfFlows):
                    for s in range(ChainLength+1):
                        for l in range(NumOfLinksFrom[i-1]):
                            ObjCoeff[count1][0]=((LinksFrom[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s)
                            ObjCoeff[count1][1]=1
                        for l in range(NumOfLinksFromCrossCore[i-1]):
                            ObjCoeff[count1][0]=(LinksFromCrossCore[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s
                            ObjCoeff[count1][1]=.25
                            count1=count1+1
                            Loc[count2]=(LinksFromCrossCore[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s
                            count2=count2+1
                        for l in range(NumOfLinksTo[i-1]):
                            ObjCoeff[count1][0]=(LinksTo[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s
                            ObjCoeff[count1][1]=1
                            count1=count1+1
                        for l in range(NumOfLinksToCrossCore[i-1]):
                            ObjCoeff[count1][0]=(LinksToCrossCore[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s
                            ObjCoeff[count1][1]=.25
                            Loc[count2]=(LinksToCrossCore[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s
                            count2=count2+1
        count=0
        LocX=(NumOfFuncNodes*NumOfFuncs*NumOfFlows+NumOfFuncNodes+NumOfFuncNodes*NumOfFuncs)*[0]
        for i in RankStart[Rank]:
            if i!=0:
                for k in range(NumOfFlows):
                    for j in range(NumOfFuncNodes):
                        if FuncNodes[j]==i:
                            if j<= NumOfFuncNodes:
                                f=FunctionChain[k]
                                LocX[count]=NumOfFuncNodes*NumOfFuncs*NumOfFlows+j
                                count=count+1
        for i in RankStart[Rank]:
            if i!=0:
                for k in range(NumOfFlows):
                    for j in range(NumOfFuncNodes):
                        if FuncNodes[j]==i:
                            if j<= NumOfFuncNodes:
                                f=FunctionChain[k]
                                LocX[count]=(j)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k
                                count=count+1
        for i in RankStart[Rank]:
            if i!=0:
                for k in range(NumOfFlows):
                    for j in range(NumOfFuncNodes):
                        if FuncNodes[j]==i:
                            if j<= NumOfFuncNodes:
                                f=FunctionChain[k]
                                LocX[count]=NumOfFuncNodes*NumOfFuncs*NumOfFlows+NumOfFuncNodes+(j)*NumOfFuncs+f-1
                                count=count+1
        I_m_0=np.zeros((NumOfR,NumOfR))
        for i in Loc:
            I_m_0[i,i]=1
        count=0
        A_m_0=np.zeros((NumOfX,NumOfX))
        for i in LocX:
            A_m_0[i,i]=1
        I_m_0M[Rank]=I_m_0
        A_m_0M[Rank]=A_m_0
        [FlowTemp,Xtemp]=LocalUpdate(I_m_0, LinkTableFrom, NumOfLinksFrom, NumOfCores, source, dest, FunctionChain, FuncNodes, NumOfFuncs, NumOfFuncNodes, beta0, beta, ObjCoeff, Rank, RecvLinkRate, DualLinkRate, FR, FX, Q, Fb, CapaVector, X, DualX, NumOfLinksFromCrossCore, NumOfLinksToCrossCore, LinkTableFromCrossCore, LinkTableToCrossCore, RankStart, SizeOfF, NumOfNodes, ChainLength, NumOfFlows, NumOfLinks, NumOfR, NumOfX)
        print FlowTemp
        for i in range(NumOfR):
            FlowTempM[Rank][i]=FlowTemp[i][0]
        for i in range(NumOfX):
            XtempM[Rank][i]=Xtemp[i][0]

    LocC=sum(NumOfLinksFromCrossCore)*(ChainLength+1)*NumOfFlows*[0]
    count=0
    for Rank in range(NumOfCores):
        for i in RankStart[Rank]:
            for k in range(NumOfFlows):
                for s in range(ChainLength+1):
                    for l in range(NumOfLinksFromCrossCore[i-1]):
                        LocC[count]=(LinksFromCrossCore[i-1][l]-1)*NumOfFlows*(ChainLength+1)+(k)*(ChainLength+1)+s
                        count=count+1

    [RecvLinkRate,X,DualLinkRateM,DualXM]=UpdateC(sigma2, Xikf, LocC, I_m_0M, LocM, LocXM, FunctionChain, FuncNodes, NumOfFuncs, NumOfFuncNodes, beta0, beta, ObjCoeff, DualLinkRate, FR, FX, Q, Fb, CapaVector, DualX, Xtemp, FlowTemp, NumOfLinksFromCrossCore, NumOfLinksToCrossCore, LinkTableFromCrossCore, LinkTableToCrossCore, RankStart, SizeOfF, NumOfNodes, ChainLength, NumOfFlows, NumOfLinks, NumOfR, NumOfX)
    sum0=0
    for uu in range(NumOfCores):
        print np.dot((np.dot(A_m_0M[uu],X)-np.dot(A_m_0M[uu],XtempM[uu])),DualXM[uu])
        sum0=sum0+np.dot((np.dot(A_m_0M[uu],X)-np.dot(A_m_0M[uu],XtempM[uu])),DualXM[uu])+rho/2*LA.norm(np.dot(A_m_0M[uu],X)-np.dot(A_m_0M[uu],XtempM[uu]))**2+np.dot(np.dot(I_m_0M[uu],RecvLinkRate)-np.dot(I_m_0M[uu],FlowTempM[uu]),DualLinkRateM[uu])+rho/2*LA.norm(np.dot(I_m_0M[uu],RecvLinkRate)-np.dot(I_m_0M[uu],FlowTempM[uu]))**2

    SimulationIte=SimulationIte+1
    for i in range(len(X)):
        Xikf[SimulationIte][i]=X[i]

    count=0
    if integer_count(X[0:NumOfFuncNodes*NumOfFuncs*NumOfFlows]) == 0:
        pen=1
    elif SimulationIte ==2:
        for k in range(NumOfFlows):
            AA=[]
            BB=[]
            for Rank in range(NumOfCores):
                for i in RankStart[Rank]:
                    if i!=0:
                        for j in range(NumOfFuncNodes):
                            if FuncNodes[j]==i:
                                if j<= NumOfFuncNodes:
                                    f=FunctionChain[k]
                                    AA.append(X[(j-1)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k])
                                    BB.append((j-1)*NumOfFuncs*NumOfFlows+(f-1)*NumOfFlows+k)
            for ut in range(len(BB)):
                X[BB[ut]]=0
            aa=findmember(AA,max(AA))
            X[BB[aa]]=1
        for ii in range(len(X)):
            Xikf[SimulationIte][ii]=X[ii]






