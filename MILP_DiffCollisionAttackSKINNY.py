import gurobipy  as gb
import numpy as np
import os
import time
import math as mt
np.set_printoptions(threshold=np.inf)

SeperateLine="-"*100

#SKINNY parameters
MC = [[1, 0, 1, 1], [1, 0, 0, 0], [0, 1, 1, 0], [1, 0, 1, 0]]
Inv_MC = [[0, 1, 0, 0], [0, 1, 1, 1], [0, 1, 0, 1], [1, 0, 0, 1]]  # inverse of MC





BPT_SKINNY=[
[1,0,0,0,0,0,0,0,		0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,		0,0,0,1,0,0,0,0],
[0,0,0,0,1,0,0,0,		0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,		0,0,0,0,0,0,0,1],
[0,1,0,0,0,0,0,0,		0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,4,0,		0,0,0,0,0,0,0,1],
[0,0,0,0,0,0,0,0,		0,4,0,1,0,0,0,0],
[0,0,8,0,0,0,4,0,		0,0,0,4,0,0,0,1],

[0,1,0,0,0,0,0,0,		0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,		0,0,4,1,0,0,0,0],
[0,0,0,0,0,1,0,0,		0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,		0,0,0,0,0,0,4,1],
[0,0,0,0,0,0,0,0,		0,0,0,0,4,1,0,0],
[0,0,0,0,0,0,0,4,		0,0,0,0,0,0,4,1],
[0,0,0,0,0,0,0,0,		8,4,0,0,4,1,0,0],
[0,0,0,0,0,8,0,4,		0,0,0,0,8,4,4,1],
]
# print(BPT_SKINNY)

def GetPattern(x,WordSize):
    Pattern=0

    MASK=2**WordSize-1
    for i in range(4):
        if (x>>(i*WordSize ))& MASK!=0:
            Pattern+=1<<i
    return Pattern
def MultiplyMatrix(x,M,WordSize):
    y=0

    MASK = 2 ** WordSize - 1
    for row in range(4):
        temp=0
        for col in range(4):
            if M[row][col]==1:
                temp^=((x>>((3-col)*WordSize))&MASK)
        y^=temp<<((3-row)*WordSize)

    return y






def GetBPT(M,WordSize=4):
    BPT=np.ones((16,16))*(-1)
    BPT[0][0]=0
    for i in range(1,16):
            NumDiffs=0
            for dx in range(2**(WordSize*4)):
                if GetPattern(dx,WordSize)==i:
                        NumDiffs+=1

                        D_y=MultiplyMatrix(dx, M,WordSize)
                        OutPattern=GetPattern(D_y,WordSize)
                        if BPT[i][OutPattern]==-1:
                            BPT[i][OutPattern]=1
                        else:
                            BPT[i][OutPattern]+=1

            print(BPT[i])
            for j in range(1,16):
                if BPT[i][j]>=0:
                    BPT[i][j]=mt.log2(BPT[i][j]/(NumDiffs))
    return BPT



class SKINNY:
    def __init__(self, WordSize=8,TweakeyMode=3 ):
        self.WordSize = WordSize
        self.n = WordSize * 2
        self.TweakeyMode = TweakeyMode
        assert (TweakeyMode>0)&(TweakeyMode<4),"Possible vaules are 1,2,3 !"

        self.MC = [[1,0,1,1],[1,0,0,0],[0,1,1,0],[1,0,1,0]]
        self.Inv_MC = [[0,1,0,0],[0,1,1,1],[0,1,0,1],[1,0,0,1]] #inverse of MC

    def column(self,A, j):
	    return [A[j], A[j+4], A[j+8], A[j+12]]


    def ShiftRow(self,A,Rightshift=True):

        if Rightshift:
            return [A[0], A[1], A[2], A[3],
                A[7], A[4], A[5], A[6],
                A[10],A[11],A[8], A[9],
                A[13],A[14],A[15],A[12]]
        else:
            return [
                A[0], A[1], A[2], A[3],
                 A[5], A[6],A[7], A[4],
                 A[10], A[11],A[8], A[9],
                A[15],A[12],A[13],A[14]
            ]



    def TD_RoundFunction(self, Model,A, B, Pr):
        #Model.addConstr(gb.quicksum(A)<=15)

        ShiftRow_A=self.ShiftRow(A)
        for i in range(4):
            TEMP=[ ShiftRow_A[(i+4*j)%16 ] for j in range(4) ]
            TEMP+= [B[ (i + 4 * j)%16] for j in range(4)]
            TEMP +=[Pr[i][0],Pr[i][1]]

            S_T = [(0, 1, 2, 1, -1, 1, -1, 0, -4, -2, 0), (0, 0, -2, -1, 1, 0, 1, 1, 2, 1, 0),
                   (0, -1, 0, 0, 1, 0, 1, -1, 2, 1, 0), (0, -1, 1, 0, 0, -1, 1, 1, 1, 0, 0),
                   (0, 0, 0, 1, -1, 0, -2, -1, -4, -2, 4), (0, 1, 1, 0, 0, 0, -1, 0, -2, 0, 0),
                   (0, 0, -2, 0, 0, 0, 1, 1, 2, 1, 0), (0, 1, -1, 0, 0, 0, 1, 0, 0, 0, 0),
                   (0, 0, -2, 0, 1, 0, 1, 0, 2, 1, 0), (0, 0, 1, 1, -1, 1, -2, 0, -4, -2, 2),
                   (0, 1, 2, 0, -2, 1, -1, -1, -4, -2, 2), (0, -1, 0, -1, 1, -1, 1, 1, 2, 1, 0)]

            # S_T = [[2, 1, 1, 1, 1, -2, -1, -1, -2, 0, 0],
            #         [-2, 0, 2, 1, -1, 2, 0, 1, -2, -2, 0],
            #         [2, 0, -2, -1, 1, -2, 1, 1, 2, 1, 0],
            #         [-2, -1, 0, 0, 1, 2, 1, -1, 2, 1, 0],
            #         [-2, 2, -2, 0, -1, 2, 1, -1, -2, -1, 3],
            #         [1, 2, 1, -1, 1, -2, -2, 2, -2, 2, 0],
            #         [2, 0, 1, 0, -2, -1, -2, -1, -2, -2, 4],
            #         [-1, 0, -1, 1, -1, 2, 0, 2, 0, 0, 0],
            #         [2, -2, 2, 0, 1, -1, 2, -2, 2, 1, 0],
            #         [-2, -2, -1, 1, 2, 2, 1, -1, 2, 1, 1],
            #         [1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 0]]
            #
            # S_T = [[1, 2, 4, 2, -2, 1, -2, 0, -8, -4, 0],
            #         [-3, 0, -8, -4, 5, 1, 5, 4, 10, 5, 0],
            #         [10, -2, -1, 1, 5, -10, 3, -5, 7, 4, 0],
            #         [-10, -2, 4, 2, -2, 7, 2, 5, 10, -1, 0],
            #         [10, 2, 3, 0, -4, -9, -2, -3, -8, -4, 6],
            #         [-9, 0, -3, -1, 1, 10, -10, 3, -8, 2, 10],
            #         [-1, 1, -1, 1, 1, -1, -1, -1, 1, -1, 5]]

             # 8p0 , 4p1
            S_T_Back = [(0, 0, -1, 0, 0, 1, 1, 0, 0, 0, 0), (1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0),
                        (0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0), (0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                        (0, 0, 0, 0, 0, -1, -1, -1, 0, -1, 3), (0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 1),
                        (-1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), (0, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0),
                        (1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0), (0, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0),
                        (0, 0, 0, -1, 0, 1, 0, 1, 0, -1, 1), (0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0),
                        (0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0), (1, 0, 0, 0, 0, -1, 0, 1, 0, -1, 1),
                        (0, 0, 0, 0, 0, 1, -1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0),
                        (0, 0, 0, 0, -1, 0, 1, 0, 1, 1, 0), (0, 0, 0, -1, 0, 0, 0, 1, 1, 1, 0),
                        (0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 2), (0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0)]

            for equ in S_T:
                TempSum=0
                for j in range(len(equ)):
                    if j==(len(equ)-1):
                        if equ[j] != 0:
                            TempSum +=equ[j]
                        break
                    if equ[j]!=0:

                        TempSum+=equ[j]*TEMP[j]
                Model.addConstr(TempSum>=0)
            Model.addConstr(ShiftRow_A[i]==B[ (i + 4 * 1)%16])



    def TD_RoundFunction_Invserse(self, Model,A, B, Pr):
        #Model.addConstr(gb.quicksum(A)<=15)

        ShiftRow_A=self.ShiftRow(A)
        for i in range(4):
            TEMP=[B[ (i + 4 * j)%16] for j in range(4)]
            TEMP+= [ ShiftRow_A[(i+4*j)%16 ] for j in range(4) ]
            TEMP +=[Pr[i][0],Pr[i][1]]

            S_T = [(0, 1, 2, 1, -1, 1, -1, 0, -4, -2, 0), (0, 0, -2, -1, 1, 0, 1, 1, 2, 1, 0),
                   (0, -1, 0, 0, 1, 0, 1, -1, 2, 1, 0), (0, -1, 1, 0, 0, -1, 1, 1, 1, 0, 0),
                   (0, 0, 0, 1, -1, 0, -2, -1, -4, -2, 4), (0, 1, 1, 0, 0, 0, -1, 0, -2, 0, 0),
                   (0, 0, -2, 0, 0, 0, 1, 1, 2, 1, 0), (0, 1, -1, 0, 0, 0, 1, 0, 0, 0, 0),
                   (0, 0, -2, 0, 1, 0, 1, 0, 2, 1, 0), (0, 0, 1, 1, -1, 1, -2, 0, -4, -2, 2),
                   (0, 1, 2, 0, -2, 1, -1, -1, -4, -2, 2), (0, -1, 0, -1, 1, -1, 1, 1, 2, 1, 0)]
            # 8p0 , 4p1

            S_T_Back = [(0, 0, -1, 0, 0, 1, 1, 0, 0, 0, 0), (1, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0),
                        (0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0), (0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                        (0, 0, 0, 0, 0, -1, -1, -1, 0, -1, 3), (0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 1),
                        (-1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0), (0, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0),
                        (1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0), (0, 0, 1, 0, 0, 1, -1, 0, 0, 0, 0),
                        (0, 0, 0, -1, 0, 1, 0, 1, 0, -1, 1), (0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0),
                        (0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0), (1, 0, 0, 0, 0, -1, 0, 1, 0, -1, 1),
                        (0, 0, 0, 0, 0, 1, -1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0),
                        (0, 0, 0, 0, -1, 0, 1, 0, 1, 1, 0), (0, 0, 0, -1, 0, 0, 0, 1, 1, 1, 0),
                        (0, 0, 0, 0, 0, -1, -1, 0, -1, 0, 2), (0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0)]

            for equ in S_T_Back:
                TempSum=0
                for j in range(len(equ)):
                    if j==(len(equ)-1):
                        if equ[j] != 0:
                            TempSum +=equ[j]
                        break
                    if equ[j]!=0:

                        TempSum+=equ[j]*TEMP[j]
                Model.addConstr(TempSum>=0)
            Model.addConstr(ShiftRow_A[i]==B[ (i + 4 * 1)%16])


    def SearchTD(self,Round,ExpectedNumSol=1,InOutDiff=[]):

        s = gb.Model()

        # Truncated difference representation
        TD_Values=[
           [
               s.addVar(vtype=gb.GRB.BINARY, name="TD_{}_R{}".format(i,k))
             for i in range(16) ]
            for k in range(Round+1)
        ]
        # Probabilistic weight representation
        Pr_Values=[
            [s.addVars(2,vtype=gb.GRB.BINARY, name="Pr_{}_R{}".format(i,k)) for i in range(4)] for k in range(Round)
        ]
        # Set the objective function
        Obj = s.addVar(lb=0, ub=self.TweakeyMode*16*self.WordSize+100, vtype=gb.GRB.INTEGER, name="obj")
        s.update()

        #2. Define the relationship between wheel functions

        for r in range(Round):
            InTD=TD_Values[r]
            OutTD=TD_Values[r+1]
            Pr=Pr_Values[r]
            self.TD_RoundFunction(s,InTD,OutTD,Pr)

        # Input is not 0
        s.addConstr(gb.quicksum(TD_Values[0])>=1)
        # Input cannot be all 1s.
        s.addConstr(gb.quicksum(TD_Values[0])<=15)

        # The following is a constraint that specifically illustrates the statistical significance of truncated difference.
        s.addConstr((gb.quicksum(TD_Values[-1] ) ) * self.WordSize + Obj <= (16 * self.WordSize - 1))

        # The following defines the truncated input and output differences.

        if len(InOutDiff)!=0:
            if InOutDiff[0]!=-1:
                for i in range(16):
                    s.addConstr(TD_Values[0][i]== ((InOutDiff[0]>>(15-i))&1))
            if InOutDiff[1]!=-1:
                for i in range(16):
                    s.addConstr(TD_Values[-1][i]== ((InOutDiff[1]>>(15-i))&1))


        ##3. Solve the problem and show the result.
        TempSum=0
        for r in range(Round):
            for i in range(4):
                TempSum+=self.n*Pr_Values[r][i][0]
        for r in range(Round):
            for i in range(4):
                TempSum+=self.WordSize*Pr_Values[r][i][1]




        s.addConstr(Obj==TempSum)
        s.setObjective(Obj, sense=gb.GRB.MINIMIZE)




        print("Start to solve the model...")
        try:
            s.setParam(gb.GRB.Param.PoolSolutions, ExpectedNumSol)
            s.setParam(gb.GRB.Param.PoolSearchMode, 2)

            s.update()

            StartTime = time.time()
            s.optimize()
            EndTime = time.time()
            print("Time used:", EndTime - StartTime)

            status = s.Status

            if status != gb.GRB.OPTIMAL:
                print('Optimization was stopped with status ' + str(status))

            nSolution = s.SolCount
            print('number of solution stored:' + str(nSolution))


            x = s.getVars()
            AllSolution = []
            # Print objective values of solutions
            for e in range(nSolution):
                tempDictionary = {}
                s.setParam(gb.GRB.Param.SolutionNumber, e)

                for i in range(len(x)):
                    tempDictionary[str(x[i].VarName)] = (x[i].Xn);
                AllSolution.append(tempDictionary)

                # Print solution values
                print(SeperateLine)
                print('Solution:\t', e)
                print("Obj:\t",Obj.x)

                print("TD_Values:")
                for r in range(len(TD_Values)):

                    for i in range(16):
                        print(int(TD_Values[r][i].x), end="")
                    print()


                print("Pr_Values:")
                for r in range(len(Pr_Values)):
                    for i in range(4):
                        print(int(Pr_Values[r][i][0].x),int(Pr_Values[r][i][1].x), end="\t\t")
                    print()



        except gb.GurobiError as e:
            print('Error reported')

    def GetDiffRoundbyRound_Upward(self,Model,A,B,MC):
        SHIFT_ROW_A=self.ShiftRow(A)

        for col in range(4):
            In=self.column(SHIFT_ROW_A,col)
            Out=self.column(B,col)


            for i in range(4):
                for j in range(4):
                    if MC[i][j]==1:
                        Model.addConstr(Out[i]<=In[j])


                TempSum=0
                for j in range(4):
                    if MC[j][i]==1:
                        TempSum+=Out[j]
                Model.addConstr(TempSum>=In[i])

    def GetDiffRoundbyRound_Up(self,Model,A,B,MC):
        SHIFT_ROW_A=self.ShiftRow(A)

        for col in range(4):
            Out=self.column(B,col)
            In=self.column(SHIFT_ROW_A,col)

            for i in range(4):
                for j in range(4):
                    if MC[i][j]==1:
                        Model.addConstr(Out[i]<=In[j])


                TempSum=0
                for j in range(4):
                    if MC[j][i]==1:
                        TempSum+=Out[j]
                Model.addConstr(TempSum>=In[i])



    def GetDiffRoundbyRound_Low(self,Model,A,B,MC):
        SHIFT_ROW_A=self.ShiftRow(A)

        for col in range(4):
            In=self.column(B,col)
            Out=self.column(SHIFT_ROW_A,col)



            for i in range(4):
                for j in range(4):
                    if MC[i][j]==1:
                        Model.addConstr(Out[i]<=In[j])

                #输入字块活跃，必须要有至少有一个影响的输出字块活跃
                TempSum=0
                for j in range(4):
                    if MC[j][i]==1:
                        TempSum+=Out[j]
                Model.addConstr(TempSum>=In[i])

    def GetRelationofK(self,Round,AdditionalKeyatRound=-1):
        PT = [9, 15, 8, 13, 10, 14, 12, 11, 0, 1, 2, 3, 4, 5, 6, 7]
        AllIndices = [[i for i in range(16)]]

        for r in range(Round - 1):
            TempIndice = AllIndices[-1]
            NewIndice = [0] * 16
            for i in range(16):
                NewIndice[i] = TempIndice[PT[i]]
            AllIndices.append(NewIndice)

        if AdditionalKeyatRound==-1:
            return AllIndices
        else:
            LastR_Index=(AdditionalKeyatRound+Round)%Round

            Reorder_AllIndices=[]
            for r in range(LastR_Index+1):
                Reorder_AllIndices.append(AllIndices[r])
            Reorder_AllIndices=AllIndices[LastR_Index+1:Round]+Reorder_AllIndices
            return Reorder_AllIndices





    def GetLeftKofEc(self,Key):
        return [
            Key[0],Key[1],Key[2],Key[3],
            Key[0], Key[1], Key[2], Key[3],
            Key[7], Key[4], Key[5], Key[6],
            Key[0], Key[1], Key[2], Key[3],


        ]

    def DisplayArray(self, Array):
        NumRow=len(Array)//4

        for i in range(NumRow):
            for j in range(4):
                print(Array[i * 4 + j], end="")
            print()

    def DisplayMatrices(self, Matrices):
        LINE = "-" * 30

        for r in range(len(Matrices[0])):
            print(LINE)
            print("Round", r)
            for i in range(len(Matrices)):
                if r <= (len(Matrices[i]) - 1):
                    print(LINE)
                    self.DisplayArray(Matrices[i][r])

    def Transform(self, X,WordSize=16):
        return [[int(X[r][i].x) for i in range(WordSize)] for r in range(len(X))]


    def AddCons(self,s,Vars,Equ):
        Temp=0
        for i in range(len(Vars)):
            if Equ[i]!=0:
                Temp+=Equ[i]*Vars[i]
        if Equ[-1]!=0:
            Temp+=Equ[-1]
        s.addConstr(Temp>=0)

    def OperationOR(self,s,List):
        Num=len(List)
        TempSum=0
        for i in range(Num-1):
            s.addConstr(List[i]  <=List[-1])
            TempSum+=List[i]
        s.addConstr(TempSum-List[-1]>=0)



    def SearchAttack(self,R_D,R_B,R_F,InOutDiffPr=[],ExpectedNumSol=1,NoImprovedTechnique=False,AdditionalKeyatRound=-1,TimeLimit=-1):

        TD_OR_Diff=(len(InOutDiffPr)<3)#True means TD-MITM attack, otherwise it is Diff-MITM attack.



        r_Additional=2

        r0=R_B
        r1=R_B+R_D
        r2=R_B+R_D+R_F
        r3=R_B+R_D+R_F+r_Additional


        s = gb.Model()



        TD_Values = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="TD_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_D+1)
        ]

        Pr_Values = [
            [s.addVars(2, vtype=gb.GRB.BINARY, name="Pr_{}_R{}".format(i, k)) for i in range(4)] for k in range(R_D)
        ]


        U_Diff=[
            [
                s.addVar(vtype=gb.GRB.BINARY, name="U_Diff_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_B + 1)
        ]

        U_P_DerivedFromP = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="U_P_DerivedFromP_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_B )
        ]

        U_P_Exhaustion = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="U_P_Exhaustion_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_B)
        ]

        U_P_Current = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="U_P_Current_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_B)
        ]




        U_Pr = [
            [s.addVars(2, vtype=gb.GRB.BINARY, name="U_Pr_{}_R{}".format(i, k)) for i in range(4)] for k in range(R_B)
        ]




        U_K = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="U_K_{}_R{}".format(i, k))
                for i in range(8)]
            for k in range(R_B-1)
        ]

        if not TD_OR_Diff:
            U_K_LastR =[
                s.addVar(vtype=gb.GRB.BINARY, name="U_K_LastR_{}_R{}".format(i,R_B-1))
                for i in range(8)]

            L_K_LastR = [
                s.addVar(vtype=gb.GRB.BINARY, name="L_K_LastR_{}_R{}".format(i, R_F - 1))
                for i in range(8)]




        L_Diff = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="L_Diff_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_F + 1)
        ]
        L_Pr = [
            [s.addVars(2, vtype=gb.GRB.BINARY, name="L_Pr_{}_R{}".format(i, k)) for i in range(4)] for k in range(R_F)
        ]

        L_P_DerivedFromP = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="L_P_DerivedFromP_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_F)
        ]

        L_P_Exhaustion = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="L_P_Exhaustion_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_F)
        ]

        L_P_Current = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="L_P_Current_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_F)
        ]


        L_K = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="L_K_{}_R{}".format(i, k))
                for i in range(8)]
            for k in range(R_F - 1)
        ]



        #The keys involved in Ec are dealt with separately below (here Ec involves 2 rounds of keys by default).
        Extended_Kin=[
            [
                s.addVar(vtype=gb.GRB.BINARY, name="Extended_Kin_{}_R{}".format(i, k))
                for i in range(8)]
            for k in range(r_Additional)
        ]

        Extended_Kout = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="Extended_Kout_{}_R{}".format(i, k))
                for i in range(8)]
            for k in range(r_Additional)
        ]

        RealKnownKey_Kin = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="RealKnownKey_Kin_{}_R{}".format(i, k))
                for i in range(8)]
            for k in range(r_Additional)
        ]

        RealKnownKey_Kout = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="RealKnownKey_Kout_{}_R{}".format(i, k))
                for i in range(8)]
            for k in range(r_Additional)
        ]






        LinearRelationship=[
            s.addVar(vtype=gb.GRB.BINARY, name="LinearRelationship_{}".format(i))
            for i in range(8)
        ]

        ReducedData = [
            s.addVar(vtype=gb.GRB.BINARY, name="ReducedData_{}".format(i))
            for i in range(8)
        ]


        #DirectLink indicates a relationship in which states are directly connected.
        DirectLink=[
            s.addVar(vtype=gb.GRB.BINARY, name="DirectLink_{}".format(i))
            for i in range(16)
        ]
        #XOR_Link indicates that the current byte is XOR with an unknown key byte.
        XOR_Link = [
            s.addVar(vtype=gb.GRB.BINARY, name="XOR_Link_{}".format(i))
            for i in range(16)
        ]

        NumFixed=s.addVar(lb=0, ub=16, vtype=gb.GRB.INTEGER, name="NumFixed")
        NumLinear=s.addVar(lb=0, ub=16, vtype=gb.GRB.INTEGER, name="NumLinear")



        ##TEMP_MK_Kin and TEMP_MK_Kout count the number of occurrences of each key block.
        TEMP_MK_Kin = [
            s.addVar(lb=0, ub=10, vtype=gb.GRB.INTEGER, name="TEMP_MK_Kin_" + str(i))
            for i in range(16)
        ]
        TEMP_MK_Kout = [
            s.addVar(lb=0, ub=10, vtype=gb.GRB.INTEGER, name="TEMP_MK_Kout_" + str(i))
            for i in range(16)
        ]


        MK_Kin=[
            s.addVar(lb=0, ub=self.TweakeyMode,vtype=gb.GRB.INTEGER, name="MK_Kin_"+str(i))
            for i in range(16)
        ]
        MK_Kout = [
            s.addVar(lb=0, ub=self.TweakeyMode, vtype=gb.GRB.INTEGER, name="MK_Kout_" + str(i))
            for i in range(16)
        ]

        MK_Known = [
            s.addVar( vtype=gb.GRB.BINARY, name="MK_Known_" + str(i))
            for i in range(16)
        ]

        MK_Kin_Known = [
            s.addVar(vtype=gb.GRB.BINARY, name="MK_Kin_Known_" + str(i))
            for i in range(16)
        ]

        MK_Kout_Known = [
            s.addVar(vtype=gb.GRB.BINARY, name="MK_Kout_Known_" + str(i))
            for i in range(16)
        ]






        Num_K_in = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="Num_K_in")
        Num_K_out = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="Num_K_out")

        IntersecionofKinAndKout=[
            s.addVar(lb=-self.TweakeyMode, ub=10, vtype=gb.GRB.INTEGER, name="IntersecionofKinAndKout_{}".format(i))
            for i in range(16)
        ]

        RealIntersecionofKinAndKout = [
            s.addVar(lb=0, ub=self.TweakeyMode, vtype=gb.GRB.INTEGER,
                     name="RealIntersecionofKinAndKout_{}".format(i))
            for i in range(16)
        ]


        NumofIntersectionofKinAndKout= s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="NumofIntersectionofKinAndKout")




        #The following are the settings of key parameters

        P_in = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER,name="P_in")
        P_out = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="P_out")
        #Exhaustive number of States
        S_in = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="S_in")
        S_out = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="S_out")
        ## Number of blocks for truncated differential input and output blocks
        Delta_in = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="Delta_in")
        Delta_out = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="Delta_out")
        DataOneTime=s.addVar(lb=0, ub=16*self.WordSize, vtype=gb.GRB.INTEGER, name="DataOneTime")


        UPBOUND=self.TweakeyMode * 16 * self.WordSize + 100

        TotalTime=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="TotalTime")
        Time1=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Time1")
        Time2=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Time2")
        Time3=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Time3")
        Pr_Obj = s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Pr_Obj")

        s.update()

        #2. Establish the relationship on the truncated difference partition.
        if TD_OR_Diff:

            for r in range(R_D):
                InTD=TD_Values[r]
                OutTD=TD_Values[r+1]
                Pr=Pr_Values[r]
                self.TD_RoundFunction(s,InTD,OutTD,Pr)

            s.addConstr(gb.quicksum(TD_Values[0])>=1)


            TempSum = 0
            for r in range(R_D):
                for i in range(4):
                    TempSum += self.n * Pr_Values[r][i][0]
            for r in range(R_D):
                for i in range(4):
                    TempSum += self.WordSize * Pr_Values[r][i][1]

            s.addConstr(Pr_Obj == TempSum)





            for i in range(16):
                if len(InOutDiffPr)>0:
                    InDiff = InOutDiffPr[0]
                    s.addConstr(TD_Values[0][i] == ((InDiff >> (15 - i)) & 1))
                if len(InOutDiffPr)>1:
                    OutDiff = InOutDiffPr[1]
                    s.addConstr(TD_Values[-1][i] == ((OutDiff >> (15 - i)) & 1))


        else:
            s.addConstr(Pr_Obj == int(InOutDiffPr[-1]))
            InDiff=InOutDiffPr[0]
            OutDiff=InOutDiffPr[1]

            for i in range(16):
                s.addConstr(TD_Values[0][i]==((InDiff>>(15-i))&1))
                s.addConstr(TD_Values[-1][i] == ((OutDiff >> (15 - i)) & 1))

            for r in range(R_D):
                InTD=TD_Values[r]
                OutTD=TD_Values[r+1]
                Pr=Pr_Values[r]
                for i in range(16):
                    if r !=0:
                        s.addConstr(InTD[i] ==0)
                    if r!=(R_D-1):
                        s.addConstr(OutTD[i]==0)
                for i in range(4):
                    for j in range(2):
                        s.addConstr(Pr[i][j]==0)








        # The following is the connection relationship between the truncated differential discriminator and both sides.
        for i in range(16):
            s.addConstr(U_Diff[-1][i]==TD_Values[0][i])
            s.addConstr(L_Diff[0][i] == TD_Values[-1][i])

        # 3. Give the definite relation of k_in.
        # 3.1 Determining the Differential Propagation Relationship
        for r in range(R_B):
            InTD=U_Diff[r]
            OutTD=U_Diff[r+1]
            Pr=U_Pr[r]
            self.TD_RoundFunction_Invserse(s,InTD,OutTD,Pr)
        # 3.2 Determine the state propagation relationship
        # establish the derivation relation of U_P_DerivedFromP
        for r in range(R_B-1):
            A=U_P_DerivedFromP[r]
            B=U_P_Current[r+1]
            self.GetDiffRoundbyRound_Up( s, A, B, self.MC)
        A=U_P_DerivedFromP[-1]
        B=U_Diff[-1]
        self.GetDiffRoundbyRound_Up(s, A, B, self.MC)



        # Establish the derivation relation of U_P_Exhaustion
        for r in range(R_B):
            A=U_P_Exhaustion[r]
            B=U_P_DerivedFromP[r]
            C=U_Diff[r]
            for i in range(16):## Ensure that if the current block does not need to know, it cannot be exhaustive.
                s.addConstr(A[i]<=B[i]+C[i])
        # establish the derivation relation of U_P_Current.
        for r in range(R_B):
            A=U_P_Current[r]
            B=U_P_Exhaustion[r]
            C=U_Diff[r]
            D=U_P_DerivedFromP[r]
            for i in range(16):
                s.addConstr(A[i]<=1-B[i])
                s.addConstr(B[i]+ (1-C[i])+ A[i]>=1)
                s.addConstr(B[i] + (1 - D[i]) + A[i] >= 1)
                s.addConstr( C[i] + D[i] +(1- A[i]) >= 1)





        for r in range(R_F):
            InTD=L_Diff[r]
            OutTD=L_Diff[r+1]
            Pr=L_Pr[r]
            self.TD_RoundFunction(s,InTD,OutTD,Pr)

        for r in range(R_F - 1):
                A = L_P_Current[r]
                B = L_P_DerivedFromP[r+1]
                self.GetDiffRoundbyRound_Low(s, A, B, self.Inv_MC)
        for i in range(16):
                s.addConstr(L_P_DerivedFromP[0][i] == 0)


        for r in range(R_F):
                A = L_P_Exhaustion[r]
                B = L_P_DerivedFromP[r]
                C = L_Diff[r]
                for i in range(16):
                    s.addConstr(A[i] <= B[i] + C[i])

        for r in range(R_F):
                A = L_P_Current[r]
                B = L_P_Exhaustion[r]
                C = L_Diff[r]
                D = L_P_DerivedFromP[r]
                for i in range(16):
                    s.addConstr(A[i] <= 1 - B[i])


                    s.addConstr(B[i] + (1 - C[i]) + A[i] >= 1)
                    s.addConstr(B[i] + (1 - D[i]) + A[i] >= 1)
                    s.addConstr(C[i] + D[i] + (1 - A[i]) >= 1)

        for r in range(R_B-1):
            for i in range(4):

                s.addConstr(U_K[r][i]>=U_P_Current[r+1][i])
                s.addConstr(U_K[r][i] >= U_P_Current[r + 1][i+4])
                s.addConstr(U_K[r][i] >= U_P_Current[r + 1][i+12])
                s.addConstr(U_P_Current[r+1][i]+U_P_Current[r+1][i+4]+U_P_Current[r+1][i+12] - U_K[r][i]>= 0)
                s.addConstr(U_K[r][i+4] ==U_P_Current[r + 1][(i+1)%4+8])

        for r in range(R_F-1):
            for i in range(8):
                s.addConstr(L_K[r][i]==L_P_Current[r+1][i])


        if not TD_OR_Diff:
            InD=TD_Values[0]
            OutD=TD_Values[-1]

            for i in range(4):

                s.addConstr(U_K_LastR[i] >= InD[i])
                s.addConstr(U_K_LastR[i] >= InD[i + 4])
                s.addConstr(U_K_LastR[i] >= InD[i + 12])
                s.addConstr(InD[i] + InD[i + 4] + InD[i + 12] - U_K_LastR[i] >= 0)
                s.addConstr(U_K_LastR[i + 4] == InD[(i + 1) % 4 + 8])

            for i in range(8):
                s.addConstr(L_K_LastR[i]==OutD[i])


        #5.Calculate the intersection of k_in and k_out.

        AllIndices=self.GetRelationofK(r3,AdditionalKeyatRound)
        print("AllIndices")
        print(AllIndices)

        for i in range(16):
            TempSum=0
            for j in range(R_B-1):
                for k in range(8):
                    if AllIndices[R_B-2-j][k]==i:
                        TempSum+=U_K[R_B-2-j][k]
            if not TD_OR_Diff:
                for k in range(8):
                    if AllIndices[R_B-1][k]==i:
                        TempSum+=U_K_LastR[k]

            for j in range(r_Additional):
                for k in range(8):
                    if AllIndices[-2+j][k]==i:
                        TempSum+=Extended_Kin[j][k]
            s.addConstr(TEMP_MK_Kin[i]==TempSum)

            #SKINNY algorithm key generation algorithm properties, determine the number of each key block.
            s.addConstr(MK_Kin[i]==gb.min_(TEMP_MK_Kin[i],self.TweakeyMode))

            TempSum = 0
            for j in range(R_F-1):
                for k in range(8):
                    if AllIndices[r1+1+j][k]==i:
                        TempSum+=L_K[j][k]
            if not TD_OR_Diff:
                for k in range(8):
                    if AllIndices[r1][k] == i:
                        TempSum += L_K_LastR[k]


            for j in range(r_Additional):
                for k in range(8):
                    if AllIndices[-2 + j][k] == i:
                        TempSum += Extended_Kout[j][k]
            s.addConstr(TEMP_MK_Kout[i] == TempSum)


            s.addConstr(MK_Kout[i] == gb.min_(TEMP_MK_Kout[i], self.TweakeyMode))
            s.addConstr(IntersecionofKinAndKout[i]==MK_Kin[i]+MK_Kout[i]-self.TweakeyMode)
            s.addConstr(RealIntersecionofKinAndKout[i]==gb.max_(IntersecionofKinAndKout[i],0))



            #Calculate whether the I-th block of the master key is completely known.
            s.addConstr(1<=self.TweakeyMode+MK_Known[i]*20-MK_Kout[i]-MK_Kin[i])
            s.addConstr(  MK_Kout[i] +MK_Kin[i]-self.TweakeyMode +(1-MK_Known[i]) *20>=0 )

            #Calculate whether the master key is known under Kin.

            s.addConstr(1 <= self.TweakeyMode - TEMP_MK_Kin[i] + MK_Kin_Known[i] * 20 )
            s.addConstr( TEMP_MK_Kin[i] - self.TweakeyMode + (1 - MK_Kin_Known[i]) * 20 >= 0)

            #Calculate whether the master key is known under Kout.
            s.addConstr(1 <= self.TweakeyMode - TEMP_MK_Kout[i] + MK_Kout_Known[i] * 20 )
            s.addConstr(TEMP_MK_Kout[i] - self.TweakeyMode + (1 - MK_Kout_Known[i]) * 20 >= 0)

        #Require that the keys Extended_Kout and Extended_Kin on Ec do not intersect.
        for r in range(r_Additional):
            for i in range(8):
                s.addConstr(Extended_Kout[r][i]+Extended_Kin[r][i] <=1)

        s.addConstr(gb.quicksum(RealIntersecionofKinAndKout)==NumofIntersectionofKinAndKout)
        s.addConstr(Num_K_in==gb.quicksum(MK_Kin))
        s.addConstr(Num_K_out==gb.quicksum(MK_Kout))



        s.addConstr(NumLinear == gb.quicksum(XOR_Link) )#
        s.addConstr(NumFixed==gb.quicksum(DirectLink)+gb.quicksum(LinearRelationship))

        #(1)To determine the relationship between DirectLink and XOR_Link, a simpler derivation is needed.

        for i in range(8):
            #低8字块  不是直接连接关系就是XOR关系
            A=XOR_Link[i+8]
            B=DirectLink[i+8]
            s.addConstr(A + B == 1)



            Index_K = self.GetLeftKofEc(AllIndices[-2])[8 + i]
            C = MK_Known[Index_K]
            D=self.GetLeftKofEc(Extended_Kout[0])[8+i]
            E=self.GetLeftKofEc(Extended_Kin[0])[8+i]
            self.OperationOR(s, [C, D, E,B])


        for i in range(8):
                Index_K = AllIndices[-2][i]
                A= MK_Kin_Known[Index_K]
                B=Extended_Kin[0][i]
                C=RealKnownKey_Kin[0][i]
                self.OperationOR(s, [A,B,C])

                D=MK_Kout_Known[Index_K]
                E=Extended_Kout[0][i]
                F=RealKnownKey_Kout[0][i]
                self.OperationOR(s, [D, E, F])

                Index_K = AllIndices[-1][i]

                A = MK_Kin_Known[Index_K]
                B = Extended_Kin[1][i]
                C = RealKnownKey_Kin[1][i]
                self.OperationOR(s, [A, B, C])

                D = MK_Kout_Known[Index_K]
                E = Extended_Kout[1][i]
                F = RealKnownKey_Kout[1][i]
                self.OperationOR(s, [D, E, F])




        #Determine the relationship between DirectLink and XOR_Link on the upper 8 nibble.
        for i in range(8):
            XOR=XOR_Link[i]
            Direct=DirectLink[i]
            Index_K_Left = self.GetLeftKofEc(AllIndices[-2])[i]
            MK_Left = MK_Known[Index_K_Left]
            Index_K_Right = AllIndices[-1][i]
            MK_Right = MK_Known[Index_K_Right]

            Kin_Left =  self.GetLeftKofEc(RealKnownKey_Kin[0]) [i]
            Kout_Left = self.GetLeftKofEc(RealKnownKey_Kout[0])[i]
            Kin_Right = RealKnownKey_Kin[1] [i]
            Kout_Right = RealKnownKey_Kout[1] [i]


            ## The relationship between XOR_Link[i] and DirectLink[i] is discussed in categories below.
            # There are two ways to write MILP programs: (1) discuss them according to conditions, and (2) directly use feasible points to generate inequalities.

            ## The following is a classified discussion. When the first situation is not satisfied, the second situation will be discussed, and so on.
            #(1) At least one of Kin_Right or Kout_Left on both sides of an S-box must be 1, otherwise the filter cannot be constructed.
            s.addConstr(Kout_Left+Kin_Right+(1-XOR)>=1)
            s.addConstr(Kout_Left + Kin_Right + (1 - Direct) >= 1)
            ## If it is not the first case, there is XOR_Link[i]+DirectLink[i]==1. After adding the following constraints, only XOR constraints will be generated later.
            s.addConstr((1-Kout_Left)+XOR+ Direct  >=1)
            s.addConstr((1 - Kout_Left) + (1-XOR) + (1-Direct) >= 1)
            s.addConstr((1 - Kin_Right) + XOR + Direct >= 1)
            s.addConstr((1 - Kin_Right) + (1 - XOR) + (1 - Direct) >= 1)

            #(2) If one side has Kin_Right or Kout_Left with 1, as long as the key of the other side is known on the master key, or K_in,K_out, a DirectLink filter can be formed.
            s.addConstr(1-Kout_Left + (1-Kout_Right) + (1 - XOR) >= 1)
            s.addConstr(1 - Kout_Left + (1 - Kin_Right) + (1 - XOR) >= 1)
            s.addConstr(1 - Kout_Left + (1 - MK_Right) + (1 - XOR) >= 1)

            s.addConstr(1 - Kin_Right + (1 - Kout_Left) + (1 - XOR) >= 1)
            s.addConstr(1 - Kin_Right + (1 - Kin_Left) + (1 - XOR) >= 1)
            s.addConstr(1 - Kin_Right + (1 - MK_Left) + (1 - XOR) >= 1)

            #(3) If both situations are not satisfied, you need to build an XOR_Link filter whose
            s.addConstr(1 - Kout_Left +Kout_Right+Kin_Right+MK_Right +  XOR >= 1)
            s.addConstr(1 - Kin_Right + Kout_Left+Kin_Left+MK_Left +  XOR >= 1)

            #(4) The derivation of LinearRelationship is given.
            Linear=LinearRelationship[i]

            # 4.1 If the current block is XOR and Kin_Right is 1, the LinearRelationship can be made 1.
            # Determine the linear relationship between the keys on Ec, for example: K24[0+i]=K24[4+i]=K24[12+i].

            s.addConstr(1- Direct - Linear >= 0)
            s.addConstr(1 - XOR - Linear >= 0)
            Reduce=ReducedData[i]
            s.addConstr(Reduce <=Linear)
            s.addConstr(Direct + XOR + 1 - Linear + Reduce >= 1)



        # The number of required LinearRelationship cannot be too large.

        s.addConstr( gb.quicksum(ReducedData)*self.WordSize+ Pr_Obj+P_in+P_out<= 16 * self.WordSize)


        #Set the time complexity as the objective function

        TempSum = 0
        for r in range(R_B):
            for i in range(4):
                TempSum += self.n * U_Pr[r][i][0]
                TempSum += self.WordSize * U_Pr[r][i][1]
        s.addConstr(P_in == TempSum)

        TempSum = 0
        for r in range(R_F):
            for i in range(4):
                TempSum += self.n * L_Pr[r][i][0]
                TempSum += self.WordSize * L_Pr[r][i][1]
        s.addConstr(P_out == TempSum)

        ## Calculate S_in and S_out
        AllSum=[ U_P_Exhaustion[j][i] for i in range(16) for j in range(R_B)]
        s.addConstr(S_in==self.WordSize *gb.quicksum(AllSum))

        AllSum = [L_P_Exhaustion[j][i] for i in range(16) for j in range(R_F)]
        s.addConstr(S_out == self.WordSize *gb.quicksum(AllSum))

        # Calculate Delta_in and Delta_out
        if TD_OR_Diff:
            s.addConstr(Delta_in==self.WordSize *gb.quicksum(TD_Values[0]))
            s.addConstr(Delta_out==self.WordSize *gb.quicksum(TD_Values[-1]))
        else:
            s.addConstr(Delta_in ==0)
            s.addConstr(Delta_out == 0)



        s.addConstr(DataOneTime==self.WordSize*(16-NumFixed))

        K_in=Num_K_in*self.WordSize
        K_out=Num_K_out*self.WordSize
        s.addConstr(Time1==Pr_Obj+P_out+S_in+K_in)
        s.addConstr(Time2==Pr_Obj+P_in+S_out+K_out+Delta_out-Delta_in+gb.quicksum(LinearRelationship)*self.WordSize )
        s.addConstr(Time3==Pr_Obj+DataOneTime+S_in+K_in+S_out+K_out+Delta_out-NumofIntersectionofKinAndKout*self.WordSize-
                    (NumLinear+NumFixed-gb.quicksum(LinearRelationship))*self.WordSize)

        s.addConstr(Time1<=TotalTime)
        s.addConstr(Time2<=TotalTime)
        s.addConstr(Time3<=TotalTime)


        s.setObjectiveN(TotalTime, index=0,priority=2)

        s.setObjectiveN(Time3, index=1,priority=1)
        s.setObjectiveN(S_in+S_out, index=2, priority=0.5)

        if NoImprovedTechnique:
            s.addConstr(P_in==0)
            s.addConstr(P_out == 0)
            s.addConstr(S_in == 0)
            s.addConstr(S_out == 0)



        if TD_OR_Diff:

            s.addConstr(
                16 * self.WordSize - 1 >= Pr_Obj + Delta_out
            )

        if not TD_OR_Diff:
            #s.addConstr((16 * self.WordSize-1)  >=(Pr_Obj+P_in+P_out) )
            s.addConstr((16 * self.WordSize - 1) >= (Pr_Obj + P_in + P_out))

        else:
            s.addConstr(16 * self.WordSize+Delta_in-1  >=(Pr_Obj+P_in+P_out))


        print("Start to solve the model...")
        try:
            s.setParam(gb.GRB.Param.PoolSolutions, ExpectedNumSol)
            s.setParam(gb.GRB.Param.PoolSearchMode, 2)

            if TimeLimit>0:
                s.setParam('TimeLimit', TimeLimit)
            s.update()

            StartTime = time.time()
            s.optimize()
            EndTime = time.time()
            print("Time used:", EndTime - StartTime)

            status = s.Status

            if status != gb.GRB.OPTIMAL:
                print('Optimization was stopped with status ' + str(status))

            nSolution = s.SolCount
            print('number of solution stored:' + str(nSolution))


            x = s.getVars()
            AllSolution = []
            # Print objective values of solutions
            for e in range(nSolution):
                tempDictionary = {}
                s.setParam(gb.GRB.Param.SolutionNumber, e)

                for i in range(len(x)):
                    tempDictionary[str(x[i].VarName)] = (x[i].Xn);
                AllSolution.append(tempDictionary)

                # Print solution values
                print(SeperateLine)
                print('Solution:\t', e)
                print("TotalTime:\t",TotalTime.x)
                print("Time1:\t\t",Time1.x)
                print("Time2:\t\t",Time2.x)
                print("Time3:\t\t",Time3.x)
                print("Pr_Obj:\t\t", Pr_Obj.x)
                print("P_in:\t\t",P_in.x)
                print("P_out:\t\t",P_out.x)
                print("S_in:\t\t",S_in.x)
                print("S_out:\t\t",S_out.x)
                print("Delta_in:\t",Delta_in.x)
                print("Delta_out:\t",Delta_out.x)
                print("Num_K_in:\t",Num_K_in.x)
                print("Num_K_out:\t",Num_K_out.x)
                print("NumofIntersectionofKinAndKout:\t",NumofIntersectionofKinAndKout.x)

                print("DataOneTime:\t",DataOneTime.x)
                print("NumFixed:\t",NumFixed.x)
                print("NumLinear:\t",NumLinear.x)

                print(SeperateLine)
                print("The number of times for TK nibbles appearing in TEMP_MK_Kin:")
                for i in range(16):
                    print(int(TEMP_MK_Kin[i].x), end=" ")
                print()
                print("The number of times for TK nibbles appearing in MK_Kin:")
                for i in range(16):
                    print(int(MK_Kin[i].x), end=" ")
                print()
                print("The TK nibbles that have been known in MK_Kin:")
                for i in range(16):
                    print(int(MK_Kin_Known[i].x), end=" ")
                print()

                print(SeperateLine)
                print("Extended_Kin:")
                self.DisplayMatrices([self.Transform(Extended_Kin, WordSize=8),
                                      ])

                print("RealKnownKey_Kin:")
                for i in range(2):
                    for j in range(len(RealKnownKey_Kin[i])):
                        print(int(RealKnownKey_Kin[i][j].x), end=" ")
                    print()







                print("The number of times for TK nibbles appearing in TEMP_MK_Kout:")
                for i in range(16):
                    print(int(TEMP_MK_Kout[i].x), end=" ")
                print()
                print("The number of times for TK nibbles appearing in MK_Kout:")
                for i in range(16):
                    print(int(MK_Kout[i].x), end=" ")
                print()

                print("The TK nibbles that have been known in MK_Kout_Known:")
                for i in range(16):
                    print(int(MK_Kout_Known[i].x), end=" ")
                print()

                print(SeperateLine)
                print("Extended_Kout:")
                self.DisplayMatrices([
                    self.Transform(Extended_Kout, WordSize=8)
                ])

                print("RealKnownKey_Kout:")
                for i in range(2):
                    for j in range(len(RealKnownKey_Kout[i])):
                        print(int(RealKnownKey_Kout[i][j].x), end=" ")
                    print()






                print("IntersecionofKinAndKout:")
                for i in range(16):
                    print(int(IntersecionofKinAndKout[i].x), end=" ")
                print()

                print("RealIntersecionofKinAndKout:")
                for i in range(16):
                    print(int(RealIntersecionofKinAndKout[i].x), end=" ")
                print()







                print("MK_Known:")
                for i in range(16):
                    print(int(MK_Known[i].x), end=" ")
                print()


                print("DirectLink:")
                for i in range(16):
                    print(int(DirectLink[i].x), end=" ")
                print()

                print("XOR_Link:")
                for i in range(16):
                    print(int(XOR_Link[i].x), end=" ")
                print()

                print("LinearRelationship:")
                for i in range(8):
                    print(int(LinearRelationship[i].x), end=" ")
                print()

                print("ReducedData:")
                for i in range(8):
                    print(int(ReducedData[i].x), end=" ")
                print()






                #Truncated differential trail
                print(SeperateLine)
                print("TD_Values:")
                for r in range(len(TD_Values)):

                    for i in range(16):
                        print(int(TD_Values[r][i].x), end="")
                    print()


                print("Pr_Values:")
                for r in range(len(Pr_Values)):
                    for i in range(4):
                        print(int(Pr_Values[r][i][0].x),int(Pr_Values[r][i][1].x), end="\t\t")
                    print()

                print(SeperateLine)
                print("U_Diff,U_P_DerivedFromP,U_P_Exhaustion,U_P_Current,U_K:")

                self.DisplayMatrices([self.Transform(U_Diff),
                                      self.Transform(U_P_DerivedFromP),
                                      self.Transform(U_P_Exhaustion),
                                      self.Transform(U_P_Current),
                                      self.Transform(U_K,WordSize=8)
                                      ])
                if not TD_OR_Diff:
                    print(SeperateLine)
                    print("U_K_LastR")
                    self.DisplayArray(self.Transform([U_K_LastR],WordSize=8)[0])


                print("U_Pr")
                for i in range(len(U_Pr)):
                    print("Round", i, ":")
                    for j in range(4):
                        print(int(U_Pr[i][j][0].x), int(U_Pr[i][j][1].x), end=" ")
                    print()

                print(SeperateLine)
                print("L_Diff,L_P_DerivedFromP,L_P_Exhaustion,L_P_Current,L_K:")

                self.DisplayMatrices([self.Transform(L_Diff),
                                      self.Transform(L_P_DerivedFromP),
                                      self.Transform(L_P_Exhaustion),
                                      self.Transform(L_P_Current),
                                      self.Transform(L_K,WordSize=8)
                                      ])

                if not TD_OR_Diff:
                    print(SeperateLine)
                    print("L_K_LastR")
                    self.DisplayArray(self.Transform([L_K_LastR], WordSize=8)[0])

                print("L_Pr")
                for i in range(len(L_Pr)):
                    print("Round", i, ":")
                    for j in range(4):
                        print(int(L_Pr[i][j][0].x), int(L_Pr[i][j][1].x), end=" ")
                    print()
                return (TotalTime.x)



        except gb.GurobiError as e:
            print('Error reported')

    def SearchAttack_Test_Different_Para(self,R_D,R_B,R_F,InOutDiffPr=[],ExpectedNumSol=1,NoImprovedTechnique=False,AdditionalKeyatRound=-1):
        Result=[]


        for r_b in range(2,R_B+R_F-1):
            r_f=R_B+R_F-r_b

            for AdditionalKeyatR in range(-1-r_b,r_f+2) :
                print(SeperateLine)
                print(R_D,r_b,r_f,AdditionalKeyatR)
                T=self.SearchAttack(R_D=R_D, R_B=r_b, R_F=r_f, InOutDiffPr=InOutDiffPr, ExpectedNumSol=ExpectedNumSol,
                                  NoImprovedTechnique=NoImprovedTechnique,AdditionalKeyatRound=AdditionalKeyatR)
                Result.append([R_D,r_b,r_f,AdditionalKeyatR,T])
                print(SeperateLine)
                print("{}-th result:".format(len(Result)))
                print(Result)

        print(Result)
        return Result







def NumToBinary(num, length):
    return [ (num>>(length-1-i))&1   for i in range(length)]


if __name__  == '__main__':

    GenerateBPT=False
    TDSearch_MILP=False
    AttackSearch_MILP=True

    if GenerateBPT:
        WordSize=4
        BPT_SKINNY = GetBPT(MC,WordSize=WordSize)
        print(SeperateLine)
        print(BPT_SKINNY)

        print(SeperateLine)
        AllPoints=[]
        for x in range(16):
            for y in range(16):
                if BPT_SKINNY[x][y]!=-1:
                    Bin_x=NumToBinary(x,4)
                    Bin_y=NumToBinary(y,4)
                    FLag=[0,0]
                    Lower=-WordSize-1
                    Upper=-WordSize+1

                    if ((BPT_SKINNY[x][y])<=Upper) & ((BPT_SKINNY[x][y])>=Lower):
                        FLag[1]=1

                    Lower = -2*WordSize - 1
                    Upper = -2*WordSize + 1

                    if ((BPT_SKINNY[x][y]) <= Upper) & ((BPT_SKINNY[x][y]) >= Lower):
                        FLag[0] = 1



                    AllPoints.append(Bin_x+Bin_y+FLag)
        print(SeperateLine)
        print(AllPoints)



    if TDSearch_MILP:
        Round=10
        WordSize=4
        #InOutDiff=[0x0001,0x4800]
        InOutDiff = []

        SKINNY_MILP=SKINNY(WordSize=WordSize)
        SKINNY_MILP.SearchTD(Round,ExpectedNumSol=1,InOutDiff=InOutDiff)

    if AttackSearch_MILP:


        R_D, R_B, R_F=15,4,4
        InOutDiffPr = [0x4fd9, 0x8098, 116.5]  

        # R_D, R_B, R_F = 13, 5, 5
        # InOutDiffPr = [0x1008, 0x2062, 105.9]  

        TimeLimit=0

        R_Additional=2
        Round = R_D + R_B + R_F+R_Additional


        WordSize = 8
        TweakeyMode=3
        print(SeperateLine)
        print("Round:", Round)
        print("WordSize:", WordSize)
        SKINNY_MILP = SKINNY(WordSize=WordSize,TweakeyMode=TweakeyMode)
        Time=SKINNY_MILP.SearchAttack( R_D, R_B, R_F,InOutDiffPr=InOutDiffPr, ExpectedNumSol=1,
                                  NoImprovedTechnique=False,AdditionalKeyatRound=-1,TimeLimit=TimeLimit)
        print(SeperateLine)
        print("Time Complexity:", Time)






