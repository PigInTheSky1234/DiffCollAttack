import gurobipy  as gb
import numpy as np
import os
import time
import math as mt
np.set_printoptions(threshold=np.inf)

SeperateLine="-"*100
NEG_INF=-1000


#CRAFT parameters
MC = [[1, 0, 1, 1], [0,1,0,1], [0, 0, 1, 0], [0, 0, 0, 1]]
Inv_MC = MC  # inverse of MC





BPT_CRAFT=[
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
# print(BPT_CRAFT)

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






def GetBPT(M,WordSize=4,UseLog2=True):
    if UseLog2:
        BPT=np.ones((16,16))*NEG_INF
    else:
        BPT = np.zeros((16, 16))
    BPT[0][0]=0
    for i in range(1,16):
            NumDiffs=0
            for dx in range(2**(WordSize*4)):
                if GetPattern(dx,WordSize)==i:
                        NumDiffs+=1

                        D_y=MultiplyMatrix(dx, M,WordSize)
                        OutPattern=GetPattern(D_y,WordSize)
                        if UseLog2:
                            if BPT[i][OutPattern]==NEG_INF:
                                BPT[i][OutPattern]=1
                            else:
                                BPT[i][OutPattern]+=1
                        else:
                            BPT[i][OutPattern] += 1


            if UseLog2:
                for j in range(1,16):
                        if BPT[i][j]>=0:
                            print(i, j, BPT[i][j], NumDiffs)
                            BPT[i][j]=mt.log2(BPT[i][j]/(NumDiffs))

    return BPT



class CRAFT:
    def __init__(self, WordSize=4,TweakeyMode=2 ):
        self.WordSize = WordSize
        self.n = WordSize * 2
        self.TweakeyMode = TweakeyMode

        self.MC =[[1, 0, 1, 1], [0,1,0,1], [0, 0, 1, 0], [0, 0, 0, 1]]
        self.Inv_MC = [[1, 0, 1, 1], [0,1,0,1], [0, 0, 1, 0], [0, 0, 0, 1]] #inverse of MC

    def column(self,A, j):
	    return [A[j], A[j+4], A[j+8], A[j+12]]
    def IsKnown(self,Model,A,B,C):
        Model.addConstr(A[2]==B[2])
        Model.addConstr(A[3]==B[3])


        Model.addConstr(A[0]>=C[0])
        Model.addConstr(A[1] >= C[1])

        Model.addConstr((1-B[0])+(1-B[2])+(1-B[3])+A[0] >= 1)
        Model.addConstr((1 - B[1])  + (1 - B[3]) + A[1] >= 1)

        Model.addConstr(C[0]+ B[0]+(1-A[0])>=1)
        Model.addConstr(C[0] + B[2] + (1 - A[0]) >= 1)
        Model.addConstr(C[0] + B[3] + (1 - A[0]) >= 1)
        Model.addConstr(C[1] + B[1] + (1 - A[1]) >= 1)
        Model.addConstr(C[1] + B[3] + (1 - A[1]) >= 1)


    def ShiftRow(self,A):

        return [A[15], A[12], A[13], A[14],
                A[10], A[9], A[8], A[11],
                A[6],A[5],A[4], A[7],
                A[1],A[2],A[3],A[0]]

    def TD_RoundFunction_UP(self, Model,A, B, Pr):
        #Model.addConstr(gb.quicksum(A)<=15)

        ShiftRow_B=self.ShiftRow(B)
        for i in range(4):

            TEMP= [ShiftRow_B[ (i + 4 * j)%16] for j in range(4)]
            TEMP+= [A[(i + 4 * j) % 16] for j in range(4)]

            TEMP +=[Pr[i][0],Pr[i][1]]

            S_T = [(1, 2, 1, 3, -1, -2, 0, 0, -5, -2, 0),
                 (-1, 0, -1, -3, 3, 2, 0, 0, 5, 2, 0),
                 (0, 1, 0, 1, -2, -1, 0, 0, -2, -2, 2),
                 (-1, -2, -1, 0, 2, 2, 0, 0, 3, 2, 0),
                (2, 1, 2, 3, -2, -1, 0, 0, -6, -2, 0),
                 (0, 0, 0, -1, 0, 1, 0, 0, 1, 1, 0),
                (0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 1),
                 (-1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0),
                (1, 0, 1, 1, -1, -2, 0, 0, -3, -2, 2)]

             # 8p0 , 4p1
            S_T_Back = S_T.copy()

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
            In2=A[(i + 4 * 2) % 16]
            In3=A[(i + 4 * 3) % 16]
            Out2=ShiftRow_B[(i + 4 * 2) % 16]
            Out3=ShiftRow_B[(i + 4 * 3) % 16]


            Model.addConstr(In2==Out2)
            Model.addConstr(In3==Out3)

    def MultiplyIEQ(self, Model,X,IEQ):
        for ieq in IEQ:
            TempSum=0
            for j in range(len(ieq)):
                if j==(len(ieq)-1):
                    if ieq[j] != 0:
                        TempSum+=ieq[j]
                        continue

                if ieq[j]!=0:
                    TempSum+=ieq[j]*X[j]
            Model.addConstr(TempSum>=0)



    def TD_RoundFunction(self, Model,A, B, Pr):
        #Model.addConstr(gb.quicksum(A)<=15)

        ShiftRow_B=self.ShiftRow(B)
        for i in range(4):
            TEMP=[ A[(i+4*j)%16 ] for j in range(4) ]
            TEMP+= [ShiftRow_B[ (i + 4 * j)%16] for j in range(4)]
            TEMP +=[Pr[i][0],Pr[i][1]]

            S_T = [(1, 2, 1, 3, -1, -2, 0, 0, -5, -2, 0),
                 (-1, 0, -1, -3, 3, 2, 0, 0, 5, 2, 0),
                 (0, 1, 0, 1, -2, -1, 0, 0, -2, -2, 2),
                 (-1, -2, -1, 0, 2, 2, 0, 0, 3, 2, 0),
                (2, 1, 2, 3, -2, -1, 0, 0, -6, -2, 0),
                 (0, 0, 0, -1, 0, 1, 0, 0, 1, 1, 0),
                (0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 1),
                 (-1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0),
                (1, 0, 1, 1, -1, -2, 0, 0, -3, -2, 2)]

             # 8p0 , 4p1
            S_T_Back = S_T.copy()

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
            In2=A[(i + 4 * 2) % 16]
            In3=A[(i + 4 * 3) % 16]
            Out2=ShiftRow_B[(i + 4 * 2) % 16]
            Out3=ShiftRow_B[(i + 4 * 3) % 16]


            Model.addConstr(In2==Out2)
            Model.addConstr(In3==Out3)





    def SearchTD(self,Round,ExpectedNumSol=1):

        s = gb.Model()


        TD_Values=[
           [
               s.addVar(vtype=gb.GRB.BINARY, name="TD_{}_R{}".format(i,k))
             for i in range(16) ]
            for k in range(Round+1)
        ]

        Pr_Values=[
            [s.addVars(2,vtype=gb.GRB.BINARY, name="Pr_{}_R{}".format(i,k)) for i in range(4)] for k in range(Round)
        ]

        Obj = s.addVar(lb=0, ub=self.TweakeyMode*16*self.WordSize+100, vtype=gb.GRB.INTEGER, name="obj")
        s.update()

        #2.定义轮函数之间的关系

        for r in range(Round):
            InTD=TD_Values[r]
            OutTD=TD_Values[r+1]
            Pr=Pr_Values[r]
            self.TD_RoundFunction(s,InTD,OutTD,Pr)

        #输入不为0
        s.addConstr(gb.quicksum(TD_Values[0])>=1)
        # 下面是特地说明截断差分具有统计意义的约束
        s.addConstr(gb.quicksum(TD_Values[-1]) * self.WordSize + Obj <= (16 * self.WordSize - 1))


        #3.进行求解，并展示结果
        TempSum=0
        for r in range(Round):
            for i in range(4):
                TempSum+=self.n*Pr_Values[r][i][0]
        for r in range(Round):
            for i in range(4):
                TempSum+=self.WordSize*Pr_Values[r][i][1]




        s.addConstr(Obj==TempSum)
        s.setObjective(Obj, sense=gb.GRB.MINIMIZE)



        #下面是求解过程
        print("Start to solve the model...")
        try:
            s.setParam(gb.GRB.Param.PoolSolutions, ExpectedNumSol)
            s.setParam(gb.GRB.Param.PoolSearchMode, 2)
            # 0：不保证其他解的质量，且是顺便求出来的，并且不能够保证解的数量已经达到了期望的PoolSolutions
            # 1：不保证其他解的质量，
            # 2:能够保证这是最有的k 个解

            # #设置最长的运行时间
            # s.setParam('TimeLimit', 5 * 60)
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

    def GetDiffRoundbyRound_Up(self, Model, A, B, MC):
        SHIFT_ROW_B = self.ShiftRow(B)

        for col in range(4):
            Out = self.column(SHIFT_ROW_B, col)
            In = self.column(A, col)

            # 输出差分块活跃，则涉及道德输入差分块应当全部活跃
            for i in range(4):
                for j in range(4):
                    if MC[i][j] == 1:
                        Model.addConstr(Out[i] <= In[j])

                # 输入字块活跃，必须要有至少有一个影响的输出字块活跃
                TempSum = 0
                for j in range(4):
                    if MC[j][i] == 1:
                        TempSum += Out[j]
                Model.addConstr(TempSum >= In[i])


    def GetDiffRoundbyRound_Low(self,Model,A,B,MC):
        SHIFT_ROW_B=self.ShiftRow(B)

        for col in range(4):
            Out=self.column(A,col)
            In=self.column(SHIFT_ROW_B,col)


            #输出差分块活跃，则涉及道德输入差分块应当全部活跃
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

    def GetRelationofK(self,Round):
        PT = [9, 15, 8, 13, 10, 14, 12, 11, 0, 1, 2, 3, 4, 5, 6, 7]
        AllIndices = [[i for i in range(16)]]

        for r in range(Round - 1):
            TempIndice = AllIndices[-1]
            NewIndice = [0] * 16
            for i in range(16):
                NewIndice[i] = TempIndice[PT[i]]
            AllIndices.append(NewIndice)
        return AllIndices

    def GetLeftKofEc(self,Key):
        return [
            Key[0],Key[1],Key[2],Key[3],
            Key[0], Key[1], Key[2], Key[3],
            Key[7], Key[4], Key[5], Key[6],
            Key[0], Key[1], Key[2], Key[3],


        ]

    def DisplayArray(self,Array):
        for i in range(4):
            for j in range(4):
                print(Array[i*4+j], end="")
            print()
    def DisplayMatrices(self,Matrices):
        LINE="-"*30

        for r  in range(len(Matrices[0])):
            print(LINE)
            print("Round",r)
            for i in range(len(Matrices)):
                if r<=(len(Matrices[i])-1):
                    print(LINE)
                    self.DisplayArray(Matrices[i][r])
    def Transform(self,X):
        return [ [int(X[r][i].x) for i in range(16) ]  for r in range(len(X))]


    def SearchAttack_AdditionalOneRound(self,R_D,R_B,R_F,ExpectedNumSol=1):

        r0=R_B
        r1=R_B+R_D
        r2=R_B+R_D+R_F
        r3=R_B+R_D+R_F+2




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
            [s.addVars(2, vtype=gb.GRB.BINARY, name="Pr_{}_R{}".format(i, k)) for i in range(4)] for k in range(R_B)
        ]


        U_K = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="U_K_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(R_B-1)
        ]






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
                for i in range(16)]
            for k in range(R_F - 1)
        ]
        L_K_LastRound = [
                s.addVar(vtype=gb.GRB.BINARY, name="L_K_LastRound_{}".format(i))
                for i in range(8)

        ]



        Extended_Kin=[
            [
                s.addVar(vtype=gb.GRB.BINARY, name="Extended_Kin_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(1)
        ]

        Extended_Kout = [
            [
                s.addVar(vtype=gb.GRB.BINARY, name="Extended_Kout_{}_R{}".format(i, k))
                for i in range(16)]
            for k in range(1)
        ]


        DirectLink=[
            s.addVar(vtype=gb.GRB.BINARY, name="DirectLink_{}".format(i))
            for i in range(16)
        ]

        XOR_Link = [
            s.addVar(vtype=gb.GRB.BINARY, name="XOR_Link_{}".format(i))
            for i in range(16)
        ]

        NumFixed=s.addVar(lb=0, ub=16, vtype=gb.GRB.INTEGER, name="NumFixed")
        NumLinear=s.addVar(lb=0, ub=16, vtype=gb.GRB.INTEGER, name="NumLinear")




        MK_Kin=[
            s.addVar(vtype=gb.GRB.BINARY, name="MK_Kin_"+str(i))
            for i in range(32)
        ]
        MK_Kout = [
            s.addVar( vtype=gb.GRB.BINARY, name="MK_Kout_" + str(i))
            for i in range(32)
        ]




        Num_K_in = s.addVar(lb=0, ub=32, vtype=gb.GRB.INTEGER, name="Num_K_in")
        Num_K_out = s.addVar(lb=0, ub=32, vtype=gb.GRB.INTEGER, name="Num_K_out")

        IntersectionofKinAndKout = [
            s.addVar(vtype=gb.GRB.BINARY, name="IntersectionofKinAndKout_" + str(i))
            for i in range(32)
        ]


        NumofIntersectionofKinAndKout= s.addVar(lb=0, ub=32, vtype=gb.GRB.INTEGER, name="NumofIntersectionofKinAndKout")



        P_in = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER,name="P_in")
        P_out = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="P_out")
        S_in = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="S_in")
        S_out = s.addVar(lb=0, ub=100, vtype=gb.GRB.INTEGER, name="S_out")
        Delta_in = s.addVar(lb=0, ub=16*self.WordSize, vtype=gb.GRB.INTEGER, name="Delta_in")
        Delta_out = s.addVar(lb=0, ub=16*self.WordSize, vtype=gb.GRB.INTEGER, name="Delta_out")
        DataOneTime=s.addVar(lb=0, ub=16*self.WordSize, vtype=gb.GRB.INTEGER, name="DataOneTime")


        UPBOUND=self.TweakeyMode * 16 * self.WordSize + 100

        TotalTime=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="TotalTime")
        Time1=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Time1")
        Time2=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Time2")
        Time3=s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Time3")
        Time4 = s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Time4")

        Pr_Obj = s.addVar(lb=0, ub=UPBOUND, vtype=gb.GRB.INTEGER, name="Pr_Obj")

        s.update()


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
            s.addConstr(U_Diff[-1][i]==TD_Values[0][i])
            s.addConstr(L_Diff[0][i] == TD_Values[-1][i])


        for r in range(R_B):
            InTD=U_Diff[r]
            OutTD=U_Diff[r+1]
            Pr=U_Pr[r]
            self.TD_RoundFunction_UP(s,InTD,OutTD,Pr)


        for r in range(R_B-1):
            A=U_P_DerivedFromP[r]
            B=U_P_Current[r+1]
            self.GetDiffRoundbyRound_Up( s, A, B, self.MC)
        for i in range(16):
            s.addConstr(U_P_DerivedFromP[-1][i]==0)


        for r in range(R_B):
            A=U_P_Exhaustion[r]
            B=U_P_DerivedFromP[r]
            C=U_Diff[r+1]
            for i in range(16):
                s.addConstr(A[i]<=B[i]+C[i])

        for r in range(R_B-1):
            A=U_P_Current[r]
            B=U_P_Exhaustion[r]
            C=U_Diff[r+1]
            D=U_P_DerivedFromP[r]
            for i in range(16):
                s.addConstr(A[i]<=1-B[i])



                s.addConstr(B[i]+ (1-C[i])+ A[i]>=1)
                s.addConstr(B[i] + (1 - D[i]) + A[i] >= 1)
                s.addConstr( C[i] + D[i] +(1- A[i]) >= 1)
        for i in range(16):
            s.addConstr(U_P_Current[-1][i]==0)


        for r in range(R_F):
            InTD=L_Diff[r]
            OutTD=L_Diff[r+1]
            Pr=L_Pr[r]
            self.TD_RoundFunction(s,InTD,OutTD,Pr)

        for r in range(R_F-1 ):
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

        for r in range(1,R_F):
                A = L_P_Current[r]
                B = L_P_Exhaustion[r]
                C = L_Diff[r]
                D = L_P_DerivedFromP[r]
                for i in range(16):
                    s.addConstr(A[i] <= 1 - B[i])


                    s.addConstr(B[i] + (1 - C[i]) + A[i] >= 1)
                    s.addConstr(B[i] + (1 - D[i]) + A[i] >= 1)
                    s.addConstr(C[i] + D[i] + (1 - A[i]) >= 1)
        for i in range(16):
            s.addConstr(L_P_Current[0][i] == 0)


        for r in range(R_B-1):
            for i in range(16):


                s.addConstr(U_K[r][i]==self.ShiftRow(U_P_Current[r ])[i])

        for r in range(R_F-2):

            A=L_P_Current[r+1]
            B=L_K[r]

            for col in range(4):
                    Out = self.column(A, col)
                    In = self.column(B, col)


                    for i in range(4):
                        for j in range(4):
                            if MC[i][j] == 1:
                                s.addConstr(Out[i] <= In[j])


                        TempSum = 0
                        for j in range(4):
                            if MC[j][i] == 1:
                                TempSum += Out[j]
                        s.addConstr(TempSum >= In[i])


        A = L_P_Current[-1]
        B = L_K[-1]
        for i in range(8):
            s.addConstr(B[i+8] >= A[i+8])
        for i in range(4):
            s.addConstr(A[8+i]+ A[i]+(1-B[8+i])>=1)
            s.addConstr(A[12 + i] + A[i] +A[4+i]+ (1 - B[12 + i]) >= 1)


        for col in range(4):
            Out = self.column(B, col)
            In = self.column(A, col)
            F1=L_K_LastRound[col]
            F2 = L_K_LastRound[col+4]

            IEQ1=[[0, -1, 0, 1, 0, 0], [0, -1, 1, 0, 0, 0], [0, 0, -1, -1, -1, 2]]
            TempVar1=[In[0],Out[0],Out[2],Out[3],F1]
            self.MultiplyIEQ(s,TempVar1,IEQ1)
            s.addConstr(TempVar1[0]-TempVar1[1]-TempVar1[4]==0)


            IEQ2 = [[0, -1, 1, 0, 0], [0, -1, 0, -1, 1]]
            TempVar2 = [In[1],  Out[1], Out[3], F2]
            self.MultiplyIEQ(s, TempVar2, IEQ2)
            s.addConstr(TempVar2[0] - TempVar2[1] - TempVar2[3] == 0)









        for i in range(16):
            TempSum0=0
            TempSum1=0
            for r in range(R_B-1):
                if r%2==0:
                    TempSum0+=U_K[r][i]
                else:
                    TempSum1+=U_K[r][i]

                s.addConstr(U_K[r][i]<=MK_Kin[i+(r%2)*16])
            for j in range(1):
                if ((R_D+R_B+R_F+j) % 2)==0:
                    TempSum0+=Extended_Kin[j][i]
                else:
                    TempSum1+=Extended_Kin[j][i]

                s.addConstr(Extended_Kin[j][i] <= MK_Kin[i + ((R_D+R_B+R_F+j) % 2) * 16])
            s.addConstr(MK_Kin[i]<=TempSum0)
            s.addConstr(MK_Kin[i+16] <= TempSum1)

        for i in range(16):
            TempSum0 = 0
            TempSum1 = 0
            for r in range(R_F-1):
                if ((r+R_D+R_B+1)%2)==0:
                    TempSum0 += L_K[r][i]
                    s.addConstr(L_K[r][i]<=MK_Kout[i])
                else:
                    TempSum1 += L_K[r][i]
                    s.addConstr(L_K[r][i]<=MK_Kout[i+16])
            for j in range(1):
                if ((R_D+R_B+R_F+j) % 2)==0:
                    TempSum0 += Extended_Kout[j][i]
                    s.addConstr(Extended_Kout[j][i]<=MK_Kout[i])
                else:
                    TempSum1 += Extended_Kout[j][i]
                    s.addConstr(Extended_Kout[j][i]<=MK_Kout[i+16])
            s.addConstr(MK_Kout[i]<=TempSum0)
            s.addConstr(MK_Kout[i+16]<=TempSum1)




        for i in range(32):
            s.addConstr((1-MK_Kin[i])+(1- MK_Kout[i])+IntersectionofKinAndKout[i]>=1)
            s.addConstr(MK_Kin[i] + 1-IntersectionofKinAndKout[i] >= 1)
            s.addConstr(MK_Kout[i] + 1 - IntersectionofKinAndKout[i] >= 1)

        s.addConstr(NumofIntersectionofKinAndKout==gb.quicksum(IntersectionofKinAndKout))
        s.addConstr(Num_K_in==gb.quicksum(MK_Kin))
        s.addConstr(Num_K_out==gb.quicksum(MK_Kout)+gb.quicksum(L_K_LastRound))



        s.addConstr(NumFixed==gb.quicksum(DirectLink))
        s.addConstr(NumLinear==gb.quicksum(XOR_Link))


        Temp=((r2)%2)*16
        Kout_R0=(MK_Kout[Temp:Temp+16  ])
        Kin_R0 = (MK_Kin[Temp:Temp+16])



        for i in range(16):
            A0=Kout_R0[i]
            B0=Kin_R0[i]
            C=DirectLink[i]
            D=XOR_Link[i]

            s.addConstr(C+D==1)
            s.addConstr((1-A0)+C>=1)
            s.addConstr((1 - B0) + C >= 1)
            s.addConstr(A0+B0 + 1-C >= 1)




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


        AllSum=[ U_P_Exhaustion[j][i] for i in range(16) for j in range(R_B)]
        s.addConstr(S_in==self.WordSize *gb.quicksum(AllSum))

        AllSum = [L_P_Exhaustion[j][i] for i in range(16) for j in range(R_F)]
        s.addConstr(S_out == self.WordSize *gb.quicksum(AllSum))


        s.addConstr(Delta_in==self.WordSize *gb.quicksum(TD_Values[0]))
        s.addConstr(Delta_out==self.WordSize *gb.quicksum(TD_Values[-1]))


        s.addConstr(DataOneTime==self.WordSize*(16-NumFixed))

        K_in=Num_K_in*self.WordSize
        K_out=Num_K_out*self.WordSize
        s.addConstr(Time1==Pr_Obj+P_out+S_in+K_in)
        s.addConstr(Time2==Pr_Obj+P_in+S_out+K_out+Delta_out-Delta_in)
        s.addConstr(Time3==Pr_Obj+DataOneTime+S_in+K_in+S_out+K_out+Delta_out-NumofIntersectionofKinAndKout*self.WordSize-
                    self.WordSize*(NumFixed+NumLinear)
                    )
        s.addConstr(Time4==128+Pr_Obj + DataOneTime + Delta_out -
            self.WordSize * (NumFixed + NumLinear)-NumLinear*self.WordSize)

        s.addConstr(
           16*self.WordSize -1 >= Pr_Obj + Delta_out
        )


        s.addConstr(Time1<=TotalTime)
        s.addConstr(Time2<=TotalTime)
        s.addConstr(Time3<=TotalTime)
        s.addConstr(Time4 <= TotalTime)
        s.setObjective(TotalTime, sense=gb.GRB.MINIMIZE)


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
                print("TotalTime:\t",TotalTime.x)
                print("Time1:\t\t",Time1.x)
                print("Time2:\t\t",Time2.x)
                print("Time3:\t\t",Time3.x)
                print("Time4:\t\t", Time4.x)
                print("Pr_Obj:\t\t",Pr_Obj.x)
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
                                      self.Transform(U_K)
                                      ])
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
                                      self.Transform(L_K)
                                      ])
                print("L_Pr")
                for i in range(len(L_Pr)):
                    print("Round", i, ":")
                    for j in range(4):
                        print(int(L_Pr[i][j][0].x), int(L_Pr[i][j][1].x), end=" ")
                    print()

                print(SeperateLine)
                print("L_K_LastRound")
                for i in range(8):
                    print(int(L_K_LastRound[i].x), end=" ")


                print(SeperateLine)
                print("Extended_Kin,Extended_Kout:")
                self.DisplayMatrices([self.Transform(Extended_Kin),
                                      self.Transform(Extended_Kout)
                                      ])
                print(SeperateLine)
                print("MK_Kin,MK_Kout:")
                self.DisplayMatrices([self.Transform([MK_Kin[:16],MK_Kin[16:]]),
                                      self.Transform([MK_Kout[:16], MK_Kout[16:]]),

                                      ])
                print(SeperateLine)


        except gb.GurobiError as e:
            print('Error reported')




def NumToBinary(num, length):
    return [ (num>>(length-1-i))&1   for i in range(length)]








if __name__  == '__main__':

    GenerateBPT=False

    TDSearch_MILP=False

    AttackSearch_AdditionalOneRound_MILP = True

    if GenerateBPT:
        WordSize=4
        BPT_CRAFT = GetBPT(MC,WordSize=WordSize)
        print(SeperateLine)
        print(BPT_CRAFT)

        print(SeperateLine)
        AllPoints=[]
        for x in range(16):
            for y in range(16):
                if BPT_CRAFT[x][y]!=-1:
                    Bin_x=NumToBinary(x,4)
                    Bin_y=NumToBinary(y,4)
                    FLag=[0,0]
                    Lower=-WordSize-1
                    Upper=-WordSize+1

                    if ((BPT_CRAFT[x][y])<=Upper) & ((BPT_CRAFT[x][y])>=Lower):
                        FLag[1]=1

                    Lower = -2*WordSize - 1
                    Upper = -2*WordSize + 1

                    if ((BPT_CRAFT[x][y]) <= Upper) & ((BPT_CRAFT[x][y]) >= Lower):
                        FLag[0] = 1



                    AllPoints.append(Bin_x+Bin_y+FLag)
        print(SeperateLine)
        print(AllPoints)

        Points_Prob1=[]
        for y in range(16):
            Bin_y = NumToBinary(y, 4)
            TempX=np.zeros(4,dtype=np.uint)
            for i in range(4):
                if Bin_y[i]==1:
                    TempX|=np.array(MC[i],dtype=np.uint)
            Points_Prob1.append(TempX.tolist()+Bin_y+[0,0])
        print(Points_Prob1)

        print(SeperateLine)
        CombinedPoints=[]









    if TDSearch_MILP:
        Round=11

        CRAFT_MILP=CRAFT()
        CRAFT_MILP.SearchTD(Round,ExpectedNumSol=1)


    if AttackSearch_AdditionalOneRound_MILP:
        R_D, R_B, R_F=11,6,5
        Round = R_D + R_B + R_F+1

        print(SeperateLine)
        print("Round:", Round)
        CRAFT_MILP = CRAFT()
        CRAFT_MILP.SearchAttack_AdditionalOneRound( R_D, R_B, R_F, ExpectedNumSol=1)











