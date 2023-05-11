__author__="Chengshun Shang"
__copyright__= "Copyright (C) 2022-present by Geo_StarLab"
__version__ = "1.0"
__maintainer__ = "Chengshun Shang"
__email__ = "chengshun.shang1996@gmail.com"
__status__ = "development"
__date__ = "Nov 25, 2022"

import numpy as np
class FemProcedure:

    def __init__(self) -> None:
        """
        Initialize the class "FemProcedure"
        """
        pass

    def TrussStaticFE(self, x, y, z, ele, Load, Constr, Scale):
        """
        The main FE calculation procedure

        Parameters
        ----------
        num : int
            the number of gauss point one want to generate
        """

        Dofs = 2 * x.shape[1] #�����ɶ���
        self.EleCount = ele.shape[0] #��Ԫ����
        self.K = np.zeros(Dofs,Dofs) #����նȾ���
        self.F = np.zeros(Dofs,1) #�����غ�����
        self.U = np.zeros(Dofs,1) #����λ������
        self.BarLength = self.BarsLength(x,y,ele)
        
        #figure('Name','Undeformed Truss')
        #RenderTruss(ele,Load,Constr,x,y,U,1,'-k',1) %�������
        #figure('Name','Deformed Truss')
        #RenderTruss(ele,Load,Constr,x,y,U,0.5,'-k',1) %���Ʊ��κ����

        #�������е�Ԫ��������Ԫ�ն���ֿ���װ������ն���
        for iEle in range(1, self.EleCount):
            #�õ�Ԫ�������ڵ�ı��
            n1=ele(iEle,2)
            n2=ele(iEle,3)
            #��������任����
            temp_BarLength = self.BarLength(iEle)
            R = self.CoordTransform([x(n1), x(n2)],[y(n1) y(n2)],temp_BarLength)
            #���㵥Ԫ�նȾ���
            ke= self.BarElementKe(ele(iEle,4),ele(iEle,5),R,temp_BarLength)
            #������Ԫ�նȷֿ���װ���ܸ���Ӧλ��
            self.K[2*n1-1:2*n1, 2*n1-1:2*n1] = self.K[2*n1-1:2*n1, 2*n1-1:2*n1] + ke[1:2, 1:2]
            self.K[2*n1-1:2*n1, 2*n2-1:2*n2] = self.K[2*n1-1:2*n1, 2*n2-1:2*n2] + ke[1:2, 3:4]
            self.K[2*n2-1:2*n2, 2*n1-1:2*n1] = self.K[2*n2-1:2*n2, 2*n1-1:2*n1] + ke[3:4, 1:2]
            self.K[2*n2-1:2*n2, 2*n2-1:2*n2] = self.K[2*n2-1:2*n2, 2*n2-1:2*n2] + ke[3:4, 3:4]

        #�γ��غ�����--1��2��3��4���ɶȸ��غ�ֵ
        for LoadNum in range(1, Load.shape[0] + 1):
            for i in range(2, 3+1):
                self.F[2*Load(LoadNum)+i-3,1] = Load[LoadNum, i]

        #ʩ��Լ��--�˴�����
        for iConstr in range(1, Constr.shape[1]+1):
            for j in range(2, 3+1):
                if Constr[iConstr,j] != None:
                    self.K[2*Constr[iConstr,1]+j-3, 2*Constr[iConstr,1]+j-3] = 1e12 * self.K[2*Constr[iConstr,1 ]+j-3, 2*Constr[iConstr,1 ]+j-3]
                    self.F[2 * Constr[iConstr,1]+j-3] = Constr[iConstr,j]* self.K[2*Constr[iConstr,1 ]+j-3, 2*Constr[iConstr,1]+j-3]

        U=pinv(K)*F; #ȫ������ϵ��λ��,�˴�ԭ��Ϊ�浫�Ǿ������죬�ʸĳ�α��

        for iEle in range(1, self.EleCount+1):
            #����˾ֲ������µ�λ��
            n1=ele(iEle,2);n2=ele(iEle,3);
            R = self.CoordTransform([x(n1) x(n2)],[y(n1) y(n2)],BarLength(iEle));
            localU = R*[U(2*n1-1:2*n1,1);U(2*n2-1:2*n2,1)];
            Strain(1, iEle)=[-1/BarLength(iEle) 1/BarLength(iEle)]*localU  #Ӧ��
            Stress(1, iEle)=ele(iEle,5)* Strain(1, iEle)  #Ӧ��
            AxialForce(1, iEle)=ele(iEle,4)* Stress(1, iEle)  #����

        print(Strain)
        '''
        #����λ�ƣ�Ӧ�����������ı��ļ�
        fp=fopen('Result.txt','a');
        str = [char(13,10)','U',' ',num2str(U'),char(13,10)','Stress',' ',...
        num2str(Stress),char(13,10)','AxialForce',' ',num2str(AxialForce)];
        fprintf(fp,str);
        fclose(fp);
        RenderTruss(ele,Load,Constr,x,y,U,1,'-.b',Scale) #���Ʊ��κ����
        '''

    def BarElementKe(self, A, E, R, Barlength):
        pass

    def BarsLength(self, x, y, z, ele):
        pass

    def CoordTransform(self, a, b, c):
        pass

if __name__ == "__main__":

    d1=30
    d2=2*d1
    d3=3*d1
    d4=4*d1
    d5=5*d1
    d6=6*d1
    d7=7*d1
    d8=8*d1
    h1=40
    h2=2*h1
    h3=3*h1
    h4=4*h1
    x=[d4 d3 d5 d2 d4 d6 d1 d3 d5 d7 d4 0 d8] # �ڵ�x�᷽������
    y=[h4 h3 h3 h2 h2 h2 h1 h1 h1 h1 0 0 0] # �ڵ�y�᷽������
    A=3
    E=2.1E005 # �����������͵���ģ��
    #��Ԫ��Ϣ����ţ��ڵ�1��ţ��ڵ�2��ţ��������������ģ��
    ele=[1 1 2 A E;2 1 3 A E;3 2 3 A E;4 2 4 A E;5 2 5 A E;6 3 5 A E;7 3 6 A E;...
        8 4 5 A E;9 5 6 A E;10 4 7 A E;11 4 8 A E;12 5 8 A E;13 5 9 A E;14 6 9 A E;...
        15 6 10 A E;16 7 8 A E;17 8 9 A E;18 9 10 A E;19 7 12 A E;20 8 12 A E;...
        21 7 11 A E;22 8 11 A E;23 9 11 A E;24 10 11 A E;25 9 13 A E;26 10 13 A E]
    #�غ���Ϣ���ڵ��ţ�x��������y��������z������
    Load=[1 200 100]
    #Լ�����ڵ��ţ�x����Լ����y����Լ����z����Լ��
    Constr=[11 0 0;12 0 0;13 0 0]
