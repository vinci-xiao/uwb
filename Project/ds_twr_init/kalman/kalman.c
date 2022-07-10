#include  "kalman.h"
#include  "stdlib.h"
#include  "math.h"


/*==============================================================================
1.���
   X(k|k-1) = A(k,k-1)*X(k-1|k-1)                 //�S���v�R0


2.���솕��tԴ������
   P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q(k)
   Q(k) = U(k)�GU(k)' 


3.����U���C���y����
   K(k) = P(k|k-1)*C(k)' / (C(k)*P(k|k-1)*C(k)' + R(k))
   R(k) = N(k)�GN(k)' 


4.�d���
   X(k|k) = X(k|k-1)+K(k)*(Y(k)-C(k)*X(k|k-1))


5.�����d��Y��t�紬����
   P(k|k) =��I-K(k)*C(k)��*P(k|k-1)


6. �d���I�b�t


A(k,k-1):     О�F�v�����
B(k,k-1):     О�Fǻ�S���v
C(k):         ��������v�����
X(k|k-1):     �H��k-1�C�Dǻ�I�b�t�k�C�Dǻ�t
X(k-1|k-1):   k-1�C�Dǻ�I�b�t
P(k|k-1):     X(k|k-1)����ǻcovariance
P(k-1|k-1):   X(k-1|k-1)����ǻcovariance
Q(k):         ���f����ǻcovariance(�ڊf�������v��tǻ�|?��H)
R(k):         ���v����ǻ�tԴ��
Y(k):         k�C�Dǻ���v�t
K(k):         �U���C���y
U(k):         k�C�D�ېF���
N(k):         k�C�D������


�?:     ���fО�FX-----------桕��y�H
          ���f����A-----------�y�H���w�v��
          О�Fǻ�S���vB-------(�o������)
          ����tY-------------�y�H����f
          �������C-----------�G�ރH-��ò�ރH
          ������Q-----------�y�H���w?��
          ���v���R-----------��f�d��

�r��保: �����H�ފf����(?�����C����u�H�ކ����t����)ǻ�f�ޔ���̒��ǻ��t,
          ���Z,�H�ފf����ǻ�f�ޔ���̒����tǻ�tԴ��;  Ն҃,�ϒ����tǻ�t
          Դ������̿U���C���y;  �I�Y,�H�������t����v�t����g?�I�b�t��?�tԴ��
          吂X,�U���C�]���ό�C�K�fǻ�ҠZ,�ݱp?�R?޷�I������
==============================================================================*/



//================================================//
//==             �I�b�tԴ���x�ʼ�               ==//
//================================================//
typedef struct  _tCovariance
{
  float PNowOpt[P_LENGTH];
  float PPreOpt[P_LENGTH];
}tCovariance;



//================================================//
//==               �I�b�t�x�ʼ�                 ==//
//================================================//
typedef struct  _tOptimal
{
  float XNowOpt[X_LENGTH];
  float XPreOpt[X_LENGTH];
}tOptimal;



tOptimal      tOpt;                                     //  ���]��B�f���v�k���a�w�^�t
tOptimal      tOpt_1;                                     //  ���]��B�f���v�k���a�w�^�t
tOptimal      tOpt_2;                                     //  ���]��B�f���v�k���a�w�^�t

tCovariance   tCov;                                     //  ���]��B�f���v�k���a�w�^�t
tCovariance   tCov_1;                                     //  ���]��B�f���v�k���a�w�^�t
tCovariance   tCov_2;                                     //  ���]��B�f���v�k���a�w�^�t

float         Y[Y_LENGTH]  = Y_VALUE;                   //  ���v�t(�������vǻ�f�ބz����?͘�f�V)
float         Y_1[Y_LENGTH]  = Y_VALUE;                   //  ���v�t(�������vǻ�f�ބz����?͘�f�V)
float         Y_2[Y_LENGTH]  = Y_VALUE;                   //  ���v�t(�������vǻ�f�ބz����?͘�f�V)
float         I[I_LENGTH]  = I_VALUE;                   //  �ȏm����
float         I_1[I_LENGTH]  = I_VALUE;                   //  �ȏm����
float         I_2[I_LENGTH]  = I_VALUE;                   //  �ȏm����
float         X[X_LENGTH]  = X_VALUE;                   //  �g?О�Fǻ�����t
float         X_1[X_LENGTH]  = X_VALUE;                   //  �g?О�Fǻ�����t
float         X_2[X_LENGTH]  = X_VALUE;                   //  �g?О�Fǻ�����t
float         P[P_LENGTH]  = P_VALUE;                   //  �g?О�Fǻ�����tǻ�tԴ��
float         P_1[P_LENGTH]  = P_VALUE;                   //  �g?О�Fǻ�����tǻ�tԴ��
float         P_2[P_LENGTH]  = P_VALUE;                   //  �g?О�Fǻ�����tǻ�tԴ��
float         K[K_LENGTH]  = K_VALUE;                   //  �U���C���y
float         K_1[K_LENGTH]  = K_VALUE;                   //  �U���C���y
float         K_2[K_LENGTH]  = K_VALUE;                   //  �U���C���y
float         Temp1[1]     = {0};                       //  �������v
float         Temp1_1[1]     = {0};                       //  �������v
float         Temp1_2[1]     = {0};                       //  �������v
//============================================================================//
//==                    �U���C�]��z�����ǻ���v                            ==//
//============================================================================//
float         A[A_LENGTH]       = A_VALUE;              //  О�F�v�����
float         A_1[A_LENGTH]       = A_VALUE;              //  О�F�v�����
float         A_2[A_LENGTH]       = A_VALUE;              //  О�F�v�����
float         B[B_LENGTH]       = B_VALUE;              //  ���f�y�f
float         B_1[B_LENGTH]       = B_VALUE;              //  ���f�y�f
float         B_2[B_LENGTH]       = B_VALUE;              //  ���f�y�f
float         Q[Q_LENGTH]       = Q_VALUE;              //  ���f����ǻ�tԴ��
float         Q_1[Q_LENGTH]       = Q_VALUE;              //  ���f����ǻ�tԴ��
float         Q_2[Q_LENGTH]       = Q_VALUE;              //  ���f����ǻ�tԴ��
float         C[C_LENGTH]       = C_VALUE;              //  ��������v�����
float         C_1[C_LENGTH]       = C_VALUE;              //  ��������v�����
float         C_2[C_LENGTH]       = C_VALUE;              //  ��������v�����
float         R[R_LENGTH]       = R_VALUE;              //  ���v����ǻ�tԴ��
float         R_1[R_LENGTH]       = R_VALUE;              //  ���v����ǻ�tԴ��
float         R_2[R_LENGTH]       = R_VALUE;              //  ���v����ǻ�tԴ��
float         Temp2[X_LENGTH]   = X_VALUE;              //  �������v, ���C����tOpt.XPreOpt[]ǻ���a�w�t
float         Temp2_1[X_LENGTH]   = X_VALUE;              //  �������v, ���C����tOpt.XPreOpt[]ǻ���a�w�t
float         Temp2_2[X_LENGTH]   = X_VALUE;              //  �������v, ���C����tOpt.XPreOpt[]ǻ���a�w�t
float         Temp22[X_LENGTH]  = X_VALUE;              //  �������v
float         Temp22_1[X_LENGTH]  = X_VALUE;              //  �������v
float         Temp22_2[X_LENGTH]  = X_VALUE;              //  �������v
float         Temp4[P_LENGTH]   = P_VALUE;              //  �������v, ���C����tCov.PPreOpt[]ǻ���a�w�t
float         Temp4_1[P_LENGTH]   = P_VALUE;              //  �������v, ���C����tCov.PPreOpt[]ǻ���a�w�t
float         Temp4_2[P_LENGTH]   = P_VALUE;              //  �������v, ���C����tCov.PPreOpt[]ǻ���a�w�t


//============================================================================//
//==                          �U���C�]��                                    ==//
//============================================================================//
//==?�O�y�f: ��                                                            ==//
//==���O�y�f: ��                                                            ==//
//==��϶�t:   ��                                                            ==//
//============================================================================//
float Watch1[N]={0};
float Watch2[N]={0};
float Watch3[N]={0};

void KalMan_Init(void)
{
	unsigned char   i;
	for (i=0; i<X_LENGTH; i++)
	{
		tOpt.XPreOpt[i] = Temp2[i];           //�^�t���a�w
	}
	for (i=0; i<P_LENGTH; i++)
	{
		tCov.PPreOpt[i] = Temp4[i];           //�^�t���a�w
	}
}

void KalMan_Init_1(void)
{
	unsigned char   i;
	for (i=0; i<X_LENGTH; i++)
	{
		tOpt_1.XPreOpt[i] = Temp2_1[i];           //�^�t���a�w
	}
	for (i=0; i<P_LENGTH; i++)
	{
		tCov_1.PPreOpt[i] = Temp4_1[i];           //�^�t���a�w
	}
}

void KalMan_Init_2(void)
{
	unsigned char   i;
	for (i=0; i<X_LENGTH; i++)
	{
		tOpt_1.XPreOpt[i] = Temp2_1[i];           //�^�t���a�w
	}
	for (i=0; i<P_LENGTH; i++)
	{
		tCov_1.PPreOpt[i] = Temp4_1[i];           //�^�t���a�w
	}
}

float KalMan(float input)
{
	unsigned char   i,j;
  //for (j=0; j<N; j++)
  {
    //Watch1[j] = 100 + j*2;
    //Watch1[j] = input;

    //Y[0] = Watch1[j] + Random1(0, 0.4);
   // Y[0] = Watch1[j] + (rand()%20)-10;
    Y[0] = input;
    //Watch2[j] = Y[0];
    MatrixMul(A, tOpt.XPreOpt, X, A_ROW, X_ROW, X_COLUMN);       //  �r�����fǻ�f��О�F�W��������О�F; X(k|k-1) = A(k,k-1)*X(k-1|k-1)
    
    MatrixCal1(A, tCov.PPreOpt, Temp4, SYS_ORDER);
    MatrixAdd(Temp4, Q, P, P_ROW, P_COLUMN);                     //  �����f��ǻ�tԴ������; P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q
    
    MatrixCal2(C, P, Temp1, C_ROW, C_COLUMN);
    MatrixAdd(Temp1, R, Temp1, R_ROW, R_COLUMN);
    Gauss_Jordan(Temp1, C_ROW);
    MatrixTrans(C, Temp2, C_ROW, C_COLUMN);
    MatrixMul(P, Temp2, Temp22, P_ROW, C_COLUMN, C_ROW);
    MatrixMul(Temp22, Temp1, K, P_ROW, C_ROW, C_ROW);            //  ����U���C���y; K(k) = P(k|k-1)*C' / (C(k)*P(k|k-1)*C(k)' + R)
    
    MatrixMul(C, X, Temp1, C_ROW, X_ROW, X_COLUMN);
    MatrixMinus(Y, Temp1, Temp1, Y_ROW, Y_COLUMN);
    MatrixMul(K, Temp1, Temp2, K_ROW, Y_ROW, Y_COLUMN);
    MatrixAdd(X, Temp2, tOpt.XNowOpt, X_ROW, X_COLUMN);          //  �H�������t����v�t����g?�I�b�t; X(k|k) = X(k|k-1)+Kg(k)*(Y(k)-C*X(k|k-1))
    
    MatrixMul(K, C, Temp4, K_ROW, C_ROW, C_COLUMN);
    MatrixMinus(I, Temp4, Temp4, I_ROW, I_COLUMN);
    MatrixMul(Temp4, P, tCov.PNowOpt, I_ROW, P_ROW, P_COLUMN);   //  �����d��Y��t�紬����; P(k|k) =��I-Kg(k)*C��*P(k|k-1)
    
    for (i=0; i<X_LENGTH; i++)
    {
      tOpt.XPreOpt[i] = tOpt.XNowOpt[i];
    }
    for (i=0; i<P_LENGTH; i++)
    {
      tCov.PPreOpt[i] = tCov.PNowOpt[i];
    }
   // Watch3[j] = tOpt.XNowOpt[0];
    
  }//end of for
  return tOpt.XNowOpt[0];
}

float KalMan_1(float input)
{
	unsigned char   i,j;
  //for (j=0; j<N; j++)
  {
    //Watch1[j] = 100 + j*2;
    //Watch1[j] = input;

    //Y[0] = Watch1[j] + Random1(0, 0.4);
   // Y[0] = Watch1[j] + (rand()%20)-10;
    Y_1[0] = input;
    //Watch2[j] = Y[0];
    MatrixMul(A_1, tOpt_1.XPreOpt, X_1, A_ROW, X_ROW, X_COLUMN);       //  �r�����fǻ�f��О�F�W��������О�F; X(k|k-1) = A(k,k-1)*X(k-1|k-1)
    
    MatrixCal1(A_1, tCov_1.PPreOpt, Temp4_1, SYS_ORDER);
    MatrixAdd(Temp4_1, Q_1, P_1, P_ROW, P_COLUMN);                     //  �����f��ǻ�tԴ������; P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q
    
    MatrixCal2(C_1, P_1, Temp1_1, C_ROW, C_COLUMN);
    MatrixAdd(Temp1_1, R_1, Temp1_1, R_ROW, R_COLUMN);
    Gauss_Jordan(Temp1_1, C_ROW);
    MatrixTrans(C_1, Temp2_1, C_ROW, C_COLUMN);
    MatrixMul(P_1, Temp2_1, Temp22_1, P_ROW, C_COLUMN, C_ROW);
    MatrixMul(Temp22_1, Temp1_1, K_1, P_ROW, C_ROW, C_ROW);            //  ����U���C���y; K(k) = P(k|k-1)*C' / (C(k)*P(k|k-1)*C(k)' + R)
    
    MatrixMul(C_1, X_1, Temp1_1, C_ROW, X_ROW, X_COLUMN);
    MatrixMinus(Y_1, Temp1_1, Temp1_1, Y_ROW, Y_COLUMN);
    MatrixMul(K_1, Temp1_1, Temp2_1, K_ROW, Y_ROW, Y_COLUMN);
    MatrixAdd(X_1, Temp2_1, tOpt_1.XNowOpt, X_ROW, X_COLUMN);          //  �H�������t����v�t����g?�I�b�t; X(k|k) = X(k|k-1)+Kg(k)*(Y(k)-C*X(k|k-1))
    
    MatrixMul(K_1, C_1, Temp4_1, K_ROW, C_ROW, C_COLUMN);
    MatrixMinus(I_1, Temp4_1, Temp4_1, I_ROW, I_COLUMN);
    MatrixMul(Temp4_1, P_1, tCov_1.PNowOpt, I_ROW, P_ROW, P_COLUMN);   //  �����d��Y��t�紬����; P(k|k) =��I-Kg(k)*C��*P(k|k-1)
    
    for (i=0; i<X_LENGTH; i++)
    {
      tOpt_1.XPreOpt[i] = tOpt_1.XNowOpt[i];
    }
    for (i=0; i<P_LENGTH; i++)
    {
      tCov_1.PPreOpt[i] = tCov_1.PNowOpt[i];
    }
   // Watch3[j] = tOpt.XNowOpt[0];
    
  }//end of for
  return tOpt_1.XNowOpt[0];
}

float KalMan_2(float input)
{
	unsigned char   i,j;
  //for (j=0; j<N; j++)
  {
    //Watch1[j] = 100 + j*2;
    //Watch1[j] = input;

    //Y[0] = Watch1[j] + Random1(0, 0.4);
   // Y[0] = Watch1[j] + (rand()%20)-10;
    Y_2[0] = input;
    //Watch2[j] = Y[0];
    MatrixMul(A_2, tOpt_2.XPreOpt, X_2, A_ROW, X_ROW, X_COLUMN);       //  �r�����fǻ�f��О�F�W��������О�F; X(k|k-1) = A(k,k-1)*X(k-1|k-1)
    
    MatrixCal1(A_2, tCov_2.PPreOpt, Temp4_2, SYS_ORDER);
    MatrixAdd(Temp4_2, Q_2, P_2, P_ROW, P_COLUMN);                     //  �����f��ǻ�tԴ������; P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q

    MatrixCal2(C_2, P_2, Temp1_2, C_ROW, C_COLUMN);
    MatrixAdd(Temp1_2, R_2, Temp1_2, R_ROW, R_COLUMN);
    Gauss_Jordan(Temp1_2, C_ROW);
    MatrixTrans(C_2, Temp2_2, C_ROW, C_COLUMN);
    MatrixMul(P_2, Temp2_2, Temp22_2, P_ROW, C_COLUMN, C_ROW);
    MatrixMul(Temp22_2, Temp1_2, K_2, P_ROW, C_ROW, C_ROW);            //  ����U���C���y; K(k) = P(k|k-1)*C' / (C(k)*P(k|k-1)*C(k)' + R)
    
    MatrixMul(C_2, X_2, Temp1_2, C_ROW, X_ROW, X_COLUMN);
    MatrixMinus(Y_2, Temp1_2, Temp1_2, Y_ROW, Y_COLUMN);
    MatrixMul(K_2, Temp1_2, Temp2_2, K_ROW, Y_ROW, Y_COLUMN);
    MatrixAdd(X_2, Temp2_2, tOpt_2.XNowOpt, X_ROW, X_COLUMN);          //  �H�������t����v�t����g?�I�b�t; X(k|k) = X(k|k-1)+Kg(k)*(Y(k)-C*X(k|k-1))
    
    MatrixMul(K_2, C_2, Temp4_2, K_ROW, C_ROW, C_COLUMN);
    MatrixMinus(I_2, Temp4_2, Temp4_2, I_ROW, I_COLUMN);
    MatrixMul(Temp4_2, P_2, tCov_2.PNowOpt, I_ROW, P_ROW, P_COLUMN);   //  �����d��Y��t�紬����; P(k|k) =��I-Kg(k)*C��*P(k|k-1)
    
    for (i=0; i<X_LENGTH; i++)
    {
      tOpt_2.XPreOpt[i] = tOpt_2.XNowOpt[i];
    }
    for (i=0; i<P_LENGTH; i++)
    {
      tCov_2.PPreOpt[i] = tCov_2.PNowOpt[i];
    }
   // Watch3[j] = tOpt.XNowOpt[0];
    
  }//end of for
  return tOpt_2.XNowOpt[0];
}
