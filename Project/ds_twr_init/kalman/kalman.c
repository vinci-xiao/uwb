#include  "kalman.h"
#include  "stdlib.h"
#include  "math.h"


/*==============================================================================
1.†•Âï”µ
   X(k|k-1) = A(k,k-1)*X(k-1|k-1)                 //ÖS¶‡ÖvR0


2.”µ…ì†•Âï”µÐtÔ´´¬“éäÇ
   P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q(k)
   Q(k) = U(k)¨GU(k)' 


3.”µ…ì¿UÏÓìC”µy“éäÇ
   K(k) = P(k|k-1)*C(k)' / (C(k)*P(k|k-1)*C(k)' + R(k))
   R(k) = N(k)¨GN(k)' 


4.ÝdÚëÂï”µ
   X(k|k) = X(k|k-1)+K(k)*(Y(k)-C(k)*X(k|k-1))


5.”µ…ìÝdÚë¾YÂï”µÐtœç´¬“éäÇ
   P(k|k) =¨×I-K(k)*C(k)¨Ø*P(k|k-1)


6. ÝdÚëàIÍb³t


A(k,k-1):     ÐžFÍvðä“éäÇ
B(k,k-1):     ÐžFÇ»ÖS¶‡Öv
C(k):         â¹ñö“éäÇÍvðä“éäÇ
X(k|k-1):     ÛHÀÞk-1ŠCÓDÇ»àIÍb³tÂï”µkŠCÓDÇ»³t
X(k-1|k-1):   k-1ŠCÓDÇ»àIÍb³t
P(k|k-1):     X(k|k-1)ÇÚÜíÇ»covariance
P(k-1|k-1):   X(k-1|k-1)ÇÚÜíÇ»covariance
Q(k):         žçÆfØ×îÇ»covariance(ÇÚŠf«”°ôñöÖvÂï”µ³tÇ»ê|?×îƒH)
R(k):         ñöÖvØ×îÇ»ÐtÔ´´¬
Y(k):         kŠCÓDÇ»ñöÖv³t
K(k):         ¿UÏÓìC”µy
U(k):         kŠCÓDÐÛFŠñšû
N(k):         kŠCÓDâ¹ñöŠñšû


î«?:     žçÆfÐžFX-----------æ¡•¿yƒH
          žçÆf“éäÇA-----------yƒH²ÜÚwÍvðä
          ÐžFÇ»ÖS¶‡ÖvB-------(»o¶¼Áƒô¬)
          â¹ñö³tY-------------yƒH”µÊò–f
          â¹ñö“éäÇC-----------’GŽÞƒH-¡÷Ã²ŽÞƒH
          Ø×îŠñšûQ-----------yƒH²ÜÚw?´¬
          ñöÖvŠñšûR-----------Êò–f•d´¬

ƒr’ìä¿: ìýçæÛHÀÞŠf«”°ô(?•ŒëÈC«”°ôŒuÛHÀÞ†•¶­³t”µ…ì)Ç»–fÀÞ”µ…ìµÌ’ì°ôÇ»Âï”µ³t,
          °¹ Z,ÛHÀÞŠf«”°ôÇ»–fÀÞ”µ…ìµÌ’ì°ôÂï”µ³tÇ»ÐtÔ´´¬;  Õ†Òƒ,²Ï’ì°ôÂï”µ³tÇ»Ðt
          Ô´´¬”µ…ìµÌ¿UÏÓìC”µy;  àI¾Y,ÛHÀÞÂïñö³tî£ñöÖv³t”µ…ì½g?àIÍb³t“´?ÐtÔ´´¬
          å‚X,¿UÏÓìCË]ÊèíÏŒëŠCàKŠfÇ»½Ò Z,ñÝ±p?R?Þ·ÇIÅø½ßÜÖ
==============================================================================*/



//================================================//
//==             àIÍb³tÔ´´¬ÙxµÊ¼«               ==//
//================================================//
typedef struct  _tCovariance
{
  float PNowOpt[P_LENGTH];
  float PPreOpt[P_LENGTH];
}tCovariance;



//================================================//
//==               àIÍb³tÙxµÊ¼«                 ==//
//================================================//
typedef struct  _tOptimal
{
  float XNowOpt[X_LENGTH];
  float XPreOpt[X_LENGTH];
}tOptimal;



tOptimal      tOpt;                                     //  ŠóË]ÊèB–f¸œÝv‚kˆöŒaÚwå^³t
tOptimal      tOpt_1;                                     //  ŠóË]ÊèB–f¸œÝv‚kˆöŒaÚwå^³t
tOptimal      tOpt_2;                                     //  ŠóË]ÊèB–f¸œÝv‚kˆöŒaÚwå^³t

tCovariance   tCov;                                     //  ŠóË]ÊèB–f¸œÝv‚kˆöŒaÚwå^³t
tCovariance   tCov_1;                                     //  ŠóË]ÊèB–f¸œÝv‚kˆöŒaÚwå^³t
tCovariance   tCov_2;                                     //  ŠóË]ÊèB–f¸œÝv‚kˆöŒaÚwå^³t

float         Y[Y_LENGTH]  = Y_VALUE;                   //  ñöÖv³t(·ª°ôñöÖvÇ»–fÀÞ„záûäÕ?Í˜–fàV)
float         Y_1[Y_LENGTH]  = Y_VALUE;                   //  ñöÖv³t(·ª°ôñöÖvÇ»–fÀÞ„záûäÕ?Í˜–fàV)
float         Y_2[Y_LENGTH]  = Y_VALUE;                   //  ñöÖv³t(·ª°ôñöÖvÇ»–fÀÞ„záûäÕ?Í˜–fàV)
float         I[I_LENGTH]  = I_VALUE;                   //  µÈm“éäÇ
float         I_1[I_LENGTH]  = I_VALUE;                   //  µÈm“éäÇ
float         I_2[I_LENGTH]  = I_VALUE;                   //  µÈm“éäÇ
float         X[X_LENGTH]  = X_VALUE;                   //  ½g?ÐžFÇ»†•ñö³t
float         X_1[X_LENGTH]  = X_VALUE;                   //  ½g?ÐžFÇ»†•ñö³t
float         X_2[X_LENGTH]  = X_VALUE;                   //  ½g?ÐžFÇ»†•ñö³t
float         P[P_LENGTH]  = P_VALUE;                   //  ½g?ÐžFÇ»†•ñö³tÇ»ÐtÔ´´¬
float         P_1[P_LENGTH]  = P_VALUE;                   //  ½g?ÐžFÇ»†•ñö³tÇ»ÐtÔ´´¬
float         P_2[P_LENGTH]  = P_VALUE;                   //  ½g?ÐžFÇ»†•ñö³tÇ»ÐtÔ´´¬
float         K[K_LENGTH]  = K_VALUE;                   //  ¿UÏÓìC”µy
float         K_1[K_LENGTH]  = K_VALUE;                   //  ¿UÏÓìC”µy
float         K_2[K_LENGTH]  = K_VALUE;                   //  ¿UÏÓìC”µy
float         Temp1[1]     = {0};                       //  ÂäÁ”²ÜÖv
float         Temp1_1[1]     = {0};                       //  ÂäÁ”²ÜÖv
float         Temp1_2[1]     = {0};                       //  ÂäÁ”²ÜÖv
//============================================================================//
//==                    ¿UÏÓìCË]Êè„záûðÀëÇ»²ÜÖv                            ==//
//============================================================================//
float         A[A_LENGTH]       = A_VALUE;              //  ÐžFÍvðä“éäÇ
float         A_1[A_LENGTH]       = A_VALUE;              //  ÐžFÍvðä“éäÇ
float         A_2[A_LENGTH]       = A_VALUE;              //  ÐžFÍvðä“éäÇ
float         B[B_LENGTH]       = B_VALUE;              //  žçÆf½y–f
float         B_1[B_LENGTH]       = B_VALUE;              //  žçÆf½y–f
float         B_2[B_LENGTH]       = B_VALUE;              //  žçÆf½y–f
float         Q[Q_LENGTH]       = Q_VALUE;              //  žçÆfØ×îÇ»ÐtÔ´´¬
float         Q_1[Q_LENGTH]       = Q_VALUE;              //  žçÆfØ×îÇ»ÐtÔ´´¬
float         Q_2[Q_LENGTH]       = Q_VALUE;              //  žçÆfØ×îÇ»ÐtÔ´´¬
float         C[C_LENGTH]       = C_VALUE;              //  â¹ñö“éäÇÍvðä“éäÇ
float         C_1[C_LENGTH]       = C_VALUE;              //  â¹ñö“éäÇÍvðä“éäÇ
float         C_2[C_LENGTH]       = C_VALUE;              //  â¹ñö“éäÇÍvðä“éäÇ
float         R[R_LENGTH]       = R_VALUE;              //  ñöÖvØ×îÇ»ÐtÔ´´¬
float         R_1[R_LENGTH]       = R_VALUE;              //  ñöÖvØ×îÇ»ÐtÔ´´¬
float         R_2[R_LENGTH]       = R_VALUE;              //  ñöÖvØ×îÇ»ÐtÔ´´¬
float         Temp2[X_LENGTH]   = X_VALUE;              //  ÂäÁ”²ÜÖv, °¹ŠCäÕtOpt.XPreOpt[]Ç»ˆöŒaÚw³t
float         Temp2_1[X_LENGTH]   = X_VALUE;              //  ÂäÁ”²ÜÖv, °¹ŠCäÕtOpt.XPreOpt[]Ç»ˆöŒaÚw³t
float         Temp2_2[X_LENGTH]   = X_VALUE;              //  ÂäÁ”²ÜÖv, °¹ŠCäÕtOpt.XPreOpt[]Ç»ˆöŒaÚw³t
float         Temp22[X_LENGTH]  = X_VALUE;              //  ÂäÁ”²ÜÖv
float         Temp22_1[X_LENGTH]  = X_VALUE;              //  ÂäÁ”²ÜÖv
float         Temp22_2[X_LENGTH]  = X_VALUE;              //  ÂäÁ”²ÜÖv
float         Temp4[P_LENGTH]   = P_VALUE;              //  ÂäÁ”²ÜÖv, °¹ŠCäÕtCov.PPreOpt[]Ç»ˆöŒaÚw³t
float         Temp4_1[P_LENGTH]   = P_VALUE;              //  ÂäÁ”²ÜÖv, °¹ŠCäÕtCov.PPreOpt[]Ç»ˆöŒaÚw³t
float         Temp4_2[P_LENGTH]   = P_VALUE;              //  ÂäÁ”²ÜÖv, °¹ŠCäÕtCov.PPreOpt[]Ç»ˆöŒaÚw³t


//============================================================================//
//==                          ¿UÏÓìCË]Êè                                    ==//
//============================================================================//
//==?ÖO½y–f: ’                                                            ==//
//==µÌÖO½y–f: ’                                                            ==//
//==µîÏ¶³t:   ’                                                            ==//
//============================================================================//
float Watch1[N]={0};
float Watch2[N]={0};
float Watch3[N]={0};

void KalMan_Init(void)
{
	unsigned char   i;
	for (i=0; i<X_LENGTH; i++)
	{
		tOpt.XPreOpt[i] = Temp2[i];           //å^³tˆöŒaÚw
	}
	for (i=0; i<P_LENGTH; i++)
	{
		tCov.PPreOpt[i] = Temp4[i];           //å^³tˆöŒaÚw
	}
}

void KalMan_Init_1(void)
{
	unsigned char   i;
	for (i=0; i<X_LENGTH; i++)
	{
		tOpt_1.XPreOpt[i] = Temp2_1[i];           //å^³tˆöŒaÚw
	}
	for (i=0; i<P_LENGTH; i++)
	{
		tCov_1.PPreOpt[i] = Temp4_1[i];           //å^³tˆöŒaÚw
	}
}

void KalMan_Init_2(void)
{
	unsigned char   i;
	for (i=0; i<X_LENGTH; i++)
	{
		tOpt_1.XPreOpt[i] = Temp2_1[i];           //å^³tˆöŒaÚw
	}
	for (i=0; i<P_LENGTH; i++)
	{
		tCov_1.PPreOpt[i] = Temp4_1[i];           //å^³tˆöŒaÚw
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
    MatrixMul(A, tOpt.XPreOpt, X, A_ROW, X_ROW, X_COLUMN);       //  ƒrôÀžçÆfÇ»Šf«”ÐžFŠW†•ñö«€ŠóÐžF; X(k|k-1) = A(k,k-1)*X(k-1|k-1)
    
    MatrixCal1(A, tCov.PPreOpt, Temp4, SYS_ORDER);
    MatrixAdd(Temp4, Q, P, P_ROW, P_COLUMN);                     //  †•ñö–fÀÞÇ»ÐtÔ´´¬“éäÇ; P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q
    
    MatrixCal2(C, P, Temp1, C_ROW, C_COLUMN);
    MatrixAdd(Temp1, R, Temp1, R_ROW, R_COLUMN);
    Gauss_Jordan(Temp1, C_ROW);
    MatrixTrans(C, Temp2, C_ROW, C_COLUMN);
    MatrixMul(P, Temp2, Temp22, P_ROW, C_COLUMN, C_ROW);
    MatrixMul(Temp22, Temp1, K, P_ROW, C_ROW, C_ROW);            //  ”µ…ì¿UÏÓìC”µy; K(k) = P(k|k-1)*C' / (C(k)*P(k|k-1)*C(k)' + R)
    
    MatrixMul(C, X, Temp1, C_ROW, X_ROW, X_COLUMN);
    MatrixMinus(Y, Temp1, Temp1, Y_ROW, Y_COLUMN);
    MatrixMul(K, Temp1, Temp2, K_ROW, Y_ROW, Y_COLUMN);
    MatrixAdd(X, Temp2, tOpt.XNowOpt, X_ROW, X_COLUMN);          //  ÛHÀÞÂïñö³tî£ñöÖv³t”µ…ì½g?àIÍb³t; X(k|k) = X(k|k-1)+Kg(k)*(Y(k)-C*X(k|k-1))
    
    MatrixMul(K, C, Temp4, K_ROW, C_ROW, C_COLUMN);
    MatrixMinus(I, Temp4, Temp4, I_ROW, I_COLUMN);
    MatrixMul(Temp4, P, tCov.PNowOpt, I_ROW, P_ROW, P_COLUMN);   //  ”µ…ìÝdÚë¾YÂï”µÐtœç´¬“éäÇ; P(k|k) =¨×I-Kg(k)*C¨Ø*P(k|k-1)
    
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
    MatrixMul(A_1, tOpt_1.XPreOpt, X_1, A_ROW, X_ROW, X_COLUMN);       //  ƒrôÀžçÆfÇ»Šf«”ÐžFŠW†•ñö«€ŠóÐžF; X(k|k-1) = A(k,k-1)*X(k-1|k-1)
    
    MatrixCal1(A_1, tCov_1.PPreOpt, Temp4_1, SYS_ORDER);
    MatrixAdd(Temp4_1, Q_1, P_1, P_ROW, P_COLUMN);                     //  †•ñö–fÀÞÇ»ÐtÔ´´¬“éäÇ; P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q
    
    MatrixCal2(C_1, P_1, Temp1_1, C_ROW, C_COLUMN);
    MatrixAdd(Temp1_1, R_1, Temp1_1, R_ROW, R_COLUMN);
    Gauss_Jordan(Temp1_1, C_ROW);
    MatrixTrans(C_1, Temp2_1, C_ROW, C_COLUMN);
    MatrixMul(P_1, Temp2_1, Temp22_1, P_ROW, C_COLUMN, C_ROW);
    MatrixMul(Temp22_1, Temp1_1, K_1, P_ROW, C_ROW, C_ROW);            //  ”µ…ì¿UÏÓìC”µy; K(k) = P(k|k-1)*C' / (C(k)*P(k|k-1)*C(k)' + R)
    
    MatrixMul(C_1, X_1, Temp1_1, C_ROW, X_ROW, X_COLUMN);
    MatrixMinus(Y_1, Temp1_1, Temp1_1, Y_ROW, Y_COLUMN);
    MatrixMul(K_1, Temp1_1, Temp2_1, K_ROW, Y_ROW, Y_COLUMN);
    MatrixAdd(X_1, Temp2_1, tOpt_1.XNowOpt, X_ROW, X_COLUMN);          //  ÛHÀÞÂïñö³tî£ñöÖv³t”µ…ì½g?àIÍb³t; X(k|k) = X(k|k-1)+Kg(k)*(Y(k)-C*X(k|k-1))
    
    MatrixMul(K_1, C_1, Temp4_1, K_ROW, C_ROW, C_COLUMN);
    MatrixMinus(I_1, Temp4_1, Temp4_1, I_ROW, I_COLUMN);
    MatrixMul(Temp4_1, P_1, tCov_1.PNowOpt, I_ROW, P_ROW, P_COLUMN);   //  ”µ…ìÝdÚë¾YÂï”µÐtœç´¬“éäÇ; P(k|k) =¨×I-Kg(k)*C¨Ø*P(k|k-1)
    
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
    MatrixMul(A_2, tOpt_2.XPreOpt, X_2, A_ROW, X_ROW, X_COLUMN);       //  ƒrôÀžçÆfÇ»Šf«”ÐžFŠW†•ñö«€ŠóÐžF; X(k|k-1) = A(k,k-1)*X(k-1|k-1)
    
    MatrixCal1(A_2, tCov_2.PPreOpt, Temp4_2, SYS_ORDER);
    MatrixAdd(Temp4_2, Q_2, P_2, P_ROW, P_COLUMN);                     //  †•ñö–fÀÞÇ»ÐtÔ´´¬“éäÇ; P(k|k-1) = A(k,k-1)*P(k-1|k-1)*A(k,k-1)'+Q

    MatrixCal2(C_2, P_2, Temp1_2, C_ROW, C_COLUMN);
    MatrixAdd(Temp1_2, R_2, Temp1_2, R_ROW, R_COLUMN);
    Gauss_Jordan(Temp1_2, C_ROW);
    MatrixTrans(C_2, Temp2_2, C_ROW, C_COLUMN);
    MatrixMul(P_2, Temp2_2, Temp22_2, P_ROW, C_COLUMN, C_ROW);
    MatrixMul(Temp22_2, Temp1_2, K_2, P_ROW, C_ROW, C_ROW);            //  ”µ…ì¿UÏÓìC”µy; K(k) = P(k|k-1)*C' / (C(k)*P(k|k-1)*C(k)' + R)
    
    MatrixMul(C_2, X_2, Temp1_2, C_ROW, X_ROW, X_COLUMN);
    MatrixMinus(Y_2, Temp1_2, Temp1_2, Y_ROW, Y_COLUMN);
    MatrixMul(K_2, Temp1_2, Temp2_2, K_ROW, Y_ROW, Y_COLUMN);
    MatrixAdd(X_2, Temp2_2, tOpt_2.XNowOpt, X_ROW, X_COLUMN);          //  ÛHÀÞÂïñö³tî£ñöÖv³t”µ…ì½g?àIÍb³t; X(k|k) = X(k|k-1)+Kg(k)*(Y(k)-C*X(k|k-1))
    
    MatrixMul(K_2, C_2, Temp4_2, K_ROW, C_ROW, C_COLUMN);
    MatrixMinus(I_2, Temp4_2, Temp4_2, I_ROW, I_COLUMN);
    MatrixMul(Temp4_2, P_2, tCov_2.PNowOpt, I_ROW, P_ROW, P_COLUMN);   //  ”µ…ìÝdÚë¾YÂï”µÐtœç´¬“éäÇ; P(k|k) =¨×I-Kg(k)*C¨Ø*P(k|k-1)
    
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
