void LegendreLobattoAndDerivative(int N, double x, double *q, double *qp, double *LN)
{
    double LNm2, LNm1, LNm2p, LNm1p, LNp, LNp1, LNp1p, LNl, o;
    
    LNm2  = 1;
    LNm1  = x;
    LNm2p = 0;
    LNm1p = 1;
    
    /* Looping through the degree */
    for (o=2; o<=N; o++){
        LNl = (2*o - 1)/o * x * LNm1 - (o - 1)/o * LNm2;
        LNp   = LNm2p + (2*o - 1)*LNm1;
        
        if (o<N){
            LNm2  = LNm1;
            LNm1  = LNl;
            
            LNm2p = LNm1p;
            LNm1p = LNp;
        }
    }
    
    /* Calculating the plus one degree*/
    
    
    LNp1  = (2*o - 1)/o * x * LNl - (o - 1)/o * LNm1;
    LNp1p = LNm1p + (2*o - 1) * LNl;
    *LN = LNl;
    *q  = LNp1  - LNm1;
    *qp = LNp1p - LNm1p;
}