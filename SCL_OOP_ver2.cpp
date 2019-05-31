/****************************************************************************************/
/* Target: Optimum Overlap (OO) partitioning with Locality Constraints                  */
/* By: Homa Esfahanizadeh (hesfahanizadeh@ucla.edu)                                     */
/* University of California, Los Angeles (UCLA)                                         */
/* Update: 05/30/2019                                                                   */
/****************************************************************************************/

#include <iostream>
#include <math.h>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <fstream>
#include <time.h>
using namespace std;

// n choose k
int comb_choose ( int n , int k );

// Finding the number of cycle-6 in H (brute force)
int cycle_6_enumeration_protograph_block ( int** H , int m , int n );

// Finding the number of cycle-4 in H (brute force)
int cycle_4_enumeration_protograph_block ( int** H , int m , int n );


// Finding the number of cycle-6 given a set of overlap parameters when all overlaps belong to one replica
int A ( int tabc , int tab , int tac , int tbc );

// Finding the number of cycle-6 given a set of overlap parameters when two overlaps belong to one replica and one overlap belongs to the other one
int B ( int tabc , int tab , int tac , int tbc );

// Finding the number of cycle-6 given a set of overlap parameters when each overlap belongs to an individual replica
int C ( int tab , int tac , int tbc );

// Finding the number of cycle-6 in H_SC using overlapping parameters
int cycle_6_enumeration_protograph_SC ( int kappa , int gamma , int L , int m , int*** OM3 , int** OM2 , int* OM1 );

// Finding a partitioning that achieve a particular set of overlapping parameters
// The output is the number of partitionig matrices that achieve the set of given overlapping parameters. But, PM is filled by one of the choices
int Partitioning_Matrix_Construction ( int kappa , int gamma_c , int gamma_l , int m , int* tdeg3 , int* tdeg2 , int* tdeg1 , int** PM );

int main ()
{
    int kappa = 13;
    int p  = 1;
    int L = 10;
    int m = 1;
    int gamma_c = 3; // coupling checks
    int gamma_l = 2; // local checks
    int gamma_a = 0; // appended local checks with no optimization
    int gamma = gamma_c + gamma_l;

    
    // Overlap Parameters
    int*** OM3 = new int** [(m+1)*gamma];
    for ( int i = 0 ; i < ((m+1)*gamma) ; i++ )
    {
        OM3[i] = new int* [(m+1)*gamma];
        for ( int j = 0 ; j < ((m+1)*gamma) ; j++ )
        {
            OM3[i][j] = new int [(m+1)*gamma];
            for ( int k = 0 ; k < ((m+1)*gamma) ; k++ )
            {
                OM3[i][j][k] = 0;
            }
        }
    }
    int** OM2 = new int* [(m+1)*gamma];
    for ( int i = 0 ; i < ((m+1)*gamma) ; i++ )
    {
        OM2[i] = new int [(m+1)*gamma];
        for ( int j = 0 ; j < ((m+1)*gamma) ; j++ )
        {
            OM2[i][j] = 0;
        }
    }
    int* OM1 = new int [(m+1)*gamma];
    for ( int i = 0 ; i < ((m+1)*gamma) ; i++ )
    {
        OM1[i] = 0;
    }
    
    int F_min = -1;
    bool first_trial = true;
    
    int* tdeg3_opt;
    int* tdeg2_opt;
    int* tdeg1_opt;
    
    if ( ( gamma_c == 3 ) && ( m == 1 ) )
    {
        tdeg3_opt = new int [1];
        tdeg2_opt = new int [3];
        tdeg1_opt = new int [3];
        
        for ( int t0 = 0 ; t0 <= kappa ; t0++ )
        {
            for ( int t01 = 0 ; t01 <= t0 ; t01++ )
            {
            for ( int t1 = t01 ; t1 <= (kappa-t0+t01) ; t1++ )
            {
                for ( int t012 = 0 ; t012 <= t01 ; t012++ )
                {
                for ( int t02 = t012 ; t02 <= (t0-t01+t012) ; t02++ )
                {
                for ( int t12 = t012 ; t12 <= (t1-t01+t012) ; t12++ )
                {
                for ( int t2 = (t02+t12-t012) ; t2 <= (kappa-t0-t1+t01+t02+t12-t012) ; t2++ )
                {
                //if ( ( ( t0 + t1 + t2 ) >= floor( kappa * gamma_c / 2 ) ) && ( ( t0 + t1 + t2 ) <= ceil( kappa * gamma_c / 2 ) ) )
                //{
                    OM3[0][1][2] = t012;
                    OM2[0][1] = t01;
                    OM2[0][2] = t02;
                    OM2[1][2] = t12;
                    OM1[0] = t0;
                    OM1[1] = t1;
                    OM1[2] = t2;
                    
                    for ( int k = gamma_c ; k < gamma ; k++ )
                    {
                        for ( int i = 0 ; i < (k-1) ; i++ )
                        {
                            for ( int j = i+1 ; j < k ; j++ )
                            {
                                if ( ( i < gamma_c ) && ( j < gamma_c ) )
                                    OM3[i][j][k] = OM2[i][j];
                                if ( ( i < gamma_c ) && ( j >= gamma_c ) )
                                    OM3[i][j][k] = OM1[i];
                                if ( ( i >= gamma_c ) && ( j >= gamma_c ) )
                                    OM3[i][j][k] = kappa;
                            }
                        }
                    }
                    
                    for ( int j = gamma_c ; j < gamma ; j++ )
                    {
                        for ( int i = 0 ; i < j ; i++ )
                        {
                            if ( i < gamma_c )
                                OM2[i][j] = OM1[i];
                            if ( i >= gamma_c )
                                OM2[i][j] = kappa;
                        }
                    }
                    
                    for ( int i = gamma_c ; i < gamma ; i++ )
                    {
                        OM1[i] = kappa;
                    }
                    
                    int F = cycle_6_enumeration_protograph_SC ( kappa , gamma , L , m , OM3 , OM2 , OM1  );
                    if ( ( first_trial ) || ( F < F_min ) )
                    {
                        F_min = F;
                        tdeg3_opt[0] = OM3[0][1][2];
                        tdeg2_opt[0] = OM2[0][1]; tdeg2_opt[1] = OM2[0][2]; tdeg2_opt[2] = OM2[1][2];
                        tdeg1_opt[0] = OM1[0]; tdeg1_opt[1] = OM1[1]; tdeg1_opt[2] = OM1[2];
                        first_trial = false;
                    }
                //}
                }
                }
                }
                }
            }
            }
        }
        cout << "Minimum Number of Cycles 6: " << F_min << endl;
        
        cout << "Optimum Set of Overlapping Parameters: " << endl;
        for ( int i = 0 ; i < 1 ; i++ ) cout << tdeg3_opt[i] << "   "; cout << endl;
        for ( int i = 0 ; i < 3 ; i++ ) cout << tdeg2_opt[i] << "   "; cout << endl;
        for ( int i = 0 ; i < 3 ; i++ ) cout << tdeg1_opt[i] << "   "; cout << endl;
        
    }

    if ( ( gamma_c == 4 ) && ( m == 1 ) )
    {
        tdeg3_opt = new int [4];
        tdeg2_opt = new int [6];
        tdeg1_opt = new int [4];
        
        for ( int t0 = 0 ; t0 <= kappa ; t0++ )
        {
            for ( int t01 = 0 ; t01 <= t0 ; t01++ )
            {
            for ( int t1 = t01 ; t1 <= (kappa-t0+t01) ; t1++ )
            {
                for ( int t012 = 0 ; t012 <= t01 ; t012++ )
                {
                for ( int t02 = t012 ; t02 <= (t0-t01+t012) ; t02++ )
                {
                for ( int t12 = t012 ; t12 <= (t1-t01+t012) ; t12++ )
                {
                for ( int t2 = (t02+t12-t012) ; t2 <= (kappa-t0-t1+t01+t02+t12-t012) ; t2++ )
                {
                    for ( int t0123 = 0 ; t0123 <= t012 ; t0123++ )
                    {
                    for ( int t013 = t0123 ; t013 <= ( t01 - t012 + t0123 ) ; t013++ )
                    {
                    for ( int t023 = t0123 ; t023 <= ( t02 - t012 + t0123 ) ; t023++ )
                    {
                    for ( int t123 = t0123 ; t123 <= ( t12 - t012 + t0123 ) ; t123++ )
                    {
                    for ( int t03 = ( t013 + t023 - t0123 ) ; t03 <= ( t0 - t01 - t02 + t012 + t013 + t023 - t0123 ) ; t03++  )
                    {
                    for ( int t13 = ( t013 + t123 - t0123 ) ; t13 <= ( t1 - t01 - t12 + t012 + t013 + t123 - t0123 ) ; t13++  )
                    {
                    for ( int t23 = ( t023 + t123 - t0123 ) ; t23 <= ( t2 - t02 - t12 + t012 + t023 + t123 - t0123 ) ; t23++  )
                    {
                    for ( int t3 = ( t03 + t13 + t23 - t013 - t023 - t123 + t0123 ) ; t3 <= ( kappa - t0 - t1 - t2 + t01 + t02 + t12 + t03 + t13 + t23 - t012 - t013 - t023 - t123 + t0123 ) ; t3++ )
                    {
                        //if ( ( t0 + t1 + t2 + t3 ) == ( kappa * gamma_c / 2 ) )
                        //{
                            OM3[0][1][2] = t012;
                            OM3[0][1][3] = t013;
                            OM3[0][2][3] = t023;
                            OM3[1][2][3] = t123;
                            OM2[0][1] = t01;
                            OM2[0][2] = t02;
                            OM2[0][3] = t03;
                            OM2[1][2] = t12;
                            OM2[1][3] = t13;
                            OM2[2][3] = t23;
                            OM1[0] = t0;
                            OM1[1] = t1;
                            OM1[2] = t2;
                            OM1[3] = t3;
                            
                            for ( int k = gamma_c ; k < gamma ; k++ )
                            {
                                for ( int i = 0 ; i < (k-1) ; i++ )
                                {
                                    for ( int j = i+1 ; j < k ; j++ )
                                    {
                                        if ( ( i < gamma_c ) && ( j < gamma_c ) )
                                            OM3[i][j][k] = OM2[i][j];
                                        if ( ( i < gamma_c ) && ( j >= gamma_c ) )
                                            OM3[i][j][k] = OM1[i];
                                        if ( ( i >= gamma_c ) && ( j >= gamma_c ) )
                                            OM3[i][j][k] = kappa;
                                    }
                                }
                            }
                            
                            for ( int j = gamma_c ; j < gamma ; j++ )
                            {
                                for ( int i = 0 ; i < j ; i++ )
                                {
                                    if ( i < gamma_c )
                                        OM2[i][j] = OM1[i];
                                    if ( i >= gamma_c )
                                        OM2[i][j] = kappa;
                                }
                            }
                            
                            for ( int i = gamma_c ; i < gamma ; i++ )
                            {
                                OM1[i] = kappa;
                            }
                            
                            int F = cycle_6_enumeration_protograph_SC ( kappa , gamma , L , m , OM3 , OM2 , OM1  );
                            if ( ( first_trial ) || ( F < F_min ) )
                            {
                                F_min = F;
                                tdeg3_opt[0] = OM3[0][1][2];
                                tdeg3_opt[1] = OM3[0][1][3];
                                tdeg3_opt[2] = OM3[0][2][3];
                                tdeg3_opt[3] = OM3[1][2][3];
                                tdeg2_opt[0] = OM2[0][1];
                                tdeg2_opt[1] = OM2[0][2];
                                tdeg2_opt[2] = OM2[0][3];
                                tdeg2_opt[3] = OM2[1][2];
                                tdeg2_opt[4] = OM2[1][3];
                                tdeg2_opt[5] = OM2[2][3];
                                tdeg1_opt[0] = OM1[0];
                                tdeg1_opt[1] = OM1[1];
                                tdeg1_opt[2] = OM1[2];
                                tdeg1_opt[3] = OM1[3];
                                
                                first_trial = false;
                            }
                        //}
                    }
                    }
                    }
                    }
                    }
                    }
                    }
                    }
                }
                }
                }
                }
            }
            }
        }
        cout << "Minimum Number of Cycles 6: " << F_min << endl;
        
        cout << "Optimum Set of Overlapping Parameters: " << endl;
        for ( int i = 0 ; i < 4 ; i++ ) cout << tdeg3_opt[i] << "   "; cout << endl;
        for ( int i = 0 ; i < 6 ; i++ ) cout << tdeg2_opt[i] << "   "; cout << endl;
        for ( int i = 0 ; i < 4 ; i++ ) cout << tdeg1_opt[i] << "   "; cout << endl;

    }
    
    if ( ( gamma_c == 3 ) && ( m == 2 ) )
    {
        tdeg3_opt = new int [8];
        tdeg2_opt = new int [12];
        tdeg1_opt = new int [6];
        
        for ( int t0 = 0 ; t0 <= kappa ; t0++ )
        {
            for ( int t01 = 0 ; t01 <= t0 ; t01++ )
            {
            for ( int t1 = t01 ; t1 <= (kappa-t0+t01) ; t1++ )
            {
                for ( int t012 = 0 ; t012 <= t01 ; t012++ )
                {
                for ( int t02 = t012 ; t02 <= (t0-t01+t012) ; t02++ )
                {
                for ( int t12 = t012 ; t12 <= (t1-t01+t012) ; t12++ )
                {
                for ( int t2 = (t02+t12-t012) ; t2 <= (kappa-t0-t1+t01+t02+t12-t012) ; t2++ )
                {
                    for ( int t123 = 0 ; t123 <= (t12-t012) ; t123++ )
                    {
                    for ( int t13 = t123 ; t13 <= (t1-t01-t12+t012+t123) ; t13++ )
                    {
                    for ( int t23 = t123 ; t23 <= (t2-t02-t12+t012+t123) ; t23++ )
                    {
                    for ( int t3 = (t13+t23-t123) ; t3 <= (kappa-t0-t1-t2+t01+t02+t12+t13+t23-t012-t123) ; t3++ )
                    {
                        for ( int t024 = 0 ; t024 <= (t02-t012) ; t024++ )
                        {
                        for ( int t234 = 0 ; t234 <= (t23-t123) ; t234++ )
                        {
                        for ( int t04 = t024 ; t04 <= (t0-t01-t02+t012+t024) ; t04++ )
                        {
                        for ( int t24 = (t024+t234) ; t24 <= (t2-t02-t12-t23+t012+t123+t024+t234) ; t24++ )
                        {
                        for ( int t34 = t234 ; t34 <= (t3-t13-t23+t123+t234) ; t34++ )
                        {
                        for ( int t4 = (t04+t24+t34-t024-t234) ; t4 <= (kappa-t0-t1-t2-t3+t01+t02+t04+t12+t13+t23+t24+t34-t012-t123-t024-t234) ; t4++ )
                        {
                            for ( int t015 = 0 ; t015 <= (t01-t012) ; t015++ )
                            {
                            for ( int t045 = 0 ; t045 <= (t04-t024) ; t045++ )
                            {
                            for ( int t135 = 0 ; t135 <= (t13-t123) ; t135++ )
                            {
                            for ( int t345 = 0 ; t345 <= (t34-t234) ; t345++ )
                            {
                            for ( int t05 = (t015+t045) ; t05 <= (t0-t01-t02-t04+t012+t024+t015+t045) ; t05++ )
                            {
                            for ( int t15 = (t015+t135) ; t15 <= (t1-t01-t12-t13+t012+t123+t015+t135) ; t15++ )
                            {
                            for ( int t35 = (t135+t345) ; t35 <= (t3-t13-t23-t34+t123+t234+t135+t345) ; t35++ )
                            {
                            for ( int t45 = (t045+t345) ; t45 <= (t4-t04-t24-t34+t024+t234+t045+t345) ; t45++ )
                            {
                            for ( int t5 = (t05+t15+t35+t45-t015-t045-t135-t345) ; t5 <= (kappa-t0-t1-t2-t3-t4+t01+t02+t04+t05+t12+t13+t15+t23+t24+t34+t35+t45-t012-t123-t024-t234-t015-t045-t135-t345) ; t5++ )
                            {
                                OM3[0][1][2] = t012;
                                OM3[0][1][5+gamma_l] = t015;
                                OM3[0][4+gamma_l][2] = t024;
                                OM3[0][4+gamma_l][5+gamma_l] = t045;
                                OM3[3+gamma_l][1][2] = t123;
                                OM3[3+gamma_l][1][5+gamma_l] = t135;
                                OM3[3+gamma_l][4+gamma_l][2] = t234;
                                OM3[3+gamma_l][4+gamma_l][5+gamma_l] = t345;
                                
                                OM2[0][1] = t01;
                                OM2[0][2] = t02;
                                OM2[0][4+gamma_l] = t04;
                                OM2[0][5+gamma_l] = t05;
                                OM2[1][2] = t12;
                                OM2[3+gamma_l][1] = t13;
                                OM2[1][5+gamma_l] = t15;
                                OM2[3+gamma_l][2] = t23;
                                OM2[4+gamma_l][2] = t24;
                                OM2[3+gamma_l][4+gamma_l] = t34;
                                OM2[3+gamma_l][5+gamma_l] = t35;
                                OM2[4+gamma_l][5+gamma_l] = t45;
                                
                                OM1[0] = t0;
                                OM1[1] = t1;
                                OM1[2] = t2;
                                OM1[3+gamma_l] = t3;
                                OM1[4+gamma_l] = t4;
                                OM1[5+gamma_l] = t5;
                                
                                for ( int k = gamma_c ; k < gamma ; k++ )
                                {
                                    for ( int i = 0 ; i < (k-1) ; i++ )
                                    {
                                        for ( int j = i+1 ; j < k ; j++ )
                                        {
                                            if ( ( i < gamma_c ) && ( j < gamma_c ) )
                                            {
                                                OM3[i][j][k] = OM2[i][j];
                                                OM3[i+gamma][j][k] = OM2[i+gamma][j];
                                                OM3[i][j+gamma][k] = OM2[i][j+gamma];
                                                OM3[i+gamma][j+gamma][k] = OM2[i+gamma][j+gamma];
                                            }
                                            if ( ( i < gamma_c ) && ( j >= gamma_c ) )
                                            {
                                                OM3[i][j][k] = OM1[i];
                                                OM3[i+gamma][j][k] = OM1[i+gamma];
                                            }
                                            if ( ( i >= gamma_c ) && ( j >= gamma_c ) )
                                            {
                                                OM3[i][j][k] = kappa;
                                            }
                                        }
                                    }
                                }
                                
                                for ( int j = gamma_c ; j < gamma ; j++ )
                                {
                                    for ( int i = 0 ; i < j ; i++ )
                                    {
                                        if ( i < gamma_c )
                                        {
                                            OM2[i][j] = OM1[i];
                                            OM2[i+gamma][j] = OM1[i+gamma];
                                        }
                                        if ( i >= gamma_c )
                                        {
                                            OM2[i][j] = kappa;
                                        }
                                    }
                                }
                                
                                for ( int i = gamma_c ; i < gamma ; i++ )
                                {
                                    OM1[i] = kappa;
                                }
                                
                                
                                int F = cycle_6_enumeration_protograph_SC ( kappa , gamma , L , m , OM3 , OM2 , OM1  );
                                if ( ( first_trial ) || ( F < F_min ) )
                                {
                                    F_min = F;
                                    tdeg3_opt[0] = t012;
                                    tdeg3_opt[1] = t015;
                                    tdeg3_opt[2] = t024;
                                    tdeg3_opt[3] = t045;
                                    tdeg3_opt[4] = t123;
                                    tdeg3_opt[5] = t135;
                                    tdeg3_opt[6] = t234;
                                    tdeg3_opt[7] = t345;
                                    tdeg2_opt[0] = t01;
                                    tdeg2_opt[1] = t02;
                                    tdeg2_opt[2] = t04;
                                    tdeg2_opt[3] = t05;
                                    tdeg2_opt[4] = t12;
                                    tdeg2_opt[5] = t13;
                                    tdeg2_opt[6] = t15;
                                    tdeg2_opt[7] = t23;
                                    tdeg2_opt[8] = t24;
                                    tdeg2_opt[9] = t34;
                                    tdeg2_opt[10] = t35;
                                    tdeg2_opt[11] = t45;
                                    tdeg1_opt[0] = t0;
                                    tdeg1_opt[1] = t1;
                                    tdeg1_opt[2] = t2;
                                    tdeg1_opt[3] = t3;
                                    tdeg1_opt[4] = t4;
                                    tdeg1_opt[5] = t5;
                                    first_trial = false;
                                }
                            
                            }
                            }
                            }
                            }
                            }
                            }
                            }
                            }
                            }
                        }
                        }
                        }
                        }
                        }
                        }
                    }
                    }
                    }
                    }
                }
                }
                }
                }
            }
            }
        }
        cout << "Minimum Number of Cycles 6: " << F_min << endl;
        
        cout << "Optimum Set of Overlapping Parameters: " << endl;
        for ( int i = 0 ; i < 8 ; i++ ) cout << tdeg3_opt[i] << "   "; cout << endl;
        for ( int i = 0 ; i < 12 ; i++ ) cout << tdeg2_opt[i] << "   "; cout << endl;
        for ( int i = 0 ; i < 6 ; i++ ) cout << tdeg1_opt[i] << "   "; cout << endl;
    }
    
    
    gamma = gamma + gamma_a;
    int** PM = new int* [ gamma ];
    for ( int i = 0 ; i < gamma ; i++ )
    {
        PM[i] = new int [ kappa ];
        for ( int j = 0 ; j < kappa ; j++ ) PM[i][j] = 0;
    }
    
    Partitioning_Matrix_Construction ( kappa , gamma_c , gamma_l , m , tdeg3_opt , tdeg2_opt , tdeg1_opt , PM );
    
    cout << "Optimum Partitioning Matrix: " << endl;
    for ( int i = 0 ; i < gamma ; i++ )
    {
        for ( int j = 0 ; j < kappa ; j++ ) cout << PM[i][j] << "  "; cout << endl;
    }
    
    int** H_SC = new int* [ (L+m) * gamma ]; for ( int i = 0 ; i < ( (L+m) * gamma ) ; i++ ) H_SC[i] = new int [ L * kappa ];
    for ( int i = 0 ; i < ( (L+m) * gamma ) ; i++ ) for ( int j = 0 ; j < ( L * kappa ) ; j++ ) H_SC[i][j] = 0;
    
    
    for ( int r = 0 ; r < L ; r++ )
    {
        for ( int i = 0 ; i < gamma ; i++ )
        {
            for ( int j = 0 ; j < kappa ; j++ )
            {
                if ( PM[i][j] == 0 )
                {
                    H_SC[r*gamma+i][r*kappa+j] = 1;
                }
                else if ( PM[i][j] == 1 )
                {
                    H_SC[(r+1)*gamma+i][r*kappa+j] = 1;
                }
                else if ( PM[i][j] == 2 )
                {
                    H_SC[(r+2)*gamma+i][r*kappa+j] = 1;
                }
            }
        }
    }
    
    cout << "C6: " << cycle_6_enumeration_protograph_block ( H_SC , (L+m) * gamma , L * kappa ) << endl;
    cout << "C4: " << cycle_4_enumeration_protograph_block ( H_SC , (L+m) * gamma , L * kappa ) << endl;
    
    for ( int i = 0 ; i < ( (L+m) * gamma ) ; i++ ) delete[] H_SC[i]; delete[] H_SC;
    
    for ( int i = 0 ; i < gamma ; i++ ) delete[] PM[i];
    delete[] PM;

    delete[] tdeg3_opt; delete[] tdeg2_opt; delete[] tdeg1_opt;
    for ( int i = 0 ; i < ((m+1)*(gamma-gamma_a)) ; i++ )
    {
        for ( int j = 0 ; j < ((m+1)*(gamma-gamma_a)) ; j++ )
        {
            delete[] OM3[i][j];
        }
        delete[] OM3[i];
    }
    delete[] OM3;
    
    for ( int i = 0 ; i < ((m+1)*(gamma-gamma_a)) ; i++ )
    {
        delete[] OM2[i];
    }
    delete[] OM2;
    
    delete[] OM1;
}

// n choose k
int comb_choose ( int n , int k )
{
    int result = 1;
    if ( ( n < 0 ) || ( k < 0 ) || ( n < k ) ) { cout << "Error!" << endl;}
    if ( ( n == 0 ) || ( k == 0 ) || ( n == k ) ) return result;
    for ( int i = 1 ; i <= k ; i++ ) result *= ( n - i + 1 );
    for ( int i = 1 ; i <= k ; i++ ) result /= i;
    return result;
}

// Finding the number of cycle-6 in H (brute force)
int cycle_6_enumeration_protograph_block ( int** H , int m , int n )
{
    //ofstream myfile;
    //myfile.open ("name.txt");
    int counter = 0;
    for ( int j0 = 0 ; j0 < n ; j0++ )
    {
        for ( int i0 = 0 ; i0 < m ; i0++ )
        {
            if ( H[i0][j0] == 1 )
            {
                for ( int j1 = (j0+1) ; j1 < n ; j1++ )
                {
                    if ( H[i0][j1] == 1 )
                    {
                        for ( int i1 = 0 ; i1 < m ; i1++ )
                        {
                            if ( ( i1 != i0 ) && ( H[i1][j1] == 1 ) )
                            {
                                for ( int j2 = (j1+1) ; j2 < n ; j2++ )
                                {
                                    if ( H[i1][j2] == 1 )
                                    {
                                        for ( int i2 = 0 ; i2 < m ; i2++ )
                                        {
                                            if ( ( i2 != i0 ) && ( i2 != i1 ) && ( H[i2][j2] == 1 ) && ( H[i2][j0] == 1 ) )
                                            {
                                                //myfile << i0 << ' ' << j0 << ' ' << i1 << ' ' << j1 << ' ' << i2 << ' ' << j2 << endl;
                                                counter++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //myfile.close();
    return counter;
}

// Finding the number of cycle-4 in H (brute force)
int cycle_4_enumeration_protograph_block ( int** H , int m , int n )
{
    //ofstream myfile;
    //myfile.open ("name.txt");
    int counter = 0;
    for ( int j0 = 0 ; j0 < n ; j0++ )
    {
        for ( int i0 = 0 ; i0 < m ; i0++ )
        {
            if ( H[i0][j0] == 1 )
            {
                for ( int j1 = (j0+1) ; j1 < n ; j1++ )
                {
                    if ( H[i0][j1] == 1 )
                    {
                        for ( int i1 = 0 ; i1 < m ; i1++ )
                        {
                            if ( ( i1 != i0 ) && ( H[i1][j1] == 1 ) && ( H[i1][j0] == 1 ) )
                            {
                                //myfile << i0 << ' ' << j0 << ' ' << i1 << ' ' << j1 << endl;
                                counter++;
                            }
                        }
                    }
                }
            }
        }
    }
    //myfile.close();
    return counter;
}

// Finding the number of cycle-6 given a set of overlap parameters when all overlaps belong to one replica
int A ( int tabc , int tab , int tac , int tbc )
{
    int result = 0;
    result += ( tabc * max ( tabc - 1 , 0 ) * max ( tbc - 2 , 0 ) );
    result += ( tabc * ( tac - tabc ) * max ( tbc - 1 , 0 ) );
    result += ( ( tab - tabc ) * tabc * max ( tbc - 1 , 0 ) );
    result += ( ( tab - tabc ) * ( tac - tabc ) * tbc );
    return result;
}

// Finding the number of cycle-6 given a set of overlap parameters when two overlaps belong to one replica and one overlap belongs to the other one
int B ( int tabc , int tab , int tac , int tbc )
{
    int result = ( tabc * max ( tac - 1 , 0 ) * tbc ) + ( ( tab - tabc ) * tac * tbc );
    return result;
}

// Finding the number of cycle-6 given a set of overlap parameters when each overlap belongs to an individual replica
int C ( int tab , int tac , int tbc )
{
    int result = tab * tac * tbc;
    return result;
}

// Finding the number of cycle-6 in H_SC using overlapping parameters
int cycle_6_enumeration_protograph_SC ( int kappa , int gamma , int L , int m , int*** OM3 , int** OM2 , int* OM1 )
{
    int F = 0;
    if ( m == 1 )
    {
        // Completing degree 3 entries
        for ( int i = 0 ; i < ( gamma - 2 ) ; i++ )
        {
            for ( int j = (i+1) ; j < ( gamma - 1 ) ; j++ )
            {
                for ( int k = (j+1) ; k < ( gamma ) ; k++ )
                {
                    OM3[i][j][k+gamma] = OM2[i][j] - OM3[i][j][k];
                    OM3[i][j+gamma][k] = OM2[i][k] - OM3[i][j][k];
                    OM3[i+gamma][j][k] = OM2[j][k] - OM3[i][j][k];
                    OM3[i][j+gamma][k+gamma] = OM1[i] - OM2[i][j] - OM2[i][k] + OM3[i][j][k];
                    OM3[i+gamma][j][k+gamma] = OM1[j] - OM2[i][j] - OM2[j][k] + OM3[i][j][k];
                    OM3[i+gamma][j+gamma][k] = OM1[k] - OM2[i][k] - OM2[j][k] + OM3[i][j][k];
                    OM3[i+gamma][j+gamma][k+gamma] = kappa - OM1[i] - OM1[j] - OM1[k] + OM2[i][j] + OM2[i][k] + OM2[j][k] - OM3[i][j][k];
                }
            }
        }
        
        // Completing degree 2 entries
        for ( int i = 0 ; i < ( gamma - 1 ) ; i++ )
        {
            for ( int j = (i+1) ; j < ( gamma ) ; j++ )
            {
                OM2[i][j+gamma] = OM1[i] - OM2[i][j];
                OM2[i+gamma][j] = OM1[j] - OM2[i][j];
                OM2[i+gamma][j+gamma] = kappa - OM1[i] - OM1[j] + OM2[i][j];
            }
        }
        
        // Enumeration
        
        int Fs = 0;
        int Fd = 0;
        for ( int i = 0 ; i < ( gamma - 2 ) ; i++ )
        {
            for ( int j = (i+1) ; j < ( gamma - 1 ) ; j++ )
            {
                for ( int k = (j+1) ; k < ( gamma ) ; k++ )
                {
                    Fs += A(OM3[i][j][k],OM2[i][j],OM2[i][k],OM2[j][k]);
                    Fs += A(OM3[i+gamma][j][k],OM2[i+gamma][j],OM2[i+gamma][k],OM2[j][k]);
                    Fs += A(OM3[i][j+gamma][k],OM2[i][j+gamma],OM2[i][k],OM2[j+gamma][k]);
                    Fs += A(OM3[i][j][k+gamma],OM2[i][j],OM2[i][k+gamma],OM2[j][k+gamma]);
                    Fs += A(OM3[i+gamma][j+gamma][k],OM2[i+gamma][j+gamma],OM2[i+gamma][k],OM2[j+gamma][k]);
                    Fs += A(OM3[i+gamma][j][k+gamma],OM2[i+gamma][j],OM2[i+gamma][k+gamma],OM2[j][k+gamma]);
                    Fs += A(OM3[i][j+gamma][k+gamma],OM2[i][j+gamma],OM2[i][k+gamma],OM2[j+gamma][k+gamma]);
                    Fs += A(OM3[i+gamma][j+gamma][k+gamma],OM2[i+gamma][j+gamma],OM2[i+gamma][k+gamma],OM2[j+gamma][k+gamma]);
                    
                    Fd += B ( OM3[i][j][k] , OM2[i][j],OM2[i][k] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i][j][k] , OM2[i][j],OM2[j][k] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i][j][k] , OM2[i][k],OM2[j][k] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i][j+gamma][k+gamma] , OM2[i][j+gamma],OM2[i][k+gamma] , OM2[j][k] );
                    Fd += B ( OM3[i+gamma][j][k+gamma] , OM2[i+gamma][j],OM2[j][k+gamma] , OM2[i][k] );
                    Fd += B ( OM3[i+gamma][j+gamma][k] , OM2[i+gamma][k],OM2[j+gamma][k] , OM2[i][j] );
                    Fd += B ( OM3[i+gamma][j][k] , OM2[i+gamma][j],OM2[i+gamma][k] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i][j+gamma][k] , OM2[i][j+gamma],OM2[j+gamma][k] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i][j][k+gamma] , OM2[i][k+gamma],OM2[j][k+gamma] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][j+gamma],OM2[i+gamma][k+gamma] , OM2[j][k] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][j+gamma],OM2[j+gamma][k+gamma] , OM2[i][k] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][k+gamma],OM2[j+gamma][k+gamma] , OM2[i][j] );
                }
            }
        }
        
        F = ( L * Fs ) + ( ( L - 1 ) * Fd );
    }
    
    if ( m == 2 )
    {
        // Completing degree 3 entries
        for ( int i = 0 ; i < ( gamma - 2 ) ; i++ )
        {
            for ( int j = (i+1) ; j < ( gamma - 1 ) ; j++ )
            {
                for ( int k = (j+1) ; k < ( gamma ) ; k++ )
                {
                    OM3[i][j][k+(2*gamma)] = OM2[i][j] - OM3[i][j][k] - OM3[i][j][k+gamma];
                    OM3[i][j+gamma][k+(2*gamma)] = OM2[i][j+gamma] - OM3[i][j+gamma][k] - OM3[i][j+gamma][k+gamma];
                    OM3[i+gamma][j][k+(2*gamma)] = OM2[i+gamma][j] - OM3[i+gamma][j][k] - OM3[i+gamma][j][k+gamma];
                    OM3[i+gamma][j+gamma][k+(2*gamma)] = OM2[i+gamma][j+gamma] - OM3[i+gamma][j+gamma][k] - OM3[i+gamma][j+gamma][k+gamma];
                    OM3[i][j+(2*gamma)][k] = OM2[i][k] - OM3[i][j][k] - OM3[i][j+gamma][k];
                    OM3[i][j+(2*gamma)][k+gamma] = OM2[i][k+gamma] - OM3[i][j][k+gamma] - OM3[i][j+gamma][k+gamma];
                    OM3[i+gamma][j+(2*gamma)][k] = OM2[i+gamma][k] - OM3[i+gamma][j][k] - OM3[i+gamma][j+gamma][k];
                    OM3[i+gamma][j+(2*gamma)][k+gamma] = OM2[i+gamma][k+gamma] - OM3[i+gamma][j][k+gamma] - OM3[i+gamma][j+gamma][k+gamma];
                    OM3[i+(2*gamma)][j][k] = OM2[j][k] - OM3[i][j][k] - OM3[i+gamma][j][k];
                    OM3[i+(2*gamma)][j+gamma][k] = OM2[j+gamma][k] - OM3[i][j+gamma][k] - OM3[i+gamma][j+gamma][k];
                    OM3[i+(2*gamma)][j][k+gamma] = OM2[j][k+gamma] - OM3[i][j][k+gamma] - OM3[i+gamma][j][k+gamma];
                    OM3[i+(2*gamma)][j+gamma][k+gamma] = OM2[j+gamma][k+gamma] - OM3[i][j+gamma][k+gamma] - OM3[i+gamma][j+gamma][k+gamma];
                    OM3[i][j+(2*gamma)][k+(2*gamma)] = OM1[i] - OM2[i][j] - OM2[i][k] - OM2[i][j+gamma] - OM2[i][k+gamma] + OM3[i][j][k] + OM3[i][j][k+gamma] + OM3[i][j+gamma][k] + OM3[i][j+gamma][k+gamma];
                    OM3[i+gamma][j+(2*gamma)][k+(2*gamma)] = OM1[i+gamma] - OM2[i+gamma][j] - OM2[i+gamma][k] - OM2[i+gamma][j+gamma] - OM2[i+gamma][k+gamma] + OM3[i+gamma][j][k] + OM3[i+gamma][j][k+gamma] + OM3[i+gamma][j+gamma][k] + OM3[i+gamma][j+gamma][k+gamma];
                    OM3[i+(2*gamma)][j][k+(2*gamma)] = OM1[j] - OM2[i][j] - OM2[j][k] - OM2[i+gamma][j] - OM2[j][k+gamma] + OM3[i][j][k] + OM3[i][j][k+gamma] + OM3[i+gamma][j][k] + OM3[i+gamma][j][k+gamma];
                    OM3[i+(2*gamma)][j+gamma][k+(2*gamma)] = OM1[j+gamma] - OM2[i][j+gamma] - OM2[j+gamma][k] - OM2[i+gamma][j+gamma] - OM2[j+gamma][k+gamma] + OM3[i][j+gamma][k] + OM3[i][j+gamma][k+gamma] + OM3[i+gamma][j+gamma][k] + OM3[i+gamma][j+gamma][k+gamma];
                    OM3[i+(2*gamma)][j+(2*gamma)][k] = OM1[k] - OM2[i][k] - OM2[j][k] - OM2[i+gamma][k] - OM2[j+gamma][k] + OM3[i][j][k] + OM3[i+gamma][j][k] + OM3[i][j+gamma][k] + OM3[i+gamma][j+gamma][k];
                    OM3[i+(2*gamma)][j+(2*gamma)][k+gamma] = OM1[k+gamma] - OM2[i][k+gamma] - OM2[j][k+gamma] - OM2[i+gamma][k+gamma] - OM2[j+gamma][k+gamma] + OM3[i][j][k+gamma] + OM3[i+gamma][j][k+gamma] + OM3[i][j+gamma][k+gamma] + OM3[i+gamma][j+gamma][k+gamma];
                    
                    OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] = kappa - OM1[i] - OM1[j] - OM1[k] - OM1[i+gamma] - OM1[j+gamma] - OM1[k+gamma] + OM2[i][j] + OM2[i+gamma][j] + OM2[i][j+gamma] + OM2[i+gamma][j+gamma] + OM2[i][k] + OM2[i+gamma][k] + OM2[i][k+gamma] + OM2[i+gamma][k+gamma] + OM2[j][k] + OM2[j+gamma][k] + OM2[j][k+gamma] + OM2[j+gamma][k+gamma] - OM3[i][j][k] - OM3[i+gamma][j][k] - OM3[i][j+gamma][k] - OM3[i][j][k+gamma] - OM3[i+gamma][j+gamma][k] - OM3[i+gamma][j][k+gamma] - OM3[i][j+gamma][k+gamma] - OM3[i+gamma][j+gamma][k+gamma];
                }
            }
        }
        
        // Completing degree 2 entries
        for ( int i = 0 ; i < ( gamma - 1 ) ; i++ )
        {
            for ( int j = (i+1) ; j < ( gamma ) ; j++ )
            {
                OM2[i][j+(2*gamma)] = OM1[i] - OM2[i][j] - OM2[i][j+gamma];
                OM2[i+gamma][j+(2*gamma)] = OM1[i+gamma] - OM2[i+gamma][j] - OM2[i+gamma][j+gamma];
                OM2[i+(2*gamma)][j] = OM1[j] - OM2[i][j] - OM2[i+gamma][j];
                OM2[i+(2*gamma)][j+gamma] = OM1[j+gamma] - OM2[i][j+gamma] - OM2[i+gamma][j+gamma];
                OM2[i+(2*gamma)][j+(2*gamma)] = kappa - OM1[i] - OM1[j] - OM1[i+gamma] - OM1[j+gamma] + OM2[i][j] + OM2[i+gamma][j] + OM2[i][j+gamma] + OM2[i+gamma][j+gamma];
            }
        }
        
        // Enumeration
        
        int Fs = 0;
        int Fd = 0;
        int Ft = 0;
        
        for ( int i = 0 ; i < ( gamma - 2 ) ; i++ )
        {
            for ( int j = (i+1) ; j < ( gamma - 1 ) ; j++ )
            {
                for ( int k = (j+1) ; k < ( gamma ) ; k++ )
                {
                    Fs += A ( OM3[i][j][k] , OM2[i][j] , OM2[i][k] , OM2[j][k] );
                    Fs += A ( OM3[i][j][k+gamma] , OM2[i][j] , OM2[i][k+gamma] , OM2[j][k+gamma] );
                    Fs += A ( OM3[i][j][k+(2*gamma)] , OM2[i][j] , OM2[i][k+(2*gamma)] , OM2[j][k+(2*gamma)] );
                    Fs += A ( OM3[i][j+gamma][k] , OM2[i][j+gamma] , OM2[i][k] , OM2[j+gamma][k] );
                    Fs += A ( OM3[i][j+gamma][k+gamma] , OM2[i][j+gamma] , OM2[i][k+gamma] , OM2[j+gamma][k+gamma] );
                    Fs += A ( OM3[i][j+gamma][k+(2*gamma)] , OM2[i][j+gamma] , OM2[i][k+(2*gamma)] , OM2[j+gamma][k+(2*gamma)] );
                    Fs += A ( OM3[i][j+(2*gamma)][k] , OM2[i][j+(2*gamma)] , OM2[i][k] , OM2[j+(2*gamma)][k] );
                    Fs += A ( OM3[i][j+(2*gamma)][k+gamma] , OM2[i][j+(2*gamma)] , OM2[i][k+gamma] , OM2[j+(2*gamma)][k+gamma] );
                    Fs += A ( OM3[i][j+(2*gamma)][k+(2*gamma)] , OM2[i][j+(2*gamma)] , OM2[i][k+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Fs += A ( OM3[i+gamma][j][k] , OM2[i+gamma][j] , OM2[i+gamma][k] , OM2[j][k] );
                    Fs += A ( OM3[i+gamma][j][k+gamma] , OM2[i+gamma][j] , OM2[i+gamma][k+gamma] , OM2[j][k+gamma] );
                    Fs += A ( OM3[i+gamma][j][k+(2*gamma)] , OM2[i+gamma][j] , OM2[i+gamma][k+(2*gamma)] , OM2[j][k+(2*gamma)] );
                    Fs += A ( OM3[i+gamma][j+gamma][k] , OM2[i+gamma][j+gamma] , OM2[i+gamma][k] , OM2[j+gamma][k] );
                    Fs += A ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][j+gamma] , OM2[i+gamma][k+gamma] , OM2[j+gamma][k+gamma] );
                    Fs += A ( OM3[i+gamma][j+gamma][k+(2*gamma)] , OM2[i+gamma][j+gamma] , OM2[i+gamma][k+(2*gamma)] , OM2[j+gamma][k+(2*gamma)] );
                    Fs += A ( OM3[i+gamma][j+(2*gamma)][k] , OM2[i+gamma][j+(2*gamma)] , OM2[i+gamma][k] , OM2[j+(2*gamma)][k] );
                    Fs += A ( OM3[i+gamma][j+(2*gamma)][k+gamma] , OM2[i+gamma][j+(2*gamma)] , OM2[i+gamma][k+gamma] , OM2[j+(2*gamma)][k+gamma] );
                    Fs += A ( OM3[i+gamma][j+(2*gamma)][k+(2*gamma)] , OM2[i+gamma][j+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Fs += A ( OM3[i+(2*gamma)][j][k] , OM2[i+(2*gamma)][j] , OM2[i+(2*gamma)][k] , OM2[j][k] );
                    Fs += A ( OM3[i+(2*gamma)][j][k+gamma] , OM2[i+(2*gamma)][j] , OM2[i+(2*gamma)][k+gamma] , OM2[j][k+gamma] );
                    Fs += A ( OM3[i+(2*gamma)][j][k+(2*gamma)] , OM2[i+(2*gamma)][j] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j][k+(2*gamma)] );
                    Fs += A ( OM3[i+(2*gamma)][j+gamma][k] , OM2[i+(2*gamma)][j+gamma] , OM2[i+(2*gamma)][k] , OM2[j+gamma][k] );
                    Fs += A ( OM3[i+(2*gamma)][j+gamma][k+gamma] , OM2[i+(2*gamma)][j+gamma] , OM2[i+(2*gamma)][k+gamma] , OM2[j+gamma][k+gamma]);
                    Fs += A ( OM3[i+(2*gamma)][j+gamma][k+(2*gamma)] , OM2[i+(2*gamma)][j+gamma] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+gamma][k+(2*gamma)] );
                    Fs += A ( OM3[i+(2*gamma)][j+(2*gamma)][k] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+(2*gamma)][k] , OM2[j+(2*gamma)][k] );
                    Fs += A ( OM3[i+(2*gamma)][j+(2*gamma)][k+gamma] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+(2*gamma)][k+gamma] , OM2[j+(2*gamma)][k+gamma] );
                    Fs += A ( OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    
                    
                    
                    Fd += B ( OM3[i][j+gamma][k+gamma] , OM2[i][j+gamma] , OM2[i][k+gamma] , OM2[j][k] );
                    Fd += B ( OM3[i][j+gamma][k+(2*gamma)] , OM2[i][j+gamma] , OM2[i][k+(2*gamma)] , OM2[j][k+gamma] );
                    Fd += B ( OM3[i][j+(2*gamma)][k+gamma] , OM2[i][j+(2*gamma)] , OM2[i][k+gamma] , OM2[j+gamma][k] );
                    Fd += B ( OM3[i][j+(2*gamma)][k+(2*gamma)] , OM2[i][j+(2*gamma)] , OM2[i][k+(2*gamma)] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][j+gamma] , OM2[i+gamma][k+gamma] , OM2[j][k] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+(2*gamma)] , OM2[i+gamma][j+gamma] , OM2[i+gamma][k+(2*gamma)] , OM2[j][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k+gamma] , OM2[i+gamma][j+(2*gamma)] , OM2[i+gamma][k+gamma] , OM2[j+gamma][k] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k+(2*gamma)] , OM2[i+gamma][j+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k+gamma] , OM2[i+(2*gamma)][j+gamma] , OM2[i+(2*gamma)][k+gamma] , OM2[j][k] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k+(2*gamma)] , OM2[i+(2*gamma)][j+gamma] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j][k+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+gamma] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+(2*gamma)][k+gamma] , OM2[j+gamma][k] );
                    Fd += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i+gamma][j][k+gamma] , OM2[i+gamma][j] , OM2[j][k+gamma] , OM2[i][k] );
                    Fd += B ( OM3[i+gamma][j][k+(2*gamma)] , OM2[i+gamma][j] , OM2[j][k+(2*gamma)] , OM2[i][k+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j][k+gamma] , OM2[i+(2*gamma)][j] , OM2[j][k+gamma] , OM2[i+gamma][k] );
                    Fd += B ( OM3[i+(2*gamma)][j][k+(2*gamma)] , OM2[i+(2*gamma)][j] , OM2[j][k+(2*gamma)] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][j+gamma] , OM2[j+gamma][k+gamma] , OM2[i][k] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+(2*gamma)] , OM2[i+gamma][j+gamma] , OM2[j+gamma][k+(2*gamma)] , OM2[i][k+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k+gamma] , OM2[i+(2*gamma)][j+gamma] , OM2[j+gamma][k+gamma] , OM2[i+gamma][k] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k+(2*gamma)] , OM2[i+(2*gamma)][j+gamma] , OM2[j+gamma][k+(2*gamma)] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k+gamma] , OM2[i+gamma][j+(2*gamma)] , OM2[j+(2*gamma)][k+gamma] , OM2[i][k] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k+(2*gamma)] , OM2[i+gamma][j+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] , OM2[i][k+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+gamma] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[j+(2*gamma)][k+gamma] , OM2[i+gamma][k] );
                    Fd += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k] , OM2[i+gamma][k] , OM2[j+gamma][k] , OM2[i][j] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k] , OM2[i+gamma][k] , OM2[j+(2*gamma)][k] , OM2[i][j+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k] , OM2[i+(2*gamma)][k] , OM2[j+gamma][k] , OM2[i+gamma][j] );
                    Fd += B ( OM3[i+(2*gamma)][j+(2*gamma)][k] , OM2[i+(2*gamma)][k] , OM2[j+(2*gamma)][k] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][k+gamma] , OM2[j+gamma][k+gamma] , OM2[i][j] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k+gamma] , OM2[i+gamma][k+gamma] , OM2[j+(2*gamma)][k+gamma] , OM2[i][j+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k+gamma] , OM2[i+(2*gamma)][k+gamma] , OM2[j+gamma][k+gamma] , OM2[i+gamma][j] );
                    Fd += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+gamma] , OM2[i+(2*gamma)][k+gamma] , OM2[j+(2*gamma)][k+gamma] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j+gamma][k+(2*gamma)] , OM2[i][j] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] , OM2[i][j+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k+(2*gamma)] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+gamma][k+(2*gamma)] , OM2[i+gamma][j] );
                    Fd += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i][j][k] , OM2[i][j] , OM2[i][k] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i][j][k+gamma] , OM2[i][j] , OM2[i][k+gamma] , OM2[j+gamma][k+(2*gamma)] );
                    Fd += B ( OM3[i][j+gamma][k] , OM2[i][j+gamma] , OM2[i][k] , OM2[j+(2*gamma)][k+gamma] );
                    Fd += B ( OM3[i][j+gamma][k+gamma] , OM2[i][j+gamma] , OM2[i][k+gamma] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j][k] , OM2[i+gamma][j] , OM2[i+gamma][k] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i+gamma][j][k+gamma] , OM2[i+gamma][j] , OM2[i+gamma][k+gamma] , OM2[j+gamma][k+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j+gamma][k] , OM2[i+gamma][j+gamma] , OM2[i+gamma][k] , OM2[j+(2*gamma)][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][j+gamma] , OM2[i+gamma][k+gamma] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Fd += B ( OM3[i+(2*gamma)][j][k] , OM2[i+(2*gamma)][j] , OM2[i+(2*gamma)][k] , OM2[j+gamma][k+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j][k+gamma] , OM2[i+(2*gamma)][j] , OM2[i+(2*gamma)][k+gamma] , OM2[j+gamma][k+(2*gamma)] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k] , OM2[i+(2*gamma)][j+gamma] , OM2[i+(2*gamma)][k] , OM2[j+(2*gamma)][k+gamma] );
                    Fd += B ( OM3[i+(2*gamma)][j+gamma][k+gamma] , OM2[i+(2*gamma)][j+gamma] , OM2[i+(2*gamma)][k+gamma] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Fd += B ( OM3[i][j][k] , OM2[i][j] , OM2[j][k] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i][j][k+gamma] , OM2[i][j] , OM2[j][k+gamma] , OM2[i+gamma][k+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j][k] , OM2[i+gamma][j] , OM2[j][k] , OM2[i+(2*gamma)][k+gamma] );
                    Fd += B ( OM3[i+gamma][j][k+gamma] , OM2[i+gamma][j] , OM2[j][k+gamma] , OM2[i+(2*gamma)][k+(2*gamma)] );
                    Fd += B ( OM3[i][j+gamma][k] , OM2[i][j+gamma] , OM2[j+gamma][k] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i][j+gamma][k+gamma] , OM2[i][j+gamma] , OM2[j+gamma][k+gamma] , OM2[i+gamma][k+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j+gamma][k] , OM2[i+gamma][j+gamma] , OM2[j+gamma][k] , OM2[i+(2*gamma)][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][j+gamma] , OM2[j+gamma][k+gamma] , OM2[i+(2*gamma)][k+(2*gamma)] );
                    Fd += B ( OM3[i][j+(2*gamma)][k] , OM2[i][j+(2*gamma)] , OM2[j+(2*gamma)][k] , OM2[i+gamma][k+gamma] );
                    Fd += B ( OM3[i][j+(2*gamma)][k+gamma] , OM2[i][j+(2*gamma)] , OM2[j+(2*gamma)][k+gamma] , OM2[i+gamma][k+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k] , OM2[i+gamma][j+(2*gamma)] , OM2[j+(2*gamma)][k] , OM2[i+(2*gamma)][k+gamma] );
                    Fd += B ( OM3[i+gamma][j+(2*gamma)][k+gamma] , OM2[i+gamma][j+(2*gamma)] , OM2[j+(2*gamma)][k+gamma] , OM2[i+(2*gamma)][k+(2*gamma)] );
                    Fd += B ( OM3[i][j][k] , OM2[i][k] , OM2[j][k] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i][j+gamma][k] , OM2[i][k] , OM2[j+gamma][k] , OM2[i+gamma][j+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j][k] , OM2[i+gamma][k] , OM2[j][k] , OM2[i+(2*gamma)][j+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k] , OM2[i+gamma][k] , OM2[j+gamma][k] , OM2[i+(2*gamma)][j+(2*gamma)] );
                    Fd += B ( OM3[i][j][k+gamma] , OM2[i][k+gamma] , OM2[j][k+gamma] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i][j+gamma][k+gamma] , OM2[i][k+gamma] , OM2[j+gamma][k+gamma] , OM2[i+gamma][j+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j][k+gamma] , OM2[i+gamma][k+gamma] , OM2[j][k+gamma] , OM2[i+(2*gamma)][j+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+gamma] , OM2[i+gamma][k+gamma] , OM2[j+gamma][k+gamma] , OM2[i+(2*gamma)][j+(2*gamma)] );
                    Fd += B ( OM3[i][j][k+(2*gamma)] , OM2[i][k+(2*gamma)] , OM2[j][k+(2*gamma)] , OM2[i+gamma][j+gamma] );
                    Fd += B ( OM3[i][j+gamma][k+(2*gamma)] , OM2[i][k+(2*gamma)] , OM2[j+gamma][k+(2*gamma)] , OM2[i+gamma][j+(2*gamma)] );
                    Fd += B ( OM3[i+gamma][j][k+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j][k+(2*gamma)] , OM2[i+(2*gamma)][j+gamma] );
                    Fd += B ( OM3[i+gamma][j+gamma][k+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j+gamma][k+(2*gamma)] , OM2[i+(2*gamma)][j+(2*gamma)] );
                    
                    
                    Ft += B ( OM3[i][j+(2*gamma)][k+(2*gamma)] , OM2[i][j+(2*gamma)] , OM2[i][k+(2*gamma)] , OM2[j][k] );
                    Ft += B ( OM3[i+gamma][j+(2*gamma)][k+(2*gamma)] , OM2[i+gamma][j+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j][k] );
                    Ft += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j][k] );
                    Ft += B ( OM3[i+(2*gamma)][j][k+(2*gamma)] , OM2[i+(2*gamma)][j] , OM2[j][k+(2*gamma)] , OM2[i][k] );
                    Ft += B ( OM3[i+(2*gamma)][j+gamma][k+(2*gamma)] , OM2[i+(2*gamma)][j+gamma] , OM2[j+gamma][k+(2*gamma)] , OM2[i][k] );
                    Ft += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] , OM2[i+(2*gamma)][j+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] , OM2[i][k] );
                    Ft += B ( OM3[i+(2*gamma)][j+(2*gamma)][k] , OM2[i+(2*gamma)][k] , OM2[j+(2*gamma)][k] , OM2[i][j] );
                    Ft += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+gamma] , OM2[i+(2*gamma)][k+gamma] , OM2[j+(2*gamma)][k+gamma] , OM2[i][j] );
                    Ft += B ( OM3[i+(2*gamma)][j+(2*gamma)][k+(2*gamma)] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+(2*gamma)][k+(2*gamma)] , OM2[i][j] );
                    Ft += B ( OM3[i][j][k] , OM2[i][j] , OM2[i][k] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Ft += B ( OM3[i+gamma][j][k] , OM2[i+gamma][j] , OM2[i+gamma][k] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Ft += B ( OM3[i+(2*gamma)][j][k] , OM2[i+(2*gamma)][j] , OM2[i+(2*gamma)][k] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Ft += B ( OM3[i][j][k] , OM2[i][j] , OM2[j][k] , OM2[i+(2*gamma)][k+(2*gamma)] );
                    Ft += B ( OM3[i][j+gamma][k] , OM2[i][j+gamma] , OM2[j+gamma][k] , OM2[i+(2*gamma)][k+(2*gamma)] );
                    Ft += B ( OM3[i][j+(2*gamma)][k] , OM2[i][j+(2*gamma)] , OM2[j+(2*gamma)][k] , OM2[i+(2*gamma)][k+(2*gamma)] );
                    Ft += B ( OM3[i][j][k] , OM2[i][k] , OM2[j][k] , OM2[i+(2*gamma)][j+(2*gamma)] );
                    Ft += B ( OM3[i][j][k+gamma] , OM2[i][k+gamma] , OM2[j][k+gamma] , OM2[i+(2*gamma)][j+(2*gamma)] );
                    Ft += B ( OM3[i][j][k+(2*gamma)] , OM2[i][k+(2*gamma)] , OM2[j][k+(2*gamma)] , OM2[i+(2*gamma)][j+(2*gamma)] );
                    
                    
                    Ft += C ( OM2[i+gamma][j+(2*gamma)] , OM2[i][k+gamma] , OM2[j][k] );
                    Ft += C ( OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+gamma][k+gamma] , OM2[j][k] );
                    Ft += C ( OM2[i+gamma][j+(2*gamma)] , OM2[i][k+(2*gamma)] , OM2[j][k+gamma] );
                    Ft += C ( OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j][k+gamma] );
                    Ft += C ( OM2[i+(2*gamma)][j+gamma] , OM2[i][k] , OM2[j][k+gamma] );
                    Ft += C ( OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i][k] , OM2[j+gamma][k+gamma] );
                    Ft += C ( OM2[i+(2*gamma)][j+gamma] , OM2[i][k+gamma] , OM2[j][k+(2*gamma)] );
                    Ft += C ( OM2[i+(2*gamma)][j+(2*gamma)] , OM2[i][k+gamma] , OM2[j+gamma][k+(2*gamma)] );
                    Ft += C ( OM2[i][j+gamma] , OM2[i+gamma][k+(2*gamma)] , OM2[j][k] );
                    Ft += C ( OM2[i+gamma][j+gamma] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j][k] );
                    Ft += C ( OM2[i][j+(2*gamma)] , OM2[i+gamma][k+(2*gamma)] , OM2[j+gamma][k] );
                    Ft += C ( OM2[i+gamma][j+(2*gamma)] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+gamma][k] );
                    Ft += C ( OM2[i][j] , OM2[i+(2*gamma)][k+gamma] , OM2[j+gamma][k] );
                    Ft += C ( OM2[i][j] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+gamma][k+gamma] );
                    Ft += C ( OM2[i][j+gamma] , OM2[i+(2*gamma)][k+gamma] , OM2[j+(2*gamma)][k] );
                    Ft += C ( OM2[i][j+gamma] , OM2[i+(2*gamma)][k+(2*gamma)] , OM2[j+(2*gamma)][k+gamma] );
                    Ft += C ( OM2[i+gamma][j] , OM2[i][k] , OM2[j+gamma][k+(2*gamma)] );
                    Ft += C ( OM2[i+gamma][j+gamma] , OM2[i][k] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Ft += C ( OM2[i+(2*gamma)][j] , OM2[i+gamma][k] , OM2[j+gamma][k+(2*gamma)] );
                    Ft += C ( OM2[i+(2*gamma)][j+gamma] , OM2[i+gamma][k] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Ft += C ( OM2[i][j] , OM2[i+gamma][k] , OM2[j+(2*gamma)][k+gamma] );
                    Ft += C ( OM2[i][j] , OM2[i+gamma][k+gamma] , OM2[j+(2*gamma)][k+(2*gamma)] );
                    Ft += C ( OM2[i+gamma][j] , OM2[i+(2*gamma)][k] , OM2[j+(2*gamma)][k+gamma] );
                    Ft += C ( OM2[i+gamma][j] , OM2[i+(2*gamma)][k+gamma] , OM2[j+(2*gamma)][k+(2*gamma)] );
                }
            }
        }
        
        F = ( L * Fs ) + ( ( L - 1 ) * Fd ) + ( ( L - 2 ) * Ft );
    }
    
    return F;
}

// Finding a partitioning that achieve a particular set of overlapping parameters
// The output is the number of partitionig matrices that achieve the set of given overlapping parameters. But, PM is filled by one of the choices
int Partitioning_Matrix_Construction ( int kappa , int gamma_c , int gamma_l , int m , int* tdeg3 , int* tdeg2 , int* tdeg1 , int** PM )
{
    
    if ( ( gamma_c == 3 ) && ( m == 2 ) )
    {
        int t012 = tdeg3[0];
        int t015 = tdeg3[1];
        int t024 = tdeg3[2];
        int t045 = tdeg3[3];
        int t123 = tdeg3[4];
        int t135 = tdeg3[5];
        int t234 = tdeg3[6];
        int t345 = tdeg3[7];
        
        int t01 = tdeg2[0];
        int t02 = tdeg2[1];
        int t04 = tdeg2[2];
        int t05 = tdeg2[3];
        int t12 = tdeg2[4];
        int t13 = tdeg2[5];
        int t15 = tdeg2[6];
        int t23 = tdeg2[7];
        int t24 = tdeg2[8];
        int t34 = tdeg2[9];
        int t35 = tdeg2[10];
        int t45 = tdeg2[11];
        
        int t0 = tdeg1[0];
        int t1 = tdeg1[1];
        int t2 = tdeg1[2];
        int t3 = tdeg1[3];
        int t4 = tdeg1[4];
        int t5 = tdeg1[5];
        
        int index = 0;
        
        for ( int i = 0 ; i < gamma_c ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) PM[i][j] = 2;
        
        for ( int n0 = 0 ; n0 < t0  ; n0++ )
        {
            index = n0 * floor ( kappa / t0 );
            PM[0][index] = 0;
        }
        
        index = 0;
        for ( int n01 = 0 ; n01 < t01  ; n01++ )
        {
            while ( PM[0][index] != 0 ) index++;
            PM[1][index] = 0; index++;
        }
        
        index = 0;
        for ( int n1 = 0 ; n1 < (t1 - t01)  ; n1++ )
        {
            while ( PM[0][index] != 2 ) index++;
            PM[1][index] = 0; index++;
        }
        
        index = 0;
        for ( int n012 = 0 ; n012 < t012 ; n012++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 0 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n02 = 0 ; n02 < ( t02 - t012 ) ; n02++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 2 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n12 = 0 ; n12 < ( t12 - t012 ) ; n12++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 0 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n2 = 0 ; n2 < ( t2 - t02 - t12 + t012 ) ; n2++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 2 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n123 = 0 ; n123 < t123 ; n123++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 0 ) ) index++;
            PM[0][index] = 1; index++;
        }
        
        index = 0;
        for ( int n13 = 0 ; n13 < ( t13 - t123 ) ; n13++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 2 ) ) index++;
            PM[0][index] = 1; index++;
        }
        
        index = 0;
        for ( int n23 = 0 ; n23 < ( t23 - t123 ) ; n23++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 0 ) ) index++;
            PM[0][index] = 1; index++;
        }
        
        index = 0;
        for ( int n3 = 0 ; n3 < ( t3 - t13 - t23 + t123 ) ; n3++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 2 ) ) index++;
            PM[0][index] = 1; index++;
        }
        
        index = 0;
        for ( int n024 = 0 ; n024 < t024 ; n024++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 0 ) ) index++;
            PM[1][index] = 1; index++;
        }
        
        index = 0;
        for ( int n234 = 0 ; n234 < t234 ; n234++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 0 ) ) index++;
            PM[1][index] = 1; index++;
        }
        
        index = 0;
        for ( int n04 = 0 ; n04 < ( t04 - t024 ) ; n04++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 2 ) ) index++;
            PM[1][index] = 1; index++;
        }
        
        index = 0;
        for ( int n24 = 0 ; n24 < ( t24 - t024 - t234 ) ; n24++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 0 ) ) index++;
            PM[1][index] = 1; index++;
        }
        
        index = 0;
        for ( int n34 = 0 ; n34 < ( t34 - t234 ) ; n34++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 2 ) ) index++;
            PM[1][index] = 1; index++;
        }
        
        index = 0;
        for ( int n4 = 0 ; n4 < ( t4 - t04 - t24 - t34 + t024 + t234 ) ; n4++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 2 ) ) index++;
            PM[1][index] = 1; index++;
        }
        
        index = 0;
        for ( int n015 = 0 ; n015 < t015 ; n015++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n045 = 0 ; n045 < t045 ; n045++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 1 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n135 = 0 ; n135 < t135 ; n135++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n345 = 0 ; n345 < t345 ; n345++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 1 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n05 = 0 ; n05 < ( t05 - t015 - t045 ) ; n05++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n15 = 0 ; n15 < ( t15 - t015 - t135 ) ; n15++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n35 = 0 ; n35 < ( t35 - t135 - t345 ) ; n35++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n45 = 0 ; n45 < ( t45 - t045 - t345 ) ; n45++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 1 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
        
        index = 0;
        for ( int n5 = 0 ; n5 < ( t5 - t05 - t15 - t35 - t45 + t015 + t045 + t135 + t345 ) ; n5++ )
        {
            while ( ( PM[0][index] != 2 ) || ( PM[1][index] != 2 ) || ( PM[2][index] != 2 ) ) index++;
            PM[2][index] = 1; index++;
        }
    }
    
    if ( ( gamma_c == 4 ) && ( m == 1 ) )
    {
        
        int t0123 = 0;
        
        int t012 = tdeg3 [0];
        int t013 = tdeg3 [1];
        int t023 = tdeg3 [2];
        int t123 = tdeg3 [3];
        
        int t01 = tdeg2 [0];
        int t02 = tdeg2 [1];
        int t03 = tdeg2 [2];
        int t12 = tdeg2 [3];
        int t13 = tdeg2 [4];
        int t23 = tdeg2 [5];
        
        int t0 = tdeg1 [0];
        int t1 = tdeg1 [1];
        int t2 = tdeg1 [2];
        int t3 = tdeg1 [3];
        
        int index = 0;
        
        for ( int i = 0 ; i < gamma_c ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) PM[i][j] = 1;
        
        for ( int n0 = 0 ; n0 < t0  ; n0++ )
        {
            index = n0 * floor ( kappa / t0 );
            PM[0][index] = 0;
        }
        
        index = 0;
        for ( int n01 = 0 ; n01 < t01  ; n01++ )
        {
            while ( PM[0][index] != 0 ) index++;
            PM[1][index] = 0; index++;
        }
        
        index = 0;
        for ( int n1 = 0 ; n1 < (t1 - t01)  ; n1++ )
        {
            while ( PM[0][index] != 1 ) index++;
            PM[1][index] = 0; index++;
        }
        
        index = 0;
        for ( int n012 = 0 ; n012 < t012 ; n012++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 0 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n02 = 0 ; n02 < ( t02 - t012 ) ; n02++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 1 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n12 = 0 ; n12 < ( t12 - t012 ) ; n12++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 0 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n2 = 0 ; n2 < ( t2 - t02 - t12 + t012 ) ; n2++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 1 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n0123 = 0 ; n0123 < t0123 ; n0123++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 0 ) ) index++;
            PM[3][index] = 0; index++;
        }
        
        index = 0;
        for ( int n013 = 0 ; n013 < ( t013 - t0123 ) ; n013++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 1 ) ) index++;
            PM[3][index] = 0; index++;
        }
        
        index = 0;
        for ( int n023 = 0 ; n023 < ( t023 - t0123 ) ; n023++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 1 ) || ( PM[2][index] != 0 ) ) index++;
            PM[3][index] = 0; index++;
        }
        
        index = 0;
        for ( int n123 = 0 ; n123 < ( t123 - t0123 ) ; n123++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 0 ) ) index++;
            PM[3][index] = 0; index++;
        }
        
        index = 0;
        for ( int n03 = 0 ; n03 < ( t03 - t013 - t023 + t0123 ) ; n03++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 1 ) || ( PM[2][index] != 1 ) ) index++;
            PM[3][index] = 0; index++;
        }
        
        index = 0;
        for ( int n13 = 0 ; n13 < ( t13 - t013 - t123 + t0123 ) ; n13++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 0 ) || ( PM[2][index] != 1 ) ) index++;
            PM[3][index] = 0; index++;
        }
        
        index = 0;
        for ( int n23 = 0 ; n23 < ( t23 - t023 - t123 + t0123 ) ; n23++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 1 ) || ( PM[2][index] != 0 ) ) index++;
            PM[3][index] = 0; index++;
        }
        
        index = 0;
        for ( int n3 = 0 ; n3 < ( t3 - t03 - t13 - t23 + t013 + t023 + t123 - t0123 ) ; n3++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 1 ) || ( PM[2][index] != 1 ) ) index++;
            PM[3][index] = 0; index++;
        }
    }
    
    if ( ( gamma_c == 3 ) && ( m == 1 ) )
    {
        int t012 = tdeg3[0];
        
        int t01 = tdeg2[0];
        int t02 = tdeg2[1];
        int t12 = tdeg2[2];
        
        int t0 = tdeg1[0];
        int t1 = tdeg1[1];
        int t2 = tdeg1[2];
        
        int index = 0;
        
        for ( int i = 0 ; i < gamma_c ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) PM[i][j] = 1;
        
        for ( int n0 = 0 ; n0 < t0  ; n0++ )
        {
            index = n0 * floor ( kappa / t0 );
            PM[0][index] = 0;
        }
        
        index = 0;
        for ( int n01 = 0 ; n01 < t01  ; n01++ )
        {
            while ( PM[0][index] != 0 ) index++;
            PM[1][index] = 0; index++;
        }
        
        index = 0;
        for ( int n1 = 0 ; n1 < (t1 - t01)  ; n1++ )
        {
            while ( PM[0][index] != 1 ) index++;
            PM[1][index] = 0; index++;
        }
        
        index = 0;
        for ( int n012 = 0 ; n012 < t012 ; n012++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 0 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n02 = 0 ; n02 < ( t02 - t012 ) ; n02++ )
        {
            while ( ( PM[0][index] != 0 ) || ( PM[1][index] != 1 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n12 = 0 ; n12 < ( t12 - t012 ) ; n12++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 0 ) ) index++;
            PM[2][index] = 0; index++;
        }
        
        index = 0;
        for ( int n2 = 0 ; n2 < ( t2 - t02 - t12 + t012 ) ; n2++ )
        {
            while ( ( PM[0][index] != 1 ) || ( PM[1][index] != 1 ) ) index++;
            PM[2][index] = 0; index++;
        }
    }
    
    for ( int i = gamma_c ; i < (gamma_c+gamma_l) ; i++ )
    {
        for (int j = 0 ; j < kappa ; j++)
        {
            PM[i][j] = 0;
        }
    }
    
    return 1;
}

