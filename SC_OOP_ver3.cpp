/****************************************************************************************/
/* Target: Enumeration of objects for SC codes with optimum overlap (OO) partitioning   */
/* By: Homa Esfahanizadeh (hesfahanizadeh@ucla.edu)                                     */
/* University of California, Los Angeles (UCLA)                                         */
/* Update: 07/31/2017                                                                   */
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

int comb_choose ( int n , int k )
{
    int result = 1;
    if ( ( n < 0 ) || ( k < 0 ) || ( n < k ) ) { cout << "Error!" << endl; while(1); }
    if ( n < k ) { cout << "Error!" << endl; while(1); }
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

// Finding the number of cycle-6 given a set of overlap parameters when all overlaps belong to one section
// The involved rows indices in the cycle are a, b, and c
// tabc is the overlap over all three rows, and the others are defined over two rows
int A ( int tabc , int tab , int tac , int tbc )
{
    int result = 0;
    result += ( tabc * max ( tabc - 1 , 0 ) * max ( tbc - 2 , 0 ) );
    result += ( tabc * ( tac - tabc ) * max ( tbc - 1 , 0 ) );
    result += ( ( tab - tabc ) * tabc * max ( tbc - 1 , 0 ) );
    result += ( ( tab - tabc ) * ( tac - tabc ) * tbc );
    return result;
}

// Finding the number of cycle-6 given a set of overlap parameters when two overlaps belong to one section and one overlap belongs to the other one
int B ( int tabc , int tab , int tac , int tbc )
{
    int result = ( tabc * max ( tac - 1 , 0 ) * tbc ) + ( ( tab - tabc ) * tac * tbc );
    return result;
}

int C ( int tab , int tac , int tbc )
{
    int result = tab * tac * tbc;
    return result;
}

// Finding the number of cycle-6 in H_SC using overlapping parameters
int cycle_6_enumeration_protograph_SC ( int kappa , int gamma , int L , int m , int* tdeg3 , int* tdeg2 , int* tdeg1 )
{
    int F;
    if ( ( gamma == 3 ) && ( m == 2 ) )
    {
        // independebt overlap parameters

        int t012 = tdeg3 [0];
        int t015 = tdeg3 [1];
        int t024 = tdeg3 [2];
        int t045 = tdeg3 [3];
        int t123 = tdeg3 [4];
        int t135 = tdeg3 [5];
        int t234 = tdeg3 [6];
        int t345 = tdeg3 [7];

        int t01 = tdeg2 [0];
        int t02 = tdeg2 [1];
        int t04 = tdeg2 [2];
        int t05 = tdeg2 [3];
        int t12 = tdeg2 [4];
        int t13 = tdeg2 [5];
        int t15 = tdeg2 [6];
        int t23 = tdeg2 [7];
        int t24 = tdeg2 [8];
        int t34 = tdeg2 [9];
        int t35 = tdeg2 [10];
        int t45 = tdeg2 [11];

        int t0 = tdeg1 [0];
        int t1 = tdeg1 [1];
        int t2 = tdeg1 [2];
        int t3 = tdeg1 [3];
        int t4 = tdeg1 [4];
        int t5 = tdeg1 [5];


        // dependent overlap parameters

        int t018 = t01 - t012 - t015;
        int t048 = t04 - t024 - t045;
        int t027 = t02 - t012 - t024;
        int t057 = t05 - t015 - t045;
        int t078 = t0 - t01 - t02 - t04 - t05 + t012 + t015 + t024 + t045;
        int t138 = t13 - t123 - t135;
        int t348 = t34 - t234 - t345;
        int t237 = t23 - t123 - t234;
        int t357 = t35 - t135 - t345;
        int t378 = t3 - t13 - t23 - t34 - t35 + t123 + t135 + t234 + t345;
        int t126 = t12 - t012 - t123;
        int t156 = t15 - t015 - t135;
        int t168 = t1 - t01 - t12 - t13 - t15 + t012 + t015 + t123 + t135;
        int t246 = t24 - t024 - t234;
        int t456 = t45 - t045 - t345;
        int t468 = t4 - t04 - t24 - t34 - t45 + t024 + t045 + t234 + t345;
        int t267 = t2 - t02 - t12 - t23 - t24 + t012 + t024 + t123 + t234;
        int t567 = t5 - t05 - t15 - t35 - t45 + t015 + t045 + t135 + t345;
        int t678 = kappa - t0 - t1 - t2 - t3 - t4 - t5 + t01 + t02 + t04 + t05 + t12 + t13 + t15 + t23 + t24 + t34 + t35 + t45 - t012 - t015 - t024 - t045 - t123 - t135 - t234 - t345;

        int t07 = t0 - t01 - t04;
        int t08 = t0 - t02 - t05;
        int t16 = t1 - t01 - t13;
        int t18 = t1 - t12 - t15;
        int t26 = t2 - t02 - t23;
        int t27 = t2 - t12 - t24;
        int t37 = t3 - t13 - t34;
        int t38 = t3 - t23 - t35;
        int t46 = t4 - t04 - t34;
        int t48 = t4 - t24 - t45;
        int t56 = t5 - t05 - t35;
        int t57 = t5 - t15 - t45;
        int t67 = kappa - t0 - t3 - t1 - t4 + t01 + t04 + t13 + t34;
        int t68 = kappa - t0 - t3 - t2 - t5 + t02 + t05 + t23 + t35;
        int t78 = kappa - t1 - t4 - t2 - t5 + t12 + t15 + t24 + t45;


        // Enumeration

        int Fs = 0; // 3 * 3 * 3 = 27 choices for abc

        Fs += A ( t012 , t01 , t02 , t12 );
        Fs += A ( t015 , t01 , t05 , t15 );
        Fs += A ( t018 , t01 , t08 , t18 );
        Fs += A ( t024 , t02 , t04 , t24 );
        Fs += A ( t045 , t45 , t04 , t05 );
        Fs += A ( t048 , t04 , t08 , t48 );
        Fs += A ( t027 , t02 , t07 , t27 );
        Fs += A ( t057 , t05 , t07 , t57 );
        Fs += A ( t078 , t78 , t07 , t08 );
        Fs += A ( t123 , t12 , t13 , t23 );
        Fs += A ( t135 , t35 , t13 , t15 );
        Fs += A ( t138 , t13 , t18 , t38 );
        Fs += A ( t234 , t34 , t23 , t24 );
        Fs += A ( t345 , t34 , t35 , t45 );
        Fs += A ( t348 , t34 , t38 , t48 );
        Fs += A ( t237 , t23 , t27 , t37 );
        Fs += A ( t357 , t35 , t37 , t57 );
        Fs += A ( t378 , t78 , t37 , t38 );
        Fs += A ( t126 , t12 , t16 , t26 );
        Fs += A ( t156 , t15 , t16 , t56 );
        Fs += A ( t168 , t68 , t16 , t18 );
        Fs += A ( t246 , t24 , t26 , t46 );
        Fs += A ( t456 , t45 , t46 , t56 );
        Fs += A ( t468 , t68 , t46 , t48 );
        Fs += A ( t267 , t67 , t26 , t27 );
        Fs += A ( t567 , t67 , t56 , t57 );
        Fs += A ( t678 , t67 , t68 , t78 );

        int Fd = 0; // 6 * 12 choices

        // ab ac|bc     3 * 2 * 2 = 12 choices
        Fd += B ( t045 , t04 , t05 , t12 );
        Fd += B ( t048 , t04 , t08 , t15 );
        Fd += B ( t057 , t05 , t07 , t24 );
        Fd += B ( t078 , t07 , t08 , t45 );
        Fd += B ( t345 , t34 , t35 , t12 );
        Fd += B ( t348 , t34 , t38 , t15 );
        Fd += B ( t357 , t35 , t37 , t24 );
        Fd += B ( t378 , t37 , t38 , t45 );
        Fd += B ( t456 , t46 , t56 , t12 );
        Fd += B ( t468 , t68 , t46 , t15 );
        Fd += B ( t567 , t67 , t56 , t24 );
        Fd += B ( t678 , t67 , t68 , t45 );

        // ab bc|ac     3 * 2 * 2 = 12 choices
        Fd += B ( t135 , t13 , t15 , t02 );
        Fd += B ( t138 , t13 , t18 , t05 );
        Fd += B ( t156 , t15 , t16 , t23 );
        Fd += B ( t168 , t16 , t18 , t35 );
        Fd += B ( t345 , t34 , t45 , t02 );
        Fd += B ( t348 , t34 , t48 , t05 );
        Fd += B ( t456 , t45 , t46 , t23 );
        Fd += B ( t468 , t46 , t48 , t35 );
        Fd += B ( t357 , t37 , t57 , t02 );
        Fd += B ( t378 , t78 , t37 , t05 );
        Fd += B ( t567 , t67 , t57 , t23 );
        Fd += B ( t678 , t67 , t78 , t35 );

        // ac bc|ab     3 * 2 * 2 = 12 choices
        Fd += B ( t234 , t23 , t24 , t01 );
        Fd += B ( t237 , t23 , t27 , t04 );
        Fd += B ( t246 , t24 , t26 , t13 );
        Fd += B ( t267 , t26 , t27 , t34 );
        Fd += B ( t345 , t35 , t45 , t01 );
        Fd += B ( t357 , t35 , t57 , t04 );
        Fd += B ( t456 , t45 , t56 , t13 );
        Fd += B ( t567 , t56 , t57 , t34 );
        Fd += B ( t348 , t38 , t48 , t67 );
        Fd += B ( t378 , t78 , t38 , t04 );
        Fd += B ( t468 , t68 , t48 , t13 );
        Fd += B ( t678 , t68 , t78 , t34 );

        // bc|ab ac     3 * 2 * 2 = 12 choices
        Fd += B ( t012 , t01 , t02 , t45 );
        Fd += B ( t015 , t01 , t05 , t48 );
        Fd += B ( t024 , t02 , t04 , t57 );
        Fd += B ( t045 , t04 , t05 , t78 );
        Fd += B ( t123 , t13 , t23 , t45 );
        Fd += B ( t135 , t35 , t13 , t48 );
        Fd += B ( t234 , t34 , t23 , t57 );
        Fd += B ( t345 , t34 , t35 , t78 );
        Fd += B ( t126 , t16 , t26 , t45 );
        Fd += B ( t156 , t16 , t56 , t48 );
        Fd += B ( t246 , t26 , t46 , t57 );
        Fd += B ( t456 , t46 , t56 , t78 );

        // ac|ab bc     3 * 2 * 2 = 12 choices
        Fd += B ( t012 , t01 , t12 , t35 );
        Fd += B ( t015 , t01 , t15 , t38 );
        Fd += B ( t123 , t12 , t13 , t56 );
        Fd += B ( t135 , t13 , t15 , t68 );
        Fd += B ( t024 , t04 , t24 , t35 );
        Fd += B ( t045 , t45 , t04 , t38 );
        Fd += B ( t234 , t34 , t24 , t56 );
        Fd += B ( t345 , t34 , t45 , t68 );
        Fd += B ( t027 , t07 , t27 , t35 );
        Fd += B ( t057 , t07 , t57 , t38 );
        Fd += B ( t237 , t27 , t37 , t56 );
        Fd += B ( t357 , t37 , t57 , t68 );

        // ab|ac bc     3 * 2 * 2 = 12 choices
        Fd += B ( t012 , t02 , t12 , t34 );
        Fd += B ( t024 , t02 , t24 , t37 );
        Fd += B ( t123 , t12 , t23 , t46 );
        Fd += B ( t234 , t23 , t24 , t67 );
        Fd += B ( t015 , t05 , t15 , t34 );
        Fd += B ( t045 , t45 , t05 , t37 );
        Fd += B ( t135 , t35 , t15 , t46 );
        Fd += B ( t345 , t35 , t45 , t67 );
        Fd += B ( t018 , t08 , t18 , t34 );
        Fd += B ( t048 , t08 , t48 , t37 );
        Fd += B ( t138 , t18 , t38 , t46 );
        Fd += B ( t348 , t38 , t48 , t01 );

        int Ft = 0; // 6 * 3 + 6 * 4 choices

        // ab ac| |bc     3 choices
        Ft += B ( t078 , t07 , t08 , t12 );
        Ft += B ( t378 , t37 , t38 , t12 );
        Ft += B ( t678 , t67 , t68 , t12 );


        // ab bc| |ac     3 choices
        Ft += B ( t168 , t16 , t18 , t02 );
        Ft += B ( t468 , t46 , t48 , t02 );
        Ft += B ( t678 , t67 , t78 , t02 );

        // bc ac| |ab     3 choices
        Ft += B ( t267 , t26 , t27 , t01 );
        Ft += B ( t567 , t56 , t57 , t01 );
        Ft += B ( t678 , t68 , t78 , t01 );

        // bc| |ab ac
        Ft += B ( t012 , t01 , t02 , t78 );
        Ft += B ( t123 , t13 , t23 , t78 );
        Ft += B ( t126 , t16 , t26 , t78 );

        // ac| |ab bc
        Ft += B ( t012 , t01 , t12 , t68 );
        Ft += B ( t024 , t04 , t24 , t68 );
        Ft += B ( t027 , t07 , t27 , t68 );

        // ab| |bc ac
        Ft += B ( t012 , t02 , t12 , t67 );
        Ft += B ( t015 , t05 , t15 , t67 );
        Ft += B ( t018 , t08 , t18 , t67 );

        // ab|ac|bc
        Ft += C ( t12 , t05 , t37 );
        Ft += C ( t12 , t35 , t67 );
        Ft += C ( t37 , t08 , t15 );
        Ft += C ( t67 , t38 , t15 );

        // ab|bc|ac
        Ft += C ( t02 , t15 , t46 );
        Ft += C ( t02 , t45 , t67 );
        Ft += C ( t46 , t18 , t05 );
        Ft += C ( t67 , t48 , t05 );

        // ac|ab|bc
        Ft += C ( t12 , t04 , t38 );
        Ft += C ( t12 , t34 , t68 );
        Ft += C ( t38 , t07 , t24 );
        Ft += C ( t68 , t24 , t37 );

        // ac|bc|ab
        Ft += C ( t01 , t24 , t56 );
        Ft += C ( t01 , t45 , t68 );
        Ft += C ( t56 , t27 , t04 );
        Ft += C ( t68 , t04 , t57 );

        // bc|ab|ac
        Ft += C ( t02 , t13 , t48 );
        Ft += C ( t02 , t34 , t78 );
        Ft += C ( t48 , t16 , t23 );
        Ft += C ( t78 , t23 , t46 );

        // bc|ac|ab
        Ft += C ( t01 , t23 , t57 );
        Ft += C ( t01 , t35 , t78 );
        Ft += C ( t57 , t26 , t13 );
        Ft += C ( t78 , t13 , t56 );

        F = ( L * Fs ) + ( ( L - 1 ) * Fd ) + ( ( L - 2 ) * Ft );

    }
    if ( ( gamma == 3 ) && ( m == 1 ) )
    {
        // independebt overlap parameters

        int t012 = tdeg3 [0];

        int t01 = tdeg2 [0];
        int t02 = tdeg2 [1];
        int t12 = tdeg2 [2];

        int t0 = tdeg1 [0];
        int t1 = tdeg1 [1];
        int t2 = tdeg1 [2];


        // dependent overlap parameters

        int t015 = t01 - t012;
        int t024 = t02 - t012;
        int t045 = t0 - t01 - t02 + t012;
        int t123 = t12 - t012;
        int t135 = t1 - t01 - t12 + t012;
        int t234 = t2 - t02 - t12 + t012;
        int t345 = kappa - t0 - t1 - t2 + t01 + t02 + t12 - t012;

        int t04 = t0 - t01;
        int t05 = t0 - t02;
        int t13 = t1 - t01;
        int t15 = t1 - t12;
        int t23 = t2 - t02;
        int t24 = t2 - t12;
        int t34 = kappa - t0 - t1 + t01;
        int t35 = kappa - t0 - t2 + t02;
        int t45 = kappa - t1 - t2 + t12;

        // Enumeration

        int Fs = 0; // 8 choices for abc

        Fs += A ( t012 , t01 , t02 , t12 );
        Fs += A ( t015 , t01 , t05 , t15 );
        Fs += A ( t024 , t02 , t04 , t24 );
        Fs += A ( t045 , t45 , t04 , t05 );
        Fs += A ( t123 , t12 , t13 , t23 );
        Fs += A ( t135 , t35 , t13 , t15 );
        Fs += A ( t234 , t34 , t23 , t24 );
        Fs += A ( t345 , t34 , t35 , t45 );

        int Fd = 0; // 12 choices

        // ab ac|bc
        Fd += B ( t045 , t04 , t05 , t12 );
        Fd += B ( t345 , t34 , t35 , t12 );

        // ab bc|ac
        Fd += B ( t135 , t13 , t15 , t02 );
        Fd += B ( t345 , t34 , t45 , t02 );

        // ac bc|ab
        Fd += B ( t234 , t23 , t24 , t01 );
        Fd += B ( t345 , t35 , t45 , t01 );

        // bc|ab ac
        Fd += B ( t012 , t01 , t02 , t45 );
        Fd += B ( t123 , t13 , t23 , t45 );

        // ac|ab bc
        Fd += B ( t012 , t01 , t12 , t35 );
        Fd += B ( t024 , t04 , t24 , t35 );

        // ab|ac bc
        Fd += B ( t012 , t02 , t12 , t34 );
        Fd += B ( t015 , t05 , t15 , t34 );

        F = ( L * Fs ) + ( ( L - 1 ) * Fd );

    }

    if ( ( gamma == 4 ) && ( m == 1 ) )
    {
        // independebt overlap parameters

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


        // dependent overlap parameters

        int t016 = t01 - t012;
        int t017 = t01 - t013;
        int t025 = t02 - t012;
        int t027 = t02 - t023;
        int t035 = t03 - t013;
        int t036 = t03 - t023;
        int t124 = t12 - t012;
        int t127 = t12 - t123;
        int t134 = t13 - t013;
        int t136 = t13 - t123;
        int t234 = t23 - t023;
        int t235 = t23 - t123;
        int t056 = t0 - t01 - t02 + t012;
        int t057 = t0 - t01 - t03 + t013;
        int t067 = t0 - t02 - t03 + t023;
        int t146 = t1 - t01 - t12 + t012;
        int t147 = t1 - t01 - t13 + t013;
        int t167 = t1 - t12 - t13 + t123;
        int t245 = t2 - t02 - t12 + t012;
        int t247 = t2 - t02 - t23 + t023;
        int t257 = t2 - t12 - t23 + t123;
        int t345 = t3 - t03 - t13 + t013;
        int t346 = t3 - t03 - t23 + t023;
        int t356 = t3 - t13 - t23 + t123;
        int t456 = kappa - t0 - t1 - t2 + t01 + t02 + t12 - t012;
        int t457 = kappa - t0 - t1 - t3 + t01 + t03 + t13 - t013;
        int t467 = kappa - t0 - t2 - t3 + t02 + t03 + t23 - t023;
        int t567 = kappa - t1 - t2 - t3 + t12 + t13 + t23 - t123;

        int t05 = t0 - t01;
        int t06 = t0 - t02;
        int t07 = t0 - t03;
        int t14 = t1 - t01;
        int t16 = t1 - t12;
        int t17 = t1 - t13;
        int t24 = t2 - t02;
        int t25 = t2 - t12;
        int t27 = t2 - t23;
        int t34 = t3 - t03;
        int t35 = t3 - t13;
        int t36 = t3 - t23;
        int t45 = kappa - t0 - t1 + t01;
        int t46 = kappa - t0 - t2 + t02;
        int t47 = kappa - t0 - t3 + t03;
        int t56 = kappa - t1 - t2 + t12;
        int t57 = kappa - t1 - t3 + t13;
        int t67 = kappa - t2 - t3 + t23;

        // Enumeration

        int Fs = 0; //  32 choices

        Fs += A ( t012 , t01 , t02 , t12 );
        Fs += A ( t013 , t01 , t03 , t13 );
        Fs += A ( t016 , t01 , t06 , t16 );
        Fs += A ( t017 , t01 , t07 , t17 );
        Fs += A ( t023 , t02 , t03 , t23 );
        Fs += A ( t025 , t02 , t05 , t25 );
        Fs += A ( t027 , t02 , t07 , t27 );
        Fs += A ( t035 , t03 , t05 , t35 );
        Fs += A ( t036 , t03 , t06 , t36 );
        Fs += A ( t056 , t05 , t06 , t56 );
        Fs += A ( t057 , t05 , t07 , t57 );
        Fs += A ( t067 , t06 , t07 , t67 );
        Fs += A ( t123 , t12 , t13 , t23 );
        Fs += A ( t124 , t12 , t14 , t24 );
        Fs += A ( t127 , t12 , t17 , t27 );
        Fs += A ( t134 , t13 , t14 , t34 );
        Fs += A ( t136 , t13 , t16 , t36 );
        Fs += A ( t146 , t14 , t16 , t46 );
        Fs += A ( t147 , t14 , t17 , t47 );
        Fs += A ( t167 , t16 , t17 , t67 );
        Fs += A ( t234 , t23 , t24 , t34 );
        Fs += A ( t235 , t23 , t25 , t35 );
        Fs += A ( t245 , t24 , t25 , t45 );
        Fs += A ( t247 , t24 , t27 , t47 );
        Fs += A ( t257 , t25 , t27 , t57 );
        Fs += A ( t345 , t34 , t35 , t45 );
        Fs += A ( t346 , t34 , t36 , t46 );
        Fs += A ( t356 , t35 , t36 , t56 );
        Fs += A ( t456 , t45 , t46 , t56 );
        Fs += A ( t457 , t45 , t47 , t57 );
        Fs += A ( t467 , t46 , t47 , t67 );
        Fs += A ( t567 , t56 , t57 , t67 );

        int Fd = 0; // 48 choices

        Fd += B ( t056 , t05 , t06 , t12 );
        Fd += B ( t057 , t05 , t07 , t13 );
        Fd += B ( t067 , t06 , t07 , t23 );
        Fd += B ( t146 , t14 , t16 , t02 );
        Fd += B ( t147 , t14 , t17 , t03 );
        Fd += B ( t167 , t16 , t17 , t23 );
        Fd += B ( t245 , t24 , t25 , t01 );
        Fd += B ( t247 , t24 , t27 , t03 );
        Fd += B ( t257 , t25 , t27 , t13 );
        Fd += B ( t345 , t34 , t35 , t01 );
        Fd += B ( t346 , t34 , t36 , t02 );
        Fd += B ( t356 , t35 , t36 , t12 );
        Fd += B ( t456 , t45 , t46 , t12 );
        Fd += B ( t456 , t45 , t56 , t02 );
        Fd += B ( t456 , t46 , t56 , t01 );
        Fd += B ( t457 , t45 , t47 , t13 );
        Fd += B ( t457 , t45 , t57 , t03 );
        Fd += B ( t457 , t47 , t57 , t01 );
        Fd += B ( t467 , t46 , t47 , t23 );
        Fd += B ( t467 , t46 , t67 , t03 );
        Fd += B ( t467 , t47 , t67 , t02 );
        Fd += B ( t567 , t56 , t57 , t23 );
        Fd += B ( t567 , t56 , t67 , t13 );
        Fd += B ( t567 , t57 , t67 , t12 );

        Fd += B ( t012 , t01 , t02 , t56 );
        Fd += B ( t012 , t01 , t12 , t46 );
        Fd += B ( t012 , t02 , t12 , t45 );
        Fd += B ( t013 , t01 , t03 , t57 );
        Fd += B ( t013 , t01 , t13 , t47 );
        Fd += B ( t013 , t03 , t13 , t45 );
        Fd += B ( t023 , t02 , t03 , t67 );
        Fd += B ( t023 , t02 , t23 , t47 );
        Fd += B ( t023 , t03 , t23 , t46 );
        Fd += B ( t123 , t12 , t13 , t67 );
        Fd += B ( t123 , t12 , t23 , t57 );
        Fd += B ( t123 , t13 , t23 , t56 );
        Fd += B ( t124 , t14 , t24 , t56 );
        Fd += B ( t134 , t14 , t34 , t57 );
        Fd += B ( t234 , t24 , t34 , t67 );
        Fd += B ( t025 , t05 , t25 , t46 );
        Fd += B ( t035 , t05 , t35 , t47 );
        Fd += B ( t235 , t25 , t35 , t67 );
        Fd += B ( t016 , t06 , t16 , t45 );
        Fd += B ( t036 , t06 , t36 , t47 );
        Fd += B ( t136 , t16 , t36 , t57 );
        Fd += B ( t017 , t07 , t17 , t45 );
        Fd += B ( t027 , t07 , t27 , t46 );
        Fd += B ( t127 , t17 , t27 , t56 );


        F = ( L * Fs ) + ( ( L - 1 ) * Fd );

    }

    return F;
}

// Finding a partitioning that achieve a particular set of overlapping parameters
// The partitioning matrix PM in the input arguments need to be filled
// The output is the number of partitionig matrices that achieve the set of given overlapping parameters. But, PM is filled by one of the choices
int Partitioning_Matrix_Construction ( int kappa , int gamma , int m , int* tdeg3 , int* tdeg2 , int* tdeg1 , int** PM )
{

    int options = 1;

    if ( ( gamma == 3 ) && ( m == 2 ) )
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

        options *= comb_choose ( kappa , t0 );

        options *= comb_choose ( t0 , t01 );
        options *= comb_choose ( kappa - t0 , t1 - t01 );

        options *= comb_choose ( t01 , t012 );
        options *= comb_choose ( t0 - t01 , t02 - t012 );
        options *= comb_choose ( t1 - t01 , t12 - t012 );
        options *= comb_choose ( kappa - t0 - t1 + t01 , t2 - t02 - t12 + t012 );

        options *= comb_choose ( t12 - t012 , t123 );
        options *= comb_choose ( t1 - t01 - t12 + t012 , t13 - t123 );
        options *= comb_choose ( t2 - t02 - t12 + t012 , t23 - t123 );
        options *= comb_choose ( kappa - t0 - t1 - t2 + t01 + t02 + t12 - t012 , t3 - t13 - t23 + t123 );

        options *= comb_choose ( t02 - t012 , t024 );
        options *= comb_choose ( t23 - t123 , t234 );
        options *= comb_choose ( t0 - t01 - t02 + t012 , t04 - t024 );
        options *= comb_choose ( t2 - t02 - t12 - t23 + t012 + t123 , t24 - t024 - t234 );
        options *= comb_choose ( t3 - t13 - t23 + t123 , t34 - t234 );
        options *= comb_choose ( kappa - t0 - t1 - t2 - t3 + t01 + t02 + t12 + t13 +t23 - t012 - t123 , t4 - t04 - t24 - t34 + t024 + t234 );

        options *= comb_choose ( t01 - t012 , t015 );
        options *= comb_choose ( t04 - t024 , t045 );
        options *= comb_choose ( t13 - t123 , t135 );
        options *= comb_choose ( t34 - t234 , t345 );
        options *= comb_choose ( t0 - t01 - t02 - t04 + t012 + t024 , t05 - t015 - t045 );
        options *= comb_choose ( t1 - t01 - t12 - t13 + t012 + t123 , t15 - t015 - t135 );
        options *= comb_choose ( t3 - t13 - t23 - t34 + t123 + t234 , t35 - t135 - t345 );
        options *= comb_choose ( t4 - t04 - t24 - t34 + t024 + t234 , t45 - t045 - t345 );
        options *= comb_choose ( kappa - t0 - t1 - t2 - t3 - t4 + t01 + t02 + t04 + t12 + t13 + t23 + t24 + t34 - t012 - t123 - t024 - t234 , t5 - t05 - t15 - t35 - t45 + t015 + t045 + t135 + t345 );

        int index = 0;

        for ( int i = 0 ; i < gamma ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) PM[i][j] = 2;

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

    if ( ( gamma == 4 ) && ( m == 1 ) )
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

        options *= comb_choose ( kappa , t0 );

        options *= comb_choose ( t0 , t01 );
        options *= comb_choose ( kappa - t0 , t1 - t01 );

        options *= comb_choose ( t01 , t012 );
        options *= comb_choose ( t0 - t01 , t02 - t012 );
        options *= comb_choose ( t1 - t01 , t12 - t012 );
        options *= comb_choose ( kappa - t0 - t1 + t01 , t2 - t02 - t12 + t012 );

        options *= comb_choose ( t012 , t0123 );
        options *= comb_choose ( t01 - t012 , t013 - t0123 );
        options *= comb_choose ( t02 - t012 , t023 - t0123 );
        options *= comb_choose ( t12 - t012 , t123 - t0123 );
        options *= comb_choose ( t0 - t01 - t02 + t012 , t03 - t013 - t023 + t0123 );
        options *= comb_choose ( t1 - t01 - t12 + t012 , t13 - t013 - t123 + t0123 );
        options *= comb_choose ( t2 - t02 - t12 + t012 , t23 - t023 - t123 + t0123 );
        options *= comb_choose ( kappa - t0 - t1 - t2 + t01 + t02 + t12 - t012 , t3 - t03 - t13 - t23 + t013 + t023 + t123 - t0123 );

        int index = 0;

        for ( int i = 0 ; i < gamma ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) PM[i][j] = 1;

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

    if ( ( gamma == 3 ) && ( m == 1 ) )
    {
        int t012 = tdeg3[0];

        int t01 = tdeg2[0];
        int t02 = tdeg2[1];
        int t12 = tdeg2[2];

        int t0 = tdeg1[0];
        int t1 = tdeg1[1];
        int t2 = tdeg1[2];

        options *= comb_choose ( kappa , t0 );

        options *= comb_choose ( t0 , t01 );
        options *= comb_choose ( kappa - t0 , t1 - t01 );

        options *= comb_choose ( t01 , t012 );
        options *= comb_choose ( t0 - t01 , t02 - t012 );
        options *= comb_choose ( t1 - t01 , t12 - t012 );
        options *= comb_choose ( kappa - t0 - t1 + t01 , t2 - t02 - t12 + t012 );

        int index = 0;

        for ( int i = 0 ; i < gamma ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) PM[i][j] = 1;

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

    return options;


}

int main ()
{
    int kappa = 33;
    int gamma = 4;
    int p  = 1;
    int L = 4;
    int m = 1;

    if ( ( gamma == 3 ) && ( m == 2 ) )
    {
        int* tdeg3 = new int [8];
        int* tdeg2 = new int [12];
        int* tdeg1 = new int [6];

        int F_min = -1;
        bool first_trial = true;

        int** tdeg3_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg3_opt_sets[i] = new int [8];
        int** tdeg2_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg2_opt_sets[i] = new int [12];
        int** tdeg1_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg1_opt_sets[i] = new int [6];

        int opt_sets_num = 0; // the number of overlapping sets that achieve the minimum number of cycles 6 in the protograph SC code

        for ( int t0 = 0 ; t0 <= kappa ; t0++ )
        {
            //cout << "--      " << t0 << endl;
            for ( int t01 = 0 ; t01 <= t0 ; t01++ )
            {
            for ( int t1 = t01 ; t1 <= (kappa-t0+t01) ; t1++ )
            {
                //if ( ( t0 + t1 ) <= kappa )
                for ( int t012 = 0 ; t012 <= t01 ; t012++ )
                {
                for ( int t02 = t012 ; t02 <= (t0-t01+t012) ; t02++ )
                {
                for ( int t12 = t012 ; t12 <= (t1-t01+t012) ; t12++ )
                {
                for ( int t2 = (t02+t12-t012) ; t2 <= (kappa-t0-t1+t01+t02+t12-t012) ; t2++ )
                {
                    //if ( ( t0 + t1 + t2 ) == kappa )
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
                            //if ( ( t3 + t4 ) <= kappa )
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
                                //if ( ( t3 + t4 + t5 ) == kappa )
                                //{
                                tdeg3[0] = t012;
                                tdeg3[1] = t015;
                                tdeg3[2] = t024;
                                tdeg3[3] = t045;
                                tdeg3[4] = t123;
                                tdeg3[5] = t135;
                                tdeg3[6] = t234;
                                tdeg3[7] = t345;

                                tdeg2[0] = t01;
                                tdeg2[1] = t02;
                                tdeg2[2] = t04;
                                tdeg2[3] = t05;
                                tdeg2[4] = t12;
                                tdeg2[5] = t13;
                                tdeg2[6] = t15;
                                tdeg2[7] = t23;
                                tdeg2[8] = t24;
                                tdeg2[9] = t34;
                                tdeg2[10] = t35;
                                tdeg2[11] = t45;

                                tdeg1[0] = t0;
                                tdeg1[1] = t1;
                                tdeg1[2] = t2;
                                tdeg1[3] = t3;
                                tdeg1[4] = t4;
                                tdeg1[5] = t5;


                                int F = cycle_6_enumeration_protograph_SC ( kappa , gamma , L , m , tdeg3 , tdeg2 , tdeg1  );
                                if ( ( first_trial ) || ( F < F_min ) )
                                {
                                    opt_sets_num = 1;
                                    F_min = F;
                                    for ( int i = 0 ; i < 8 ; i++ ) tdeg3_opt_sets[0][i] = tdeg3[i];
                                    for ( int i = 0 ; i < 12 ; i++ ) tdeg2_opt_sets[0][i] = tdeg2[i];
                                    for ( int i = 0 ; i < 6 ; i++ ) tdeg1_opt_sets[0][i] = tdeg1[i];
                                    first_trial = false;
                                    cout << F << endl;
                                }
                                else if ( F == F_min )
                                {
                                    for ( int i = 0 ; i < 8 ; i++ ) tdeg3_opt_sets[opt_sets_num][i] = tdeg3[i];
                                    for ( int i = 0 ; i < 12 ; i++ ) tdeg2_opt_sets[opt_sets_num][i] = tdeg2[i];
                                    for ( int i = 0 ; i < 6 ; i++ ) tdeg1_opt_sets[opt_sets_num][i] = tdeg1[i];
                                    opt_sets_num ++;
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
        for ( int i = 0 ; i < 8 ; i++ ) cout << tdeg3_opt_sets[0][i] << "   "; cout << endl;
        for ( int i = 0 ; i < 12 ; i++ ) cout << tdeg2_opt_sets[0][i] << "   "; cout << endl;
        for ( int i = 0 ; i < 6 ; i++ ) cout << tdeg1_opt_sets[0][i] << "   "; cout << endl;

        int** PM = new int* [ gamma ];
        for ( int i = 0 ; i < gamma ; i++ ) PM[i] = new int [ kappa ];

        Partitioning_Matrix_Construction ( kappa , gamma , m , tdeg3_opt_sets[0] , tdeg2_opt_sets[0] , tdeg1_opt_sets[0] , PM );

        cout << "Optimum Partitioning Matrix: " << endl;
        for ( int i = 0 ; i < gamma ; i++ )
        {
            for ( int j = 0 ; j < kappa ; j++ ) cout << PM[i][j] << "  "; cout << endl;
        }

        for ( int i = 0 ; i < gamma ; i++ ) delete[] PM[i];
        delete[] PM;
        delete[] tdeg3; delete[] tdeg2; delete[] tdeg1;
        for ( int i = 0 ; i < 1000 ; i++ ) { delete[] tdeg3_opt_sets[i]; delete[] tdeg2_opt_sets[i]; delete[] tdeg1_opt_sets[i]; }
        delete[] tdeg3_opt_sets; delete[] tdeg2_opt_sets; delete[] tdeg1_opt_sets;
    }

    if ( ( gamma == 3 ) && ( m == 1 ) )
    {
        int* tdeg3 = new int [1];
        int* tdeg2 = new int [3];
        int* tdeg1 = new int [3];

        int F_min = -1;
        bool first_trial = true;

        int** tdeg3_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg3_opt_sets[i] = new int [1];
        int** tdeg2_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg2_opt_sets[i] = new int [3];
        int** tdeg1_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg1_opt_sets[i] = new int [3];

        int opt_sets_num = 0; // the number of overlapping sets that achieve the minimum number of cycles 6 in the protograph SC code

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
                if ( ( ( t0 + t1 + t2 ) >= floor( kappa * gamma / 2 ) ) && ( ( t0 + t1 + t2 ) <= ceil( kappa * gamma / 2 ) ) )
                {
                    tdeg3[0] = t012;
                    tdeg2[0] = t01;
                    tdeg2[1] = t02;
                    tdeg2[2] = t12;
                    tdeg1[0] = t0;
                    tdeg1[1] = t1;
                    tdeg1[2] = t2;
                    int F = cycle_6_enumeration_protograph_SC ( kappa , gamma , L , m , tdeg3 , tdeg2 , tdeg1  );
                    if ( ( first_trial ) || ( F < F_min ) )
                    {
                        opt_sets_num = 1;
                        F_min = F;
                        for ( int i = 0 ; i < 1 ; i++ ) tdeg3_opt_sets[0][i] = tdeg3[i];
                        for ( int i = 0 ; i < 3 ; i++ ) tdeg2_opt_sets[0][i] = tdeg2[i];
                        for ( int i = 0 ; i < 3 ; i++ ) tdeg1_opt_sets[0][i] = tdeg1[i];
                        first_trial = false;
                    }
                    else if ( F == F_min )
                    {
                        for ( int i = 0 ; i < 1 ; i++ ) tdeg3_opt_sets[opt_sets_num][i] = tdeg3[i];
                        for ( int i = 0 ; i < 3 ; i++ ) tdeg2_opt_sets[opt_sets_num][i] = tdeg2[i];
                        for ( int i = 0 ; i < 3 ; i++ ) tdeg1_opt_sets[opt_sets_num][i] = tdeg1[i];
                        opt_sets_num ++;
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
        for ( int i = 0 ; i < 1 ; i++ ) cout << tdeg3_opt_sets[0][i] << "   "; cout << endl;
        for ( int i = 0 ; i < 3 ; i++ ) cout << tdeg2_opt_sets[0][i] << "   "; cout << endl;
        for ( int i = 0 ; i < 3 ; i++ ) cout << tdeg1_opt_sets[0][i] << "   "; cout << endl;

        int** PM = new int* [ gamma ];
        for ( int i = 0 ; i < gamma ; i++ ) PM[i] = new int [ kappa ];

        Partitioning_Matrix_Construction ( kappa , gamma , m , tdeg3_opt_sets[0] , tdeg2_opt_sets[0] , tdeg1_opt_sets[0] , PM );

        cout << "Optimum Partitioning Matrix: " << endl;
        for ( int i = 0 ; i < gamma ; i++ )
        {
            for ( int j = 0 ; j < kappa ; j++ ) cout << PM[i][j] << "  "; cout << endl;
        }

        for ( int i = 0 ; i < gamma ; i++ ) delete[] PM[i];
        delete[] PM;
        delete[] tdeg3; delete[] tdeg2; delete[] tdeg1;
        for ( int i = 0 ; i < 1000 ; i++ ) { delete[] tdeg3_opt_sets[i]; delete[] tdeg2_opt_sets[i]; delete[] tdeg1_opt_sets[i]; }
        delete[] tdeg3_opt_sets; delete[] tdeg2_opt_sets; delete[] tdeg1_opt_sets;
    }

    if ( ( gamma == 4 ) && ( m == 1 ) )
    {
        int* tdeg3 = new int [4];
        int* tdeg2 = new int [6];
        int* tdeg1 = new int [4];

        int F_min = -1;
        bool first_trial = true;

        int** tdeg3_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg3_opt_sets[i] = new int [4];
        int** tdeg2_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg2_opt_sets[i] = new int [6];
        int** tdeg1_opt_sets = new int* [1000]; for ( int i = 0 ; i < 1000 ; i++ ) tdeg1_opt_sets[i] = new int [4];

        int opt_sets_num = 0; // the number of overlapping sets that achieve the minimum number of cycles 6 in the protograph SC code

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
                        if ( ( t0 + t1 + t2 + t3 ) == ( kappa * gamma / 2 ) )
                        {
                            tdeg3 [0] = t012; tdeg3 [1] = t013; tdeg3 [2] = t023; tdeg3 [3] = t123;
                            tdeg2 [0] = t01; tdeg2 [1] = t02; tdeg2 [2] = t03; tdeg2 [3] = t12; tdeg2 [4] = t13; tdeg2 [5] = t23;
                            tdeg1 [0] = t0; tdeg1 [1] = t1; tdeg1 [2] = t2; tdeg1 [3] = t3;

                            int F = cycle_6_enumeration_protograph_SC ( kappa , gamma , L , m , tdeg3 , tdeg2 , tdeg1  );
                            if ( ( first_trial ) || ( F < F_min ) )
                            {
                                opt_sets_num = 1;
                                F_min = F;
                                for ( int i = 0 ; i < 4 ; i++ ) tdeg3_opt_sets[0][i] = tdeg3[i];
                                for ( int i = 0 ; i < 6 ; i++ ) tdeg2_opt_sets[0][i] = tdeg2[i];
                                for ( int i = 0 ; i < 6 ; i++ ) tdeg1_opt_sets[0][i] = tdeg1[i];
                                first_trial = false;
                            }
                            else if ( F == F_min )
                            {
                                for ( int i = 0 ; i < 4 ; i++ ) tdeg3_opt_sets[opt_sets_num][i] = tdeg3[i];
                                for ( int i = 0 ; i < 6 ; i++ ) tdeg2_opt_sets[opt_sets_num][i] = tdeg2[i];
                                for ( int i = 0 ; i < 4 ; i++ ) tdeg1_opt_sets[opt_sets_num][i] = tdeg1[i];
                                opt_sets_num ++;
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
        for ( int i = 0 ; i < 4 ; i++ ) cout << tdeg3_opt_sets[0][i] << "   "; cout << endl;
        for ( int i = 0 ; i < 6 ; i++ ) cout << tdeg2_opt_sets[0][i] << "   "; cout << endl;
        for ( int i = 0 ; i < 4 ; i++ ) cout << tdeg1_opt_sets[0][i] << "   "; cout << endl;

        int** PM = new int* [ gamma ];
        for ( int i = 0 ; i < gamma ; i++ ) PM[i] = new int [ kappa ];

        Partitioning_Matrix_Construction ( kappa , gamma , m , tdeg3_opt_sets[0] , tdeg2_opt_sets[0] , tdeg1_opt_sets[0] , PM );

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
                        //H_SC[((r+1)*gamma+i)%(L*gamma)][r*kappa+j] = 1;
                        H_SC[(r+1)*gamma+i][r*kappa+j] = 1;
                    }
                }
            }
        }

        cout << cycle_6_enumeration_protograph_block ( H_SC , (L+m) * gamma , L * kappa ) << endl;
        cout << cycle_4_enumeration_protograph_block ( H_SC , (L+m) * gamma , L * kappa ) << endl;

        for ( int i = 0 ; i < ( (L+m) * gamma ) ; i++ ) delete[] H_SC[i]; delete[] H_SC;


        for ( int i = 0 ; i < gamma ; i++ ) delete[] PM[i];
        delete[] PM;

        delete[] tdeg3; delete[] tdeg2; delete[] tdeg1;
        for ( int i = 0 ; i < 1000 ; i++ ) { delete[] tdeg3_opt_sets[i]; delete[] tdeg2_opt_sets[i]; delete[] tdeg1_opt_sets[i]; }
        delete[] tdeg3_opt_sets; delete[] tdeg2_opt_sets; delete[] tdeg1_opt_sets;

    }
}
