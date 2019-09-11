/****************************************************************************************/
/* Target: Enumeration of objects for SC codes with minimum overlap (MO) partitioning   */
/* By: Homa Esfahanizadeh (hesfahanizadeh@ucla.edu)                                     */
/* University of California, Los Angeles (UCLA)                                         */
/* Update: 09/06/2019                                                                   */
/* Comments: Repoting cycles that span different number of replicas separately          */
/*           Irregularity is added                                                      */
/****************************************************************************************/

#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;

int cycle_4_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP, int R1_limit );

int cycle_6_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP, int R1_limit );

int cycle_8_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP, int R1_limit );

int object_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP , int a );

int main ()
{
    // input parameters
    int p = 13; // Circulant size
    int kappa = 13; // Row weight
    int gamma = 6; // Column weight
    int L = 10; // Coupling length
    int m = 1;

    // PM is the Patitioning Matrix
    int** PM = new int* [ gamma ];
    for ( int i = 0 ; i < gamma ; i++ ) PM[i] = new int [ kappa ];

    // CP is the Circulant Powers
    int** CP = new int* [ gamma ];
    for ( int i = 0 ; i < gamma ; i++ ) CP[i] = new int [ kappa ];

    //for ( int i = 0 ; i < gamma ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) PM[i][j]=0;
    for ( int i = 0 ; i < gamma ; i++ ) for ( int j = 0 ; j < kappa ; j++ ) CP[i][j]=i*j;


    // constrained OO approach
    //PM[0][0]=0; PM[0][1]=1; PM[0][2]=0; PM[0][3]=1; PM[0][4]=0; PM[0][5]=1; PM[0][6]=0; PM[0][7]=1; PM[0][8]=0; PM[0][9]=1; PM[0][10]=0; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=1; PM[1][1]=0; PM[1][2]=1; PM[1][3]=0; PM[1][4]=1; PM[1][5]=0; PM[1][6]=1; PM[1][7]=0; PM[1][8]=1; PM[1][9]=0; PM[1][10]=1; PM[1][11]=0; PM[1][12]=0;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=1; PM[2][7]=1; PM[2][8]=1; PM[2][9]=1; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;

    // appended OO approach
    //PM[0][0]=0; PM[0][1]=1; PM[0][2]=1; PM[0][3]=1; PM[0][4]=1; PM[0][5]=1; PM[0][6]=1; PM[0][7]=1; PM[0][8]=1; PM[0][9]=1; PM[0][10]=1; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=1; PM[1][1]=0; PM[1][2]=0; PM[1][3]=0; PM[1][4]=0; PM[1][5]=0; PM[1][6]=0; PM[1][7]=1; PM[1][8]=1; PM[1][9]=1; PM[1][10]=1; PM[1][11]=1; PM[1][12]=1;
    //PM[2][0]=0; PM[2][1]=1; PM[2][2]=1; PM[2][3]=1; PM[2][4]=1; PM[2][5]=1; PM[2][6]=1; PM[2][7]=0; PM[2][8]=0; PM[2][9]=0; PM[2][10]=0; PM[2][11]=0; PM[2][12]=0;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;

    // appended CV approach
    //PM[0][0]=1; PM[0][1]=1; PM[0][2]=1; PM[0][3]=0; PM[0][4]=0; PM[0][5]=0; PM[0][6]=0; PM[0][7]=0; PM[0][8]=0; PM[0][9]=0; PM[0][10]=0; PM[0][11]=0; PM[0][12]=0;
    //PM[1][0]=1; PM[1][1]=1; PM[1][2]=1; PM[1][3]=1; PM[1][4]=1; PM[1][5]=1; PM[1][6]=0; PM[1][7]=0; PM[1][8]=0; PM[1][9]=0; PM[1][10]=0; PM[1][11]=0; PM[1][12]=0;
    //PM[2][0]=1; PM[2][1]=1; PM[2][2]=1; PM[2][3]=1; PM[2][4]=1; PM[2][5]=1; PM[2][6]=1; PM[2][7]=1; PM[2][8]=1; PM[2][9]=0; PM[2][10]=0; PM[2][11]=0; PM[2][12]=0;
    //PM[3][0]=1; PM[3][1]=1; PM[3][2]=1; PM[3][3]=1; PM[3][4]=1; PM[3][5]=1; PM[3][6]=1; PM[3][7]=1; PM[3][8]=1; PM[3][9]=1; PM[3][10]=1; PM[3][11]=1; PM[3][12]=1;
    //PM[4][0]=1; PM[4][1]=1; PM[4][2]=1; PM[4][3]=1; PM[4][4]=1; PM[4][5]=1; PM[4][6]=1; PM[4][7]=1; PM[4][8]=1; PM[4][9]=1; PM[4][10]=1; PM[4][11]=1; PM[4][12]=1;

    //PM[0][0]=0; PM[0][1]=0; PM[0][2]=0; PM[0][3]=1; PM[0][4]=1; PM[0][5]=1; PM[0][6]=1; PM[0][7]=1; PM[0][8]=1; PM[0][9]=1; PM[0][10]=1; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=0; PM[1][1]=0; PM[1][2]=0; PM[1][3]=0; PM[1][4]=0; PM[1][5]=0; PM[1][6]=1; PM[1][7]=1; PM[1][8]=1; PM[1][9]=1; PM[1][10]=1; PM[1][11]=1; PM[1][12]=1;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=0; PM[2][7]=0; PM[2][8]=0; PM[2][9]=1; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;

    // OO-CPO
    // PM[0][0]=0; PM[0][1]=1; PM[0][2]=0; PM[0][3]=1; PM[0][4]=0; PM[0][5]=1; PM[0][6]=0; PM[0][7]=1; PM[0][8]=0; PM[0][9]=1; PM[0][10]=0; PM[0][11]=1; PM[0][12]=1;
    // PM[1][0]=1; PM[1][1]=0; PM[1][2]=1; PM[1][3]=0; PM[1][4]=1; PM[1][5]=0; PM[1][6]=1; PM[1][7]=0; PM[1][8]=1; PM[1][9]=0; PM[1][10]=1; PM[1][11]=0; PM[1][12]=0;
    // PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=1; PM[2][7]=1; PM[2][8]=1; PM[2][9]=1; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;

    // CP[0][0]=0; CP[0][1]=11; CP[0][2]=11; CP[0][3]=4; CP[0][4]=11; CP[0][5]=0; CP[0][6]=10; CP[0][7]=3; CP[0][8]=9; CP[0][9]=0; CP[0][10]=0; CP[0][11]=0; CP[0][12]=4;
    // CP[1][0]=0; CP[1][1]=1; CP[1][2]=2; CP[1][3]=3; CP[1][4]=4; CP[1][5]=5; CP[1][6]=6; CP[1][7]=7; CP[1][8]=8; CP[1][9]=9; CP[1][10]=10;CP[1][11]=11; CP[1][12]=12;
    // CP[2][0]=8; CP[2][1]=2; CP[2][2]=4; CP[2][3]=6; CP[2][4]=8; CP[2][5]=10; CP[2][6]=12; CP[2][7]=1; CP[2][8]=3; CP[2][9]=5; CP[2][10]=7; CP[2][11]=0; CP[2][12]=0;


    // constrained OO approach
    //PM[0][0]=0; PM[0][1]=1; PM[0][2]=0; PM[0][3]=1; PM[0][4]=0; PM[0][5]=1; PM[0][6]=0; PM[0][7]=1; PM[0][8]=0; PM[0][9]=1; PM[0][10]=0; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=1; PM[1][1]=0; PM[1][2]=1; PM[1][3]=0; PM[1][4]=1; PM[1][5]=0; PM[1][6]=1; PM[1][7]=0; PM[1][8]=1; PM[1][9]=0; PM[1][10]=1; PM[1][11]=0; PM[1][12]=0;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=1; PM[2][3]=0; PM[2][4]=1; PM[2][5]=0; PM[2][6]=1; PM[2][7]=0; PM[2][8]=1; PM[2][9]=0; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0;

    // constrained OO approach and some circulants are removed
    //PM[0][0]=0; PM[0][1]=1; PM[0][2]=0; PM[0][3]=1; PM[0][4]=0; PM[0][5]=1; PM[0][6]=0; PM[0][7]=1; PM[0][8]=0; PM[0][9]=1; PM[0][10]=0; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=1; PM[1][1]=0; PM[1][2]=1; PM[1][3]=0; PM[1][4]=1; PM[1][5]=0; PM[1][6]=1; PM[1][7]=0; PM[1][8]=1; PM[1][9]=0; PM[1][10]=1; PM[1][11]=0; PM[1][12]=0;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=1; PM[2][3]=0; PM[2][4]=1; PM[2][5]=0; PM[2][6]=1; PM[2][7]=0; PM[2][8]=1; PM[2][9]=0; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=-1; PM[3][9]=-1; PM[3][10]=-1; PM[3][11]=-1; PM[3][12]=-1;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0;


    // appended OO approach
    PM[0][0]=0; PM[0][1]=1; PM[0][2]=1; PM[0][3]=1; PM[0][4]=1; PM[0][5]=1; PM[0][6]=1; PM[0][7]=1; PM[0][8]=1; PM[0][9]=1; PM[0][10]=1; PM[0][11]=1; PM[0][12]=1;
    PM[1][0]=1; PM[1][1]=0; PM[1][2]=0; PM[1][3]=0; PM[1][4]=0; PM[1][5]=0; PM[1][6]=0; PM[1][7]=1; PM[1][8]=1; PM[1][9]=1; PM[1][10]=1; PM[1][11]=1; PM[1][12]=1;
    PM[2][0]=0; PM[2][1]=1; PM[2][2]=1; PM[2][3]=1; PM[2][4]=1; PM[2][5]=1; PM[2][6]=1; PM[2][7]=0; PM[2][8]=0; PM[2][9]=0; PM[2][10]=0; PM[2][11]=0; PM[2][12]=0;
    PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0;
    PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;
    PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0;

    // appended CV approach
    //PM[0][0]=1; PM[0][1]=1; PM[0][2]=1; PM[0][3]=0; PM[0][4]=0; PM[0][5]=0; PM[0][6]=0; PM[0][7]=0; PM[0][8]=0; PM[0][9]=0; PM[0][10]=0; PM[0][11]=0; PM[0][12]=0;
    //PM[1][0]=1; PM[1][1]=1; PM[1][2]=1; PM[1][3]=1; PM[1][4]=1; PM[1][5]=1; PM[1][6]=0; PM[1][7]=0; PM[1][8]=0; PM[1][9]=0; PM[1][10]=0; PM[1][11]=0; PM[1][12]=0;
    //PM[2][0]=1; PM[2][1]=1; PM[2][2]=1; PM[2][3]=1; PM[2][4]=1; PM[2][5]=1; PM[2][6]=1; PM[2][7]=1; PM[2][8]=1; PM[2][9]=0; PM[2][10]=0; PM[2][11]=0; PM[2][12]=0;
    //PM[3][0]=1; PM[3][1]=1; PM[3][2]=1; PM[3][3]=1; PM[3][4]=1; PM[3][5]=1; PM[3][6]=1; PM[3][7]=1; PM[3][8]=1; PM[3][9]=1; PM[3][10]=1; PM[3][11]=1; PM[3][12]=1;
    //PM[4][0]=1; PM[4][1]=1; PM[4][2]=1; PM[4][3]=1; PM[4][4]=1; PM[4][5]=1; PM[4][6]=1; PM[4][7]=1; PM[4][8]=1; PM[4][9]=1; PM[4][10]=1; PM[4][11]=1; PM[4][12]=1;
    //PM[5][0]=1; PM[5][1]=1; PM[5][2]=1; PM[5][3]=1; PM[5][4]=1; PM[5][5]=1; PM[5][6]=1; PM[5][7]=1; PM[5][8]=1; PM[5][9]=1; PM[5][10]=1; PM[5][11]=1; PM[5][12]=1;

    //PM[0][0]=0; PM[0][1]=0; PM[0][2]=0; PM[0][3]=1; PM[0][4]=1; PM[0][5]=1; PM[0][6]=1; PM[0][7]=1; PM[0][8]=1; PM[0][9]=1; PM[0][10]=1; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=0; PM[1][1]=0; PM[1][2]=0; PM[1][3]=0; PM[1][4]=0; PM[1][5]=0; PM[1][6]=1; PM[1][7]=1; PM[1][8]=1; PM[1][9]=1; PM[1][10]=1; PM[1][11]=1; PM[1][12]=1;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=0; PM[2][7]=0; PM[2][8]=0; PM[2][9]=1; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0;


    // constrained OO approach and appended OO approach
    //PM[0][0]=0; PM[0][1]=1; PM[0][2]=0; PM[0][3]=1; PM[0][4]=0; PM[0][5]=1; PM[0][6]=0; PM[0][7]=1; PM[0][8]=0; PM[0][9]=1; PM[0][10]=0; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=1; PM[1][1]=0; PM[1][2]=1; PM[1][3]=0; PM[1][4]=1; PM[1][5]=0; PM[1][6]=1; PM[1][7]=0; PM[1][8]=1; PM[1][9]=0; PM[1][10]=1; PM[1][11]=0; PM[1][12]=0;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=1; PM[2][7]=1; PM[2][8]=1; PM[2][9]=1; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;
    //PM[3][0]=1; PM[3][1]=1; PM[3][2]=1; PM[3][3]=1; PM[3][4]=1; PM[3][5]=1; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0;

    // appended CV approach
    //PM[0][0]=0; PM[0][1]=0; PM[0][2]=0; PM[0][3]=1; PM[0][4]=1; PM[0][5]=1; PM[0][6]=1; PM[0][7]=1; PM[0][8]=1; PM[0][9]=1; PM[0][10]=1; PM[0][11]=1; PM[0][12]=1;
    //PM[1][0]=0; PM[1][1]=0; PM[1][2]=0; PM[1][3]=0; PM[1][4]=0; PM[1][5]=1; PM[1][6]=1; PM[1][7]=1; PM[1][8]=1; PM[1][9]=1; PM[1][10]=1; PM[1][11]=1; PM[1][12]=1;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=0; PM[2][7]=0; PM[2][8]=1; PM[2][9]=1; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=1; PM[3][12]=1;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0;

    // constrained OO approach
    //PM[0][0]=0; PM[0][1]=1; PM[0][2]=0; PM[0][3]=1; PM[0][4]=0; PM[0][5]=1; PM[0][6]=0; PM[0][7]=1; PM[0][8]=0; PM[0][9]=1; PM[0][10]=0; PM[0][11]=1; PM[0][12]=0; PM[0][13]=1; PM[0][14]=0; PM[0][15]=1; PM[0][16]=1;
    //PM[1][0]=1; PM[1][1]=0; PM[1][2]=1; PM[1][3]=0; PM[1][4]=1; PM[1][5]=0; PM[1][6]=1; PM[1][7]=0; PM[1][8]=1; PM[1][9]=0; PM[1][10]=1; PM[1][11]=0; PM[1][12]=1; PM[1][13]=0; PM[1][14]=1; PM[1][15]=0; PM[1][16]=0;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=1; PM[2][7]=0; PM[2][8]=1; PM[2][9]=0; PM[2][10]=1; PM[2][11]=1; PM[2][12]=1; PM[2][13]=1; PM[2][14]=1; PM[2][15]=1; PM[2][16]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0; PM[3][13]=0; PM[3][14]=0; PM[3][15]=0; PM[3][16]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0; PM[4][13]=0; PM[4][14]=0; PM[4][15]=0; PM[4][16]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0; PM[5][13]=0; PM[5][14]=0; PM[5][15]=0; PM[5][16]=0;

    // appended OO approach
    //PM[0][0]=0; PM[0][1]=1; PM[0][2]=1; PM[0][3]=1; PM[0][4]=1; PM[0][5]=1; PM[0][6]=1; PM[0][7]=1; PM[0][8]=1; PM[0][9]=1; PM[0][10]=1; PM[0][11]=1; PM[0][12]=1; PM[0][13]=1; PM[0][14]=1; PM[0][15]=1; PM[0][16]=1;
    //PM[1][0]=1; PM[1][1]=0; PM[1][2]=0; PM[1][3]=0; PM[1][4]=0; PM[1][5]=0; PM[1][6]=0; PM[1][7]=0; PM[1][8]=0; PM[1][9]=1; PM[1][10]=1; PM[1][11]=1; PM[1][12]=1; PM[1][13]=1; PM[1][14]=1; PM[1][15]=1; PM[1][16]=1;
    //PM[2][0]=0; PM[2][1]=1; PM[2][2]=1; PM[2][3]=1; PM[2][4]=1; PM[2][5]=1; PM[2][6]=1; PM[2][7]=1; PM[2][8]=1; PM[2][9]=0; PM[2][10]=0; PM[2][11]=0; PM[2][12]=0; PM[2][13]=0; PM[2][14]=0; PM[2][15]=0; PM[2][16]=0;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0; PM[3][13]=0; PM[3][14]=0; PM[3][15]=0; PM[3][16]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0; PM[4][13]=0; PM[4][14]=0; PM[4][15]=0; PM[4][16]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0; PM[5][13]=0; PM[5][14]=0; PM[5][15]=0; PM[5][16]=0;

    // appended CV approach
    //PM[0][0]=0; PM[0][1]=0; PM[0][2]=0; PM[0][3]=0; PM[0][4]=1; PM[0][5]=1; PM[0][6]=1; PM[0][7]=1; PM[0][8]=1; PM[0][9]=1; PM[0][10]=1; PM[0][11]=1; PM[0][12]=1; PM[0][13]=1; PM[0][14]=1; PM[0][15]=1; PM[0][16]=1;
    //PM[1][0]=0; PM[1][1]=0; PM[1][2]=0; PM[1][3]=0; PM[1][4]=0; PM[1][5]=0; PM[1][6]=0; PM[1][7]=0; PM[1][8]=1; PM[1][9]=1; PM[1][10]=1; PM[1][11]=1; PM[1][12]=1; PM[1][13]=1; PM[1][14]=1; PM[1][15]=1; PM[1][16]=1;
    //PM[2][0]=0; PM[2][1]=0; PM[2][2]=0; PM[2][3]=0; PM[2][4]=0; PM[2][5]=0; PM[2][6]=0; PM[2][7]=0; PM[2][8]=0; PM[2][9]=0; PM[2][10]=0; PM[2][11]=0; PM[2][12]=1; PM[2][13]=1; PM[2][14]=1; PM[2][15]=1; PM[2][16]=1;
    //PM[3][0]=0; PM[3][1]=0; PM[3][2]=0; PM[3][3]=0; PM[3][4]=0; PM[3][5]=0; PM[3][6]=0; PM[3][7]=0; PM[3][8]=0; PM[3][9]=0; PM[3][10]=0; PM[3][11]=0; PM[3][12]=0; PM[3][13]=0; PM[3][14]=0; PM[3][15]=0; PM[3][16]=0;
    //PM[4][0]=0; PM[4][1]=0; PM[4][2]=0; PM[4][3]=0; PM[4][4]=0; PM[4][5]=0; PM[4][6]=0; PM[4][7]=0; PM[4][8]=0; PM[4][9]=0; PM[4][10]=0; PM[4][11]=0; PM[4][12]=0; PM[4][13]=0; PM[4][14]=0; PM[4][15]=0; PM[4][16]=0;
    //PM[5][0]=0; PM[5][1]=0; PM[5][2]=0; PM[5][3]=0; PM[5][4]=0; PM[5][5]=0; PM[5][6]=0; PM[5][7]=0; PM[5][8]=0; PM[5][9]=0; PM[5][10]=0; PM[5][11]=0; PM[5][12]=0; PM[5][13]=0; PM[5][14]=0; PM[5][15]=0; PM[5][16]=0;

    for (int i = 0 ; i < gamma ; i++)
    {
        for ( int j = 0 ; j < kappa ; j++ )
        {
          if ( PM[i][j] == -1 ) cout << "X  ";
          else cout << PM[i][j] << "  ";
        }
        cout << endl;
    }

    int counter_8 = object_enumeration_SC ( kappa , gamma , p , L , PM , CP , 4 );
    int counter_6 = object_enumeration_SC ( kappa , gamma , p , L , PM , CP , 3 );
    int counter_4 = object_enumeration_SC ( kappa , gamma , p , L , PM , CP , 2 );
    cout << "C4: " << counter_4 << endl << "C6: " << counter_6 << endl << "C8: " << counter_8 << endl;

    for ( int i = 0 ; i < gamma ; i++ ) delete[] PM[i]; delete[] PM;
    for ( int i = 0 ; i < gamma ; i++ ) delete[] CP[i]; delete[] CP;

    return 0;
}

int cycle_4_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP , int R1_limit )
{
    int counter = 0;

    for ( int R1 = 0 ; R1 < R1_limit ; R1++ )
    {
    for ( int j1 = 0 ; j1 < kappa ; j1++ )
    {
    int k1 = 0;
        for ( int i1 = 0 ; i1 < gamma ; i1++ )
        {
            if ( PM[i1][j1] > -1 )
            {
            for ( int R2 = R1 ; R2 < L ; R2++ )
            {
            int init_j2 = ( R2 == R1 ) ? ( j1 + 1 ) : 0;
            for ( int j2 = init_j2 ; j2 < kappa ; j2++ )
            {
            if ( ( j2 != j1 ) && ( PM[i1][j2] > -1 ) )
            {
            int k2 = ( (gamma * p ) + CP[i1][j1] - CP[i1][j2] + k1  ) % p;
            if ( ( ( R1 + PM[i1][j1] ) == ( R2 + PM[i1][j2] ) ) )
            {
                for ( int i2 = 0 ; i2 < gamma ; i2++ )
                {
                if ( ( i2 != i1 ) && ( PM[i2][j1] > -1 ) && ( PM[i2][j2] > -1 ) )
                {
                if ( ( ( R1 + PM[i2][j1] ) == ( R2 + PM[i2][j2] ) ) && ( ( ( CP[i2][j2] + k2 ) % p ) == ( ( CP[i2][j1] + k1 ) % p ) ) )
                {
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

    counter = counter * p; // considering different k1\in{0,1,...,p}

    return counter;
}

int cycle_6_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP , int R1_limit )
{
    int counter = 0;

    for ( int R1 = 0 ; R1 < R1_limit ; R1++ )
    {
    for ( int j1 = 0 ; j1 < kappa ; j1++ )
    {
    int k1 = 0;
    for ( int i1 = 0 ; i1 < gamma ; i1++ )
    {
        if ( PM[i1][j1] > -1 )
        {
        for ( int R2 = R1 ; R2 < L ; R2++ )
        {
        int init_j2 = ( R2 == R1 ) ? ( j1 + 1 ) : 0;
        for ( int j2 = init_j2 ; j2 < kappa ; j2++ )
        {
        if ( ( j2 != j1 ) && ( PM[i1][j2] > -1 ) )
        {
        int k2 = ( (gamma * p ) + CP[i1][j1] - CP[i1][j2] + k1  ) % p;
        if ( ( ( R1 + PM[i1][j1] ) == ( R2 + PM[i1][j2] ) ) )
        {
        for ( int i2 = 0 ; i2 < gamma ; i2++ )
        {
        if ( ( i2 != i1 ) && ( PM[i2][j2] > -1 ) )
        {
            for ( int R3 = R2 ; R3 < L ; R3++ )
            {
            int init_j3 = ( R3 == R2 ) ? ( j2 + 1 ) : 0;
            for ( int j3 = init_j3 ; j3 < kappa ; j3++ )
            {
            if ( ( j3 != j2 ) && ( j3 != j1 ) && ( PM[i2][j3] > -1 ) )
            {
            int k3 = ( (gamma * p ) + CP[i2][j2] - CP[i2][j3] + k2  ) % p;
            if ( ( ( R2 + PM[i2][j2] ) == ( R3 + PM[i2][j3] ) ) )
            {
            for ( int i3 = 0 ; i3 < gamma ; i3++ )
            {
            if ( ( i3 != i2 ) && ( i3 != i1 ) && ( PM[i3][j3] > -1 ) && ( PM[i3][j1] > -1 ) )
            {
            if ( ( ( R1 + PM[i3][j1] ) == ( R3 + PM[i3][j3] ) ) && ( ( ( CP[i3][j3] + k3 ) % p ) == ( ( CP[i3][j1] + k1 ) % p ) ) )
            {
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
        }
        }
        }
    }
    }
    }

    counter = counter * p; // considering different k1\in{0,1,...,p}

    return counter;
}


int cycle_8_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP , int R1_limit )
{
    int counter = 0;

    for ( int R1 = 0 ; R1 < R1_limit ; R1++ )
    {
    for ( int j1 = 0 ; j1 < kappa ; j1++ )
    {
    for ( int k1 = 0 ; k1 < p ; k1++ )
    {
    int v1_temp = ( j1 * p ) + k1;
    for ( int i1 = 0 ; i1 < gamma ; i1++ )
    {
    if ( PM[i1][j1] > -1 )
    {
    int c1_temp = ( i1 * p ) + ( ( CP[i1][j1] + k1 ) % p );
    int c1 = ( ( R1 + PM[i1][j1] ) * gamma * p ) + c1_temp;
    int v1 = ( R1 * kappa * p ) + v1_temp;
        for ( int R2 = R1 ; R2 < L ; R2++ )
        {
        int init_j2 = ( R2 == R1 ) ? ( j1 + 1 ) : 0;
        for ( int j2 = init_j2 ; j2 < kappa ; j2++ )
        {
        if ( ( j2 != j1 ) && ( PM[i1][j2] > -1 ) )
        {
        int k2 = ( (gamma * p ) + CP[i1][j1] - CP[i1][j2] + k1  ) % p;
        int v2_temp = ( j2 * p ) + k2;
        if ( ( ( R1 + PM[i1][j1] ) == ( R2 + PM[i1][j2] ) ) )
        {
        for ( int i2 = 0 ; i2 < gamma ; i2++ )
        {
        if ( ( i2 != i1 ) && ( PM[i2][j2] > -1 ) )
        {
        int c2_temp = ( i2 * p ) + ( ( CP[i2][j2] + k2 ) % p );
        int c2 = ( ( R2 + PM[i2][j2] ) * gamma * p ) + c2_temp;
        int v2 = ( R2 * kappa * p ) + v2_temp;
            for ( int R3 = R1 ; R3 < L ; R3++ )
            {
            int init_j3 = ( R3 == R1 ) ? j1 : 0;
            for ( int j3 = init_j3 ; j3 < kappa ; j3++ )
            {
            if ( PM[i2][j3] > -1 )
            {
            int k3 = ( (gamma * p ) + CP[i2][j2] - CP[i2][j3] + k2  ) % p;
            int v3_temp = ( j3 * p ) + k3;
            int v3 = ( R3 * kappa * p ) + v3_temp;
            if ( ( j3 != j2 ) && ( v3 > v1 ) )
            {
            if ( ( ( R2 + PM[i2][j2] ) == ( R3 + PM[i2][j3] ) ) )
            {
            for ( int i3 = 0 ; i3 < gamma ; i3++ )
            {
            if ( PM[i3][j3] > -1 )
            {
            int c3_temp = ( i3 * p ) + ( ( CP[i3][j3] + k3 ) % p );
            int c3 = ( ( R3 + PM[i3][j3] ) * gamma * p ) + c3_temp;
            if ( ( i3 != i2 ) && ( c3 != c1 ) )
            {
                for ( int R4 = R2 ; R4 < L ; R4++ )
                {
                int init_j4 = ( R4 == R2 ) ? j2 : 0;
                for ( int j4 = init_j4 ; j4 < kappa ; j4++ )
                {
                if ( PM[i3][j4] > -1 )
                {
                int k4 = ( (gamma * p ) + CP[i3][j3] - CP[i3][j4] + k3  ) % p;
                int v4_temp = ( j4 * p ) + k4;
                int v4 = ( R4 * kappa * p ) + v4_temp;
                if ( ( j4 != j3 ) && ( j4 != j1 ) && ( v4 > v2 ) )
                {
                if ( ( ( R3 + PM[i3][j3] ) == ( R4 + PM[i3][j4] ) ) )
                {
                for ( int i4 = 0 ; i4 < gamma ; i4++ )
                {
                if ( ( PM[i4][j4] > -1 ) && ( PM[i4][j1] > -1 ) )
                {
                int c4_temp = ( i4 * p ) + ( ( CP[i4][j4] + k4 ) % p );
                int c4 = ( ( R4 + PM[i4][j4] ) * gamma * p ) + c4_temp;
                if ( ( i4 != i3 ) && ( i4 != i1 ) && ( c4 != c2 ) )
                {
                if ( ( ( R1 + PM[i4][j1] ) == ( R4 + PM[i4][j4] ) ) && ( ( ( CP[i4][j4] + k4 ) % p ) == ( ( CP[i4][j1]+ k1 ) % p ) ) )
                {
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
    return counter;
}

int object_enumeration_SC ( int kappa , int gamma , int p , int L , int** PM , int** CP , int a )
{

    int span = 3;
    int counter; // total number of objects
    int counter1; // number of objects that span 1 replica
    int counter2; // number of objects that span 2 replicas
    int counter3; // number of objects that span 3 replicas

    if ( a == 4 )
    {
      counter1 = L * cycle_8_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter2 = (L-1) * cycle_8_enumeration_SC ( kappa , gamma , p , 2 , PM , CP , 1 );
      counter2 = counter2 + cycle_8_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter2 = counter2 - counter1;
      counter3 = (L-2) * cycle_8_enumeration_SC ( kappa , gamma , p , 3 , PM , CP , 1 );
      counter3 = counter3 + cycle_8_enumeration_SC ( kappa , gamma , p , 2 , PM , CP , 1 );
      counter3 = counter3 + cycle_8_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter3 = counter3 - counter2 - counter1;
      counter = (L-(span-1)) * cycle_8_enumeration_SC ( kappa , gamma , p , span , PM , CP , 1 );
      for ( int i = 1 ; i < span ; i++ )
      {
          counter = counter + cycle_8_enumeration_SC ( kappa , gamma , p , span - i , PM , CP , 1 );
      }
    }

    if ( a == 3 )
    {
      counter1 = L * cycle_6_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter2 = (L-1) * cycle_6_enumeration_SC ( kappa , gamma , p , 2 , PM , CP , 1 );
      counter2 = counter2 + cycle_6_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter2 = counter2 - counter1;
      counter3 = (L-2) * cycle_6_enumeration_SC ( kappa , gamma , p , 3 , PM , CP , 1 );
      counter3 = counter3 + cycle_6_enumeration_SC ( kappa , gamma , p , 2 , PM , CP , 1 );
      counter3 = counter3 + cycle_6_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter3 = counter3 - counter2 - counter1;
      counter = (L-(span-1)) * cycle_6_enumeration_SC ( kappa , gamma , p , span , PM , CP , 1 );
      for ( int i = 1 ; i < span ; i++ )
      {
        counter = counter + cycle_6_enumeration_SC ( kappa , gamma , p , span - i , PM , CP , 1 );
      }
    }

    if ( a == 2 )
    {
      counter1 = L * cycle_4_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter2 = (L-1) * cycle_4_enumeration_SC ( kappa , gamma , p , 2 , PM , CP , 1 );
      counter2 = counter2 + cycle_4_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter2 = counter2 - counter1;
      counter3 = (L-2) * cycle_4_enumeration_SC ( kappa , gamma , p , 3 , PM , CP , 1 );
      counter3 = counter3 + cycle_4_enumeration_SC ( kappa , gamma , p , 2 , PM , CP , 1 );
      counter3 = counter3 + cycle_4_enumeration_SC ( kappa , gamma , p , 1 , PM , CP , 1 );
      counter3 = counter3 - counter2 - counter1;
      counter = (L-(span-1)) * cycle_4_enumeration_SC ( kappa , gamma , p , span , PM , CP , 1 );
      for ( int i = 1 ; i < span ; i++ )
      {
        counter = counter + cycle_4_enumeration_SC ( kappa , gamma , p , span - i , PM , CP , 1 );
      }
    }
    cout << "Type 1 cycles = " << counter1 << ", Type 2 cycles = " << counter2 << ", Type 3 cycles = " << counter3 << endl;
    return counter;
}
