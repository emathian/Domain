/*
 * File:    atom_array.c
 * Purpose: Read PDB atom records into an array of "atom" structures.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ATOMS       10000
#define LINE_LENGTH     81


double MSV = 17.05;
int MSSE = 10;
int MSSM= 30;

typedef struct {
    double  x, y, z;
} Point;


typedef struct Atom{
    int serial;
    char    atomName[5];
    char    altLoc[2];
    char    resName[4];
    char    chainID[2];
    int resSeq;
    char    iCode[2];
    Point   centre;
} Atom;


char* substring(const char* str, size_t begin, size_t len);

double dist_atom2(struct Atom *atom1, struct Atom *atom2);
 
float score_split(struct Atom *seq, int x1, int y1, int x2, int y2, int len, double sc);


void  DOMAK_2Seg(struct Atom *seq, int *x1_MAXF, int *y1_MAXF, int *x2_MAXF, int *y2_MAXF,  int len);

/*
 * Declare an array to hold data read from the ATOM records of a PDB file.
 */

Atom    atom[MAX_ATOMS+1];


int read_data(filename)
    char    *filename;
{
    FILE    *stream;
    char    line[LINE_LENGTH];
        char    s_serial[6];
        char    s_name[5];      /* Alpha-carbon is " CA "; calcium is "CA  " */
        char    ss_name[5]; 
        char    s_altLoc[2];    /* Usually " " */
        char    s_resName[4];
        char    s_chainID[2];
        char    s_resSeq[5];
        char    s_iCode[2];     /* Usually " " */
        char    s_x[9];
        char    s_y[9];
        char    s_z[9];

        int     serial;
        int     resSeq;
        double  x;
        double  y;
        double  z;

    int i=0;

        if ( (stream = fopen(filename, "r")) == NULL ) {
                (void) fprintf(stderr, "Unable to open %s\n", filename);
                exit(0);
        }

        while ( fgets(line, LINE_LENGTH, stream) ) {
                if ( strncmp(line, "ATOM  ", 6) == 0 ){
            
                         strncpy(ss_name,    &line[12], 4); ss_name[4]    = '\0';
    
                        int t1 = 0;
                        for (int i = 1; i<3; i++){
                            if (i == 1 ){
                                char b = 'C';
                                if ( ss_name[i] ==  b){
                                 t1++;
                               }
                            }

                            if (i == 2 ){
                                char b = 'A';
                                if ( ss_name[i] ==  b){
                                 t1++;
                               }
                            }
                        }
                        if (t1 ==2){
                   
                        if ( strncmp(ss_name, "C", 2) == 0 ){
                            printf("OK \n");

                        }

        
                         /* Split the line into its constituent fields.
                         * We are only interested in columns 1-54.
                         */

                        strncpy(s_serial,  &line[6],  5); s_serial[5]  = '\0';
                        strncpy(s_name,    &line[12], 4); s_name[4]    = '\0';
                        strncpy(s_altLoc,  &line[16], 1); s_altLoc[1]  = '\0';
                        strncpy(s_resName, &line[17], 3); s_resName[3] = '\0';
                        strncpy(s_chainID, &line[21], 1); s_chainID[1] = '\0';
                        strncpy(s_resSeq,  &line[22], 4); s_resSeq[4]  = '\0';
                        strncpy(s_iCode,   &line[26], 1); s_iCode[1]   = '\0';
                        strncpy(s_x,       &line[30], 8); s_x[8]       = '\0';
                        strncpy(s_y,       &line[38], 8); s_y[8]       = '\0';
                        strncpy(s_z,       &line[46], 8); s_z[8]       = '\0';

                        /*
                         * Convert the numeric fields to integers or doubles.
                         * The library functions atoi() and atof() are
                         * described in the UNIX manual pages ('man atoi' and
                         * 'man atof').
                         */

                        serial = atoi(s_serial);
                        resSeq = atoi(s_resSeq);
                        x      = atof(s_x);
                        y      = atof(s_y);
                        z      = atof(s_z);

            /*
             * Copy values to the next element in the atom array.
             */

            if ( ++i > MAX_ATOMS ) {
                (void) fprintf(stderr, "Too many atoms read\n");
                exit(0);
            }
            atom[i].serial = serial;
            strcpy(atom[i].atomName, s_name);
            strcpy(atom[i].altLoc, s_altLoc);
            strcpy(atom[i].resName, s_resName);
            strcpy(atom[i].chainID, s_chainID);
            atom[i].resSeq = resSeq;
            strcpy(atom[i].iCode, s_iCode);
            atom[i].centre.x = x;
            atom[i].centre.y = y;
            atom[i].centre.z = z;
        }
      }
    }
    return i;
}


void write_pdb_atom(serial, s_name, s_altLoc, s_resName, s_chainID,
        resSeq, s_iCode, centre)
    int serial;
    char    *s_name;
    char    *s_altLoc;
    char    *s_resName;
    char    *s_chainID;
    int resSeq;
    char    *s_iCode;
    Point   centre;
{
    printf("ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f\n",
                                serial,
                                s_name,
                                s_altLoc,
                                s_resName,
                                s_chainID,
                                resSeq,
                                s_iCode,
                                centre.x,
                                centre.y,
                                centre.z);
}


int main(argc, argv)
    int argc;
    char    **argv;
{
    int numAtoms;
    //int i;
    //int j;
    double threshold = 8;
    //double d;

        if ( argc<2 ) {
                (void) fprintf(stderr, "usage: atom_array file.pdb\n");
                exit(0);
        }

    numAtoms = read_data(argv[1]);
    printf("%d\n", numAtoms );
    /*for (i=1; i<=numAtoms; ++i) {
        write_pdb_atom(
            atom[i].serial,
            atom[i].atomName,
            atom[i].altLoc,
             atom[i].resName,
            atom[i].chainID,
            atom[i].resSeq,
            atom[i].iCode,
            atom[i].centre);
    }
    */
    /*
    double dd;
    dd= dist_atom2(&atom[3], &atom[3]);
    printf("%f\n", dd);
      for (i=1; i<=numAtoms; ++i) {
        for (j=1; j<=numAtoms; ++j){
            d= dist_atom2(&atom[i], &atom[j]);
            if (d <= threshold){
                printf("%d %d \n", atom[i].serial, atom[j].serial );
            }
        }
        }
    */
        double sc1 = 10;
        double new_sc;
        int x1 = 0; int y1 = 50; int x2 =51; int y2 = 120;
        new_sc  = score_split(&atom[0], x1, y1, x2, y2, numAtoms, sc1);        
        printf("NEW SC %f\n",new_sc );

         int maxy2 = numAtoms; int maxx2 = numAtoms - MSSE + 1;  int maxy1 = numAtoms - MSSE; int maxx1 = numAtoms - MSSE - MSSM +1;
        int minx1 = 0; int miny1 = MSSE -1 ; int minx2 = MSSE; int miny2 = MSSE + MSSM -1;
        printf("X1minF %d , y1minF %d , x2minF  %d , y2_MinF %d\n",minx1 , miny1, minx2 , miny2 );

        printf("X1maxF %d , y1maxF %d , x2maxF  %d , y2_MAXF %d\n", maxx1, maxy1, maxx2 , maxy2 );
        DOMAK_2Seg(&atom[0], &minx1, &miny1, &minx2, &miny2, numAtoms);

        return 0;
}


char* substring(const char* str, size_t begin, size_t len) 
{ 
  if (str == 0 || strlen(str) == 0 || strlen(str) < begin || strlen(str) < (begin+len)) 
    return 0; 

  return strndup(str + begin, len); 
} 


double dist_atom2(struct Atom *atom1, struct Atom *atom2) 
{ 
  double d  = sqrt(pow(atom1->centre.x - atom2->centre.x,2) + pow(atom1->centre.y - atom2->centre.y,2) + pow(atom1->centre.z - atom2->centre.z,2));

  return d;
} 

float score_split(struct Atom *seq, int x1, int y1, int x2, int y2, int len, double sc){
    int i, j;
    int intA = 0; int intB = 0; int extAB = 0; int intA1 = 0; int intA2 = 0; int extA1A2 = 0;
    double D, corA1A2, score; 

    /* Interaction A1*/
    for(i=x1; i<=y1; i++){
        for(j=x1; j<=y1; j++){
            if(i != j){
                D= dist_atom2(&seq[i], &seq[j]);
                if (D < 5){
                    intA1++;
                }
            }
        }
    }

     /* Interaction A2*/
    for(i=x2; i<=y2; i++){
        for(j=x2; j<=y2; j++){
            if(i != j){
                D= dist_atom2(&seq[i], &seq[j]);
                if (D < 5){
                    intA2++;
                }
            }
        }
    }

     /* Interaction A1_A2*/
    for(i=x1; i<=y1; i++){
        for(j=x2; j<=y2; j++){
            if(i != j){
                D= dist_atom2(&seq[i], &seq[j]);
                if (D < 5){
                }
            }
        }
    }
    if (extA1A2 == 0){
        extA1A2 = 1;
    }

    corA1A2 = (intA1 / extA1A2) * (intA2 / extA1A2);
    //printf("corA1A2 %f\n", corA1A2);
    //if (corA1A2 < MSV){ // A1 and A2 are correlated

        if(x1 !=0){
            /* Interaction B1*/
            for(i=0; i<=x1; i++){
                for(j=0; j<=x1; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                             intB++;
                        }   
                    }
                }
            }

            /* Interaction B1A1*/
            for(i=0; i<=x1; i++){
                for(j=x1; j<=y1; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                            extAB++;
                        }   
                    }
                }
            }


            /* Interaction B1A2*/
            for(i=0; i<=x1; i++){
                for(j=x2; j<=y2; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                            extAB++;
                        }   
                    }
                }
            }
        }

        if(x2-y1 !=1){
            /* Interaction B2*/
            for(i=y1; i<=x2; i++){
                for(j=y1; j<=x2; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                             intB++;
                        }   
                    }
                }
            }



            /* Interaction B2A1*/
            for(i=y1; i<=x2; i++){
                for(j=x1; j<=y1; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                            extAB++;
                        }   
                    }
                }
            }

            /* Interaction B2A2*/
            for(i=y1; i<=x2; i++){
                for(j=x2; j<=y2; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                            extAB++;
                        }   
                    }
                }
            }
        }

        if(y2 != len){
            /* Interaction B3*/
            for(i=y2; i<=len; i++){
                for(j=y2; j<=len; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                             intB++;
                        }   
                    }
                }
            }

            /* Interaction B3A1*/
            for(i=y2; i<=len; i++){
                for(j=x1; j<=y1; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                            extAB++;
                        }   
                    }
                }
            }

            /* Interaction B2A2*/
            for(i=y2; i<=len; i++){
                for(j=x2; j<=y2; j++){
                    if(i != j){
                         D= dist_atom2(&seq[i], &seq[j]);
                        if (D < 5){
                            extAB++;
                        }   
                    }
                }
            }
        }
    intA = intA1 + intA2;
    if (extAB == 0){
        extAB = -1;
    }
    
    score = (intA/extAB) * (intB/extAB);
   // printf("intA %d , extAB %d, intB %d  \n", intA, extAB , intB);
    if (score > sc){
        
        return score;
    }
    else{
        return sc;
    }
  //  }
   // else{
        //printf("A1 A2 distinct \n");
    //    return sc;
    //}
}


void  DOMAK_2Seg(struct Atom *seq, int *x1_MAXF, int *y1_MAXF, int *x2_MAXF, int *y2_MAXF,  int len){
    int maxy2 = len; int maxx2 = len - MSSE + 1;  int maxy1 = len - MSSE; int maxx1 = len - MSSE - MSSM +1;
    int minx1 = 0; int miny1 = MSSE -1 ; int minx2 = MSSE; int miny2 = MSSE + MSSM -1;
    int x1_MAX; int y1_MAX; int x2_MAX; int y2_MAX;
    int x1; int y1; int x2; int y2;
    int sA1; int sA2; int sB1; int sB2; int sB3;
    int t = -1 ; int tmax = 1;
    double c_score;
    double last_score = -1;
    int Bool = 0;

    for(x1 =minx1; x1 <= maxx1; x1 = x1 + 3 ){
        printf("X1 %d\n", x1 );
        for(y1 = miny1; y1 <= maxy1; y1 = y1 + 3){
            printf("y1 %d\n", y1 );
            for(x2 = minx2; x2 <= maxx2; x2 = x2  +3  ){
                for(y2 = miny2; y2 <=maxy2; y2 = y2 +3){
             
                    sA1 = y1 - x1 +1;
                    sA2 = y2 - x2 +1;
                    sB1  = x1;
                    sB2  = x2 - y1 -1;
                    sB3 = len - y2;
                    if(x1 < y1 &&  y1 < x2 && x2 < y2 && sA1 >= MSSE && sA2 >= MSSE){

                       // printf(" x1 %d , y1 %d ,  x2 %d , y2 %d \n",  x1, y1 ,x2, y2);
                        if(x1==0){
                                
                            if(x2 == y1+1){ // Case n°1
                                if(sA1 >= MSSE && sA2 >= MSSM && sB3 >= MSSE ){
                                    t = 3;
                                    tmax = 3;
                                   
                                }
                            }
                            else{
                                if (y2 == len){ //  Case n°6
                                    if(sA1 >= MSSE && sA2 >= MSSE && sB2 >= MSSM){
                                        t =3;
                                        tmax =3;
                                   
                                    }
                                }
                                else{ // Case n°2
                                  if( sA1 >= MSSE && sB2 >= MSSM && sA2 >= MSSM && sB3 >= MSSE){
                                     t= 4;
                                     tmax = 4;
                                     
                                  }
                                }
                            }
                        }

                        else{
                            if( x2 == y1+1){
                                if(x1 == maxx1){
                                    if(sA1 >= MSSM && sA2 >= MSSE && sB1 >= MSSE){ // Case n°5
                                      
                                         t=3;
                                         tmax =3;
                                    }
                                }
                                else{
                                    if (sA1 >= MSSM && sA2 >= MSSM && sB1 >= MSSE && sB3 >= MSSE){ // Case n°7
                                         t= 4;
                                         tmax = 4;
                                    }
                                }
                            }
                            else{
                                if(y2 == len){
                                    if(sA1 >= MSSM && sA2 >= MSSE && sB1 >= MSSE && sB2 >= MSSM){ // Case n°4
                                         // printf("C4\n");
                                        t =4;
                                        tmax = 4;
                                    }
                                }
                                else{
                                    if (sA1 >= MSSM && sA2 >= MSSM && sB1 >= MSSE && sB2 >= MSSM && sB3 >= MSSE){ // Case n°3
                                        t = 5;
                                        tmax = 5;
                                    }
                                }                                
                            }
                        }


                        if(t/tmax ==1){
                            c_score = score_split(&atom[0], x1, y1, x2, y2, len, last_score);
                             //printf("c_score (%f) >= last_score(%f) \n", c_score, last_score); 
                            if (c_score > last_score){
                                x1_MAX = x1;
                                y1_MAX = y1;
                                x2_MAX = x2;
                                y2_MAX = y2;
                            
                              //printf(" last_score %f\n", last_score );
                                if (c_score > last_score){
                                    printf("c_score (%f) >= last_score(%f): X1_MAX %d , y1_MAX %d , x2_MAX %d, y2_MAX %d \n", c_score, last_score,x1_MAX, y1_MAX, x2_MAX, y2_MAX);
                                }
                                last_score =  c_score;
                            }
                        }
                    }
                }
            }
        }
    }
    if (last_score >= MSV){ // The domain are not correlated
        *x1_MAXF = x1_MAX;
        *y1_MAXF = y1_MAX;
        *x2_MAXF = x2_MAX;
        *y2_MAXF = y2_MAX; 
    }

}

