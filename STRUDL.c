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
    int class; // u = 0 , v =1
    int flag; // not flag  = 0, flag = 1
    Point   centre;
} Atom;


char* substring(const char* str, size_t begin, size_t len);

double dist_atom2(struct Atom *atom1, struct Atom *atom2);
 
int score(struct Atom *seq, const int numAtoms );

int opt_change(struct Atom *seq, const int numAtoms, const int eval_class, const int move_class, const int current_score);

 void STRUDL(struct Atom *seq, const int numAtoms);


int find_min_score(const int score_array[], const int k);

int find_min_scored(const double score_array[], const int k);
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
            atom[i].class = 0; //v
            atom[i].flag = 0; // Unflaged
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
    int i;
    int j;
    double threshold = 8;
    double d;

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
    double dd;
    dd= dist_atom2(&atom[3], &atom[3]);
    printf("%f\n", dd);
      /*for (i=1; i<=numAtoms; ++i) {
        for (j=1; j<=numAtoms; ++j){
            d= dist_atom2(&atom[i], &atom[j]);
            if (d <= threshold){
                printf("%d %d \n", atom[i].serial, atom[j].serial );
            }
        }
        }
        */
        int sc;
        for(int i=0; i <= numAtoms; i++){
        atom[i].class = 1; // v
    }
        printf("% d\n", numAtoms);
        //sc =score(&atom[0],  numAtoms );
       STRUDL(&atom[0], numAtoms);
        
        printf("U class \n");
        for (int i = 0; i < numAtoms; ++i)
        {
            if(atom[i].class == 0){
                printf("atom  : %d \n", i);
            }
        }

        printf("V class \n");
        for (int i = 0; i < numAtoms; ++i)
        {
            if(atom[i].class == 1){
                printf("atom  : %d \n", i);
            }
        }


    
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

int score(struct Atom *seq, const int numAtoms ){
    int c= 0;
   // printf("IN SCORE\n");
    for(int i=0; i< numAtoms -2; i++){
        for(int j=0; j<  numAtoms -2; j++){
           // printf(" i %d ,  j %d , class i % d , class j %d\n", i,j,seq[i].class  ,  seq[j].class);
            if(seq[i].class != seq[j].class){
              //printf("dist_atom2(&seq[i],  &seq[j] ) %f\n",dist_atom2(&seq[i],  &seq[j] ) );
                if (dist_atom2(&seq[i],  &seq[j] ) < 5){

                    c++;
                }
            }
        }
    }
    // printf("contact %d \n", c);
    return c;
}


int opt_change(struct Atom *seq, const int numAtoms, const int eval_class, const int move_class, const int current_score){

    int aug = numAtoms*100000000; 
    int pos = -1;
    int i = 0;
    int old_class;
     int c_aug;
    for(i=0; i< numAtoms; i++){

        old_class =seq[i].class;


        if (seq[i].class == eval_class && seq[i].flag == 0 )
        {
           seq[i].class  = move_class; 
           int new_score = score(&seq[0],  numAtoms );
            seq[i].class = old_class;
             c_aug = new_score - current_score;
            //printf("c_aug %d \n", c_aug);
           if (c_aug  < aug ){
            aug =  c_aug;
            pos= i;
           }
        
      }
    }
    if(pos == -1){

        printf("c_aug %d \n",c_aug );
    }
    return pos;
}



void STRUDL_knowK(struct Atom *seq, const int numAtoms, const int k_size){
    struct Atom cp_atom[MAX_ATOMS+1];
    int current_score;
    int N2 = numAtoms /2;
    int old_part[numAtoms];
    int min_score;
    for (int i=0; i < numAtoms; i++){
            cp_atom[i] = seq[i];
        }

    for (int i=0; i < numAtoms; i++){
        
         old_part[i] = cp_atom[i].class;
    }
    
    int Best_part[N2][numAtoms];
    int Best_score[N2];
    int kturn =0;
    for (int k=1; k<6; k++){   
        printf("k %d\n", k);   
        // Step 1:
 
        int start_pos = rand() % (numAtoms-1); // random start
        cp_atom[start_pos].class = 0; // u
        int nb_in_u = 1;
        while(nb_in_u < k_size){
            current_score= score(&cp_atom[0],  numAtoms );
            int new_pos_u =  opt_change(&cp_atom[0], numAtoms, 1, 0, current_score ); // Change from v to u
            cp_atom[new_pos_u].class  = 0;
            nb_in_u ++;
        }
        current_score = score(&cp_atom[0],  numAtoms );
        printf("Init score %d \n", current_score );

        // Step 2: 
        int c = 0;
        int score_array[numAtoms - k];
        int stop = 0;
        int n_flag =0;
        int stop2 = 0;
        while (stop == 0){
        int cu2 =0;


        int UV_mat[numAtoms - k][numAtoms];
        for (int i = 0; i < k; ++i)
        {
          for (int j = 0; j < numAtoms; ++j){
            UV_mat[i][j]  = 1;
            if(cp_atom[j].class == 0){
                UV_mat[i][j]  = 0;
            }
          }
        }

     
         while(stop2 ==0){  
            //printf("C %d\n ", c);
            //printf("k %d  \n", k);
            int pos_v_to_u =  opt_change(&cp_atom[0], numAtoms, 1, 0, current_score );
            int pos_u_to_v =  opt_change(&cp_atom[0], numAtoms, 0, 1, current_score );
            //printf("pos_u_to_v %d,  pos_v_to_u %d \n", pos_u_to_v, pos_v_to_u );
            for (int i = 0; i < numAtoms; ++i){
                UV_mat[c][i] = cp_atom[i].class ;
            }



            //printf("pos_u_to_v %d pos_v_to_u = %d \n", pos_u_to_v, pos_v_to_u );
            if(pos_v_to_u != -1){
                UV_mat[c][pos_v_to_u] =0;
                cp_atom[pos_v_to_u].class  = 0;
            }
            if(pos_u_to_v != -1 ){
                UV_mat[c][pos_u_to_v] =1;
                cp_atom[pos_u_to_v].class  = 1;
                cp_atom[pos_u_to_v].flag  = 1;
            }
            /*
            for (int i=0; i < numAtoms; i++){
                    if(cp_atom[i].flag ==1){
                        printf("i %d\n", i);
                    }
            }
            */
            score_array[c] = score(&cp_atom[0],  numAtoms );
        
            int count_V = 0;
            n_flag = 0;
            for (int i = 0; i < numAtoms; ++i)
             {
                if(cp_atom[i].class ==1 && cp_atom[i].flag ==1){
                    n_flag++;
                }
                if(cp_atom[i].class ==1 ){
                    count_V++;
                }
             }
             if(count_V==n_flag){
                stop2 =1;
             }
         c++;

          }
        
      
        for (int i=0; i < numAtoms; i++){
            cp_atom[i].flag = 0;
        }
        printf("END STEP 2\n \n \n");

        printf("CLASS IN WHILE WHILe\n");
        int cv =0;
        int cu = 0;
         for (int j = 0; j < numAtoms; ++j){
            if(cp_atom[j].class == 1){
                 cv++;
            }
            else{
                cu++;
            }  
            
         }
         printf("cu %d \n", cu);
         printf("cv %d \n", cv );

        int min_s = find_min_score(score_array,numAtoms - k );
  
        printf("\n");
        // Step 4
        printf("score_array[min_s] %d \n",score_array[min_s] );
        min_score  =score_array[min_s];

        if (min_score < current_score){
            for (int ii=0; ii < numAtoms; ii++){
                cp_atom[ii].class = UV_mat[min_s][ii];
            }
        
            current_score =  min_score;
            stop = 0;
        }
        else{
            stop =1;
         }
        }
       // printf("HERE\n");
        for (int i=0; i < numAtoms; i++){
            Best_part[kturn][i]= cp_atom[i].class;
            Best_score[kturn] = current_score;
        }
        //printf("CURRENT SCORE %d \n", current_score);



        for (int i=0; i < numAtoms; i++){
            cp_atom[i].class = old_part[i];
        }
        current_score = score(&cp_atom[0],  numAtoms );
        /*printf("Best_part\n" );
        for (int i = 0; i < numAtoms; ++i)
        {
            printf("%d",Best_part[kturn][i] );
        }*/

        kturn++;
    } 
    double best_best_score[N2];
    printf("\n\n\n\n");
    printf("Best_part \n");
    for (int i = 0; i < N2-5; ++i)
    {
        for (int j = 0; j < numAtoms; ++j)
        {
            printf(" % d ", Best_part[i][j]);
        }
        printf("\n");
          printf("\n");
    }
    for (int i = 0; i < N2-5; ++i)
    {
        for (int j=0; j < numAtoms; j++){
           cp_atom[j].class = Best_part[i][j];
        }
        best_best_score[i] =  score(&cp_atom[0],  numAtoms ) / i;
    }
    int min_min_s = find_min_scored(best_best_score, N2-5);
    for (int j=0; j < numAtoms; j++){
            cp_atom[j].class =  Best_part[min_min_s ][j];
    }

      for (int j=0; j < numAtoms; j++){
            seq[j].class = cp_atom[j].class;
            printf("j %d , class , %d \n", j, cp_atom[j].class);
    }
    
}


void STRUDL(struct Atom *seq, const int numAtoms){
    struct Atom cp_atom[MAX_ATOMS+1];
    int current_score;
    int N2 = numAtoms /2;
    int old_part[numAtoms];
    int min_score;
    for (int i=0; i < numAtoms; i++){
            cp_atom[i] = seq[i];
        }

    for (int i=0; i < numAtoms; i++){
        
         old_part[i] = cp_atom[i].class;
    }
    
    int Best_part[N2-5][numAtoms];
    int Best_score[N2];
    int kturn =0;
    for (int k=5; k<N2; k++){   
        printf("k %d\n", k);   
        // Step 1:
 
        int start_pos = rand() % (numAtoms-1); // random start
        cp_atom[start_pos].class = 0; // u
        int nb_in_u = 1;
        while(nb_in_u < k){
            current_score= score(&cp_atom[0],  numAtoms );
            int new_pos_u =  opt_change(&cp_atom[0], numAtoms, 1, 0, current_score ); // Change from v to u
            cp_atom[new_pos_u].class  = 0;
            nb_in_u ++;
        }
        current_score = score(&cp_atom[0],  numAtoms );


        // Step 2: 
        int c = 0;
        int score_array[numAtoms - k];
        int stop = 0;
        int n_flag =0;
        int stop2 = 0;
        while (stop == 0){
        int cu2 =0;


        int UV_mat[numAtoms - k][numAtoms];
        for (int i = 0; i < k; ++i)
        {
          for (int j = 0; j < numAtoms; ++j){
            UV_mat[i][j]  = 1;
            if(cp_atom[j].class == 0){
                UV_mat[i][j]  = 0;
            }
          }
        }

     
         while(stop2 ==0){  
       
            int pos_v_to_u =  opt_change(&cp_atom[0], numAtoms, 1, 0, current_score );
            int pos_u_to_v =  opt_change(&cp_atom[0], numAtoms, 0, 1, current_score );
            //printf("pos_u_to_v %d,  pos_v_to_u %d \n", pos_u_to_v, pos_v_to_u );
            for (int i = 0; i < numAtoms; ++i){
                UV_mat[c][i] = cp_atom[i].class ;
            }

            if(pos_v_to_u != -1){
                UV_mat[c][pos_v_to_u] =0;
                cp_atom[pos_v_to_u].class  = 0;
            }
            if(pos_u_to_v != -1 ){
                UV_mat[c][pos_u_to_v] =1;
                cp_atom[pos_u_to_v].class  = 1;
                cp_atom[pos_u_to_v].flag  = 1;
            }

            score_array[c] = score(&cp_atom[0],  numAtoms );
        
            int count_V = 0;
            n_flag = 0;
            for (int i = 0; i < numAtoms; ++i)
             {
                if(cp_atom[i].class ==1 && cp_atom[i].flag ==1){
                    n_flag++;
                }
                if(cp_atom[i].class ==1 ){
                    count_V++;
                }
             }
             if(count_V==n_flag){
                stop2 =1;
             }
         c++;

          }
        
      
        for (int i=0; i < numAtoms; i++){
            cp_atom[i].flag = 0;
        }

        int cv =0;
        int cu = 0;
         for (int j = 0; j < numAtoms; ++j){
            if(cp_atom[j].class == 1){
                 cv++;
            }
            else{
                cu++;
            }  
            
         }

        int min_s = find_min_score(score_array,numAtoms - k );
  
        printf("\n");
        // Step 4
        printf("score_array[min_s] %d \n",score_array[min_s] );
        min_score  =score_array[min_s];

        if (min_score < current_score){
            for (int ii=0; ii < numAtoms; ii++){
                cp_atom[ii].class = UV_mat[min_s][ii];
            }
        
            current_score =  min_score;
            stop = 0;
        }
        else{
            stop =1;
         }
        }

        for (int i=0; i < numAtoms; i++){
            Best_part[kturn][i]= cp_atom[i].class;
            Best_score[kturn] = current_score;
        }

        for (int i=0; i < numAtoms; i++){
            cp_atom[i].class = old_part[i];
        }
        current_score = score(&cp_atom[0],  numAtoms );

        kturn++;
    } 
    double best_best_score[N2-5];
    

    for (int i = 0; i < N2-5; ++i)
    {
        printf("i  %d \n", i );
        for (int j=0; j < numAtoms; j++){
            printf("HERE\n");
           cp_atom[j].class = Best_part[i][j];
        }
        best_best_score[i] =  score(&cp_atom[0],  numAtoms ) / i;
        printf("HEREEEE\n");
    }
    int min_min_s = find_min_scored(best_best_score, N2);
    printf("min_min_s %d \n", min_min_s);
    for (int j=0; j < numAtoms; j++){
            cp_atom[j].class =  Best_part[min_min_s ][j];
    }

      for (int j=0; j < numAtoms; j++){
            seq[j].class = cp_atom[j].class;
            printf("j %d , class , %d \n", j, cp_atom[j].class);
    }
    
}


int find_min_score(const int score_array[], const int k){
    int min = 10000;
    int min_pos = -1;
    for (int i = 0; i < k; ++i)
    {
        if ( score_array[i] < min)
        {
            min = score_array[i];
            min_pos = i;
        }
    }
    return min_pos;
}


int find_min_scored(const double score_array[], const int k){
    int min = 10000;
    int min_pos = -1;
    for (int i = 0; i < k; ++i)
    {
        if ( score_array[i] < min)
        {
            min = score_array[i];
            min_pos = i;
        }
    }
    return min_pos;

   } 