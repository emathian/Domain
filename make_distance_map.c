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
    Point   centre;
} Atom;


char* substring(const char* str, size_t begin, size_t len);

double dist_atom2(struct Atom *atom1, struct Atom *atom2);
 

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
      for (i=1; i<=numAtoms; ++i) {
        for (j=1; j<=numAtoms; ++j){
            d= dist_atom2(&atom[i], &atom[j]);
            if (d <= threshold){
                printf("%d %d \n", atom[i].serial, atom[j].serial );
            }
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
