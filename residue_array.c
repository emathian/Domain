/*
 * File:	residue_array.c
 * Purpose:	Read PDB atom records into an array of "residue" structures.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_RESIDUES		1000
#define MAX_ATOMS_PER_RESIDUE	14
#define LINE_LENGTH     81


typedef struct {
	double	x, y, z;
} Point;


typedef struct {
	int	serial;
	char	atomName[5];
	char	altLoc[2];
	Point	centre;
} Atom;


typedef struct {
	int	numAtoms;
        char    resName[4];
        char    chainID[2];
        int     resSeq;
        char    iCode[2];
	Atom	atom[MAX_ATOMS_PER_RESIDUE+1];
} Residue;


/*
 * Declare an array to hold data read from the ATOM records of a PDB file.
 */

Residue	residue[MAX_RESIDUES+1];


int read_data(filename)
	char	*filename;
{
	FILE	*stream;
	char	line[LINE_LENGTH];
        char    s_serial[6];
        char    s_name[5];      /* Alpha-carbon is " CA "; calcium is "CA  " */
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

	int	i=0;
	int	j;
	int	previousSeq=0;
	char	previousID[2];
	char	previousICode[2];

        if ( (stream = fopen(filename, "r")) == NULL ) {
                (void) fprintf(stderr, "Unable to open %s\n", filename);
                exit(0);
        }

	strcpy(previousID, "");
	strcpy(previousICode, "");
        while ( fgets(line, LINE_LENGTH, stream) ) {
                if ( strncmp(line, "ATOM  ", 6) == 0 ) {

                        /*
                         * Split the line into its constituent fields.
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
			 * Copy values to the next element in the array.
			 */
			if ( ( resSeq != previousSeq ) ||
			     ( strcmp(s_chainID, previousID) != 0 ) ||
			     ( strcmp(s_iCode, previousICode) != 0 ) ) {
				if ( ++i > MAX_RESIDUES ) {
					(void) fprintf(stderr,
						"Too many residues\n");
					exit (0);
				}
				previousSeq = resSeq;
				strcpy(previousID, s_chainID);
				strcpy(previousICode, s_iCode);
				residue[i].numAtoms = 0;
				strcpy(residue[i].resName, s_resName);
				strcpy(residue[i].chainID, s_chainID);
				residue[i].resSeq = resSeq;
				strcpy(residue[i].iCode, s_iCode);
			}
			if ( ++residue[i].numAtoms > MAX_ATOMS_PER_RESIDUE ) {
				(void) fprintf(stderr,
					"Too many atoms in residue %d\n", i);
				exit (0);
			}
			j = residue[i].numAtoms;
			residue[i].atom[j].serial = serial;
			strcpy(residue[i].atom[j].atomName, s_name);
			strcpy(residue[i].atom[j].altLoc, s_altLoc);
			residue[i].atom[j].centre.x = x;
			residue[i].atom[j].centre.y = y;
			residue[i].atom[j].centre.z = z;
		}
	}
	return i;
}


void write_pdb_atom(serial, s_name, s_altLoc, s_resName, s_chainID,
		resSeq, s_iCode, centre)
	int	serial;
	char	*s_name;
	char	*s_altLoc;
	char	*s_resName;
	char	*s_chainID;
	int	resSeq;
	char	*s_iCode;
	Point	centre;
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
	int	argc;
	char	**argv;
{
	int	numResidues;
	int	i, j;

        if ( argc<2 ) {
                (void) fprintf(stderr, "usage: residue_array file.pdb\n");
                exit(0);
        }

	numResidues = read_data(argv[1]);
	for (i=1; i<=numResidues; ++i) {
		for (j=1; j<=residue[i].numAtoms; ++j) {
			write_pdb_atom(
				residue[i].atom[j].serial,
				residue[i].atom[j].atomName,
				residue[i].atom[j].altLoc,
				residue[i].resName,
				residue[i].chainID,
				residue[i].resSeq,
				residue[i].iCode,
				residue[i].atom[j].centre);
		}
	}

        return 0;
}
