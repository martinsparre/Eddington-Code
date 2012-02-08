#include "gadget2conv.h"
#include<stdio.h>
#include<stdlib.h>
	
struct io_header header;

void WriteGadget2(struct Particle *P1, int Ndm, char *Name)
{
	int i;

	FILE *OutputFile;
	
	OutputFile = fopen(Name,"wb");

	float *Mass;
	unsigned int *ID;
	
	 float (*Pos)[3] =
		( float (*)[3])malloc(Ndm * sizeof(*Pos));

	float (*Vel)[3] =
		( float (*)[3])malloc(Ndm * sizeof(*Vel));	
		
	Mass=malloc(Ndm*sizeof(float));
	ID=malloc(Ndm*sizeof(unsigned int));

	for(i=0;i<Ndm;i++)
	{
		Pos[i][0] = P1[i].Pos[0];
		Pos[i][1] = P1[i].Pos[1];
		Pos[i][2] = P1[i].Pos[2];
		
		Vel[i][0] = P1[i].Vel[0];
		Vel[i][1] = P1[i].Vel[1];
		Vel[i][2] = P1[i].Vel[2];
		Mass[i] = P1[i].Mass;

	}
	
	for(i=0;i<Ndm;i++)
	{

		ID[i] = P1[i].ID;	
	}
	
	int blklen; // Block length
	#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, OutputFile);
	
	//Ndm=20000;//atoi(argv[3]);
	
	// INITIAL DEFINITIONS
	header.npart[0] = 0; //SPH
	header.npart[1] = Ndm; //Halo
	header.npart[2] = 0; //Disk
	header.npart[3] = 0; //Bulge
	header.npart[4] = 0; //Stars
	header.npart[5] = 0; //Boundary
	
	for (i=0;i<6;i++) {
		header.mass[i] = 0;
		header.npartTotal[i] = header.npart[i];
	}
	
	header.time = 0;
	header.flag_sfr = 0;
	header.redshift = 0;
	header.flag_feedback = 0;
	header.HubbleParam = 0;
	
	//header.num_files = 1;
	

	
	blklen = sizeof(header);
	BLKLEN;
	fwrite(&header, sizeof(header),1, OutputFile);
	BLKLEN;
	
	blklen = sizeof(float)*3*Ndm;
	BLKLEN;
	fwrite(Pos, sizeof(float),3*Ndm, OutputFile);
	BLKLEN;

	BLKLEN;
	fwrite(Vel, sizeof(float),3*Ndm, OutputFile);
	BLKLEN;
	
	blklen = sizeof(unsigned int)*Ndm;
	BLKLEN;
	fwrite(ID, sizeof(unsigned int),Ndm, OutputFile);
	BLKLEN;
	
	blklen = sizeof(float)*Ndm;
	BLKLEN;
	fwrite(Mass, sizeof(float),Ndm, OutputFile);
	BLKLEN;
	
	fclose(OutputFile);
}
	







