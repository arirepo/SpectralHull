
#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <stdio.h>
#include <stdlib.h>
#include "triangle.h"

/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/

void report(io, markers, reporttriangles, reportneighbors, reportsegments,
            reportedges, reportnorms)
struct triangulateio *io;
int markers;
int reporttriangles;
int reportneighbors;
int reportsegments;
int reportedges;
int reportnorms;
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}

/*****************************************************************************/
/*                                                                           */
/*  A wrapper for Triangle library                                           */
/*  Adopted by : ghasemi.arash@gmail.com                                     */
/*****************************************************************************/

void custom_Tri_wrapper(cmd_arg_string, numberofpoints, pointlist, /* <--- inputs */
		       numberofsegments, segmentlist,
		       segmentmarkerlist, numberofholes,
		       holelist, report_before, report_after,
		       default_numberofpoints, default_numberoftriangles, default_numberofsegments, 
     /*outputs --> */  new_numberofpoints, new_pointlist,
		       numberoftriangles, numberofcorners, trianglelist, neighborlist,
		       new_numberofsegments, new_segmentlist, new_segmentmarkerlist)

/* inputs */
char *cmd_arg_string;
int numberofpoints;
REAL *pointlist;
int numberofsegments;
int *segmentlist;
int *segmentmarkerlist;
int numberofholes;
REAL * holelist;
int report_before;
int report_after;
int default_numberofpoints;
int default_numberoftriangles;
int default_numberofsegments;

/* outputs */
int *new_numberofpoints;
REAL *new_pointlist;
int *numberoftriangles;
int *numberofcorners;
int *trianglelist;
int *neighborlist;
int *new_numberofsegments;
int *new_segmentlist;
int *new_segmentmarkerlist;

{
  /* local vars  */
  int i = 0;
  struct triangulateio in, mid, vorout;

  /* Define input points. */
  
  in.numberofpoints = numberofpoints;
  in.numberofpointattributes = 0;
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

  for( i = 0; i < (2*in.numberofpoints); i++)
       in.pointlist[i] = pointlist[i];

  in.pointattributelist = (REAL *) NULL;
  in.pointmarkerlist = (int *) NULL;

  in.numberofsegments = numberofsegments;
  in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  for ( i = 0; i < (in.numberofsegments * 2); i++)
       in.segmentlist[i] = segmentlist[i];

  in.segmentmarkerlist = (int *) malloc(in.numberofsegments * 1 * sizeof(int));
  for ( i = 0; i < (in.numberofsegments * 1); i++)
       in.segmentmarkerlist[i] = segmentmarkerlist[i];

  in.numberofholes = numberofholes;
  in.holelist = (REAL *) malloc(in.numberofholes * 2 * sizeof(REAL));
  for ( i = 0; i < (in.numberofholes * 2); i++)
       in.holelist[i] = holelist[i];
  
  in.numberofregions = 0;
  in.regionlist = (REAL *) NULL;

  /* /\* Needed only if -n switch used. *\/ */
  /* mid.neighborlist = (int *) malloc(default_numberoftriangles*3 *sizeof(int));          */
  mid.neighborlist = (int *)NULL;
  
  if ( report_before ){
       printf("the input cmdargs is ---> %s", cmd_arg_string);
       report(&in, 0, 0, 0, 1, 0, 0);
  }

  
  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */
  mid.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  mid.pointattributelist = (REAL *) NULL;
  mid.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  mid.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  mid.triangleattributelist = (REAL *) NULL;
  /* Needed only if segments are output (-p or -c) and -P not used: */
  mid.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  mid.segmentmarkerlist = (int *) NULL;
  mid.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  mid.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

  vorout.pointlist = (REAL *) NULL;        /* Needed only if -v switch used. */
  /* Needed only if -v switch used and number of attributes is not zero: */
  vorout.pointattributelist = (REAL *) NULL;
  vorout.edgelist = (int *) NULL;          /* Needed only if -v switch used. */
  vorout.normlist = (REAL *) NULL;         /* Needed only if -v switch used. */

  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), preserve the convex hull (c), number everything from  */
  /*   zero (z), assign a regional attribute to each element (A), and  */
  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  /*   neighbor list (n).                                              */

  /* triangulate("pczAevn", &in, &mid, &vorout); */
  triangulate(cmd_arg_string, &in, &mid, &vorout);  

  if ( report_after )
       report(&mid, 0, 1, 0, 1, 0, 0);

  /* check the consistency of the generated trimesh */
  /*      and any possible memory leak and segfaults */
  if ( default_numberofpoints < mid.numberofpoints ||
       (default_numberoftriangles * 3) < (mid.numberoftriangles * mid.numberofcorners) ||
       default_numberofsegments < mid.numberofsegments )
  {
       printf("\nFATAL error : the generated trimesh doesn't fit in the memory provided! increase default values for default_numberofpoints, default_numberoftriangles and default_numberofsegments. exit(0)\n");
       exit(0);
  }


  
  /* start extracting the generated triangles and boundary segments */
  /* and returning them */
  *new_numberofpoints = mid.numberofpoints;
  for(i = 0; i < (2 * mid.numberofpoints); i++)
       new_pointlist[i] = mid.pointlist[i];
  *numberoftriangles = mid.numberoftriangles;
  *numberofcorners = mid.numberofcorners;
  printf("\n\nnumber of corners : %d \n\n" ,mid.numberofcorners);
  
  for( i = 0; i < (mid.numberofcorners * mid.numberoftriangles); i++)
       trianglelist[i] = mid.trianglelist[i];
  /* printf("\n\n\n __________________________I'm here. OKKKK\n\n\n");       */
  for( i = 0; i < (3 * mid.numberoftriangles) ; i++)
       neighborlist[i] = mid.neighborlist[i];
  /* if( mid.neighborlist == (int *)NULL) printf("\nOh, trash ...\n"); */
  
  *new_numberofsegments = mid.numberofsegments;
  for ( i = 0; i < (mid.numberofsegments*2) ; i++)
       new_segmentlist[i] = mid.segmentlist[i];
  for ( i = 0; i < mid.numberofsegments ; i++)  
       new_segmentmarkerlist[i] = mid.segmentmarkerlist[i];


  
  /* Free all allocated arrays, including those allocated by Triangle. */
  /* free(in.pointlist); */
  /* free(in.pointattributelist); */
  /* free(in.pointmarkerlist); */
  /* free(in.regionlist); */
  /* free(mid.pointlist); */
  /* free(mid.pointattributelist); */
  /* free(mid.pointmarkerlist); */
  /* free(mid.trianglelist); */
  /* free(mid.triangleattributelist); */
  /* free(mid.trianglearealist); */
  /* free(mid.neighborlist); */
  /* free(mid.segmentlist); */
  /* free(mid.segmentmarkerlist); */
  /* free(mid.edgelist); */
  /* free(mid.edgemarkerlist); */
  /* free(vorout.pointlist); */
  /* free(vorout.pointattributelist); */
  /* free(vorout.edgelist); */
  /* free(vorout.normlist); */
  /* free(out.pointlist); */
  /* free(out.pointattributelist); */
  /* free(out.trianglelist); */
  /* free(out.triangleattributelist); */

}
