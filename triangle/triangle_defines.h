#pragma once

/*
 * Because the code for triunsuitable got too large to comfortably fit into the main triangle code,
 * as well as having multiple versions of the tolerance, this file now contains all the structure information
 * that can be included into multiple files
 */


/* For single precision (which will save some memory and reduce paging),     */
/*   define the symbol SINGLE by using the -DSINGLE compiler switch or by    */
/*   writing "#define SINGLE" below.                                         */
/*                                                                           */
/* For double precision (which will allow you to refine meshes to a smaller  */
/*   edge length), leave SINGLE undefined.                                   */
/*                                                                           */
/* Double precision uses more memory, but improves the resolution of the     */
/*   meshes you can generate with Triangle.  It also reduces the likelihood  */
/*   of a floating exception due to overflow.  Finally, it is much faster    */
/*   than single precision on 64-bit architectures like the DEC Alpha.  I    */
/*   recommend double precision unless you want to generate a mesh for which */
/*   you do not have enough memory.                                          */

/* #define SINGLE */

#include <gdal.h>
#include "ogr_srs_api.h"
#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

/* If yours is not a Unix system, define the NO_TIMER compiler switch to     */
/*   remove the Unix-specific timing code.                                   */

/* #define NO_TIMER */

/* To insert lots of self-checks for internal errors, define the SELF_CHECK  */
/*   symbol.  This will slow down the program significantly.  It is best to  */
/*   define the symbol using the -DSELF_CHECK compiler switch, but you could */
/*   write "#define SELF_CHECK" below.  If you are modifying this code, I    */
/*   recommend you turn self-checks on until your work is debugged.          */

/* #define SELF_CHECK */

/* To compile Triangle as a callable object library (triangle.o), define the */
/*   TRILIBRARY symbol.  Read the file triangle.h for details on how to call */
/*   the procedure triangulate() that results.                               */

/* #define TRILIBRARY */

/* It is possible to generate a smaller version of Triangle using one or     */
/*   both of the following symbols.  Define the REDUCED symbol to eliminate  */
/*   all features that are primarily of research interest; specifically, the */
/*   -i, -F, -s, and -C switches.  Define the CDT_ONLY symbol to eliminate   */
/*   all meshing algorithms above and beyond constrained Delaunay            */
/*   triangulation; specifically, the -r, -q, -a, -u, -D, -S, and -s         */
/*   switches.  These reductions are most likely to be useful when           */
/*   generating an object library (triangle.o) by defining the TRILIBRARY    */
/*   symbol.                                                                 */

/* #define REDUCED */
/* #define CDT_ONLY */

/* On some machines, my exact arithmetic routines might be defeated by the   */
/*   use of internal extended precision floating-point registers.  The best  */
/*   way to solve this problem is to set the floating-point registers to use */
/*   single or double precision internally.  On 80x86 processors, this may   */
/*   be accomplished by setting the CPU86 symbol for the Microsoft C         */
/*   compiler, or the LINUX symbol for the gcc compiler running on Linux.    */
/*                                                                           */
/* An inferior solution is to declare certain values as `volatile', thus     */
/*   forcing them to be stored to memory and rounded off.  Unfortunately,    */
/*   this solution might slow Triangle down quite a bit.  To use volatile    */
/*   values, write "#define INEXACT volatile" below.  Normally, however,     */
/*   INEXACT should be defined to be nothing.  ("#define INEXACT".)          */
/*                                                                           */
/* For more discussion, see http://www.cs.cmu.edu/~quake/robust.pc.html .    */
/*   For yet more discussion, see Section 5 of my paper, "Adaptive Precision */
/*   Floating-Point Arithmetic and Fast Robust Geometric Predicates" (also   */
/*   available as Section 6.6 of my dissertation).                           */

/* #define CPU86 */
/* #define LINUX */

#define INEXACT /* Nothing */
/* #define INEXACT volatile */

/* Maximum number of characters in a file name (including the null).         */

#define FILENAMESIZE 2048

/* Maximum number of characters in a line read from a file (including the    */
/*   null).                                                                  */

#define INPUTLINESIZE 1024

/* For efficiency, a variety of data structures are allocated in bulk.  The  */
/*   following constants determine how many of each structure is allocated   */
/*   at once.                                                                */

#define TRIPERBLOCK 4092           /* Number of triangles allocated at once. */
#define SUBSEGPERBLOCK 508       /* Number of subsegments allocated at once. */
#define VERTEXPERBLOCK 4092         /* Number of vertices allocated at once. */
#define VIRUSPERBLOCK 1020   /* Number of virus triangles allocated at once. */
/* Number of encroached subsegments allocated at once. */
#define BADSUBSEGPERBLOCK 252
/* Number of skinny triangles allocated at once. */
#define BADTRIPERBLOCK 4092
/* Number of flipped triangles allocated at once. */
#define FLIPSTACKERPERBLOCK 252
/* Number of splay tree nodes allocated at once. */
#define SPLAYNODEPERBLOCK 508

/* The vertex types.   A DEADVERTEX has been deleted entirely.  An           */
/*   UNDEADVERTEX is not part of the mesh, but is written to the output      */
/*   .node file and affects the node indexing in the other output files.     */

#define INPUTVERTEX 0
#define SEGMENTVERTEX 1
#define FREEVERTEX 2
#define DEADVERTEX -32768
#define UNDEADVERTEX -32767

/* The next line is used to outsmart some very stupid compilers.  If your    */
/*   compiler is smarter, feel free to replace the "int" with "void".        */
/*   Not that it matters.                                                    */

#define VOID int

/* Two constants for algorithms based on random sampling.  Both constants    */
/*   have been chosen empirically to optimize their respective algorithms.   */

/* Used for the point location scheme of Mucke, Saias, and Zhu, to decide    */
/*   how large a random sample of triangles to inspect.                      */

#define SAMPLEFACTOR 11

/* Used in Fortune's sweepline Delaunay algorithm to determine what fraction */
/*   of boundary edges should be maintained in the splay tree for point      */
/*   location on the front.                                                  */

#define SAMPLERATE 10

/* A number that speaks for itself, every kissable digit.                    */

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

/* Another fave.                                                             */

#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732

/* And here's one for those of you who are intimidated by math.              */

#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef NO_TIMER

#include <sys/time.h>
#include <ogr_core.h>

#endif /* not NO_TIMER */
#ifdef CPU86
#include <float.h>
#endif /* CPU86 */
#ifdef LINUX
//#include <fpu_control.h>
#endif /* LINUX */
#ifdef TRILIBRARY
#include "triangle.h"
#endif /* TRILIBRARY */

/* A few forward declarations.                                               */

#ifndef TRILIBRARY

char *readline();

char *findfield();

#endif /* not TRILIBRARY */

/* Labels that signify the result of point location.  The result of a        */
/*   search indicates that the point falls in the interior of a triangle, on */
/*   an edge, on a vertex, or outside the mesh.                              */

enum locateresult
{
    INTRIANGLE, ONEDGE, ONVERTEX, OUTSIDE
};

/* Labels that signify the result of vertex insertion.  The result indicates */
/*   that the vertex was inserted with complete success, was inserted but    */
/*   encroaches upon a subsegment, was not inserted because it lies on a     */
/*   segment, or was not inserted because another vertex occupies the same   */
/*   location.                                                               */

enum insertvertexresult
{
    SUCCESSFULVERTEX, ENCROACHINGVERTEX, VIOLATINGVERTEX,
    DUPLICATEVERTEX
};

/* Labels that signify the result of direction finding.  The result          */
/*   indicates that a segment connecting the two query points falls within   */
/*   the direction triangle, along the left edge of the direction triangle,  */
/*   or along the right edge of the direction triangle.                      */

enum finddirectionresult
{
    WITHIN, LEFTCOLLINEAR, RIGHTCOLLINEAR
};

/*****************************************************************************/
/*                                                                           */
/*  The basic mesh data structures                                           */
/*                                                                           */
/*  There are three:  vertices, triangles, and subsegments (abbreviated      */
/*  `subseg').  These three data structures, linked by pointers, comprise    */
/*  the mesh.  A vertex simply represents a mesh vertex and its properties.  */
/*  A triangle is a triangle.  A subsegment is a special data structure used */
/*  to represent an impenetrable edge of the mesh (perhaps on the outer      */
/*  boundary, on the boundary of a hole, or part of an internal boundary     */
/*  separating two triangulated regions).  Subsegments represent boundaries, */
/*  defined by the user, that triangles may not lie across.                  */
/*                                                                           */
/*  A triangle consists of a list of three vertices, a list of three         */
/*  adjoining triangles, a list of three adjoining subsegments (when         */
/*  segments exist), an arbitrary number of optional user-defined            */
/*  floating-point attributes, and an optional area constraint.  The latter  */
/*  is an upper bound on the permissible area of each triangle in a region,  */
/*  used for mesh refinement.                                                */
/*                                                                           */
/*  For a triangle on a boundary of the mesh, some or all of the neighboring */
/*  triangles may not be present.  For a triangle in the interior of the     */
/*  mesh, often no neighboring subsegments are present.  Such absent         */
/*  triangles and subsegments are never represented by NULL pointers; they   */
/*  are represented by two special records:  `dummytri', the triangle that   */
/*  fills "outer space", and `dummysub', the omnipresent subsegment.         */
/*  `dummytri' and `dummysub' are used for several reasons; for instance,    */
/*  they can be dereferenced and their contents examined without violating   */
/*  protected memory.                                                        */
/*                                                                           */
/*  However, it is important to understand that a triangle includes other    */
/*  information as well.  The pointers to adjoining vertices, triangles, and */
/*  subsegments are ordered in a way that indicates their geometric relation */
/*  to each other.  Furthermore, each of these pointers contains orientation */
/*  information.  Each pointer to an adjoining triangle indicates which face */
/*  of that triangle is contacted.  Similarly, each pointer to an adjoining  */
/*  subsegment indicates which side of that subsegment is contacted, and how */
/*  the subsegment is oriented relative to the triangle.                     */
/*                                                                           */
/*  The data structure representing a subsegment may be thought to be        */
/*  abutting the edge of one or two triangle data structures:  either        */
/*  sandwiched between two triangles, or resting against one triangle on an  */
/*  exterior boundary or hole boundary.                                      */
/*                                                                           */
/*  A subsegment consists of a list of four vertices--the vertices of the    */
/*  subsegment, and the vertices of the segment it is a part of--a list of   */
/*  two adjoining subsegments, and a list of two adjoining triangles.  One   */
/*  of the two adjoining triangles may not be present (though there should   */
/*  always be one), and neighboring subsegments might not be present.        */
/*  Subsegments also store a user-defined integer "boundary marker".         */
/*  Typically, this integer is used to indicate what boundary conditions are */
/*  to be applied at that location in a finite element simulation.           */
/*                                                                           */
/*  Like triangles, subsegments maintain information about the relative      */
/*  orientation of neighboring objects.                                      */
/*                                                                           */
/*  Vertices are relatively simple.  A vertex is a list of floating-point    */
/*  numbers, starting with the x, and y coordinates, followed by an          */
/*  arbitrary number of optional user-defined floating-point attributes,     */
/*  followed by an integer boundary marker.  During the segment insertion    */
/*  phase, there is also a pointer from each vertex to a triangle that may   */
/*  contain it.  Each pointer is not always correct, but when one is, it     */
/*  speeds up segment insertion.  These pointers are assigned values once    */
/*  at the beginning of the segment insertion phase, and are not used or     */
/*  updated except during this phase.  Edge flipping during segment          */
/*  insertion will render some of them incorrect.  Hence, don't rely upon    */
/*  them for anything.                                                       */
/*                                                                           */
/*  Other than the exception mentioned above, vertices have no information   */
/*  about what triangles, subfacets, or subsegments they are linked to.      */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  Handles                                                                  */
/*                                                                           */
/*  The oriented triangle (`otri') and oriented subsegment (`osub') data     */
/*  structures defined below do not themselves store any part of the mesh.   */
/*  The mesh itself is made of `triangle's, `subseg's, and `vertex's.        */
/*                                                                           */
/*  Oriented triangles and oriented subsegments will usually be referred to  */
/*  as "handles."  A handle is essentially a pointer into the mesh; it       */
/*  allows you to "hold" one particular part of the mesh.  Handles are used  */
/*  to specify the regions in which one is traversing and modifying the mesh.*/
/*  A single `triangle' may be held by many handles, or none at all.  (The   */
/*  latter case is not a memory leak, because the triangle is still          */
/*  connected to other triangles in the mesh.)                               */
/*                                                                           */
/*  An `otri' is a handle that holds a triangle.  It holds a specific edge   */
/*  of the triangle.  An `osub' is a handle that holds a subsegment.  It     */
/*  holds either the left or right side of the subsegment.                   */
/*                                                                           */
/*  Navigation about the mesh is accomplished through a set of mesh          */
/*  manipulation primitives, further below.  Many of these primitives take   */
/*  a handle and produce a new handle that holds the mesh near the first     */
/*  handle.  Other primitives take two handles and glue the corresponding    */
/*  parts of the mesh together.  The orientation of the handles is           */
/*  important.  For instance, when two triangles are glued together by the   */
/*  bond() primitive, they are glued at the edges on which the handles lie.  */
/*                                                                           */
/*  Because vertices have no information about which triangles they are      */
/*  attached to, I commonly represent a vertex by use of a handle whose      */
/*  origin is the vertex.  A single handle can simultaneously represent a    */
/*  triangle, an edge, and a vertex.                                         */
/*                                                                           */
/*****************************************************************************/

/* The triangle data structure.  Each triangle contains three pointers to    */
/*   adjoining triangles, plus three pointers to vertices, plus three        */
/*   pointers to subsegments (declared below; these pointers are usually     */
/*   `dummysub').  It may or may not also contain user-defined attributes    */
/*   and/or a floating-point "area constraint."  It may also contain extra   */
/*   pointers for nodes, when the user asks for high-order elements.         */
/*   Because the size and structure of a `triangle' is not decided until     */
/*   runtime, I haven't simply declared the type `triangle' as a struct.     */

typedef REAL **triangle;            /* Really:  typedef triangle *triangle   */

/* An oriented triangle:  includes a pointer to a triangle and orientation.  */
/*   The orientation denotes an edge of the triangle.  Hence, there are      */
/*   three possible orientations.  By convention, each edge always points    */
/*   counterclockwise about the corresponding triangle.                      */

struct otri
{
    triangle *tri;
    int orient;                                         /* Ranges from 0 to 2. */
};

/* The subsegment data structure.  Each subsegment contains two pointers to  */
/*   adjoining subsegments, plus four pointers to vertices, plus two         */
/*   pointers to adjoining triangles, plus one boundary marker, plus one     */
/*   segment number.                                                         */

typedef REAL **subseg;                  /* Really:  typedef subseg *subseg   */

/* An oriented subsegment:  includes a pointer to a subsegment and an        */
/*   orientation.  The orientation denotes a side of the edge.  Hence, there */
/*   are two possible orientations.  By convention, the edge is always       */
/*   directed so that the "side" denoted is the right side of the edge.      */

struct osub
{
    subseg *ss;
    int ssorient;                                       /* Ranges from 0 to 1. */
};

/* The vertex data structure.  Each vertex is actually an array of REALs.    */
/*   The number of REALs is unknown until runtime.  An integer boundary      */
/*   marker, and sometimes a pointer to a triangle, is appended after the    */
/*   REALs.                                                                  */

typedef REAL *vertex;

/* A queue used to store encroached subsegments.  Each subsegment's vertices */
/*   are stored so that we can check whether a subsegment is still the same. */

struct badsubseg
{
    subseg encsubseg;                             /* An encroached subsegment. */
    vertex subsegorg, subsegdest;                         /* Its two vertices. */
};

/* A queue used to store bad triangles.  The key is the square of the cosine */
/*   of the smallest angle of the triangle.  Each triangle's vertices are    */
/*   stored so that one can check whether a triangle is still the same.      */

struct badtriang
{
    triangle poortri;                       /* A skinny or too-large triangle. */
    REAL key;                             /* cos^2 of smallest (apical) angle. */
    vertex triangorg, triangdest, triangapex;           /* Its three vertices. */
    struct badtriang *nexttriang;             /* Pointer to next bad triangle. */
};

/* A stack of triangles flipped during the most recent vertex insertion.     */
/*   The stack is used to undo the vertex insertion if the vertex encroaches */
/*   upon a subsegment.                                                      */

struct flipstacker
{
    triangle flippedtri;                       /* A recently flipped triangle. */
    struct flipstacker *prevflip;               /* Previous flip in the stack. */
};

/* A node in a heap used to store events for the sweepline Delaunay          */
/*   algorithm.  Nodes do not point directly to their parents or children in */
/*   the heap.  Instead, each node knows its position in the heap, and can   */
/*   look up its parent and children in a separate array.  The `eventptr'    */
/*   points either to a `vertex' or to a triangle (in encoded format, so     */
/*   that an orientation is included).  In the latter case, the origin of    */
/*   the oriented triangle is the apex of a "circle event" of the sweepline  */
/*   algorithm.  To distinguish site events from circle events, all circle   */
/*   events are given an invalid (smaller than `xmin') x-coordinate `xkey'.  */

struct event
{
    REAL xkey, ykey;                              /* Coordinates of the event. */
    VOID *eventptr;      /* Can be a vertex or the location of a circle event. */
    int heapposition;              /* Marks this event's position in the heap. */
};

/* A node in the splay tree.  Each node holds an oriented ghost triangle     */
/*   that represents a boundary edge of the growing triangulation.  When a   */
/*   circle event covers two boundary edges with a triangle, so that they    */
/*   are no longer boundary edges, those edges are not immediately deleted   */
/*   from the tree; rather, they are lazily deleted when they are next       */
/*   encountered.  (Since only a random sample of boundary edges are kept    */
/*   in the tree, lazy deletion is faster.)  `keydest' is used to verify     */
/*   that a triangle is still the same as when it entered the splay tree; if */
/*   it has been rotated (due to a circle event), it no longer represents a  */
/*   boundary edge and should be deleted.                                    */

struct splaynode
{
    struct otri keyedge;                     /* Lprev of an edge on the front. */
    vertex keydest;           /* Used to verify that splay node is still live. */
    struct splaynode *lchild, *rchild;              /* Children in splay tree. */
};

/* A type used to allocate memory.  firstblock is the first block of items.  */
/*   nowblock is the block from which items are currently being allocated.   */
/*   nextitem points to the next slab of free memory for an item.            */
/*   deaditemstack is the head of a linked list (stack) of deallocated items */
/*   that can be recycled.  unallocateditems is the number of items that     */
/*   remain to be allocated from nowblock.                                   */
/*                                                                           */
/* Traversal is the process of walking through the entire list of items, and */
/*   is separate from allocation.  Note that a traversal will visit items on */
/*   the "deaditemstack" stack as well as live items.  pathblock points to   */
/*   the block currently being traversed.  pathitem points to the next item  */
/*   to be traversed.  pathitemsleft is the number of items that remain to   */
/*   be traversed in pathblock.                                              */
/*                                                                           */
/* alignbytes determines how new records should be aligned in memory.        */
/*   itembytes is the length of a record in bytes (after rounding up).       */
/*   itemsperblock is the number of items allocated at once in a single      */
/*   block.  itemsfirstblock is the number of items in the first block,      */
/*   which can vary from the others.  items is the number of currently       */
/*   allocated items.  maxitems is the maximum number of items that have     */
/*   been allocated at once; it is the current number of items plus the      */
/*   number of records kept on deaditemstack.                                */

struct memorypool
{
    VOID **firstblock, **nowblock;
    VOID *nextitem;
    VOID *deaditemstack;
    VOID **pathblock;
    VOID *pathitem;
    int alignbytes;
    int itembytes;
    int itemsperblock;
    int itemsfirstblock;
    long items, maxitems;
    int unallocateditems;
    int pathitemsleft;
};


/* Global constants.                                                         */

REAL splitter;       /* Used to split REAL factors for exact multiplication. */
REAL epsilon;                             /* Floating-point machine epsilon. */
REAL resulterrbound;
REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
REAL iccerrboundA, iccerrboundB, iccerrboundC;
REAL o3derrboundA, o3derrboundB, o3derrboundC;

/* Random number seed is not constant, but I've made it global anyway.       */

unsigned long randomseed;                     /* Current random number seed. */


/* Mesh data structure.  Triangle operates on only one mesh, but the mesh    */
/*   structure is used (instead of global variables) to allow reentrancy.    */

struct mesh
{

/* Variables used to allocate memory for triangles, subsegments, vertices,   */
/*   viri (triangles being eaten), encroached segments, bad (skinny or too   */
/*   large) triangles, and splay tree nodes.                                 */

    struct memorypool triangles;
    struct memorypool subsegs;
    struct memorypool vertices;
    struct memorypool viri;
    struct memorypool badsubsegs;
    struct memorypool badtriangles;
    struct memorypool flipstackers;
    struct memorypool splaynodes;

/* Variables that maintain the bad triangle queues.  The queues are          */
/*   ordered from 4095 (highest priority) to 0 (lowest priority).            */

    struct badtriang *queuefront[4096];
    struct badtriang *queuetail[4096];
    int nextnonemptyq[4096];
    int firstnonemptyq;

/* Variable that maintains the stack of recently flipped triangles.          */

    struct flipstacker *lastflip;

/* Other variables. */

    REAL xmin, xmax, ymin, ymax;                            /* x and y bounds. */
    REAL xminextreme;      /* Nonexistent x value used as a flag in sweepline. */
    int invertices;                               /* Number of input vertices. */
    int inelements;                              /* Number of input triangles. */
    int insegments;                               /* Number of input segments. */
    int holes;                                       /* Number of input holes. */
    int regions;                                   /* Number of input regions. */
    int undeads;    /* Number of input vertices that don't appear in the mesh. */
    long edges;                                     /* Number of output edges. */
    int mesh_dim;                                /* Dimension (ought to be 2). */
    int nextras;                           /* Number of attributes per vertex. */
    int eextras;                         /* Number of attributes per triangle. */
    long hullsize;                          /* Number of edges in convex hull. */
    int steinerleft;                 /* Number of Steiner points not yet used. */
    int vertexmarkindex;         /* Index to find boundary marker of a vertex. */
    int vertex2triindex;     /* Index to find a triangle adjacent to a vertex. */
    int highorderindex;  /* Index to find extra nodes for high-order elements. */
    int elemattribindex;            /* Index to find attributes of a triangle. */
    int areaboundindex;             /* Index to find area bound of a triangle. */
    int checksegments;         /* Are there segments in the triangulation yet? */
    int checkquality;                  /* Has quality triangulation begun yet? */
    int readnodefile;                           /* Has a .node file been read? */
    long samples;              /* Number of random samples for point location. */

    long incirclecount;                 /* Number of incircle tests performed. */
    long counterclockcount;     /* Number of counterclockwise tests performed. */
    long orient3dcount;           /* Number of 3D orientation tests performed. */
    long hyperbolacount;      /* Number of right-of-hyperbola tests performed. */
    long circumcentercount;  /* Number of circumcenter calculations performed. */
    long circletopcount;       /* Number of circle top calculations performed. */

/* Triangular bounding box vertices.                                         */

    vertex infvertex1, infvertex2, infvertex3;

/* Pointer to the `triangle' that occupies all of "outer space."             */

    triangle *dummytri;
    triangle *dummytribase;    /* Keep base address so we can free() it later. */

/* Pointer to the omnipresent subsegment.  Referenced by any triangle or     */
/*   subsegment that isn't really connected to a subsegment at that          */
/*   location.                                                               */

    subseg *dummysub;
    subseg *dummysubbase;      /* Keep base address so we can free() it later. */

/* Pointer to a recently visited triangle.  Improves point location if       */
/*   proximate vertices are inserted sequentially.                           */

    struct otri recenttri;

};                                                  /* End of `struct mesh'. */


/* Data structure for command line switches and file names.  This structure  */
/*   is used (instead of global variables) to allow reentrancy.              */

struct behavior
{

/* Switches for the triangulator.                                            */
/*   poly: -p switch.  refine: -r switch.                                    */
/*   quality: -q switch.                                                     */
/*     minangle: minimum angle bound, specified after -q switch.             */
/*     goodangle: cosine squared of minangle.                                */
/*     offconstant: constant used to place off-center Steiner points.        */
/*   vararea: -a switch without number.                                      */
/*   fixedarea: -a switch with number.                                       */
/*     maxarea: maximum area bound, specified after -a switch.               */
/*   usertest: -u switch.                                                    */
/*   regionattrib: -A switch.  convex: -c switch.                            */
/*   weighted: 1 for -w switch, 2 for -W switch.  jettison: -j switch        */
/*   firstnumber: inverse of -z switch.  All items are numbered starting     */
/*     from `firstnumber'.                                                   */
/*   edgesout: -e switch.  voronoi: -v switch.                               */
/*   neighbors: -n switch.  geomview: -g switch.                             */
/*   nobound: -B switch.  nopolywritten: -P switch.                          */
/*   nonodewritten: -N switch.  noelewritten: -E switch.                     */
/*   noiterationnum: -I switch.  noholes: -O switch.                         */
/*   noexact: -X switch.                                                     */
/*   order: element order, specified after -o switch.                        */
/*   nobisect: count of how often -Y switch is selected.                     */
/*   steiner: maximum number of Steiner points, specified after -S switch.   */
/*   incremental: -i switch.  sweepline: -F switch.                          */
/*   dwyer: inverse of -l switch.                                            */
/*   splitseg: -s switch.                                                    */
/*   conformdel: -D switch.  docheck: -C switch.                             */
/*   quiet: -Q switch.  verbose: count of how often -V switch is selected.   */
/*   usesegments: -p, -r, -q, or -c switch; determines whether segments are  */
/*     used at all.                                                          */
/*                                                                           */
/* Read the instructions to find out the meaning of these switches.          */

    int poly, refine, quality, vararea, fixedarea, usertest;
    int regionattrib, convex, weighted, jettison;
    int firstnumber;
    int edgesout, voronoi, neighbors, geomview;
    int nobound, nopolywritten, nonodewritten, noelewritten, noiterationnum;
    int noholes, noexact, conformdel;
    int incremental, sweepline, dwyer;
    int splitseg;
    int docheck;
    int quiet, verbose;
    int usesegments;
    int order;
    int nobisect;
    int steiner;
    REAL minangle, goodangle, offconstant;
    REAL minarea; // used in our custom triunsuitable to ensure we aren't refining past the point of no return
    REAL maxarea;
    REAL maxtolerance;
    int tol_method; //tolerance method to use
    GDALDatasetH  hDataset;


/* Variables for file names.                                                 */

#ifndef TRILIBRARY
    char innodefilename[FILENAMESIZE];
    char inelefilename[FILENAMESIZE];
    char inpolyfilename[FILENAMESIZE];
    char areafilename[FILENAMESIZE];
    char outnodefilename[FILENAMESIZE];
    char outelefilename[FILENAMESIZE];
    char outpolyfilename[FILENAMESIZE];
    char edgefilename[FILENAMESIZE];
    char vnodefilename[FILENAMESIZE];
    char vedgefilename[FILENAMESIZE];
    char neighborfilename[FILENAMESIZE];
    char offfilename[FILENAMESIZE];
    char demfilename[FILENAMESIZE];
#endif /* not TRILIBRARY */

};                                              /* End of `struct behavior'. */


/*****************************************************************************/
/*                                                                           */
/*  Mesh manipulation primitives.  Each triangle contains three pointers to  */
/*  other triangles, with orientations.  Each pointer points not to the      */
/*  first byte of a triangle, but to one of the first three bytes of a       */
/*  triangle.  It is necessary to extract both the triangle itself and the   */
/*  orientation.  To save memory, I keep both pieces of information in one   */
/*  pointer.  To make this possible, I assume that all triangles are aligned */
/*  to four-byte boundaries.  The decode() routine below decodes a pointer,  */
/*  extracting an orientation (in the range 0 to 2) and a pointer to the     */
/*  beginning of a triangle.  The encode() routine compresses a pointer to a */
/*  triangle and an orientation into a single pointer.  My assumptions that  */
/*  triangles are four-byte-aligned and that the `unsigned long' type is     */
/*  long enough to hold a pointer are two of the few kludges in this program.*/
/*                                                                           */
/*  Subsegments are manipulated similarly.  A pointer to a subsegment        */
/*  carries both an address and an orientation in the range 0 to 1.          */
/*                                                                           */
/*  The other primitives take an oriented triangle or oriented subsegment,   */
/*  and return an oriented triangle or oriented subsegment or vertex; or     */
/*  they change the connections in the data structure.                       */
/*                                                                           */
/*  Below, triangles and subsegments are denoted by their vertices.  The     */
/*  triangle abc has origin (org) a, destination (dest) b, and apex (apex)   */
/*  c.  These vertices occur in counterclockwise order about the triangle.   */
/*  The handle abc may simultaneously denote vertex a, edge ab, and triangle */
/*  abc.                                                                     */
/*                                                                           */
/*  Similarly, the subsegment ab has origin (sorg) a and destination (sdest) */
/*  b.  If ab is thought to be directed upward (with b directly above a),    */
/*  then the handle ab is thought to grasp the right side of ab, and may     */
/*  simultaneously denote vertex a and edge ab.                              */
/*                                                                           */
/*  An asterisk (*) denotes a vertex whose identity is unknown.              */
/*                                                                           */
/*  Given this notation, a partial list of mesh manipulation primitives      */
/*  follows.                                                                 */
/*                                                                           */
/*                                                                           */
/*  For triangles:                                                           */
/*                                                                           */
/*  sym:  Find the abutting triangle; same edge.                             */
/*  sym(abc) -> ba*                                                          */
/*                                                                           */
/*  lnext:  Find the next edge (counterclockwise) of a triangle.             */
/*  lnext(abc) -> bca                                                        */
/*                                                                           */
/*  lprev:  Find the previous edge (clockwise) of a triangle.                */
/*  lprev(abc) -> cab                                                        */
/*                                                                           */
/*  onext:  Find the next edge counterclockwise with the same origin.        */
/*  onext(abc) -> ac*                                                        */
/*                                                                           */
/*  oprev:  Find the next edge clockwise with the same origin.               */
/*  oprev(abc) -> a*b                                                        */
/*                                                                           */
/*  dnext:  Find the next edge counterclockwise with the same destination.   */
/*  dnext(abc) -> *ba                                                        */
/*                                                                           */
/*  dprev:  Find the next edge clockwise with the same destination.          */
/*  dprev(abc) -> cb*                                                        */
/*                                                                           */
/*  rnext:  Find the next edge (counterclockwise) of the adjacent triangle.  */
/*  rnext(abc) -> *a*                                                        */
/*                                                                           */
/*  rprev:  Find the previous edge (clockwise) of the adjacent triangle.     */
/*  rprev(abc) -> b**                                                        */
/*                                                                           */
/*  org:  Origin          dest:  Destination          apex:  Apex            */
/*  org(abc) -> a         dest(abc) -> b              apex(abc) -> c         */
/*                                                                           */
/*  bond:  Bond two triangles together at the resepective handles.           */
/*  bond(abc, bad)                                                           */
/*                                                                           */
/*                                                                           */
/*  For subsegments:                                                         */
/*                                                                           */
/*  ssym:  Reverse the orientation of a subsegment.                          */
/*  ssym(ab) -> ba                                                           */
/*                                                                           */
/*  spivot:  Find adjoining subsegment with the same origin.                 */
/*  spivot(ab) -> a*                                                         */
/*                                                                           */
/*  snext:  Find next subsegment in sequence.                                */
/*  snext(ab) -> b*                                                          */
/*                                                                           */
/*  sorg:  Origin                      sdest:  Destination                   */
/*  sorg(ab) -> a                      sdest(ab) -> b                        */
/*                                                                           */
/*  sbond:  Bond two subsegments together at the respective origins.         */
/*  sbond(ab, ac)                                                            */
/*                                                                           */
/*                                                                           */
/*  For interacting tetrahedra and subfacets:                                */
/*                                                                           */
/*  tspivot:  Find a subsegment abutting a triangle.                         */
/*  tspivot(abc) -> ba                                                       */
/*                                                                           */
/*  stpivot:  Find a triangle abutting a subsegment.                         */
/*  stpivot(ab) -> ba*                                                       */
/*                                                                           */
/*  tsbond:  Bond a triangle to a subsegment.                                */
/*  tsbond(abc, ba)                                                          */
/*                                                                           */
/*****************************************************************************/

