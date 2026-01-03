#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#define DOUBLE_INFINITY (double)INFINITY

struct point3D
{
    float x;
    float y;
    float z;
};

struct MOBB3D
{
    struct point3D corners[8];
};

inline double dot3D(struct point3D P1, struct point3D P2)
{
    return P1.x*P2.x + P1.y*P2.y + P1.z*P2.z;
}

inline struct point3D cross3D(struct point3D P1, struct point3D P2)
{
    struct point3D result;
    result.x = P1.y*P2.z - P1.z*P2.y;
    result.y = P1.z*P2.x - P1.x*P2.z;
    result.z = P1.x*P2.y - P1.y*P2.x;
    return result;
}

// Brute Forces the 3D MOBB by considering all possible orientations (and a bunch extra)
// Runs in time O(n^7)
void MOBB3D_Brute_Force(struct point3D *new_vals, size_t num_points, struct MOBB3D *output)
{
    double mobbVolume = DOUBLE_INFINITY;
    struct point3D dirVec = {DOUBLE_INFINITY, DOUBLE_INFINITY, DOUBLE_INFINITY}; // this is to cache a non-zero vector to help if the points are all collinear
    // Loop through triples of all pairs of points (brute force doesn't consider anything special about the edges)
    for ( size_t i = 0; i < num_points; ++i )
    {
        for ( size_t j = i + 1; j < num_points; ++j )
        {
            struct point3D Pij = {new_vals[i].x-new_vals[j].x, new_vals[i].y-new_vals[j].y, new_vals[i].z-new_vals[j].z};
            // Early out if this vector is the zero vector as all cross products that use it will be zero
            if ( Pij.x != 0 || Pij.y != 0 || Pij.z != 0 )
            {
                dirVec = Pij;
                for ( size_t k = 0; k < num_points; ++k )
                {
                    for ( size_t l = k + 1; l < num_points; ++l )
                    {
                        struct point3D Pkl = {new_vals[k].x-new_vals[l].x, new_vals[k].y-new_vals[l].y, new_vals[k].z-new_vals[l].z};
                        struct point3D PijklCross = cross3D(Pij, Pkl);
                        // Skip this pair of (i,j) and (k,l) if it is the zero vector since the volume will be undefined
                        if ( PijklCross.x != 0 || PijklCross.y != 0 || PijklCross.z != 0 )
                        {
                            for ( size_t m = 0; m < num_points; ++m )
                            {
                                for ( size_t n = m + 1; n < num_points; ++n )
                                {
                                  struct point3D Pmn = {new_vals[m].x-new_vals[n].x, new_vals[m].y-new_vals[n].y, new_vals[m].z-new_vals[n].z};
                                  struct point3D PijklmnCross = cross3D(PijklCross, Pmn);
                                  // Skip this triple of (i,j), (k,l), and (m,n) if it is the zero vector since the volume will be undefined
                                  if ( PijklmnCross.x != 0 || PijklmnCross.y != 0 || PijklmnCross.z != 0 )
                                  {
                                        // We got a set of vectors we can finally work with. Compute the MOBB with this set
                                        struct point3D PijklmnDoubleCross = cross3D(PijklCross, PijklmnCross);
                                        double max_k_1 = dot3D(new_vals[0], PijklCross);
                                        double min_k_1 = max_k_1;
                                        double max_k_2 = dot3D(new_vals[0], PijklmnCross);
                                        double min_k_2 = max_k_2;
                                        double max_k_3 = dot3D(new_vals[0], PijklmnDoubleCross);
                                        double min_k_3 = max_k_3;
                                        // find the extents for the given triple
                                        for ( size_t p = 1; p < num_points; ++p )
                                        {
                                            double pijklcross_dot = dot3D(new_vals[p], PijklCross);
                                            double pijklmncross_dot = dot3D(new_vals[p], PijklmnCross);
                                            double pijklmndoublecross_dot = dot3D(new_vals[p], PijklmnDoubleCross);
                                            if (pijklcross_dot > max_k_1)
                                            {
                                                max_k_1 = pijklcross_dot;
                                            }
                                            if (pijklcross_dot < min_k_1)
                                            {
                                                min_k_1 = pijklcross_dot;
                                            }
                                            if (pijklmncross_dot > max_k_2)
                                            {
                                                max_k_2 = pijklmncross_dot;
                                            }
                                            if (pijklmncross_dot < min_k_2)
                                            {
                                                min_k_2 = pijklmncross_dot;
                                            }
                                            if (pijklmndoublecross_dot > max_k_3)
                                            {
                                                max_k_3 = pijklmndoublecross_dot;
                                            }
                                            if (pijklmndoublecross_dot < min_k_3)
                                            {
                                                min_k_3 = pijklmndoublecross_dot;
                                            }
                                        }
                                        double dir_normsquared1 = dot3D(PijklCross,PijklCross);
                                        double dir_normsquared2 = dot3D(PijklmnCross,PijklmnCross);
                                        double volume_ijklmn = (max_k_1-min_k_1)*(max_k_2-min_k_2)*(max_k_3-min_k_3)/dir_normsquared1/dir_normsquared2;
                                        if (volume_ijklmn < mobbVolume)
                                        {
                                            mobbVolume = volume_ijklmn;
                                            // NOTE: We divide by norm squared here since the vectors are not normed and the dot products calculated were not dotted with normed vectors either
                                            output->corners[0].x = (PijklCross.x*min_k_1/dir_normsquared1 + PijklmnCross.x*min_k_2/dir_normsquared2 + PijklmnDoubleCross.x*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[0].y = (PijklCross.y*min_k_1/dir_normsquared1 + PijklmnCross.y*min_k_2/dir_normsquared2 + PijklmnDoubleCross.y*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[0].z = (PijklCross.z*min_k_1/dir_normsquared1 + PijklmnCross.z*min_k_2/dir_normsquared2 + PijklmnDoubleCross.z*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[1].x = (PijklCross.x*min_k_1/dir_normsquared1 + PijklmnCross.x*min_k_2/dir_normsquared2 + PijklmnDoubleCross.x*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[1].y = (PijklCross.y*min_k_1/dir_normsquared1 + PijklmnCross.y*min_k_2/dir_normsquared2 + PijklmnDoubleCross.y*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[1].z = (PijklCross.z*min_k_1/dir_normsquared1 + PijklmnCross.z*min_k_2/dir_normsquared2 + PijklmnDoubleCross.z*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[2].x = (PijklCross.x*min_k_1/dir_normsquared1 + PijklmnCross.x*max_k_2/dir_normsquared2 + PijklmnDoubleCross.x*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[2].y = (PijklCross.y*min_k_1/dir_normsquared1 + PijklmnCross.y*max_k_2/dir_normsquared2 + PijklmnDoubleCross.y*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[2].z = (PijklCross.z*min_k_1/dir_normsquared1 + PijklmnCross.z*max_k_2/dir_normsquared2 + PijklmnDoubleCross.z*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[3].x = (PijklCross.x*min_k_1/dir_normsquared1 + PijklmnCross.x*max_k_2/dir_normsquared2 + PijklmnDoubleCross.x*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[3].y = (PijklCross.y*min_k_1/dir_normsquared1 + PijklmnCross.y*max_k_2/dir_normsquared2 + PijklmnDoubleCross.y*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[3].z = (PijklCross.z*min_k_1/dir_normsquared1 + PijklmnCross.z*max_k_2/dir_normsquared2 + PijklmnDoubleCross.z*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[4].x = (PijklCross.x*max_k_1/dir_normsquared1 + PijklmnCross.x*min_k_2/dir_normsquared2 + PijklmnDoubleCross.x*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[4].y = (PijklCross.y*max_k_1/dir_normsquared1 + PijklmnCross.y*min_k_2/dir_normsquared2 + PijklmnDoubleCross.y*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[4].z = (PijklCross.z*max_k_1/dir_normsquared1 + PijklmnCross.z*min_k_2/dir_normsquared2 + PijklmnDoubleCross.z*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[5].x = (PijklCross.x*max_k_1/dir_normsquared1 + PijklmnCross.x*min_k_2/dir_normsquared2 + PijklmnDoubleCross.x*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[5].y = (PijklCross.y*max_k_1/dir_normsquared1 + PijklmnCross.y*min_k_2/dir_normsquared2 + PijklmnDoubleCross.y*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[5].z = (PijklCross.z*max_k_1/dir_normsquared1 + PijklmnCross.z*min_k_2/dir_normsquared2 + PijklmnDoubleCross.z*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[6].x = (PijklCross.x*max_k_1/dir_normsquared1 + PijklmnCross.x*max_k_2/dir_normsquared2 + PijklmnDoubleCross.x*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[6].y = (PijklCross.y*max_k_1/dir_normsquared1 + PijklmnCross.y*max_k_2/dir_normsquared2 + PijklmnDoubleCross.y*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[6].z = (PijklCross.z*max_k_1/dir_normsquared1 + PijklmnCross.z*max_k_2/dir_normsquared2 + PijklmnDoubleCross.z*min_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[7].x = (PijklCross.x*max_k_1/dir_normsquared1 + PijklmnCross.x*max_k_2/dir_normsquared2 + PijklmnDoubleCross.x*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[7].y = (PijklCross.y*max_k_1/dir_normsquared1 + PijklmnCross.y*max_k_2/dir_normsquared2 + PijklmnDoubleCross.y*max_k_3/(dir_normsquared1*dir_normsquared2));
                                            output->corners[7].z = (PijklCross.z*max_k_1/dir_normsquared1 + PijklmnCross.z*max_k_2/dir_normsquared2 + PijklmnDoubleCross.z*max_k_3/(dir_normsquared1*dir_normsquared2));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // TODO: Do we need to consider the unique angle case as well?
        }
    }
    // This is to check to see that the MOBB even has volume (i.e., that the most inner loop actually ran once)
    if (!isfinite(mobbVolume))
    {
        if ( num_points > 0 )
        {
            // Collinear case: The cross product calculation will always come out as the zero in the main computation due to collinearity (which causes us to never run the innermost loop).
            // (dirVec would have been set if at least two points were unique)
            if ( !isfinite(dirVec.x) )
            {
                double max_k = dot3D(new_vals[0], dirVec);
                double min_k = max_k;
                struct point3D maxPoint = new_vals[0];
                struct point3D minPoint = new_vals[0];
                for ( size_t p = 1; p < num_points; ++p )
                {
                    double dir_dot = dot3D(new_vals[p], dirVec);
                    if (dir_dot > max_k)
                    {
                        max_k = dir_dot;
                        maxPoint = new_vals[p];
                    }
                    if (dir_dot < min_k)
                    {
                        min_k = dir_dot;
                        minPoint = new_vals[p];
                    }
                }
                // We now have the two most extreme points of the MOBB. Copy them into the output
                for ( size_t i = 0; i < 4; ++i )
                {
                    output->corners[i] = maxPoint;
                    output->corners[i+4] = minPoint;
                }
            }
            // All points are the same case
            else
            {
                // If we have duplicate points only, then the MOBB is just that singular point, so copy it into all eight corners
                for ( size_t i = 0; i < 8; ++i )
                {
                    output->corners[i] = new_vals[0];
                }
            }
        }
        else
        {
            printf("bruh, please at least calculate the MOBB with points. you gave me nothing\n");
        }
    }
}

// TODO: MOBB that uses Convex Hull Algorithms like O'Rourke's algorithm
