#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

// This computes the Max OBB for a 2D point set the brute force way
// Runs in time O(n^5)
void MaxOBB2D_Brute_Force(struct point2D *new_vals, size_t num_points,struct MOBB2D *output)
{
    double mobbArea = -1;
    // Loop through all pairs of points (brute force doesn't consider anything special about the edges)
    for ( int i = 0; i < num_points; ++i )
    {
        for ( int j = i + 1; j < num_points; ++j )
        {
            struct point2D Pij = {new_vals[i].x-new_vals[j].x, new_vals[i].y-new_vals[j].y};
            // Make sure this is not the zero vector or the critical points is no longer defined
            if ( Pij.x != 0 || Pij.y != 0 )
            {
                for ( int k = 0; k < num_points; ++k )
                {
                    for ( int l = k + 1; l < num_points; ++l )
                    {
                        struct point2D Pkl = {new_vals[k].x-new_vals[l].x, new_vals[k].y-new_vals[l].y};
                        // Make sure this is not the zero vector or the critical points is no longer defined
                        if ( Pkl.x != 0 || Pkl.y != 0 )
                        {
                            double Pij_norm = sqrt(Pij.x*Pij.x + Pij.y+Pij.y);
                            double Pkl_norm = sqrt(Pkl.x*Pkl.x + Pkl.y+Pkl.y);
                            double Mixed_term = fabs(Pij.x*Pkl.y + Pkl.x*Pij.y);
                            double sign_term = ((Pij.x*Pkl.x-Pij.y*Pkl.y)/(Pij.x*Pkl.y+Pkl.x*Pij.y)>0?1.0:-1.0);
                            struct point2D Qijkl = {sqrt(Pij_norm*Pkl_norm+Mixed_term), -sign_term*sqrt(Pij_norm*Pkl_norm-Mixed_term)};
                            struct point2D Qijklperp = {Qijkl.y, -Qijkl.x};
                            double max_k = dot2D(new_vals[0], Qijkl);
                            double min_k = max_k;
                            double max_k_perp = dot2D(new_vals[0], Qijklperp);
                            double min_k_perp = max_k_perp;
                            // find the extents for the given Qijkl
                            for ( int k = 1; k < num_points; ++k )
                            {
                                double qijkl_dot = dot2D(new_vals[k], Qijkl);
                                double qijkl_perp_dot = dot2D(new_vals[k], Qijklperp);
                                if (qijkl_dot > max_k)
                                {
                                    max_k = qijkl_dot;
                                }
                                if (qijkl_dot < min_k)
                                {
                                    min_k = qijkl_dot;
                                }
                                if (qijkl_perp_dot > max_k_perp)
                                {
                                    max_k_perp = qijkl_perp_dot;
                                }
                                if (qijkl_perp_dot < min_k_perp)
                                {
                                    min_k_perp = qijkl_perp_dot;
                                }
                            }
                            double dir_norm = Qijkl.x*Qijkl.x+Qijkl.y*Qijkl.y; // squared since dot product was not calculated with normed vectors
                            double area_ij = (max_k-min_k)*(max_k_perp-min_k_perp)/dir_norm;
                            if (area_ij > mobbArea)
                            {
                                mobbArea = area_ij;
                                output->corners[0].x = (Qijkl.x*min_k + Qijklperp.x*min_k_perp)/dir_norm;
                                output->corners[0].y = (Qijkl.y*min_k + Qijklperp.y*min_k_perp)/dir_norm;
                                output->corners[1].x = (Qijkl.x*min_k + Qijklperp.x*max_k_perp)/dir_norm;
                                output->corners[1].y = (Qijkl.y*min_k + Qijklperp.y*max_k_perp)/dir_norm;
                                output->corners[2].x = (Qijkl.x*max_k + Qijklperp.x*min_k_perp)/dir_norm;
                                output->corners[2].y = (Qijkl.y*max_k + Qijklperp.y*min_k_perp)/dir_norm;
                                output->corners[3].x = (Qijkl.x*max_k + Qijklperp.x*max_k_perp)/dir_norm;
                                output->corners[3].y = (Qijkl.y*max_k + Qijklperp.y*max_k_perp)/dir_norm;
                            }
                        }
                    }
                }
            }
        }
    }
    // This is to check to see that the MOBB even has area (i.e., that the most inner loop actually ran once)
    if ( mobbArea < 0 )
    {
        if ( num_points > 0 )
        {
            // If we have duplicate points only, then the MOBB is just that singular point, so copy it into all four corners
            for ( size_t i = 0; i < 4; ++i )
            {
                output->corners[i] = new_vals[0];
            }
        }
        else
        {
            printf("bruh, please at least calculate the Max OBB with points. you gave me nothing\n");
        }
    }
}


// This computes the Max OBB for a 2D point set the brute force way, first culling the interior points as the function describing the MOBB does not rely on them
// runs in time O(n*log(n)+h^5), where h is the number of convex hull vertices. In the worst case when all points are convex hull points, this degenerates to O(n^5)
void MaxOBB2D_Convex_Hull_Brute_Force(struct point2D *new_vals, size_t num_points,struct MOBB2D *output)
{
    double mobbArea = -1;
    struct point2D convex_hull[CONVEX_HULL_WORKSPACE];

    // Assert that the program has enough memory to do this task. Crash if not
    assert(CONVEX_HULL_WORKSPACE >= num_points);

    // Cull non-hull points
    size_t convex_hull_num_points = Convex_Hull_2D(new_vals, num_points, convex_hull);

    // Loop through all pairs of points (brute force doesn't consider anything special about the edges)
    for ( int i = 0; i < convex_hull_num_points; ++i )
    {
        for ( int j = i + 1; j < convex_hull_num_points; ++j )
        {
            struct point2D Pij = {convex_hull[i].x-convex_hull[j].x, convex_hull[i].y-convex_hull[j].y};
            // Make sure this is not the zero vector or the critical points is no longer defined
            if ( Pij.x != 0 || Pij.y != 0 )
            {
                for ( int k = 0; k < convex_hull_num_points; ++k )
                {
                    for ( int l = k + 1; l < convex_hull_num_points; ++l )
                    {
                        struct point2D Pkl = {convex_hull[k].x-convex_hull[l].x, convex_hull[k].y-convex_hull[l].y};
                        // Make sure this is not the zero vector or the critical points is no longer defined
                        if ( Pkl.x != 0 || Pkl.y != 0 )
                        {
                            double Pij_norm = sqrt(Pij.x*Pij.x + Pij.y+Pij.y);
                            double Pkl_norm = sqrt(Pkl.x*Pkl.x + Pkl.y+Pkl.y);
                            double Mixed_term = fabs(Pij.x*Pkl.y + Pkl.x*Pij.y);
                            double sign_term = ((Pij.x*Pkl.x-Pij.y*Pkl.y)/(Pij.x*Pkl.y+Pkl.x*Pij.y)>0?1.0:-1.0);
                            struct point2D Qijkl = {sqrt(Pij_norm*Pkl_norm+Mixed_term), -sign_term*sqrt(Pij_norm*Pkl_norm-Mixed_term)};
                            struct point2D Qijklperp = {Qijkl.y, -Qijkl.x};
                            double max_k = dot2D(convex_hull[0], Qijkl);
                            double min_k = max_k;
                            double max_k_perp = dot2D(convex_hull[0], Qijklperp);
                            double min_k_perp = max_k_perp;
                            // find the extents for the given Qijkl
                            for ( int k = 1; k < convex_hull_num_points; ++k )
                            {
                                double qijkl_dot = dot2D(convex_hull[k], Qijkl);
                                double qijkl_perp_dot = dot2D(convex_hull[k], Qijklperp);
                                if (qijkl_dot > max_k)
                                {
                                    max_k = qijkl_dot;
                                }
                                if (qijkl_dot < min_k)
                                {
                                    min_k = qijkl_dot;
                                }
                                if (qijkl_perp_dot > max_k_perp)
                                {
                                    max_k_perp = qijkl_perp_dot;
                                }
                                if (qijkl_perp_dot < min_k_perp)
                                {
                                    min_k_perp = qijkl_perp_dot;
                                }
                            }
                            double dir_norm = Qijkl.x*Qijkl.x+Qijkl.y*Qijkl.y; // squared since dot product was not calculated with normed vectors
                            double area_ij = (max_k-min_k)*(max_k_perp-min_k_perp)/dir_norm;
                            if (area_ij > mobbArea)
                            {
                                mobbArea = area_ij;
                                output->corners[0].x = (Qijkl.x*min_k + Qijklperp.x*min_k_perp)/dir_norm;
                                output->corners[0].y = (Qijkl.y*min_k + Qijklperp.y*min_k_perp)/dir_norm;
                                output->corners[1].x = (Qijkl.x*min_k + Qijklperp.x*max_k_perp)/dir_norm;
                                output->corners[1].y = (Qijkl.y*min_k + Qijklperp.y*max_k_perp)/dir_norm;
                                output->corners[2].x = (Qijkl.x*max_k + Qijklperp.x*min_k_perp)/dir_norm;
                                output->corners[2].y = (Qijkl.y*max_k + Qijklperp.y*min_k_perp)/dir_norm;
                                output->corners[3].x = (Qijkl.x*max_k + Qijklperp.x*max_k_perp)/dir_norm;
                                output->corners[3].y = (Qijkl.y*max_k + Qijklperp.y*max_k_perp)/dir_norm;
                            }
                        }
                    }
                }
            }
        }
    }
    // This is to check to see that the MOBB even has area (i.e., that the most inner loop actually ran once)
    if ( mobbArea < 0 )
    {
        if ( num_points > 0 )
        {
            // If we have duplicate points only, then the MOBB is just that singular point, so copy it into all four corners
            for ( size_t i = 0; i < 4; ++i )
            {
                output->corners[i] = new_vals[0];
            }
        }
        else
        {
            printf("bruh, please at least calculate the Max OBB with points. you gave me nothing\n");
        }
    }
}
