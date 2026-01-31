#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <algorithm> // for std::sort

// Workspace size for computing convex hull
#define CONVEX_HULL_WORKSPACE 4096

struct point2D
{
    float x;
    float y;
};

struct MOBB2D
{
    struct point2D corners[4];
};

inline double dot2D(struct point2D P1, struct point2D P2)
{
    return P1.x*P2.x + P1.y*P2.y;
}


// This computes the MOBB for a 2D point set the brute force way
// Runs in time O(n^3)
void MOBB2D_Brute_Force(struct point2D *new_vals, size_t num_points,struct MOBB2D *output)
{
    double mobbArea = DOUBLE_INFINITY;
    // Loop through all pairs of points (brute force doesn't consider anything special about the edges)
    for ( int i = 0; i < num_points; ++i )
    {
        for ( int j = i + 1; j < num_points; ++j )
        {
            struct point2D Pij = {new_vals[i].x-new_vals[j].x, new_vals[i].y-new_vals[j].y};
            // Make sure this is not the zero vector or the critical points is no longer defined
            if ( Pij.x != 0 || Pij.y != 0 )
            {
                struct point2D Pij_perp = {-Pij.y, Pij.x};
                double max_k = dot2D(new_vals[0], Pij);
                double min_k = max_k;
                double max_k_perp = dot2D(new_vals[0], Pij_perp);
                double min_k_perp = max_k_perp;
                // find the extents for the given Pij
                for ( int k = 1; k < num_points; ++k )
                {
                    double pij_dot = dot2D(new_vals[k], Pij);
                    double pij_perp_dot = dot2D(new_vals[k], Pij_perp);
                    if (pij_dot > max_k)
                    {
                        max_k = pij_dot;
                    }
                    if (pij_dot < min_k)
                    {
                        min_k = pij_dot;
                    }
                    if (pij_perp_dot > max_k_perp)
                    {
                        max_k_perp = pij_perp_dot;
                    }
                    if (pij_perp_dot < min_k_perp)
                    {
                        min_k_perp = pij_perp_dot;
                    }
                }
                double area_ij = (max_k-min_k)*(max_k_perp-min_k_perp)/(Pij.x*Pij.x+Pij.y*Pij.y);
                if (area_ij < mobbArea)
                {
                    mobbArea = area_ij;
                    double dir_norm = Pij.x*Pij.x+Pij.y*Pij.y; // squared since dot product was not calculated with normed vectors
                    output->corners[0].x = (Pij.x*min_k + Pij_perp.x*min_k_perp)/dir_norm;
                    output->corners[0].y = (Pij.y*min_k + Pij_perp.y*min_k_perp)/dir_norm;
                    output->corners[1].x = (Pij.x*min_k + Pij_perp.x*max_k_perp)/dir_norm;
                    output->corners[1].y = (Pij.y*min_k + Pij_perp.y*max_k_perp)/dir_norm;
                    output->corners[2].x = (Pij.x*max_k + Pij_perp.x*min_k_perp)/dir_norm;
                    output->corners[2].y = (Pij.y*max_k + Pij_perp.y*min_k_perp)/dir_norm;
                    output->corners[3].x = (Pij.x*max_k + Pij_perp.x*max_k_perp)/dir_norm;
                    output->corners[3].y = (Pij.y*max_k + Pij_perp.y*max_k_perp)/dir_norm;
                }
            }
        }
    }
    // This is to check to see that the MOBB even has area (i.e., that the most inner loop actually ran once)
    if (!isfinite(mobbArea))
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
            printf("bruh, please at least calculate the MOBB with points. you gave me nothing\n");
        }
    }
}

#include "ConvexHull2D.cpp"

// This computes the MOBB for a 2D point set the brute force way (it culls the interior convex hull points only, but still considers interior edges between the vertices)
// Runs in time O(n log(n)+h^3) where h is the number of convex hull vertices of the point set
// Worst case, all points are convex hull vertices and it runs in the same time complexity as the pure brute force approach O(n^3)
void MOBB2D_Convex_Hull_Culled_Brute_Force(struct point2D *new_vals, size_t num_points,struct MOBB2D *output)
{
    double mobbArea = DOUBLE_INFINITY;
    struct point2D convex_hull[CONVEX_HULL_WORKSPACE];

    // Assert that the program has enough memory to do this task. Crash if not
    assert(CONVEX_HULL_WORKSPACE >= num_points);

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
                struct point2D Pij_perp = {-Pij.y, Pij.x};
                double max_k = dot2D(convex_hull[0], Pij);
                double min_k = max_k;
                double max_k_perp = dot2D(convex_hull[0], Pij_perp);
                double min_k_perp = max_k_perp;
                // find the extents for the given Pij
                for ( int k = 1; k < convex_hull_num_points; ++k )
                {
                    double pij_dot = dot2D(convex_hull[k], Pij);
                    double pij_perp_dot = dot2D(convex_hull[k], Pij_perp);
                    if (pij_dot > max_k)
                    {
                        max_k = pij_dot;
                    }
                    if (pij_dot < min_k)
                    {
                        min_k = pij_dot;
                    }
                    if (pij_perp_dot > max_k_perp)
                    {
                        max_k_perp = pij_perp_dot;
                    }
                    if (pij_perp_dot < min_k_perp)
                    {
                        min_k_perp = pij_perp_dot;
                    }
                }
                double dir_norm = Pij.x*Pij.x+Pij.y*Pij.y; // squared since dot product was not calculated with normed vectors
                double area_ij = (max_k-min_k)*(max_k_perp-min_k_perp)/dir_norm;
                if (area_ij < mobbArea)
                {
                    mobbArea = area_ij;
                    output->corners[0].x = (Pij.x*min_k + Pij_perp.x*min_k_perp)/dir_norm;
                    output->corners[0].y = (Pij.y*min_k + Pij_perp.y*min_k_perp)/dir_norm;
                    output->corners[1].x = (Pij.x*min_k + Pij_perp.x*max_k_perp)/dir_norm;
                    output->corners[1].y = (Pij.y*min_k + Pij_perp.y*max_k_perp)/dir_norm;
                    output->corners[2].x = (Pij.x*max_k + Pij_perp.x*min_k_perp)/dir_norm;
                    output->corners[2].y = (Pij.y*max_k + Pij_perp.y*min_k_perp)/dir_norm;
                    output->corners[3].x = (Pij.x*max_k + Pij_perp.x*max_k_perp)/dir_norm;
                    output->corners[3].y = (Pij.y*max_k + Pij_perp.y*max_k_perp)/dir_norm;
                }
            }
        }
    }
    // This is to check to see that the MOBB even has area (i.e., that the most inner loop actually ran once)
    if (!isfinite(mobbArea))
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
            printf("bruh, please at least calculate the MOBB with points. you gave me nothing\n");
        }
    }
}

// This computes the MOBB for a 2D point set the brute force way
// It runs an O(h) calculation on top of the h Convex hull edges. As the convex hull calculation takes at most O(n log(n)) steps,
// This runs in time O(n log(n) + h^2).
// Worst case, all points are convex hull vertices, giving the total algorithm runtime to be O(n^2)
// This was discovered independently by the author through algrebraic means and later found to have already been written
// and implemented by the CGAL library. This is included here for comparison.
void MOBB2D_Convex_Hull_Brute_Force(struct point2D *new_vals, size_t num_points,struct MOBB2D *output)
{
    double mobbArea = DOUBLE_INFINITY;
    struct point2D convex_hull[CONVEX_HULL_WORKSPACE];

    // Assert that the program has enough memory to do this task. Crash if not
    assert(CONVEX_HULL_WORKSPACE >= num_points);

    size_t convex_hull_num_points = Convex_Hull_2D(new_vals, num_points, convex_hull);

    // Loop through all convex hull edges
    for ( int i = 0; i < convex_hull_num_points; ++i )
    {
        int j = (i + 1) % convex_hull_num_points; // the adjacent edge
        struct point2D Pij = {convex_hull[i].x-convex_hull[j].x, convex_hull[i].y-convex_hull[j].y};
        // Make sure this is not the zero vector or the critical points is no longer defined
        if ( Pij.x != 0 || Pij.y != 0 )
        {
            struct point2D Pij_perp = {-Pij.y, Pij.x};
            double max_k = dot2D(convex_hull[0], Pij);
            double min_k = max_k;
            double max_k_perp = dot2D(convex_hull[0], Pij_perp);
            double min_k_perp = max_k_perp;
            // find the extents for the given Pij
            for ( int k = 1; k < convex_hull_num_points; ++k )
            {
                double pij_dot = dot2D(convex_hull[k], Pij);
                double pij_perp_dot = dot2D(convex_hull[k], Pij_perp);
                if (pij_dot > max_k)
                {
                    max_k = pij_dot;
                }
                if (pij_dot < min_k)
                {
                    min_k = pij_dot;
                }
                if (pij_perp_dot > max_k_perp)
                {
                    max_k_perp = pij_perp_dot;
                }
                if (pij_perp_dot < min_k_perp)
                {
                    min_k_perp = pij_perp_dot;
                }
            }
            double area_ij = (max_k-min_k)*(max_k_perp-min_k_perp)/(Pij.x*Pij.x+Pij.y*Pij.y);
            if (area_ij < mobbArea)
            {
                mobbArea = area_ij;
                double dir_norm = Pij.x*Pij.x+Pij.y*Pij.y; // squared since dot product was not calculated with normed vectors
                output->corners[0].x = (Pij.x*min_k + Pij_perp.x*min_k_perp)/dir_norm;
                output->corners[0].y = (Pij.y*min_k + Pij_perp.y*min_k_perp)/dir_norm;
                output->corners[1].x = (Pij.x*min_k + Pij_perp.x*max_k_perp)/dir_norm;
                output->corners[1].y = (Pij.y*min_k + Pij_perp.y*max_k_perp)/dir_norm;
                output->corners[2].x = (Pij.x*max_k + Pij_perp.x*min_k_perp)/dir_norm;
                output->corners[2].y = (Pij.y*max_k + Pij_perp.y*min_k_perp)/dir_norm;
                output->corners[3].x = (Pij.x*max_k + Pij_perp.x*max_k_perp)/dir_norm;
                output->corners[3].y = (Pij.y*max_k + Pij_perp.y*max_k_perp)/dir_norm;
            }
        }
        else
        {
            printf("bruh, your convex hull algo is broken. it gives me duplicate points in the convex hull\n");
        }
    }
    // This is to check to see that the MOBB even has area (i.e., that the most inner loop actually ran once)
    if (!isfinite(mobbArea))
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
            printf("bruh, please at least calculate the MOBB with some points. you gave me nothing\n");
        }
    }
}

// This computes the MOBB for a 2D point set the brute force way (same as the above, but with a small optimization that allows the removal of an if statement in the inner most loop to a single computation outside the loop)
// It runs an O(h) calculation on top of the h Convex hull edges. As the convex hull calculation takes at most O(n log(n)) steps,
// This runs in time O(n log(n) + h^2).
// Worst case, all points are convex hull vertices, giving the total algorithm runtime to be O(n^2)
// This was discovered independently by the author through algrebraic means and later found to have already been written
// and implemented by the CGAL library. This is included here for comparison.
void MOBB2D_Convex_Hull_Directed_Brute_Force(struct point2D *new_vals, size_t num_points,struct MOBB2D *output)
{
    double mobbArea = DOUBLE_INFINITY;
    struct point2D convex_hull[CONVEX_HULL_WORKSPACE];

    // Assert that the program has enough memory to do this task. Crash if not
    assert(CONVEX_HULL_WORKSPACE >= num_points);

    size_t convex_hull_num_points = Convex_Hull_2D(new_vals, num_points, convex_hull);

    // Loop through all convex hull edges
    for ( int i = 0; i < convex_hull_num_points; ++i )
    {
        int j = (i + 1) % convex_hull_num_points; // the adjacent edge
        struct point2D Pij = {convex_hull[i].x-convex_hull[j].x, convex_hull[i].y-convex_hull[j].y};
        // Make sure this is not the zero vector or the critical points is no longer defined
        if ( Pij.x != 0 || Pij.y != 0 )
        {
            struct point2D Pij_perp = {-Pij.y, Pij.x};
            double max_k = dot2D(convex_hull[0], Pij);
            double min_k = max_k;
            double max_k_perp = dot2D(convex_hull[0], Pij_perp);
            //Optimization: We know Pij_perp points inwards due to ordering of convex hull points and how we compute it. Since Pi and Pj form a convex hull edge, they must be the min value in the Pij_perp direction
            //              Also, choice of using Pi vs Pj doesn't matter since dot2D(Pi, Pij_perp) == dot2D(Pj, Pij_perp)
            double min_k_perp = dot2D(convex_hull[i], Pij_perp);
            // find the extents for the given Pij
            for ( int k = 1; k < convex_hull_num_points; ++k )
            {
                double pij_dot = dot2D(convex_hull[k], Pij);
                double pij_perp_dot = dot2D(convex_hull[k], Pij_perp);
                if (pij_dot > max_k)
                {
                    max_k = pij_dot;
                }
                if (pij_dot < min_k)
                {
                    min_k = pij_dot;
                }
                // Only need to compute max since we know the min.
                if (pij_perp_dot > max_k_perp)
                {
                    max_k_perp = pij_perp_dot;
                }
            }
            double area_ij = (max_k-min_k)*(max_k_perp-min_k_perp)/(Pij.x*Pij.x+Pij.y*Pij.y);
            if (area_ij < mobbArea)
            {
                mobbArea = area_ij;
                double dir_norm = Pij.x*Pij.x+Pij.y*Pij.y; // squared since dot product was not calculated with normed vectors
                output->corners[0].x = (Pij.x*min_k + Pij_perp.x*min_k_perp)/dir_norm;
                output->corners[0].y = (Pij.y*min_k + Pij_perp.y*min_k_perp)/dir_norm;
                output->corners[1].x = (Pij.x*min_k + Pij_perp.x*max_k_perp)/dir_norm;
                output->corners[1].y = (Pij.y*min_k + Pij_perp.y*max_k_perp)/dir_norm;
                output->corners[2].x = (Pij.x*max_k + Pij_perp.x*min_k_perp)/dir_norm;
                output->corners[2].y = (Pij.y*max_k + Pij_perp.y*min_k_perp)/dir_norm;
                output->corners[3].x = (Pij.x*max_k + Pij_perp.x*max_k_perp)/dir_norm;
                output->corners[3].y = (Pij.y*max_k + Pij_perp.y*max_k_perp)/dir_norm;
            }
        }
        else
        {
            printf("bruh, your convex hull algo is broken. it gives me duplicate points in the convex hull\n");
        }
    }
    // This is to check to see that the MOBB even has area (i.e., that the most inner loop actually ran once)
    if (!isfinite(mobbArea))
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
            printf("bruh, please at least calculate the MOBB with some points. you gave me nothing\n");
        }
    }
}

// This computes the MOBB for a 2D point set the brute force way (same as the above, but with the optimization more stable)
// It runs an O(h) calculation on top of the h Convex hull edges. As the convex hull calculation takes at most O(n log(n)) steps,
// This runs in time O(n log(n) + h^2).
// Worst case, all points are convex hull vertices, giving the total algorithm runtime to be O(n^2)
// This was discovered independently by the author through algrebraic means and later found to have already been written
// and implemented by the CGAL library. This is included here for comparison.
void MOBB2D_Convex_Hull_Directed_Stable_Brute_Force(struct point2D *new_vals, size_t num_points,struct MOBB2D *output)
{
    double mobbArea = DOUBLE_INFINITY;
    struct point2D convex_hull[CONVEX_HULL_WORKSPACE];

    // Assert that the program has enough memory to do this task. Crash if not
    assert(CONVEX_HULL_WORKSPACE >= num_points);

    size_t convex_hull_num_points = Convex_Hull_2D(new_vals, num_points, convex_hull);

    // Loop through all convex hull edges
    for ( int i = 0; i < convex_hull_num_points; ++i )
    {
        int j = (i + 1) % convex_hull_num_points; // the adjacent edge
        struct point2D Pij = {convex_hull[i].x-convex_hull[j].x, convex_hull[i].y-convex_hull[j].y};
        // Make sure this is not the zero vector or the critical points is no longer defined
        if ( Pij.x != 0 || Pij.y != 0 )
        {
            struct point2D Pij_perp = {-Pij.y, Pij.x};
            double max_k = dot2D(convex_hull[0], Pij);
            double min_k = max_k;
            double max_k_perp = dot2D(convex_hull[0], Pij_perp);
            //Optimization: We know Pij_perp points inwards due to ordering of convex hull points and how we compute it. Since Pi and Pj form a convex hull edge, they must be the min value in the Pij_perp direction
            double min_k_perp = convex_hull[i].x*convex_hull[j].y-convex_hull[j].x*convex_hull[i].y;
            // find the extents for the given Pij
            for ( int k = 1; k < convex_hull_num_points; ++k )
            {
                double pij_dot = dot2D(convex_hull[k], Pij);
                double pij_perp_dot = dot2D(convex_hull[k], Pij_perp);
                if (pij_dot > max_k)
                {
                    max_k = pij_dot;
                }
                if (pij_dot < min_k)
                {
                    min_k = pij_dot;
                }
                // Only need to compute max since we know the min.
                if (pij_perp_dot > max_k_perp)
                {
                    max_k_perp = pij_perp_dot;
                }
            }
            double area_ij = (max_k-min_k)*(max_k_perp-min_k_perp)/(Pij.x*Pij.x+Pij.y*Pij.y);
            if (area_ij < mobbArea)
            {
                mobbArea = area_ij;
                double dir_norm = Pij.x*Pij.x+Pij.y*Pij.y; // squared since dot product was not calculated with normed vectors
                output->corners[0].x = (Pij.x*min_k + Pij_perp.x*min_k_perp)/dir_norm;
                output->corners[0].y = (Pij.y*min_k + Pij_perp.y*min_k_perp)/dir_norm;
                output->corners[1].x = (Pij.x*min_k + Pij_perp.x*max_k_perp)/dir_norm;
                output->corners[1].y = (Pij.y*min_k + Pij_perp.y*max_k_perp)/dir_norm;
                output->corners[2].x = (Pij.x*max_k + Pij_perp.x*min_k_perp)/dir_norm;
                output->corners[2].y = (Pij.y*max_k + Pij_perp.y*min_k_perp)/dir_norm;
                output->corners[3].x = (Pij.x*max_k + Pij_perp.x*max_k_perp)/dir_norm;
                output->corners[3].y = (Pij.y*max_k + Pij_perp.y*max_k_perp)/dir_norm;
            }
        }
        else
        {
            printf("bruh, your convex hull algo is broken. it gives me duplicate points in the convex hull\n");
        }
    }
    // This is to check to see that the MOBB even has area (i.e., that the most inner loop actually ran once)
    if (!isfinite(mobbArea))
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
            printf("bruh, please at least calculate the MOBB with some points. you gave me nothing\n");
        }
    }
}

// TODO: Add an implementation of Toussaint's rotating calliper method here
