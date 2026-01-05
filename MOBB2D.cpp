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


// HACK: Remove this in future as we do not want to use stl calls so this can be a pure C implementation
bool comparePoints(const struct point2D p1, const struct point2D p2) {
    if (p1.x == p2.x) {
        return p1.y < p2.y;
    }
    return p1.x < p2.x;
}

// This is an implementation of graham scan for calculating the convex hull
// Graham Scan runs in time O(n log(n)) due to sorting
size_t Convex_Hull_2D(struct point2D *new_vals, size_t num_points, struct point2D *convex_hull)
{
    // Assert that the program has enough memory to do this task. Crash if not
    assert( CONVEX_HULL_WORKSPACE >= num_points );

    size_t convex_hull_num_points = 0;
    if (num_points == 0)
    {
        return 0;
    }

    // sort the points
    struct point2D sorted_points[CONVEX_HULL_WORKSPACE];

    // HACK: Improve this in the future by not memcopying and std sorting
    memcpy(sorted_points, new_vals, num_points*sizeof(struct point2D));
    std::sort(sorted_points, sorted_points+num_points, comparePoints);

    // intialize the hull with the first point
    convex_hull[0] = sorted_points[0];
    ++convex_hull_num_points;

    // find the next point that is different than the first point
    int i = 1;
    for ( ; i < num_points; ++i )
    {
        if (convex_hull[0].x != sorted_points[i].x || convex_hull[0].y != sorted_points[i].y)
        {
            convex_hull[1] = sorted_points[i];
            ++convex_hull_num_points;
            ++i;
            break;
        }
    }

    // Check to see if there was at least one second different point (break early if none)
    if ( convex_hull_num_points == 1 )
    {
        return 1;
    }

    // find the upper hull
    for ( ; i < num_points; ++i )
    {
        // right turn test
        struct point2D pq = {convex_hull[convex_hull_num_points-1].x-convex_hull[convex_hull_num_points-2].x,convex_hull[convex_hull_num_points-1].y-convex_hull[convex_hull_num_points-2].y};
        struct point2D qr = {sorted_points[i].x-convex_hull[convex_hull_num_points-1].x,sorted_points[i].y-convex_hull[convex_hull_num_points-1].y};
        // if it makes a right turn, then add it to the convex hull
        if ( pq.x*qr.y < pq.y*qr.x )
        {
            // this pushes onto the "stack"
            convex_hull[convex_hull_num_points] = sorted_points[i];
            ++convex_hull_num_points;
        }
        else
        {
            // if there are now only two points upon removing one point, then just replace the previous point and move on
            if ( convex_hull_num_points < 3 )
            {
                convex_hull[convex_hull_num_points-1] = sorted_points[i];
            }
            else
            {
                // pop the stack and check again
                --convex_hull_num_points;
                --i;
            }
        }
    }


    // find the lower hull
    struct point2D lower_hull[CONVEX_HULL_WORKSPACE];
    size_t lower_hull_points = 0;

    lower_hull[0] = sorted_points[num_points-1];
    ++lower_hull_points;

    // find the next point that is different than the first point
    int j = num_points;
    for ( ; j >= 0; --j )
    {
        if (lower_hull[0].x != sorted_points[j].x || lower_hull[0].y != sorted_points[j].y)
        {
            lower_hull[1] = sorted_points[j];
            ++lower_hull_points;
            --j;
            break;
        }
    }

    for ( ; j >= 0; --j )
    {
        // right turn test
        struct point2D pq = {lower_hull[lower_hull_points-1].x-lower_hull[lower_hull_points-2].x,lower_hull[lower_hull_points-1].y-lower_hull[lower_hull_points-2].y};
        struct point2D qr = {sorted_points[j].x-lower_hull[lower_hull_points-1].x,sorted_points[j].y-lower_hull[lower_hull_points-1].y};
        // if it makes a right turn, then add it to the convex hull
        if ( pq.x*qr.y < pq.y*qr.x )
        {
            // this pushes onto the "stack"
            lower_hull[lower_hull_points] = sorted_points[j];
            ++lower_hull_points;
        }
        else
        {
            // if there are now only two points upon removing one point, then just replace the previous point and move on
            if ( lower_hull_points < 3 )
            {
                lower_hull[lower_hull_points-1] = sorted_points[j];
            }
            else
            {
                // pop the stack and check again
                --lower_hull_points;
                ++j;
            }
        }
    }
    // merge append the lower hull to the upper hull
    memcpy(&convex_hull[convex_hull_num_points], &lower_hull[1], (lower_hull_points-2)*sizeof(point2D));
    convex_hull_num_points += lower_hull_points-2; // get rid of the duplicated last point
    return convex_hull_num_points;
}

// This computes the MOBB for a 2D point set the brute force way (it culls the interior convex hull points only, but still considers interior edges between the vertices)
// Worst case, all points are convex hull vertices and it runs in the same time complexity as the pure brute force approach
// Runs in time O(n log(n)+h^3) where h is the number of convex hull vertices of the point set
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
// It runs an O(n) calculation on top of the <= n Convex hull edges. As the convex hull calculation takes at most O(n log(n)) steps,
// This runs in time O(n^2)
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

// TODO: Add an implementation of Toussaint's rotating calliper method here
