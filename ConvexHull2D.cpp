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
