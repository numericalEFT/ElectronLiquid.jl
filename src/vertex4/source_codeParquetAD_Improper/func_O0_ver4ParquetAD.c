#include <math.h>

void eval_graph000(double *root, double *leafVal)
{
    double  g18630, g18632, g18641, g18642;
    g18630 = leafVal[0];
    g18632 = leafVal[1];
    g18641 = (g18630 * -1.0 + g18632);
    root[0] = g18641;
    g18642 = (g18630 * -1.0);
    root[1] = g18642;
}