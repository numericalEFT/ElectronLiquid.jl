
void eval_graph100(double *root, double *leafVal)
{
    double  g1, g4, g9, g10, g13;
    g1 = leafVal[0] * -1.0;
    g4 = leafVal[1];
    g9 = (g1 + g4);
    g10 = (g9);
    root[0] = g10;
    g13 = (g1);
    root[1] = g13;
}