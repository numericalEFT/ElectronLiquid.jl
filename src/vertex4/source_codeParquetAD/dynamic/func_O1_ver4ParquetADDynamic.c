
void eval_graph100(double *root, double *leafVal)
{
    double  g1, g4, g9, g10, g11, g13, g14, g16, g19, g22;
    g1 = leafVal[0] * -1.0;
    g4 = leafVal[1];
    g9 = (g1 + g4);
    g10 = (g9);
    root[0] = g10;
    g11 = leafVal[2] * -1.0;
    g13 = (g11);
    root[1] = g13;
    g14 = leafVal[3];
    g16 = (g14);
    root[2] = g16;
    g19 = (g1);
    root[3] = g19;
    g22 = (g11);
    root[4] = g22;
}