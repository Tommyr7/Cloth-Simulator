#include "tensormax.hpp"
#include "util.hpp"
using namespace std;

struct Disk {
    Vec2 c;
    double r;
    Disk (): c(Vec2(0)), r(0) {}
    Disk (const Vec2 &c, double r): c(c), r(r) {}
};

Disk minidisk (const vector<Disk> &P, const vector<Disk> &R) {
    if (P.empty() || R.size() == 3)
        return md(R);

    Disk p = P.front();//choose Disk p in P

    vector<Disk> delta_P(P.size()-1);
    for (int i = 0; i < delta_P.size(); i++)
        depta_P[i] = P[i+1];

    Disk D = minidisk(delta_P,R);

    if (norm(p.c-D.c) + p.r <= D.r + 1e-6)
        return D;
    else {
        vector<Disk> p_with_R(R.size()+1);
        p_with_R[0] = p;
        for (int i = 1; i < p_with_R.size(); i++)
            p_with_R[i] = R[i-1];
        return minidisk(delta_P, p_with_R);
    }
}

Disk md (const vector<Disk> &R) {
    if (R.empty())
        return Disk();
 
    if (R.size() == 1)
        return R.front();
 
    if (R.size() == 2) {
        double d = norm(R[0].c - R[1].c);
        double r = (R[0].r + d + R[1].r)/2;
        double t = (r - R[0].r)/d;
        
        if(debug){
            cout << R[0].c << "," << t << "," << R[1].c << "," << r << endl;
        }
        
        return Disk(R[0].c + t*(R[1].c - R[0].c), r);
    }

//http://rosettacode.org/mw/index.php?title=Problem_of_Apollonius&oldid=88212
    double x1 = R[0].c[0], y1 = R[0].c[1], r1 = R[0].r;
    double x2 = R[1].c[0], y2 = R[1].c[1], r2 = R[1].r;
    double x3 = R[2].c[0], y3 = R[2].c[1], r3 = R[2].r;

    int s1 = 1, s2 = 1, s3 = 1;
    double v11 = 2*x2 - 2*x1, v12 = 2*y2 - 2*y1;
    double v13 = x1*x1 - x2*x2 + y1*y1 - y2*y2 - r1*r1 + r2*r2;
    double v14 = 2*s2*r2 - 2*s1*r1;

    double v21 = 2*x3 - 2*x2, v22 = 2*y3 - 2*y2;
    double v23 = x2*x2 - x3*x3 + y2*y2 - y3*y3 - r2*r2 + r3*r3;
    double v24 = 2*s3*r3 - 2*s2*r2;

    double w12 = v12/v11, w13 = v13/v11, w14 = v14/v11;
    double w22 = v22/v21-w12, w23 = v23/v21-w13, w24 = v24/v21-w14;

    double P = -w23/w22, Q = w24/w22, M = -w12*P-w13, N = w14 - w12*Q;
    double a = N*N + Q*Q - 1, b = 2*M*N - 2*N*x1 + 2*P*Q - 2*Q*y1 + 2*s1*r1;
    double c = x1*x1 + M*M - 2*M*x1 + P*P + y1*y1 - 2*P*y1 - r1*r1;
    double D = b*b-4*a*c;
    double rs = (-b-sqrt(D))/(2*a), xs = M+N*rs, ys = P+Q*rs;

    if(debug){
        cout << xs << "," << ys << "," << rs << endl;
    }

    return Disk(Vec2(xs,ys), rs);
}

Mat2x2 tensor_max (const vector<Mat2x2> &Ms) {

    vector<Disk> disks;

    for (int i = 0; i < (int)Ms.size(); ++i) {
        const Mat2x2 &M = Ms[i];
        if (trace(M) == 0)
            continue;
        disks.push_back(Disk(Vec2((M(0,0)-M(1,1))/2, (M(0,1)+M(1,0))/2),
                             (M(0,0)+M(1,1))/2));
    }
    
    Disk disk = minidisk(disks,vector<Disk>(0));

    return disk.c[0]*Mat2x2(Vec2(1,0),Vec2(0,-1))
         + disk.c[1]*Mat2x2(Vec2(0,1),Vec2(1,0))
         + disk.r*Mat2x2(Vec2(1,0),Vec2(0,1));
}
