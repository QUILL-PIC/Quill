#ifndef THINPARTICLE_H_
#define THINPARTICLE_H_

struct thinparticle
{
    double x, y, z, ux, uy, uz, g, q, cmr, chi; // cmr - charge to mass ratio
    int trn; // track name
    thinparticle(): x(0), y(0), z(0), ux(0), uy(0), uz(0), g(1), q(0), cmr(-1), chi(0), trn(0) {}
};

#endif /* THINPARTICLE_H_ */

