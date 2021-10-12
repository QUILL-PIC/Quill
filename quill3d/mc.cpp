#include <cmath>
#include "compilation_defines.h"
#include "main.h"
#include "containers.h"

#ifdef QUILL_NOQED
const bool qed_enabled = false;
#else
extern bool qed_enabled;
#endif

double spatial_region::get_rand()
{
    // Метод Фибоначчи с запаздываниями
    double a;
    /* для неограниченного массива было бы
       a = random[n_random-55] - random[n_random-24]; */
    a = random[n_random%55] - random[(n_random+31)%55]; // 55-24=31
    if (a<0) a = a + 1;
    random[n_random%55] = a;
    n_random = n_random < 55 ? n_random + 1 : 0;
    return a;
}

#ifndef QUILL_NOQED
double spatial_region::chi(vector3d& e, vector3d& b, vector3d& u, double& g)
{
    double a1,a2;
    a1 = g*e.x + u.y*b.z - u.z*b.y;
    a2 = a1*a1;
    a1 = g*e.y - u.x*b.z + u.z*b.x;
    a2 += a1*a1;
    a1 = g*e.z + u.x*b.y - u.y*b.x;
    a2 += a1*a1;
    a1 = u.x*e.x + u.y*e.y + u.z*e.z;
    a2 -= a1*a1;
    return sqrt(a2)/e_s;
}

double spatial_region::w(const double& g, double& r, const double& chi)
{ // \int_0^g w \, dg = W, время нормировано на c/\omega
    /* Используются следующие аппроксимации: Ai'(x) \approx -exp(
     * -2/3*x^(3/2) )/(2*sqrt(pi)) * ( x + 0.70861*(1-0.65*x/(1+x^2))
     * )^(1/4); \int_x^\infty Ai(y) dy \approx exp( -2/3*x^(3/2)
     * )/(sqrt(2*pi))/( x + 0.80049 )^(3/4) */
    if ((1-r)*g<1||r==0||chi==0)
    {
        return 0;
    }
    else if (chi>0.133)
    { /* 0.133 соответствует обрезке высоких частот в спектре (для
         классики; здесь же, наоборот, частоты не обрезаются) при
         e^{-5} и переходу к квантовому случаю, если обрезка по
         частотам переваливает за 1/2 всего частотного интервала */
        double d = r/((1-r)*chi);
        double a = pow(d,(double) 1/3);
        double kappa = a*a;
        // 1/( 137*2*sqrt(pi) ) = 0.0020591
        return -0.0020591*e_s/(g*g)*exp(-2*d/3)*( 1/sqrt(sqrt((kappa+0.80049)*(kappa+0.80049)*(kappa+0.80049))) - (2*(1-r)+r*r)*chi*a/r*sqrt(sqrt(kappa+0.70861*(1-0.65*kappa/(1+kappa*kappa)))) );
    }
    else
    {
        double r_max; /* излучение на частотах, больших g*r_max/hbar,
                         не вычисляется */
        r_max = 1/(1+2/(3*5*chi));
        r = r_max*r;
        double d = r/((1-r)*chi);
        double a = pow(d,(double) 1/3);
        double kappa = a*a;
        // 1/( 137*2*sqrt(pi) ) = 0.0020591
        // *rmax, поскольку частотный интервал уменьшен
        return -r_max*0.0020591*e_s/(g*g)*exp(-2*d/3)*( 1/sqrt(sqrt((kappa+0.80049)*(kappa+0.80049)*(kappa+0.80049))) - (2*(1-r)+r*r)*chi*a/r*sqrt(sqrt(kappa+0.70861*(1-0.65*kappa/(1+kappa*kappa)))) );
    }
}

double spatial_region::tilde_w(double& g, double& r, double& chi)
{ // \int_0^g w \, dg = W, время нормировано на c/\omega
    /* Используются следующие аппроксимации: Ai'(x) \approx -exp(
     * -2/3*x^(3/2) )/(2*sqrt(pi)) * ( x + 0.70861*(1-0.65*x/(1+x^2))
     * )^(1/4); \int_x^\infty Ai(y) dy \approx exp( -2/3*x^(3/2)
     * )/(sqrt(2*pi))/( x + 0.80049 )^(3/4) */
    if (chi==0||r*g<1||(1-r)*g<1)
    {
        return 0;
    }
    else
    {
        double a;
        double d;
        double kappa;
        d = 1/(r*(1-r)*chi);
        a = pow(d,(double)1/3);
        kappa = a*a;
        // 1/( 137*2*sqrt(pi) ) = 0.0020591
        return 0.0020591*e_s/(g*g)*exp(-2*d/3)*(1/sqrt(sqrt((kappa+0.80049)*(kappa+0.80049)*(kappa+0.80049))) - (2*r*(1-r)-1)*chi*a*sqrt(sqrt(kappa+0.70861*(1-0.65*kappa/(1+kappa*kappa)))));
    }
}

double spatial_region::mathcal_W(vector3d& e, vector3d& b)
{ /* вероятность рождения электрон-позитронной пары из вакуума в
     единичном объёме в единицу времени, нормированная на \omega^4/c^3
   */
    double a,f,g,xi,theta;
    f = 0.5*( e.x*e.x + e.y*e.y + e.z*e.z - b.x*b.x - b.y*b.y - b.z*b.z );
    g = e.x*b.x + e.y*b.y + e.z*b.z;
    a = sqrt(f*f+g*g);
    xi = sqrt(a+f);
    theta = sqrt(a-f);
    if (xi!=0)
        if (theta!=0)
            return e_s*e_s/(4*PI)*xi*theta/tanh(PI*theta/xi)*exp(-PI*e_s/xi);
        else
            return e_s*e_s/(4*PI)*xi*xi/PI*exp(-PI*e_s/xi);
    else
        return 0;
}

particle* spatial_region::bear_particle(double cmr, vector3d& position, vector3d& direction, double g, double chi, double q)
{
    particle* p;
    double a;
    p = new_particle();
    p->cmr = cmr;
    p->q = q;
    p->g = g;
    p->chi = chi;
    p->x = position.x;
    p->y = position.y;
    p->z = position.z;
    if (cmr==0)
        a = g/sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);
    else
        a = sqrt((g*g-1)/(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z));
    p->ux = a*direction.x;
    p->uy = a*direction.y;
    p->uz = a*direction.z;
    p->trn = 0;
    return p;
}
#endif

void spatial_region::birth_from_vacuum(double q)
{
    if (!qed_enabled)
        return;
    #ifndef QUILL_NOQED
    vector3d e,b;
    vector3d position;
    vector3d direction;
    direction.x = 1;
    direction.y = 0;
    direction.z = 0;
    double a = -dt/2;
    particle* p1;
    particle* p2;
    // границы взяты из ndfx для поля e
    for(int i=2;i<nx-2;i++)
    {
        for(int j=2;j<ny-2;j++)
        {
            for(int k=2;k<nz-2;k++)
            {
                /* 0.5 - чтобы избежать пересечений, вызванных очень
                 * незначительными смещениями */
                position.x = i + 0.5;
                position.y = j + 0.5;
                position.z = k + 0.5;
                e = e_to_particle(position.x,position.y,position.z);
                b = b_to_particle(position.x,position.y,position.z);
                if (get_rand()<mathcal_W(e,b)*dt*dx*dy*dz)
                {
                    /* начальное направление движения безразлично при
                     * g = 1 */
                    p1 = bear_particle(-1,position,direction,1,0,-q);
                    advance_momentum(*p1, e, b, a);
                    p2 = bear_particle(1,position,direction,1,0,q);
                    advance_momentum(*p2, e, b, a);
                    p1->next = p2;
                    p2->previous = p1;
                    if (cp[i][j][k].pl.head==0)
                    {
                        cp[i][j][k].pl.head = p1;
                    }
                    else
                    {
                        p2->next = cp[i][j][k].pl.head;
                        cp[i][j][k].pl.head->previous = p2;
                        cp[i][j][k].pl.head = p1;
                    }
                    cp[i][j][k].pl.start = cp[i][j][k].pl.head;
                }
            }
        }
    }
    #endif
}

double spatial_region::_ppd(double q, double q0) {
    return 1/(1+pow(2*q/q0,4));
}

void spatial_region::pmerging(double* ppd, string pmerging)
{
    /* Из-за флуктуаций отношение удаляемого заряда к полному не равно
     * ppd, поэтому для сохранения нейтральности приходится считать
     * это отношение для каждого сорта частиц отдельно и использовать
     * именно его при при увеличении заряда */
    particle* current;
    particle* tmp;
    double* qs = new double[3+n_ion_populations];
    double* a = new double[3+n_ion_populations];
    for (int i=0;i<3+n_ion_populations;i++) {
        qs[i] = 0;
        a[i] = 0;
    }
    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                current = cp[i][j][k].pl.head;
                while(current!=0) {
                    if (current->cmr==-1)
                        a[0] -= current->q;
                    else if (current->cmr==1)
                        a[1] += current->q;
                    else if (current->cmr==0)
                        a[2] += current->q;
                    else {
                        for (int ii=0;ii<n_ion_populations;ii++) {
                            if (current->cmr==icmr[ii])
                                a[3+ii] += current->q;
                        }
                    }
                    current = current->next;
                }
            }
        }
    }
    if (pmerging=="ti") {
        for(int i=0;i<nx;i++) {
            for(int j=0;j<ny;j++) {
                for(int k=0;k<nz;k++) {
                    current = cp[i][j][k].pl.head;
                    while(current!=0) {
                        if (get_rand()<ppd[0]) {
                            tmp = current->next;
                            if (current->cmr==-1)
                                qs[0] -= current->q;
                            else if (current->cmr==1)
                                qs[1] += current->q;
                            else if (current->cmr==0)
                                qs[2] += current->q;
                            else {
                                for (int ii=0;ii<n_ion_populations;ii++) {
                                    if (current->cmr==icmr[ii])
                                        qs[3+ii] += current->q;
                                }
                            }
                            //
                            if (current->next!=0)
                                (current->next)->previous = current->previous;
                            //
                            if (current->previous!=0)
                                (current->previous)->next = current->next;
                            else {
                                cp[i][j][k].pl.head = current->next;
                                cp[i][j][k].pl.start = cp[i][j][k].pl.head;
                            }
                            delete_particle(current);
                            //
                            current = tmp;
                        }
                        else
                            current = current->next;
                    }
                }
            }
        }
        //
        if (a[0]-qs[0]!=0)
            a[0] = a[0]/(a[0]-qs[0]);
        if (a[1]-qs[1]!=0)
            a[1] = a[1]/(a[1]-qs[1]);
        if (a[2]-qs[2]!=0)
            a[2] = a[2]/(a[2]-qs[2]);
        for (int i=0;i<n_ion_populations;i++) {
            if (a[3+i]-qs[3+i]!=0)
                a[3+i] = a[3+i]/(a[3+i]-qs[3+i]);
        }
        for(int i=0;i<nx;i++) {
            for(int j=0;j<ny;j++) {
                for(int k=0;k<nz;k++) {
                    current = cp[i][j][k].pl.head;
                    while(current!=0) {
                        if (current->cmr==-1)
                            current->q = a[0]*current->q;
                        else if (current->cmr==1)
                            current->q = a[1]*current->q;
                        else if (current->cmr==0)
                            current->q = a[2]*current->q;
                        else {
                            for (int ii=0;ii<n_ion_populations;ii++) {
                                if (current->cmr==icmr[ii])
                                    current->q = a[3+ii]*current->q;
                            }
                        }
                        current = current->next;
                    }
                }
            }
        }
    } else if (pmerging=="nl") {
        double* av_q = new double[3+n_ion_populations];
        for(int i=0;i<3+n_ion_populations;i++)
            av_q[i] = 0;
        for(int i=0;i<nx;i++) {
            for(int j=0;j<ny;j++) {
                for(int k=0;k<nz;k++) {
                    current = cp[i][j][k].pl.head;
                    while(current!=0) {
                        if (current->cmr==-1)
                            av_q[0] += 1;
                        else if (current->cmr==1)
                            av_q[1] += 1;
                        else if (current->cmr==0)
                            av_q[2] += 1;
                        else {
                            for (int ii=0;ii<n_ion_populations;ii++) {
                                if (current->cmr==icmr[ii])
                                    av_q[3+ii] += 1;
                            }
                        }
                        current = current->next;
                    }
                }
            }
        }
        for (int i=0;i<3+n_ion_populations;i++)
            av_q[i] = a[i]/av_q[i];
        for(int i=0;i<nx;i++) {
            for(int j=0;j<ny;j++) {
                for(int k=0;k<nz;k++) {
                    current = cp[i][j][k].pl.head;
                    while(current!=0) {
                        tmp = current->next;
                        if (current->cmr==-1) {
                            if (get_rand()<ppd[0]*_ppd(-current->q,av_q[0])) {
                                if (current->next!=0)
                                    (current->next)->previous = current->previous;
                                if (current->previous!=0)
                                    (current->previous)->next = current->next;
                                else {
                                    cp[i][j][k].pl.head = current->next;
                                    cp[i][j][k].pl.start = current->next;
                                }
                                delete_particle(current);
                            } 
                        } else if (current->cmr==1) {
                            if (get_rand()<ppd[1]*_ppd(current->q,av_q[1])) {
                                if (current->next!=0)
                                    (current->next)->previous = current->previous;
                                if (current->previous!=0)
                                    (current->previous)->next = current->next;
                                else {
                                    cp[i][j][k].pl.head = current->next;
                                    cp[i][j][k].pl.start = current->next;
                                }
                                delete_particle(current);
                            } 
                        } else if (current->cmr==0) {
                            if (get_rand()<ppd[2]*_ppd(current->q,av_q[2])) {
                                if (current->next!=0)
                                    (current->next)->previous = current->previous;
                                if (current->previous!=0)
                                    (current->previous)->next = current->next;
                                else {
                                    cp[i][j][k].pl.head = current->next;
                                    cp[i][j][k].pl.start = current->next;
                                }
                                delete_particle(current);
                            } 
                        } else {
                            for (int ii=0;ii<n_ion_populations;ii++) {
                                if (current->cmr==icmr[ii]) {
                                    if (get_rand()<ppd[3+ii]*_ppd(current->q,av_q[3+ii])) {
                                        if (current->next!=0)
                                            (current->next)->previous = current->previous;
                                        if (current->previous!=0)
                                            (current->previous)->next = current->next;
                                        else {
                                            cp[i][j][k].pl.head = current->next;
                                            cp[i][j][k].pl.start = current->next;
                                        }
                                        delete_particle(current);
                                    } 
                                }
                            }
                        }
                        current = tmp;
                    }
                }
            }
        }
        for(int i=0;i<nx;i++) {
            for(int j=0;j<ny;j++) {
                for(int k=0;k<nz;k++) {
                    current = cp[i][j][k].pl.head;
                    while(current!=0) {
                        if (current->cmr==-1) {
                            current->q = current->q/(1-ppd[0]*_ppd(-current->q,av_q[0]));
                            qs[0] -= current->q;
                        } else if (current->cmr==1) {
                            current->q = current->q/(1-ppd[1]*_ppd(current->q,av_q[1]));
                            qs[1] += current->q;
                        } else if (current->cmr==0) {
                            current->q = current->q/(1-ppd[2]*_ppd(current->q,av_q[2]));
                            qs[2] += current->q;
                        } else {
                            for (int ii=0;ii<n_ion_populations;ii++) {
                                if (current->cmr==icmr[ii]) {
                                    current->q = current->q/(1-ppd[3+ii]*_ppd(current->q,av_q[3+ii]));
                                    qs[3+ii] += current->q;
                                }
                            }
                        }
                        current = current->next;
                    }
                }
            }
        }
        for(int i=0;i<3+n_ion_populations;i++) {
            if (qs[i]!=0)
                a[i] = a[i]/qs[i];
            else
                a[i] = 0;
        }
        for(int i=0;i<nx;i++) {
            for(int j=0;j<ny;j++) {
                for(int k=0;k<nz;k++) {
                    current = cp[i][j][k].pl.head;
                    while(current!=0) {
                        if (current->cmr==-1)
                            current->q = a[0]*current->q;
                        else if (current->cmr==1)
                            current->q = a[1]*current->q;
                        else if (current->cmr==0)
                            current->q = a[2]*current->q;
                        else {
                            for (int ii=0;ii<n_ion_populations;ii++) {
                                if (current->cmr==icmr[ii])
                                    current->q = a[3+ii]*current->q;
                            }
                        }
                        current = current->next;
                    }
                }
            }
        }
        delete[] av_q;
    }
    delete[] qs;
    delete[] a;
}
