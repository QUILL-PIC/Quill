#include <cmath>
#include <memory>
#include "main.h"
#include "maxwell.h"

extern bool qed_enabled;
extern void (*pusher)(spatial_region::plist::particle* p, vector3d& e_field, vector3d& b_field, double& dt);
extern std::shared_ptr<maxwell_solver> solver;

void spatial_region::fadvance()
{
    solver->fadvance(ce, cb, cj, dt, dx, dy, dz, nx, ny, nz);

    // b in locations of e for magnetic field interpolation
    for(int i=2;i<nx-1;i++)
    {
        for(int j=2;j<ny-1;j++)
        {
            for(int k=2;k<nz-1;k++)
            {
                cbe[i][j][k].bex = 0.125*(cb[i][j][k].bx + cb[i+1][j][k].bx + cb[i][j][k-1].bx + cb[i+1][j][k-1].bx + cb[i][j-1][k].bx + cb[i+1][j-1][k].bx + cb[i][j-1][k-1].bx + cb[i+1][j-1][k-1].bx);
                cbe[i][j][k].bey = 0.125*(cb[i][j][k].by + cb[i-1][j][k].by + cb[i][j][k-1].by + cb[i-1][j][k-1].by + cb[i][j+1][k].by + cb[i-1][j+1][k].by + cb[i][j+1][k-1].by + cb[i-1][j+1][k-1].by);
                cbe[i][j][k].bez = 0.125*(cb[i][j][k].bz + cb[i-1][j][k].bz + cb[i][j][k+1].bz + cb[i-1][j][k+1].bz + cb[i][j-1][k].bz + cb[i-1][j-1][k].bz + cb[i][j-1][k+1].bz + cb[i-1][j-1][k+1].bz);
            }
        }
    }
}

void spatial_region::f_zeroing_on_boundaries()
{
    solver->fzeroing_on_boundaries(ce, cb, cj, dt, dx, dy, dz, nx, ny, nz);
}

void spatial_region::padvance(bool freezing)
{
    spatial_region::plist::particle* current;
    spatial_region::plist::particle* tmp;
    vector3d displ;
    vector3d e;
    vector3d b;
    vector3d u_interim;
    vector3d u_prev;
    double g_prev,chi_prev,g_interim;
    double a, r;
    bool calc_qed;
    spatial_region::plist::particle* born;
    int_vector3d position;
    //
    p_boundary();
    //
    N_freezed = 0;
    //
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                cj[i][j][k].jx = 0;
                cj[i][j][k].jy = 0;
                cj[i][j][k].jz = 0;
            }
        }
    }
    //
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                current = cp[i][j][k].pl.start;
                while(current!=0)
                {
                    e = e_to_particle(current->x,current->y,current->z);
                    b = b_to_particle(current->x,current->y,current->z);
                    /* Теперь определяется импульс частицы в момент
                     * времени t_{n+1} (для электронов и позитронов),
                     * а координата в момент времени t_{n+3/2} ещё не
                     * определена. Зная u_{n+1} и u_n, можно вычислить
                     * поперечную силу при t=t_{n+1/2}, за счёт
                     * хранения поперечной силы на предыдущем шаге - и
                     * поперечную силу в момент времени t_n. Таким
                     * образом, можно найти w_n и произвести
                     * излучение/рождение пар в момент времени t_n.
                     * Импульсы частиц, участвовавших в квантовых
                     * процессах, приходится заново продвигать из
                     * точки t_n в точку t_{n+1} (а также заново
                     * пересчитывать chi_{n+1/2}). Рождённым частицам
                     * присваиваются такие координаты (при t_{n+3/2},
                     * поскольку данный цикл по ним не пробежит),
                     * будто при t=t_n их координаты совпадали с
                     * координатами родившей их частицы.
                     */
                    if (current->cmr==0) // *current - фотон
                    {
                        if (qed_enabled) {
                            u_interim.x = current->ux;
                            u_interim.y = current->uy;
                            u_interim.z = current->uz;
                            g_prev = current->g;
                            chi_prev = current->chi;
                            current->chi = chi(e, b, u_interim, current->g); // chi_{n+1/2}
                            a = 0.5 * (current->chi + chi_prev);
                            r = get_rand();
                            calc_qed = !(get_rand()>=tilde_w(g_prev,r,a)*g_prev*dt);
                        } else {
                            calc_qed = false;
                        }

                        if (!calc_qed)
                        {
                            displ = current->get_displacement(dt);
                            displ.x=displ.x/dx;
                            displ.y=displ.y/dy;
                            displ.z=displ.z/dz;
                            current->coordinate_advance(displ);
                            tmp = current;
                            current = current->next;
                            place(*tmp,i,j,k);
                        }
                        else
                        { // pair production
                            /* Координаты заряженных частиц не
                             * пересчитываются в точку n и обратно в
                             * n+1/2. Пара рождается с координатами
                             * фотона x^{n+1/2} для сохранения
                             * нейтральности и уменьшения
                             * электромагнитных шумов */
                            // electron
                            displ.x = current->x;
                            displ.y = current->y;
                            displ.z = current->z;
                            born = bear_particle(-1,displ,u_interim,r*g_prev,r*current->chi,-current->q);
                            u_prev.x = born->ux;
                            u_prev.y = born->uy;
                            u_prev.z = born->uz;
                            g_prev = born->g;
                            chi_prev = born->chi;
                            born->momentum_advance(e,b,dt); // p_n -> p_{n+1}
                            u_interim.x = 0.5*(born->ux+u_prev.x);
                            u_interim.y = 0.5*(born->uy+u_prev.y);
                            u_interim.z = 0.5*(born->uz+u_prev.z);
                            g_interim = 0.5*(born->g+g_prev);
                            born->chi = chi(e,b,u_interim,g_interim); // chi_{n+1/2}
                            displ = born->get_displacement(dt);
                            displ.x=displ.x/dx;
                            displ.y=displ.y/dy;
                            displ.z=displ.z/dz;
                            jdeposition(*born,displ);
                            place(*born);
                            // positron
                            displ.x = current->x;
                            displ.y = current->y;
                            displ.z = current->z;
                            u_interim.x = current->ux;
                            u_interim.y = current->uy;
                            u_interim.z = current->uz;
                            born = bear_particle(1,displ,u_interim,(1-r)*current->g,(1-r)*current->chi,current->q);
                            u_prev.x = born->ux;
                            u_prev.y = born->uy;
                            u_prev.z = born->uz;
                            g_prev = born->g;
                            chi_prev = born->chi;
                            born->momentum_advance(e,b,dt); // p_n -> p_{n+1}
                            u_interim.x = 0.5*(born->ux+u_prev.x);
                            u_interim.y = 0.5*(born->uy+u_prev.y);
                            u_interim.z = 0.5*(born->uz+u_prev.z);
                            g_interim = 0.5*(born->g+g_prev);
                            born->chi = chi(e,b,u_interim,g_interim); // chi_{n+1/2}
                            displ = born->get_displacement(dt);
                            displ.x=displ.x/dx;
                            displ.y=displ.y/dy;
                            displ.z=displ.z/dz;
                            jdeposition(*born,displ);
                            place(*born);
                            // photon destruction
                            tmp = current->next;
                            if (current->previous!=0)
                                (current->previous)->next = current->next;
                            else
                                cp[i][j][k].pl.head = current->next;
                            if (current->next!=0)
                                (current->next)->previous = current->previous;
                            delete_particle(current);
                            current = tmp;
                        }
                    }
                    else if (fabs(current->cmr)==1)
                    { // *current - электрон или позитрон
                        if (freezing==1&&(b.x*b.x+b.y*b.y+b.z*b.z-e.x*e.x-e.y*e.y-e.z*e.z)/(current->g*current->g)>0.1/(dt*dt))
                        {
                            current->ux = 0;
                            current->uy = 0;
                            current->uz = 0;
                            current->g = 1;
                            N_freezed += fabs(current->q);
                            tmp = current;
                            current = current->next;
                        }
                        else
                        {
                            if (qed_enabled) {
                                u_prev.x = current->ux;
                                u_prev.y = current->uy;
                                u_prev.z = current->uz;
                                g_prev = current->g;
                                chi_prev = current->chi;
                                current->momentum_advance(e, b, dt); // p_n -> p_{n+1}
                                u_interim.x = 0.5 * (current->ux + u_prev.x);
                                u_interim.y = 0.5 * (current->uy + u_prev.y);
                                u_interim.z = 0.5 * (current->uz + u_prev.z);
                                g_interim = 0.5 * (current->g + g_prev);
                                current->chi = chi(e, b, u_interim, g_interim); // chi_{n+1/2}
                                a = 0.5 * (current->chi + chi_prev); // chi_n
                                r = get_rand();
                                calc_qed = !(get_rand()>=w(g_prev,r,a)*g_prev*dt);
                            } else {
                                current->momentum_advance(e, b, dt); // p_n -> p_{n+1}
                                calc_qed = false;
                            }

                            if (!calc_qed)
                            {
                                displ = current->get_displacement(dt);
                                displ.x=displ.x/dx;
                                displ.y=displ.y/dy;
                                displ.z=displ.z/dz;
                                jdeposition(*current,displ);
                                tmp = current;
                                current = current->next;
                                place(*tmp,i,j,k);
                            }
                            else
                            { // photon emission
                                a = 0.5*dt/g_prev;
                                displ.x = current->x - a*u_prev.x/dx;
                                displ.y = current->y - a*u_prev.y/dy;
                                displ.z = current->z - a*u_prev.z/dz;
                                born = bear_particle(0,displ,u_prev,r*g_prev,r*current->chi,fabs(current->q));
                                a = 3*a*g_prev/born->g;
                                born->x += a*born->ux/dx;
                                born->y += a*born->uy/dy;
                                born->z += a*born->uz/dz;
                                place(*born);
                                //
                                current->chi = (1-r)*current->chi;
                                current->g = (1-r)*g_prev;
                                a = sqrt((current->g*current->g-1)/(u_prev.x*u_prev.x + u_prev.y*u_prev.y + u_prev.z*u_prev.z));
                                current->ux = u_prev.x*a;
                                current->uy = u_prev.y*a;
                                current->uz = u_prev.z*a;
                                a = 0.5*dt/g_prev;
                                current->x -= a*u_prev.x/dx;
                                current->y -= a*u_prev.y/dy;
                                current->z -= a*u_prev.z/dz;
                                a = a*g_prev/current->g;
                                current->x += a*current->ux/dx;
                                current->y += a*current->uy/dy;
                                current->z += a*current->uz/dz;
                                //
                                tmp = current->next;
                                //
                                if (is_inside(current->x,current->y,current->z))
                                {
                                    place(*current,i,j,k);
                                    //
                                    position.i = current->x; // where particle was; for place(...)
                                    position.j = current->y;
                                    position.k = current->z;
                                    //
                                    e = e_to_particle(current->x,current->y,current->z);
                                    b = b_to_particle(current->x,current->y,current->z);
                                    current->momentum_advance(e,b,dt); // p_n -> p_{n+1}
                                    displ = current->get_displacement(dt);
                                    displ.x=displ.x/dx;
                                    displ.y=displ.y/dy;
                                    displ.z=displ.z/dz;
                                    jdeposition(*current,displ);
                                    //
                                    place(*current,position.i,position.j,position.k);
                                }
                                else
                                {
                                    if (current->previous!=0)
                                        (current->previous)->next = current->next;
                                    else
                                        cp[i][j][k].pl.head = current->next;
                                    if (current->next!=0)
                                        (current->next)->previous = current->previous;
                                    delete_particle(current);
                                }
                                current = tmp;
                            }
                        }
                    }
                    else
                    { // *current - ион
                        current->momentum_advance(e,b,dt); // p_n -> p_{n+1}
                    displ = current->get_displacement(dt);
                    displ.x=displ.x/dx;
                    displ.y=displ.y/dy;
                    displ.z=displ.z/dz;
                    jdeposition(*current,displ);
                    tmp = current;
                    current = current->next;
                    place(*tmp,i,j,k);
                    }
                }
            }
        }
    }
    // adjustment of the start pointer
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                cp[i][j][k].pl.start = cp[i][j][k].pl.head;
            }
        }
    }
}

void spatial_region::moving_window(int l, int nmw, double mwspeed)
{
    spatial_region::plist::particle* current;
    //
    if((l+1)*dt*mwspeed>nmw*dx)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                erase(cp[0][j][k].pl);
            }
        }
        //
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx-1;i++)
            {
                for(int k=0;k<nz;k++)
                {
                    ce[i][j][k].ex = ce[i+1][j][k].ex;
                    ce[i][j][k].ey = ce[i+1][j][k].ey;
                    ce[i][j][k].ez = ce[i+1][j][k].ez;
                    cb[i][j][k].bx = cb[i+1][j][k].bx;
                    cb[i][j][k].by = cb[i+1][j][k].by;
                    cb[i][j][k].bz = cb[i+1][j][k].bz;
                    cbe[i][j][k].bex = cbe[i+1][j][k].bex;
                    cbe[i][j][k].bey = cbe[i+1][j][k].bey;
                    cbe[i][j][k].bez = cbe[i+1][j][k].bez;
                    cp[i][j][k].pl = cp[i+1][j][k].pl;
                }
            }
        }
        //
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                ce[nx-1][j][k].ex = 0;
                ce[nx-1][j][k].ey = 0;
                ce[nx-1][j][k].ez = 0;
                cb[nx-1][j][k].bx = 0;
                cb[nx-1][j][k].by = 0;
                cb[nx-1][j][k].bz = 0;
                cbe[nx-1][j][k].bex = 0;
                cbe[nx-1][j][k].bey = 0;
                cbe[nx-1][j][k].bez = 0;
                cp[nx-1][j][k].pl.head = 0;
                cp[nx-1][j][k].pl.start = 0;
            }
        }
        //
        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                for(int k=0;k<nz;k++)
                {
                    current = cp[i][j][k].pl.start;
                    while(current!=0)
                    {
                        current->x = current->x - 1;
                        current = current->next;
                    }
                }
            }
        }
        //
    }
}

void spatial_region::plist::xplus(double x_adjunct)
{
    spatial_region::plist::particle* current;
    current = head;
    while(current!=0)
    {
        current->x = current->x + x_adjunct;
        current = current->next;
    }
}

vector3d spatial_region::plist::particle::get_displacement(double& dt)
{
    vector3d displacement;
    double a;
    a = dt/g;
    displacement.x = ux*a;
    displacement.y = uy*a;
    displacement.z = uz*a;
    return displacement;
}

void spatial_region::plist::particle::coordinate_advance(vector3d& a)
{
    x = x + a.x;
    y = y + a.y;
    z = z + a.z;
}

void spatial_region::plist::particle::momentum_advance(vector3d& e_field, vector3d& b_field, double& dt)
{
    pusher(this, e_field, b_field, dt);
}
