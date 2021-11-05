#include "main.h"
#include "containers.h"

extern int n_sr;
extern int nm;

void spatial_region::place(particle& p, int& ci, int& cj, int& ck)
{
    // ci, cj, ck - cell there particle was
    if(((int)p.x)!=ci||((int)p.y)!=cj||((int)p.z)!=ck)
    {
        if (p.previous!=0)
            p.previous->next = p.next;
        else
            cp[ci][cj][ck].pl.head = p.next;
        if (p.next!=0)
            p.next->previous = p.previous;
        place(p);
    }
}

void spatial_region::place(particle& p)
{ /* добавляет частицу с previous = next = 0, не привязанную к
     какой-либо ячейке (ни один из указателей head не равен &p) */
    if (is_inside(p.x,p.y,p.z))
    {
        p.previous = 0; /* new_particle не обязательно возвращает
                           память, занятую нулями, поэтому при
                           размещении таких частиц p.previous следует
                           обнулять */
        p.next = cp[(int)p.x][(int)p.y][(int)p.z].pl.head;
        if(p.next!=0)
            (p.next)->previous = &p;
        cp[(int)p.x][(int)p.y][(int)p.z].pl.head = &p;
    }
    else
    {
        delete_particle(&p);
    }
}

bool spatial_region::is_inside(int i,int j, int k)
{
    if (i>0&&i<nx-2&&j>0&&j<ny-2&&k>0&&k<nz-2)
        return 1;
    else
        return 0;
}

// Indicates if the point is inside the whole box (not only the current spatial region)
bool spatial_region::is_inside_global(int i, int j, int k)
{
    if (sr_id == 0 && i<=0)
        return 0;
    if (sr_id == n_sr-1 && i>=nx-2)
        return 0;
        
    if (j>0 && j<ny-2 && k>0 && k<nz-2)
        return 1;
    else
        return 0;    
}

bool spatial_region::is_in_exchange_area(int i)
{
    if (sr_id > 0 && i >= 0 && i < nm)
        return true;
    if (sr_id < n_sr-1 && i >= nx-nm && i < nx)
        return true;
    return false;
}

void spatial_region::p_boundary()
{
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            erase(cp[i][j][0].pl);
            erase(cp[i][j][nz-1].pl);
            erase(cp[i][j][nz-2].pl);
        }
    }
    for(int i=0;i<nx;i++)
    {
        for(int k=0;k<nz;k++)
        {
            erase(cp[i][0][k].pl);
            erase(cp[i][ny-1][k].pl);
            erase(cp[i][ny-2][k].pl);
        }
    }
    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            erase(cp[0][j][k].pl);
            erase(cp[nx-1][j][k].pl);
            erase(cp[nx-2][j][k].pl);
        }
    }
}

void spatial_region::erase(plist& a)
{
    particle* current;
    particle* tmp;
    current = a.head;
    while(current!=0)
    {
        tmp = current->next;
        delete_particle(current);
        current = tmp;
    }
    a.head = 0;
    a.start = 0;
}
