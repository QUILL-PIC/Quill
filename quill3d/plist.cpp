#include "main.h"

void spatial_region::place(spatial_region::plist::particle& p, int& ci, int& cj, int& ck)
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

void spatial_region::place(spatial_region::plist::particle& p)
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

void spatial_region::copy(plist& b, plist& a)
{
    erase(a);
    //
    plist::particle* current;
    plist::particle* tmp;
    current = b.head;
    //
    a.start = 0;
    while(current!=0)
    {
	tmp = new_particle();
	*tmp = *current;
	tmp->previous = a.start;
	a.start = tmp;
	current = current->next;
    }
    a.head = 0;
    while(a.start!=0)
    {
	a.start->next = a.head;
	a.head = a.start;
	a.start = a.start->previous;
    }
    a.start = a.head;
}

void spatial_region::erase(plist& a)
{
    plist::particle* current;
    plist::particle* tmp;
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

void copy(spatial_region& a, int a1, int a2, int a3, spatial_region& b, int b1, int b2, int b3)
{ // copy fields from a to b
    b.ce[b1][b2][b3].ex = a.ce[a1][a2][a3].ex;
    b.ce[b1][b2][b3].ey = a.ce[a1][a2][a3].ey;
    b.ce[b1][b2][b3].ez = a.ce[a1][a2][a3].ez;
    b.cb[b1][b2][b3].bx = a.cb[a1][a2][a3].bx;
    b.cb[b1][b2][b3].by = a.cb[a1][a2][a3].by;
    b.cb[b1][b2][b3].bz = a.cb[a1][a2][a3].bz;
    b.cj[b1][b2][b3].jx = a.cj[a1][a2][a3].jx;
    b.cj[b1][b2][b3].jy = a.cj[a1][a2][a3].jy;
    b.cj[b1][b2][b3].jz = a.cj[a1][a2][a3].jz;
    b.cbe[b1][b2][b3].bex = a.cbe[a1][a2][a3].bex;
    b.cbe[b1][b2][b3].bey = a.cbe[a1][a2][a3].bey;
    b.cbe[b1][b2][b3].bez = a.cbe[a1][a2][a3].bez;
}
