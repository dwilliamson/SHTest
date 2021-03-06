#ifndef lint
static const char	RCSid[] = "$Id: calfunc.c,v 2.10 2003/03/05 16:16:52 greg Exp $";
#endif
/*
 *  calfunc.c - routines for calcomp using functions.
 *
 *      If VARIABLE is not set, only library functions
 *  can be accessed.
 *
 *  2/19/03	Eliminated conditional compiles in favor of esupport extern.
 */

#include "copyright.h"

#include  <stdio.h>

#include  <errno.h>

#include  <math.h>

#include  "calcomp.h"

				/* bits in argument flag (better be right!) */
#define  AFLAGSIZ	(8*sizeof(unsigned long))
#define  ALISTSIZ	6	/* maximum saved argument list */

typedef struct activation {
    char  *name;		/* function name */
    struct activation  *prev;	/* previous activation */
    double  *ap;		/* argument list */
    unsigned long  an;		/* computed argument flags */
    EPNODE  *fun;		/* argument function */
}  ACTIVATION;		/* an activation record */

static ACTIVATION  *curact = NULL;

static double  libfunc();

#ifndef  MAXLIB
#define  MAXLIB		64	/* maximum number of library functions */
#endif

static double  l_if(), l_select(), l_rand();
static double  l_floor(), l_ceil();
static double  l_sqrt();
static double  l_sin(), l_cos(), l_tan();
static double  l_asin(), l_acos(), l_atan(), l_atan2();
static double  l_exp(), l_log(), l_log10();

			/* functions must be listed alphabetically */
static LIBR  library[MAXLIB] = {
    { "acos", 1, ':', l_acos },
    { "asin", 1, ':', l_asin },
    { "atan", 1, ':', l_atan },
    { "atan2", 2, ':', l_atan2 },
    { "ceil", 1, ':', l_ceil },
    { "cos", 1, ':', l_cos },
    { "exp", 1, ':', l_exp },
    { "floor", 1, ':', l_floor },
    { "if", 3, ':', l_if },
    { "log", 1, ':', l_log },
    { "log10", 1, ':', l_log10 },
    { "rand", 1, ':', l_rand },
    { "select", 1, ':', l_select },
    { "sin", 1, ':', l_sin },
    { "sqrt", 1, ':', l_sqrt },
    { "tan", 1, ':', l_tan },
};

static int  libsize = 16;

#define  resolve(ep)	((ep)->type==VAR?(ep)->v.ln:argf((ep)->v.chan))


int
fundefined(fname)		/* return # of arguments for function */
char  *fname;
{
    register LIBR  *lp;
    register VARDEF  *vp;

    if ((vp = varlookup(fname)) != NULL && vp->def != NULL
		&& vp->def->v.kid->type == FUNC)
	return(nekids(vp->def->v.kid) - 1);
    lp = vp != NULL ? vp->lib : liblookup(fname);
    if (lp == NULL)
	return(0);
    return(lp->nargs);
}


double
funvalue(fname, n, a)		/* return a function value to the user */
char  *fname;
int  n;
double  *a;
{
    ACTIVATION  act;
    register VARDEF  *vp;
    double  rval;
					/* push environment */
    act.name = fname;
    act.prev = curact;
    act.ap = a;
    if (n >= AFLAGSIZ)
	act.an = ~0;
    else
	act.an = (1L<<n)-1;
    act.fun = NULL;
    curact = &act;

    if ((vp = varlookup(fname)) == NULL || vp->def == NULL
		|| vp->def->v.kid->type != FUNC)
	rval = libfunc(fname, vp);
    else
	rval = evalue(vp->def->v.kid->sibling);

    curact = act.prev;			/* pop environment */
    return(rval);
}


void
funset(fname, nargs, assign, fptr)	/* set a library function */
char  *fname;
int  nargs;
int  assign;
double  (*fptr)();
{
    int  oldlibsize = libsize;
    char *cp;
    register LIBR  *lp;
						/* check for context */
    for (cp = fname; *cp; cp++)
	;
    if (cp == fname)
	return;
    if (cp[-1] == CNTXMARK)
	*--cp = '\0';
    if ((lp = liblookup(fname)) == NULL) {	/* insert */
	if (libsize >= MAXLIB) {
	    eputs("Too many library functons!\n");
	    quit(1);
	}
	for (lp = &library[libsize]; lp > library; lp--)
	    if (strcmp(lp[-1].fname, fname) > 0) {
		lp[0].fname = lp[-1].fname;
		lp[0].nargs = lp[-1].nargs;
		lp[0].atyp = lp[-1].atyp;
		lp[0].f = lp[-1].f;
	    } else
		break;
	libsize++;
    }
    if (fptr == NULL) {				/* delete */
	while (lp < &library[libsize-1]) {
	    lp[0].fname = lp[1].fname;
	    lp[0].nargs = lp[1].nargs;
	    lp[0].atyp = lp[1].atyp;
	    lp[0].f = lp[1].f;
	    lp++;
	}
	libsize--;
    } else {					/* or assign */
	lp[0].fname = fname;		/* string must be static! */
	lp[0].nargs = nargs;
	lp[0].atyp = assign;
	lp[0].f = fptr;
    }
    if (libsize != oldlibsize)
	libupdate(fname);			/* relink library */
}


int
nargum()			/* return number of available arguments */
{
    register int  n;

    if (curact == NULL)
	return(0);
    if (curact->fun == NULL) {
	for (n = 0; (1L<<n) & curact->an; n++)
	    ;
	return(n);
    }
    return(nekids(curact->fun) - 1);
}


double
argument(n)			/* return nth argument for active function */
register int  n;
{
    register ACTIVATION  *actp = curact;
    register EPNODE  *ep;
    double  aval;

    if (actp == NULL || --n < 0) {
	eputs("Bad call to argument!\n");
	quit(1);
    }
						/* already computed? */
    if (n < AFLAGSIZ && 1L<<n & actp->an)
	return(actp->ap[n]);

    if (actp->fun == NULL || (ep = ekid(actp->fun, n+1)) == NULL) {
	eputs(actp->name);
	eputs(": too few arguments\n");
	quit(1);
    }
    curact = actp->prev;			/* pop environment */
    aval = evalue(ep);				/* compute argument */
    curact = actp;				/* push back environment */
    if (n < ALISTSIZ) {				/* save value */
	actp->ap[n] = aval;
	actp->an |= 1L<<n;
    }
    return(aval);
}


VARDEF *
argf(n)				/* return function def for nth argument */
int  n;
{
    register ACTIVATION  *actp;
    register EPNODE  *ep;

    for (actp = curact; actp != NULL; actp = actp->prev) {

	if (n <= 0)
	    break;

	if (actp->fun == NULL)
	    goto badarg;

	if ((ep = ekid(actp->fun, n)) == NULL) {
	    eputs(actp->name);
	    eputs(": too few arguments\n");
	    quit(1);
	}
	if (ep->type == VAR)
	    return(ep->v.ln);			/* found it */

	if (ep->type != ARG)
	    goto badarg;

	n = ep->v.chan;				/* try previous context */
    }
    eputs("Bad call to argf!\n");
    quit(1);

badarg:
    eputs(actp->name);
    eputs(": argument not a function\n");
    quit(1);
}


char *
argfun(n)			/* return function name for nth argument */
int  n;
{
    return(argf(n)->name);
}


double
efunc(ep)				/* evaluate a function */
register EPNODE  *ep;
{
    ACTIVATION  act;
    double  alist[ALISTSIZ];
    double  rval;
    register VARDEF  *dp;
					/* push environment */
    dp = resolve(ep->v.kid);
    act.name = dp->name;
    act.prev = curact;
    act.ap = alist;
    act.an = 0;
    act.fun = ep;
    curact = &act;

    if (dp->def == NULL || dp->def->v.kid->type != FUNC)
	rval = libfunc(act.name, dp);
    else
	rval = evalue(dp->def->v.kid->sibling);
    
    curact = act.prev;			/* pop environment */
    return(rval);
}


LIBR *
liblookup(fname)		/* look up a library function */
char  *fname;
{
    int  upper, lower;
    register int  cm, i;

    lower = 0;
    upper = cm = libsize;

    while ((i = (lower + upper) >> 1) != cm) {
	cm = strcmp(fname, library[i].fname);
	if (cm > 0)
	    lower = i;
	else if (cm < 0)
	    upper = i;
	else
	    return(&library[i]);
	cm = i;
    }
    return(NULL);
}


/*
 *  The following routines are for internal use:
 */


static double
libfunc(fname, vp)			/* execute library function */
char  *fname;
VARDEF  *vp;
{
    register LIBR  *lp;
    double  d;
    int  lasterrno;

    if (vp != NULL)
	lp = vp->lib;
    else
	lp = liblookup(fname);
    if (lp == NULL) {
	eputs(fname);
	eputs(": undefined function\n");
	quit(1);
    }
    lasterrno = errno;
    errno = 0;
    d = (*lp->f)(lp->fname);
#ifdef  IEEE
    if (errno == 0)
	if (isnan(d))
	    errno = EDOM;
	else if (isinf(d))
	    errno = ERANGE;
#endif
    if (errno == EDOM || errno == ERANGE) {
	wputs(fname);
	if (errno == EDOM)
		wputs(": domain error\n");
	else if (errno == ERANGE)
		wputs(": range error\n");
	else
		wputs(": error in call\n");
	return(0.0);
    }
    errno = lasterrno;
    return(d);
}


/*
 *  Library functions:
 */


static double
l_if()			/* if(cond, then, else) conditional expression */
			/* cond evaluates true if greater than zero */
{
    if (argument(1) > 0.0)
	return(argument(2));
    else
	return(argument(3));
}


static double
l_select()		/* return argument #(A1+1) */
{
	register int  n;

	n = argument(1) + .5;
	if (n == 0)
		return(nargum()-1);
	if (n < 1 || n > nargum()-1) {
		errno = EDOM;
		return(0.0);
	}
	return(argument(n+1));
}


static double
l_rand()		/* random function between 0 and 1 */
{
    double  x;

    x = argument(1);
    x *= 1.0/(1.0 + x*x) + 2.71828182845904;
    x += .785398163397447 - floor(x);
    x = 1e5 / x;
    return(x - floor(x));
}


static double
l_floor()		/* return largest integer not greater than arg1 */
{
    return(floor(argument(1)));
}


static double
l_ceil()		/* return smallest integer not less than arg1 */
{
    return(ceil(argument(1)));
}


static double
l_sqrt()
{
    return(sqrt(argument(1)));
}


static double
l_sin()
{
    return(sin(argument(1)));
}


static double
l_cos()
{
    return(cos(argument(1)));
}


static double
l_tan()
{
    return(tan(argument(1)));
}


static double
l_asin()
{
    return(asin(argument(1)));
}


static double
l_acos()
{
    return(acos(argument(1)));
}


static double
l_atan()
{
    return(atan(argument(1)));
}


static double
l_atan2()
{
    return(atan2(argument(1), argument(2)));
}


static double
l_exp()
{
    return(exp(argument(1)));
}


static double
l_log()
{
    return(log(argument(1)));
}


static double
l_log10()
{
    return(log10(argument(1)));
}
