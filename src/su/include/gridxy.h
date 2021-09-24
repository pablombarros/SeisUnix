
void gridset(double *gvals, int *errwarn) ; 
void gridrawxycdpic(double *gvals,double dx,double dy,int *icdp,int *igi,int *igc) ;
void gridicrawxy(double *gvals,int igi,int igc,double *dx,double *dy) ;
void gridicgridxy(double *gvals,int igi,int igc,double *dx,double *dy) ;
void gridiccdp(double *gvals,int igi,int igc,int *icdp) ;
void gridcdpic(double *gvals,int icdp,int *igi,int *igc) ;
void gridrawxygridxy(double *gvals,double dx,double dy,double *tx,double *ty) ;
void gridgridxyrawxy(double *gvals,double dx,double dy,double *tx,double *ty) ;
void gridcheck(double *gvals, int icheck, int *errwarn) ; 

