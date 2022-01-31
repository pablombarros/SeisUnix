
void binterpfind(int kigi, int *mgi, int mgi_tot, int mgiextr,
                 int kigc, int *mgc, int mgc_tot, int mgcextr,
                 int *mgix, int *mgcx, double *wi, double *wc) ;

void binterpapply(double *lwitims, double *hiitims, int mgi_tot, double wi,  
                  double *lwctims, double *hictims, int mgc_tot, double wc,
                  int lwinto, double *valsout) ;

int bhighi(int *all, int last, int iguy) ;
