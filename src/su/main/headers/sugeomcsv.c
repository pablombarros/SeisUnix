/* Copyright (c) Colorado School of Mines, 2021.*/
/* All rights reserved.                       */

/* SUGEOMCSV: $Revision: 1.2 $ ; $Date: 2023/04/11 16:00:01 $        */

#include <stdio.h>
#include <string.h>
#include "math.h"

/* #include <math.h> */

#include "su.h"
#include "segy.h"
#include "header.h"


/* Structure for point information */
struct  PointInfo {
     double *dfield;
     long long int *lfield;
};

segy tr;

struct PointInfo *RecInfo; /* For Storage area for all record values */
struct PointInfo  guy;     /* For storage of one record's values.    */
int num_to_sort_by;        /* needed within compSort                 */
int num_of_others;         /* needed within compOther                */

int GetCase(char* cbuf) ;
int countRec(FILE *ufile, char *textraw, int maxtext, 
             char *Rid, int lenid, int nicerecord) ;
void getCSV(char *textraw, char *textbeg, int maxtext, char rdel, 
            double *dfield, int *nspot, int numcases,   
            int ncount, int *comerr,int *morerr,int *numerr,int *nblank);
long long int longt (double dvalue, double dtolh, double dtol) ;
int compSort (const void * q1, const void * q2) ;
int compOther (const void * q1, const void * q2) ;
int bhigh(struct PointInfo * all, int last, struct PointInfo* guy) ;
double fromhead(segy tr, int k) ; 
void tohead(segy *tr, int k, double dvalue) ;
void tparse(char *tbuf, char d, char **fields, int *numfields) ; 


/*********************** self documentation **********************/
char *sdoc[] = {
"                          ",
" SUGEOMCSV - 3D/2D land/OBotm geometry from Comma Separated or fixed text.",
"                                                                          ",
" sugeomcsv <intraces >outtraces rfile=in.txt   (other parameters)         ",
"                                                                          ",
" Parameter overview:                                                      ",
"                                                                          ",
"       rfile= input text values file                                      ",
"       ufile= text file containing C_SU options (default is same rfile)   ", 
"       rtype= type of records in input text file (comma separated, fixed) ",
"       names= assign SU names to text values (with SPS2 and SPS1 options) ",
"       setid= accept text records based on first characters (X,S,R etc.)  ",
"       match= trace SU key values to be used to find correct text records ",
"      rdelim= rfile delimiter when rtype is csv (default is comma).       ",
"  nicerecord= allows skipping bad records at start of some text files     ",
"  maxrecords= maximum number of records to allocate memory for            ", 
"    unrepeat= help when duplicate fldr values exist (default is off)      ",
"      scalco= set scaler for sx,sy,gx,gy coordinate values (default 10)   ",
"      scalel= set scaler for elevation and related values (default 10)    ",
"     missing= options when trace match= values not found in text file     ", 
"      action= how to update trace header values (default is to replace)   ", 
"      create= create output traces with no input trace file               ", 
"      spikes= put spikes at specified times in the created output traces  ", 
"                                                                          ",
"  **********************************************************              ",
"   To output this documentation:  sugeomcsv 2> geomdoc.txt                ",
"  **********************************************************              ",
"                                                                          ",
" This program updates headers from general csv or fixed format text files.",
" (You cannot use UKOOA streamer-marine R-records as input here).          ",
"                                                                          ",
" This program can update full geometry using the 3 SPS fixed format files.",
" Parameter setting for SPS2 and SPS1 is so complex that this program has  ",
" explicit option names=sps2 and names=sps1 for those files. But the SPS   ",
" files often have numbers that are too large to be stored in some SU keys.",
" Or, the SPS files may be unusable for other reasons. In these cases, you ",
" need to use the SUTOOLCSV program to repair the SPS files or to convert  ",
" those files to comma-separated files for input to SpreadSheet programs   ",
" in order to make other repairs before running this program.              ",
"                                                                          ",
" Refer to the documentation of SUTOOLCSV for SPS files since it has       ",
" extensive explanations and examples.                                     ",
"                                                                          ",
"                                                                          ",
" For SPS2 files, 3 runs of this program are needed. Initially, I suggest  ",
" using create=all option for the X-file. This reads no input traces but   ",
" creates traces as defined by the X-file records (see create= below).     ",
"                                                                          ",
" sugeomcsv rfile=X.txt setid=x create=all >dummy1.su names=sps2 match=sps2",
"                                                                          ",
" sugeomcsv rfile=S.txt setid=s <dummy1.su >dummy2.su names=sps2 match=sps2",
"                                                                          ",
" sugeomcsv rfile=R.txt setid=r <dummy2.su >dummy3.su names=sps2 match=sps2",
"                                                                          ",
"                                                                          ",
" If the SPS2 files have been processed by SUTOOLCSV then these 3 runs     ",
" can be done instead.                                                     ",
"                                                                          ",
" sugeomcsv rfile=X.csv create=all >dummy1.su                              ",
" sugeomcsv rfile=S.csv <dummy1.su >dummy2.su                              ",
" sugeomcsv rfile=R.csv <dummy2.su >dummy3.su                              ",
"                                                                          ",
" ----------------------------------------------------------------------   ",
"                                                                          ",
" The documentation herein covers simple text files, but mostly concerns   ",
" itself with other aspects of the use of SUGEOMCSV.                       ",
"                                                                          ",
" In other words, the documentation for SUGEOMCSV is complex enough that   ",
" I do not want to repeat all the documentation contained in SUTOOLCSV.    ",
"                                                                          ",
"                                                                          ",
" -------------                                                            ",
"                                                                          ",
" Consider a comma separated text file with 3 values per record.           ",
" The 3 values are: record id, energy source number, and elevation.        ",
"                                                                          ",
"       S,11,343                                                           ",
"       S,42,342                                                           ",
"       S,25,340                                                           ",
"       C,some comment,                                                    ",
"       S,45,347                                                           ",
"                                                                          ",
"  To update shot elevation to trace headers requires this:                ",
"                                                                          ",
"  sugeomcsv <intraces >outtraces rfile=values3.csv setid=S match=es       ",
"            names=C_su_id,es,selev                                        ",
"                                                                          ",
"  setid=S says only read values from text records that start with S.      ",
"  match=es says get the value of es from input trace header.              ",
" names=C_su_id,es,selev assigns names to the 3 values in the_text_file    ",
" Other details for the above parameters will be discussed later.          ",
"                                                                          ",
" So, the value of es is gotten from input trace, then the text record     ",
" with same es value is found, and the selev is updated to output trace.   ",
" Note that every trace with the same es value is updated with same selev. ",
"   You can use multiple matches to find text records (match=sx,sy).       ",
" If the text records contain a fourth value with depth then you can       ",
" update it at the same time using (names=C_su_id,es,selev,sdepth).        ",
" If you want to update just sdepth then (names=C_su_id,es,null,sdepth)    ",
"                                                                          ",
" Note that trace input order does not matter. Neither does text record    ",
" order. The text file is read and all non-null and non-numb name values   ",
" are stored and sorted by match=. But there is no interpolation, so the   ",
" match= numbers must be found in text records (within a small tolerance). ",
"                                                                          ",
" Now consider the same 3 values but in a fixed format text file:          ",
"                                                                          ",
"       S       11  343                                                    ",
"       S       42  342                                                    ",
"       S       25  340                                                    ",
"       C  some comment                                                    ",
"       S       45  347                                                    ",
"                                                                          ",
"  To update shot elevation to trace headers requires this:                ",
"                                                                          ",
"  sugeomcsv <intraces >outtraces rfile=values3.txt setid=S match=es       ",
"            names=C_su_id,2_es_10,11_selev_15                             ",
"                                                                          ",
"  The first difference is that the rfile now has .txt rather than .csv    ",
"  The other difference is in names=C_su_id,2_es_10,11_selev_15 where      ",
"  both es and selev have leading and trailing numbers indicating the      ",
"  character range that they occupy in the fixed format records of rfile.  ",
"  Note that C_su_id does not have the leading and trailing numbers since  ",
"  it is special in several ways (detailed later).                         ",
"                                                                          ",
"  When parsing fixed format files it is not actually needed to specify    ", 
"  fields that you decide to null. But it is easier to make changes in     ", 
"  the future if you retain all known fields in the names= lists. This is  ", 
"  why SUTOOLCSV produces C_SU records that contain all officially defined ", 
"  fields for SPS2 and SPS1 formats (even though it nulls some names).     ", 
"                                                                          ",
"  --------------                                                          ",
"                                                                          ",
"  Alternately, instead of putting parameters on the command line,         ",
"  you can add similar parameter records into the rfile.                   ",
"  Details follow later, but here is an example:                           ",
"                                                                          ",
"  sugeomcsv <intraces >outtraces rfile=values3a.txt                       ",
"                                                                          ",
"       C_SU_MATCH,es                                                      ", 
"       C_SU_SETID,S                                                       ", 
"       C_SU_NAMES                                                         ", 
"       C_su_id,2_es_10,11_selev_15                                        ",
"       S       11  343                                                    ",
"       S       42  342                                                    ",
"       S       25  340                                                    ",
"       C  some comment                                                    ",
"       S       45  347                                                    ",
"                                                                          ",
"  Note that the default is to look for C_SU_ parameter records within     ",
"  the same rfile as the values records. And match= setid= names= are      ",
"  optional (but their C_SU_ records must be found if they are not on the  ",
"  command line). If they are specified on the command line and also exist ",
"  in the rfile, their command line options override rfile.                ",
"                                                                          ",
"  Note that commas are still specified in the C_SU_ parameter records     ",
"  even when the rest of the rfile is fixed format.                        ",
"                                                                          ",
"  --------------                                                          ",
"                                                                          ",
"  For comma separated rfiles, it looks like this:                         ",
"                                                                          ",
"  sugeomcsv <intraces >outtraces rfile=values3a.csv                       ",
"                                                                          ",
"       C_SU_MATCH,es,                                                     ", 
"       C_SU_SETID,S,                                                      ", 
"       C_SU_NAMES,,                                                       ", 
"       C_su_id,2_es_10,11_selev_15                                        ",
"       S,11,343                                                           ",
"       S,42,342                                                           ",
"       S,25,340                                                           ",
"       C,some comment,                                                    ",
"       S,45,347                                                           ",
"                                                                          ",
"  Notice that the C_SU_NAMES are allowed to have leading and trailing     ",
"  integers even though the data records have commas. When rtype=csv they  ",
"  are ignored in this program (but you should retain leading and trailing ",
"  numbers because SUTOOLCSV uses them to convert csv back to fixed files).",
"                                                                          ",
"  Notice also that C_SU_MATCH and C_SU_SETID have extra commas. This may  ",
"  occur for comma-separated files output by SpreadSheet programs. These   ",
"  extra commas are also permitted for C_SU_MATCH and C_SU_SETID records   ",
"  in fixed format files. (The C_SU_MATCH read-in stops at the first null, ",
"  all-blank, or sequential commas. And C_SU_SETID only reads first id).   ",
"                                                                          ",
"  If you are wondering why the actual names are in the row following      ",
"  C_SU_NAMES, consider how this will be output from a SpreadSheet program.",
"  (The 3 actual names need to be in the same columns with the data that   ",
"  they refer to).                                                         ",
"                                                                          ",
"  --------------                                                          ",
"                                                                          ",
"  For simple text files, using the command line to specify everything may ",
"  be more convenient. But specification for SPS and similar files gets    ",
"  quite intricate and is usually better done via C_SU_ parameter records. ",
"  For this reason, the SUTOOLCSV program actually outputs the high-level  ",
"  SPS2 and SPS1 options in expanded version into its output text files.   ",
"  This give you the opportunity to see and modify them as needed.         ",
"                                                                          ",
"  You can also specify ufile= to read the parameter records from another  ",
"  file. This allows you to use the same ufile for different surveys.      ",
"                                                                          ",
" ------------------------------------------------------------------------ ",
"                                                                          ",
" For full initial geometry this program usually must be run multiple times",
" For instance, if you have 3 SPS files, this program must be run 3 times. ",
" The SPS X-file must be used in first run, then S then R (or R then S).   ", 
"                                                                          ", 
" The first run of this program inputs the SPS X-file and then uses the    ", 
" trace header fldr and tracf numbers to find the correct X-record inside  ",
" the X-file. The computed gaps numbers is updated to the output trace     ",
" (gaps is the SU key name I decided to use, it is usually called          ",
" receiver point or receiver number in SPS and other documents).           ",
" Other values from the X-record are also updated (these are usually       ",
" receiver line number and shot line number and shot point number.         ",
" For 3D surveys, in the next 2 runs of SUGEOMCSV:                         ", 
"      Two of the values output from the X-file are used                   ", 
"      to find the correct records within the S-file.                      ", 
"      Two of the values output from the X-file are used                   ", 
"      to find the correct records within the R-file.                      ", 
" (For 2D surveys, in the next 2 runs of SUGEOMCSV only 1 value may be     ", 
"  needed to find the correct records for S-file and 1 for R-file).        ", 
"                                                                          ", 
" NOTE: The SPS discussions above assume that fldr contains field record   ", 
"       number, and tracf contains channel number. This is not necessarily ", 
"       the case. In particular, the channel number might not be in tracf. ", 
"       It might be in some other SU key name, or it might not exist       ", 
"       at all, in which case you are going to have to figure out how to   ", 
"       alter some input trace key name to correspond to the channel       ", 
"       ranges on the X-records.                                           ", 
"                                                                          ",
" NOTE: Sometimes SPS X-files deliberately specify receiver line numbers   ", 
"       or point ranges that do not exist in the R-files. For instance,    ", 
"       with a 10 line by 240 point pattern layout, it is sometimes easier ", 
"       for recording instruments to simply pretend that the line numbering", 
"       and point numbering continues infinitely. The traces corresponding ", 
"       to this are simply allowed to record noise. Othertimes the traces  ", 
"       are not recorded at all. For this situaton you need to specify     ", 
"       missing=delete or pass (or pre-delete the traces).                 ", 
"                                                                          ",
" ------------------------------------------------------------------------ ",
"                                                                          ",
" Required parameters:                                                     ",
"                                                                          ",
"       rfile=  text values file (fixed or comma separated values)         ", 
"                                                                          ", 
" Optional parameters:                                                     ", 
"                                                                          ", 
"       rtype= type of records in rfile text file. Default is csv if the   ", 
"              file name ends in csv, otherwise defaults to fixed.         ", 
"            =csv     comma separated values                               ", 
"            =fixed   This option is required for standard SPS files       ", 
"                     that have not been converted to csv by SUTOOLCSV     ", 
"                     - along with other specifications in names= list     ", 
"                     (see extensive examples below).                      ", 
"                                                                          ", 
"       ufile=  file to search for C_SU parameter records. The default is  ", 
"               to search in the rfile text file.                          ", 
"            =command   means ignore any C_SU records in rfile text file.  ", 
"            All required parameters must be specified on the command line.", 
"                                                                          ", 
" The following 3 parameters must be found on the command line or in       ", 
" their corresponding C_SU records in the rfile (or ufile) text file.      ", 
"                                                                          ", 
"       match= any number of input trace header keys needed to find exact  ",
"             record in the rfile text file. Must also be listed in names= ",
"             Example on command line:                                     ", 
"             match=fldr,tracf                                             ", 
"       match=SPS2 This is just a standard way to specify the match= list  ", 
"                  for SPS Revison 2 files (see examples in SUTOOLCSV).    ", 
"                  The setid= option must also be X,S, or R.               ", 
"       match=SPS1 This is just a standard way to specify the match= list  ", 
"                  for SPS Revison 1 files (see examples in SUTOOLCSV).    ", 
"                  The setid= option must also be X,S, or R.               ", 
"                                                                          ", 
"       names= is used to assign names to values in rfile text file (you   ",
"              are telling this program the names of values in text file). ",
"             For files with comma-separated values, a name must be listed ",
"             sequentially for each field in the rfile text file.          ", 
"             The names must also include the match= SU keys above.        ", 
"             Note C_su_id means this is field used for record acceptance. ", 
"        ***  Read the note   c_su_id IS SPECIAL   later. ***              ", 
"             Special names: null and numb with integer (such as null1)    ", 
"                           Means do not read/output this field. (You can  ", 
"                           also put nothing between sequential commas).   ", 
"               Example on command line:                                   ", 
"               names=C_su_id,cdp,null3,cx,cy,,ce                          ", 
"       names=SPS2 This is just a standard way to specify the names= list  ", 
"                  for SPS Revison 2 files (see examples in SUTOOLCSV).    ", 
"                  The setid= option must also be X,S, or R.               ", 
"       names=SPS1 This is just a standard way to specify the names= list  ", 
"                  for SPS Revison 1 files (see examples in SUTOOLCSV).    ", 
"                  The setid= option must also be X,S, or R.               ", 
"       names=SPS2ALL  Same as SPS2. The SUTOOLCSV program has this option ", 
"                      to output fields with no SU key (using numb names), ", 
"                      but numb names are treated like null names herein.  ", 
"       names=SPS1ALL  Same as SPS1.                                       ", 
"                                                                          ", 
"      setid= is used to accept text records based on their first field.   ", 
"        setid=S     means accept text records if their first field is S   ", 
"                          (any characters allowed, such as R,X,cdp,FRED)  ", 
"                    Note: this value is automatically upper-cased unless  ", 
"                          you surround it by double-quotes.               ", 
"                          So s becomes S unless you use double-quotes.    ", 
"        setid=ANY   means read all records (except those starting C_SU)   ", 
"                    and those records have an id field at front.          ", 
"        setid=NONE  means read all records (except those starting C_SU)   ", 
"                    but those records do not have an id field at front.   ", 
"                    (For csv files this means the field before the first  ", 
"                     comma is a value, not an identifier).                ", 
"             Example on command line:                                     ", 
"             setid=S                                                      ", 
"                                                                          ", 
" -----------------                                                        ", 
"                                                                          ", 
"       If match= is not specified on command line, this program searches  ",
"       for a text record starting with C_SU_MATCH and reads keys from it. ",
"             Example within the text file:                                ", 
"             C_SU_MATCH,fldr,tracf                                        ", 
"                                                                          ", 
"       If names= is not specified on command line, this program searches  ",
"       for a text record starting with C_SU_NAMES and reads names from    ",
"       the record after the C_SU_NAMES record.                            ",
"             Example within the text file:                                ", 
"             C_SU_NAMES                                                   ", 
"             C_su_id,cdp,null,cx,cy,null,ce                               ", 
"                                                                          ", 
"       If setid= is not specified on command line, this program searches  ",
"       for a text record starting with C_SU_SETID and reads id from it.   ",
"             Example within the text file:                                ", 
"             C_SU_SETID,S                                                 ", 
"                                                                          ", 
" Note these C_SU_ parameter records can be in any order within text file. ", 
"                                                                          ", 
" -----------------                                                        ", 
"                                                                          ", 
"      rdelim= If rtype=fixed you cannot specify this parameter.           ", 
"              If rtype=csv the default is comma. You can specify any      ",
"              single character here either by itself or surrounded by     ", 
"              double-quotes (some characters such as semi-colon may have  ", 
"              trouble getting through the command line).                  ", 
"     **Note** Specifying a blank here usually will not give good results  ", 
"              unless your input rfile has exactly 1 blank between fields. ", 
"                                                                          ", 
" -----------------                                                        ", 
"                                                                          ", 
"      scalco= multiply coordinates by this power of 10 (1,10,100...)      ",
"              before putting them in traces. Default is 10 which means    ",
"              that sx,sy,gx,gy from the text file are multiplied by 10.   ",
"              The actual value of scalco in the traces is therefore       ",
"              set to -10 (meaning these values need to be divided by 10). ",
"          *** You cannot specify scalco= if you are not updating any      ",
"          *** of the 4 coordinate values.                                 ",
"         **** The SU key offset value is recomputed if any of sx,sy,gx,gy ",
"         **** is in the list of names (and offset itself is NOT in list). ",
"            * If you are confident your text files contain coordinates    ",
"            * with only whole numbers, you can set this to 1. But you     ",
"            * cannot change this just for some of these 4 coordinates.    ",
"           ** This SEEMS like a problem about size, but it is actually a  ",
"           ** problem about decimal digits. It you use 1 here you will    ",
"           ** find that values like 544444.6 get rounded to 544445        ",
"           ** when SUGEOMCSV updates them to traces.                      ",
"          *** If scalco is in the list of names= then that value is set   ",
"              in output header. But only the value specified or defaulted ",
"              HERE is actually applied. So do not put it in names= unless ",
"              you are repairing some odd situation.                       ",
"                                                                          ", 
"      scalel= multiply elevation and other related values by this power   ",
"              of 10 (1,10,100...) before putting them in traces. Default  ",
"              is 10 which means gelev,selev,sdepth,gdel,sdel,swdep,gwdep  ",
"              from the text file are multiplied by 10.                    ",
"              The actual value of scalel in the traces is therefore       ",
"              set to -10 (meaning these values need to be divided by 10). ",
"          *** You cannot specify scalel= if you are not updating any      ",
"          *** of the 7 elevation related values.                          ",
"            * If you are confident your text files contain these values   ",
"            * with only whole numbers, you can set this to 1. But you     ",
"            * cannot change this just for some of these 7 values.         ",
"           ** This SEEMS like a problem about size, but it is actually a  ",
"           ** problem about decimal digits. It you use 1 here you will    ",
"           ** find that values like 3333.6 get rounded to 3334            ",
"           ** when SUGEOMCSV updates them to traces.                      ",
"          *** If scalel is in the list of names= then that value is set   ",
"              in output header. But only the value specified or defaulted ",
"              HERE is actually applied. So do not put it in names= unless ",
"              you are repairing some odd situation.                       ",
"                                                                          ", 
"     unrepeat= The default is not to enable this option.                  ", 
"               This option is general but most likely usefull for X-files ", 
"               where the field record number (fldr) increases but then    ", 
"               re-starts at a lower number. Such as 1->7800 then 5->4000. ", 
"               Normally, the finding-logic within this program would not  ", 
"               be able to distinguish the first fldr 5 from the second 5  ", 
"               and so on. (For that situation if you do not use this      ", 
"               option you will most likely get a channel-range error as   ",
"               that X-file is read-in and stored).                        ", 
"      unrepeat=1 Read the text file and generate an integer from 1 and    ", 
"                 increment by 1 every time the first match= reverses.     ", 
"                 Typically, the first match= is fldr for X-files so this  ", 
"                 increments +1 when fldr is increasing and then decreases,",
"                 and also increments +1 if fldr is decreasing and then    ",
"                 increases. The comparison is done using order of records ",
"                 as they exist in the text file (before sorting herein).  ",
"                 Another incrementing integer is generated the same way   ",
"                 except using the order of the input traces. These two    ",
"                 integers are used to match which (fldr) value in the     ",
"                 traces belongs to which (fldr) value from X-file.        ",
"             **  Note that this option is not fool-proof and should be    ",
"             **  used with caution and copious checking of results. This  ",
"             **  option is better than a pure sequential re-numbering     ",
"             **  because it still primarily uses the field record numbers,",
"             **  but it is not guaranteed.                                ",
"      unrepeat=n (any integer). The input text record incrementing integer", 
"                 will still start at 1 and increment +1. But the trace    ", 
"                 incrementing integer starts at this number. This allows  ", 
"                 you to input reduced sets of fldr records without need   ", 
"                 to edit the X-file.                                      ", 
"                                                                          ", 
" -----------------                                                        ", 
" -----------------                                                        ", 
"                                                                          ", 
"     nicerecord= record number to start trying to read data values from   ", 
"                 the rfile text file (default is 1). The beginning records", 
"                 of some text files are badly composed (comments or info.)", 
"                 When the setid= option is not able to reject them,       ", 
"                 specify a record number here where setid= will work.     ", 
"                 (This program also always knows that records starting    ", 
"                 with C_SU are not data records, and will not try to      ", 
"                 read data values from them even when setid=ALL). But     ", 
"                 it will read C_SU parameter records even if they are     ", 
"                 previous to this nicerecord number.                      ", 
"                                                                          ", 
"     maxrecords= maximum number of records to allocate memory for.        ", 
"                 If not specified, this program reads through the records ", 
"                 once and allocates for the number found. Then reads them ", 
"                 again. This double reading takes more time. If you want  ", 
"                 to avoid this, specify a maximum via this parameter.     ", 
"                                                                          ", 
"    missing= specifies what to do if an input trace match= values cannot  ", 
"             be found in the rfile text file. Default is to error-halt    ", 
"             at the first trace where this occurs.                        ", 
"     missing=delete   means delete all traces for which this occurs.      ", 
"     missing=pass     means pass those traces to output without update.   ", 
"                      This allows you to, for instance, change statics    ", 
"                      values for just some shots or receivers without     ", 
"                      having every shot or receiver in the rfile text file", 
"                                                                          ", 
"    action= specifies how to update the trace header values               ", 
"            Default is to set the headers to the text record values.      ", 
"            The option here is not applied to the values of               ", 
"            names on the match= list.                                     ", 
"     action=add       Add the text values to the trace header values.     ", 
"     action=subtract  Subtract text values from the trace header values.  ", 
"                                                                          ", 
" -----------------                                                        ", 
" -----------------                                                        ", 
"                                                                          ", 
"    create=    create output traces and update from the rfile text file.  ", 
"          ***  Attempting to also input a trace file will error-halt      ", 
"          ***  unless spikes=n).                                          ", 
"               Values of all non-null names on names= list are updated to ", 
"               traces (including those on match= list). The behavior of   ",
"               this program is the same as when traces are input. This    ", 
"               includes sorting of the text records into match= order and ", 
"               therefore created traces are output in that order. Most    ", 
"               text files result in one trace per record, but SPS X-files ",
"               (and similar) usually create many traces for each record.  ",
"  create=all   create traces for all data records in the rfile text file  ", 
"  create=n     create n or all (whichever is fewer)                       ", 
"                                                                          ", 
"  spikes= ***  This option will error-halt if create= is not specified.   ", 
"               By default, create= outputs only 1 zero-sample at 4 ms.    ", 
"               This parameter accepts a list of time,amplitude pairs.     ", 
"                                                                          ", 
"                 Example: spikes=4,0.001,800,1000,1500,-2000,3000,0       ", 
"                                                                          ", 
"               The maximum time sets the trace length (3000 ms. above).   ", 
"               The first pair in the list is special, it specifies the    ", 
"               sample interval (in ms.) and base-amplitude.               ", 
"               In this example, the sample interval is 4, and the         ", 
"               base-amplitude is 0.001. The base-amplitude is the value   ", 
"               samples are set to if they are not spiked by remainder of  ", 
"               the list. (For geometry checking and other QC and survey   ", 
"               design reasons, base-amplitude of 0.0 is not a nice value. ", 
"               For instance, 0.0 makes it difficult to check mutes).      ", 
"               All times should be whole multiples of the sample interval ", 
"               but will be rounded to nearest sample time if not.         ", 
"               Times do not have to increase but cannot be negative.      ", 
"     spikes=n  Where n is the sequential number of an input trace. The    ", 
"               sample values of this trace will be copied to all created  ", 
"               output traces. Values of all non-null names on names= list ", 
"               are updated to traces (including those on the match= list).",
"               The intention here is to allow you to input a reflectivity ", 
"               or well-log type of trace, duplicate the samples and change", 
"               some header values via the text file (for testing purposes)", 
"                                                                          ", 
" -----------------                                                        ", 
" -----------------                                                        ", 
"                                                                          ", 
" ADVICE: Run create=all for your first relation file (SPS X-file) setup.  ", 
"         Then check the values in the output trace headers.               ", 
"         Then run SPS S-file setup to update those headers (not create).  ",
"         Then check the values in the output trace headers.               ", 
"         Then run SPS R-file setup to update those headers (not create).  ",
"         Then check the values in the output trace headers.               ", 
"         This will allow you to check the mechanics of your setups.       ", 
"                                                                          ", 
" -----------------                                                        ", 
"                                                                          ", 
"  MORE ADVICE: Use spike options to test your setups of other programs.   ", 
"               For instance, after running the sequence advised above,    ", 
"               put positive and negative spikes on the traces.            ", 
"               Then run a simple bandpass filter. At this point the       ", 
"               results will begin to look like post-NMO shot gathers.     ", 
"               Apply inverse-NMO and offset-Mute to make your traces look ", 
"               like pre-NMO shot gathers.                                 ", 
"               Then you can run some interesting tests, like:             ", 
"           1.  Apply surface-consistent residual statics to some shots and", 
"               receivers. Surface-consistant statics analysis programs    ", 
"               should easily be able to find the statics you applied.     ", 
"           2.  Put surface-consistant gains on some shots and receivers.  ", 
"               Surface-consistant amplitude analysis programs should      ", 
"               easily find the gains you applied.                         ", 
"           3.  Put surface-consistent signatures on some shots and        ", 
"               receivers. Surface-consistent deconvolution programs       ",
"               should easily find the signatures you applied.             ", 
"                                                                          ", 
"   CAUTION: Just because you apply something and another program computes ", 
"            the correct inverse, it does not mean the program is setup    ", 
"            well for actual seismic situations. These kinds of tests are  ", 
"            too simple to be considered forward-modelling.                ", 
"                                                                          ", 
"   OF COURSE: You might do this kind of thing without having any actual   ", 
"              seismic data, just to test a few survey-design ideas.       ", 
"                                                                          ", 
"                                                                          ", 
" ------------------------------------------------------------------------ ",
" ------------------------------------------------------------------------ ",
"                                                                          ", 
" C_su_id IS SPECIAL ********************                                  ", 
"                                                                          ", 
"  (a) Its character range is taken from the length of the value specified ", 
"      for setid (for instance S has range 1 to 1, FRED has range 1 to 4). ", 
"  (b) The value specified for id is case-sensitive (r is not R). The id   ", 
"      value is the only thing case-sensitive in this program except for   ", 
"      the file names. So this program does not care if you use parameter  ", 
"      records starting with C_SU or c_su, but you want other programs to  ", 
"      ignore C_SU records, so use a capital C, not a lower case c.        ", 
"                                                                          ", 
" -----------------------------------------------------------------        ",
"                                                                          ",
NULL};

/* Credits:                                                       */
/* Andre Latour                                                   */ 
/*                                                                */
/* Started from sugeom by Fernando Roxo <petro@roxo.org>          */
/* This program is now almost completely different from SUGEOM,   */
/* but it did save me much effort getting familiar with SU.       */
/*                                                                */
/* Trace header fields accessed: Potentially any, or all.         */
/*                                                                */
/**************** end self doc ************************************/

int main(int argc, char **argv) {

   cwp_String Rname=NULL;  /* text file name for values            */
   cwp_String Uname=NULL;  /* text file name for C_SU_ parameters  */
   cwp_String Rid  =NULL;  /* rejection id for records             */
   FILE *fpR=NULL;         /* file pointer for Rname file          */
   FILE *fpU=NULL;         /* file pointer for Uname file          */

   int idxR;     

   int nproct = 0; 

   int i=0;
   int j=0;
   int k=0;
   int n=0;
   int m=0;

   int ineed=0;

/* Most of following are about 10 times bigger than will be used.  */    
/* But difficult to dynamically allocate since they often will be  */    
/* read-in from C_SU_ records in the input text file.              */

   cwp_String match[99];
   cwp_String names[999];   
   cwp_String namex[999];   
   int maxtext = 10001;
   char  textraw[10001]; /* fgets puts a \0 after contents */
   char  textbeg[10001]; /* so this wastes memory but not time */
   char  textfront[10];    

   int *ilead = NULL;
   int *itrail = NULL;
   int *ncase = NULL;
   int *nspot = NULL;
   double *valmx = NULL;
   double *dvals = NULL;
   int *kcase = NULL;
   int *ktol = NULL;
   int *klocn = NULL;
   int mapx[10];
        
   /* Initialize */
   initargs(argc, argv);
   requestdoc(1);

   int nicerecord; 
   if (!getparint("nicerecord", &nicerecord)) nicerecord = 1;
   if (nicerecord<1) err("**** Error: nicerecord= cannot be less than 1");

   int numR;          
   if (!getparint("maxrecords", &numR)) numR = 0;
   if (numR<0) err("**** Error: maxrecords= cannot be less than 0");

   float ftol; /* note: deliberately undocumented. Users should not set it. */ 
   if (!getparfloat("tolr", &ftol)) ftol = 0.01;
   if (ftol<0.000000001) err("**** Error: tolr= must be larger."); /* see longt */
   double dtol = ftol; 
   double dtolh = dtol / 2.;
   dtol = 1.0/dtol;

   int unrepeat; /* conceivable users will want to start traces at -2 or -3 */
   if (!getparint("unrepeat", &unrepeat)) unrepeat = -2147483645;

   int iscalel;
   int iscaleldef = 0;
   if (!getparint("scalel", &iscalel)) {
     iscalel = 10;
     iscaleldef = 1;
   }
   if(iscalel==-1) iscalel = 1;
   int nscalel = abs(iscalel);
   if(nscalel!=1 && nscalel!=10 && nscalel!=100 && nscalel!=1000 && nscalel!=10000 
     && nscalel!=100000 && nscalel!=1000000 && nscalel!=10000000) {
     err("**** Error: scalel= must be signed powers of 10 (1,10,100...-10,-100,...)");
   }
   double dscalel;
   if(iscalel>0) dscalel = (double)(iscalel);
   else          dscalel = -1./(double)(iscalel);
   iscalel *= -1; /* to output in header */

   int iscalco;
   int iscalcodef = 0;
   if (!getparint("scalco", &iscalco)) {
     iscalco = 10;
     iscalcodef = 1;
   }
   if(iscalco==-1) iscalco = 1;
   int nscalco = abs(iscalco);
   if(nscalco!=1 && nscalco!=10 && nscalco!=100 && nscalco!=1000 && nscalco!=10000 
     && nscalco!=100000 && nscalco!=1000000 && nscalco!=10000000) {
     err("**** Error: scalco= must be signed powers of 10 (1,10,100...-10,-100,...)");
   }
   double dscalco;
   if(iscalco>0) dscalco = (double)(iscalco);
   else          dscalco = -1./(double)(iscalco);
   iscalco *= -1; /* to output in header */

   int missing = -1; /* default error-halt at first trace with no match= in text file */
   char *tmiss;
   getparstring("missing", &tmiss);
   if(tmiss != NULL) {
     for (n=0; n<strlen(tmiss); n++) tmiss[n] = tolower(tmiss[n]);
     if(strlen(tmiss)==6 && strcmp(tmiss,"delete") == 0) missing = 0;
     else if(strlen(tmiss)==4 && strcmp(tmiss,"pass") == 0) missing = 1;
     else err("**** Error: missing= option not recognized.");
   }

   int iaction = 0; /* default is to set trace header values to text record values */
   char *tact;
   getparstring("action", &tact);
   if(tact != NULL) {
     for (n=0; n<strlen(tact); n++) tact[n] = tolower(tact[n]);
     if(strlen(tact)==3 && strcmp(tact,"add") == 0) iaction = 1;
     else if(strlen(tact)==8 && strcmp(tact,"subtract") == 0) iaction = -1;
     else err("**** Error: action= option not recognized.");
   }

   int ncreate = 0;
   char *tall;
   getparstring("create", &tall); /* create traces and populate their headers with rfile values? */
   if (tall != NULL) {
     if(strlen(tall)==3) {
       tall[0] = tolower(tall[0]);
       tall[1] = tolower(tall[1]);
       tall[2] = tolower(tall[2]);
       if(strcmp(tall,"all") == 0) ncreate = 2147483645;
     }
     if(ncreate==0) getparint("create", &ncreate);
     if(ncreate<0) err("**** Error create= cannot be negative.");
   }

   float *spikes;
   spikes = calloc(1,sizeof(float)); /* eliminates annoying compilor warning when unused for some option */ 
   int ireplicate = 0;
   int nspikes = countparval("spikes");
   if(nspikes>0) {
     if(ncreate<1) err("**** Error: spikes= can only be specified on command line if create= is also there.");
     if(nspikes>1) {
       if(nspikes%2 == 1) err("**** Error: spikes= must have an even number of values (time amplitude pairs).");
       spikes = calloc(nspikes,sizeof(float)); /* there are actually half this number of spikes */
       getparfloat("spikes",spikes);
     }
     else {
       if (!getparint("spikes", &ireplicate)) ireplicate = 0;
       if(ireplicate<1) err("**** Error: For spikes=n the n must be +integer (the trace to input).");
     }
   }

   if(ncreate>0) { /* creating traces */ 
     if(isatty(STDIN_FILENO)!=1) { /* have input trace file */
       if(ireplicate<1) err("**** Error: Input trace file not allowed when create= (unless spikes=n).");
     }
     else {
       if(ireplicate>0) err("**** Error: Input trace file must be specified when spikes=n.");
     }
     if (unrepeat > -2147483645 && unrepeat != 1) {
       err("**** Error: If specified, only unrepeat=1 allowed when creating traces.");
     }
   }
   else { /* not creating traces */
     if(isatty(STDIN_FILENO)==1) err("**** Error: Input trace file must be specified when not create= ");
   }

   getparstring("rfile", &Rname);
   if ( Rname == NULL) err("**** Error: rfile= text file name must be specified.");

   char rdel = ',';
   char *tdel;
   getparstring("rdelim", &tdel);
   if (tdel != NULL) {
     if(strlen(tdel)==1) rdel = tdel[0];
     else if(strlen(tdel)==3 && tdel[0]=='"' && tdel[2]=='"') rdel = tdel[1];
     else err("**** Error: rdelim= specification not recognized.");
   }

   int irtype = 1;
   char *rtype;
   if(countparval("rtype") > 0) {
     getparstring("rtype", &rtype);
     for (n=0; n<strlen(rtype); n++) rtype[n] = tolower(rtype[n]);
     if(strlen(rtype)==3 && strcmp(rtype,"csv") == 0) irtype = 1;
     else if(strlen(rtype)==5 && strcmp(rtype,"fixed") == 0) irtype = 0;
     else err("**** Error: rtype= option not recognized.");
   }
   else {
     if(strlen(Rname)>2) {
       if(strcmp(Rname+strlen(Rname)-3,"csv") == 0) irtype = 1;
       else irtype = 0;
     }
   }

   if(irtype==1) {
     if (tdel == NULL) rdel = ',';
   }
   else if(irtype==0) {
     if (tdel != NULL) err("**** Error: you cannot specify rdelim= if rtype=fixed.");
   }

   fpR = fopen(Rname, "r");
   if (fpR == NULL) err("**** Error opening the rfile text file.");

   int iread_c_su = 1;
   fpU = fpR;
   getparstring("ufile", &Uname);
   if (Uname == NULL) { /* default is to (try) to read C_SU parameter records from rfile */
     iread_c_su = 1;
     fpU = fpR;
   }
   else {
     if(strcmp(Uname,"command") == 0) { /* required parameters must be on command line */
       iread_c_su = 0;
       fpU = 0;
     }
     else if(strcmp(Uname,Rname) == 0) { /* user specified the same text file as rfile */
       iread_c_su = 1;
       fpU = fpR;
     }
     else {           /* (try) to read C_SU parameter records from different text file */
       iread_c_su = 2;
       fpU = fopen(Uname, "r");
       if (fpU == NULL) err("**** Error opening the ufile");
     }
   }
   
   int nerr= 0;
   num_to_sort_by = countparval("match");
   if(num_to_sort_by>0) {
     getparstringarray("match",match);
   }
   else {
     if(iread_c_su==0) {
       nerr = 1;
       warn("**** Must have match= on command line when ufile=command.");
     }
   }

   int lenid = 0;  /* when id is set, this is length (S,R,X have length 1) */
   if(countparval("setid") > 0) {
     getparstring("setid", &Rid);
     lenid = strlen(Rid);
   }

   int num_names = countparval("names");
   if(num_names>0) {
     getparstringarray("names",names);
     for (n=0; n<num_names; n++) {
       for (m=0; m<strlen(names[n]); m++) {
         names[n][m] = tolower(names[n][m]);
       }
     }
   }
   else {
     if(iread_c_su==0) {
       nerr = 1;
       warn("**** Must have names= on command line when ufile=command.");
     }
   }
  
   if(nerr>0) {
     if(iread_c_su == 0) {
       warn("**** If trying to read C_SU records in text file, cannot have ufile=command.");
     }
     err("**** A required parameter has not been specified. See previous messages.");
   }

/* Cycle over ufile records to get some C_SU_ parameters? */ 

   int names_more = 0;

   if(num_names==0 || num_to_sort_by<1 || lenid==0) { /* no check iread_c_su==0, already err */ 

     int in_num_names      = num_names;
     int in_num_to_sort_by = num_to_sort_by;
     int in_lenid          = lenid;
     int num_c_su_names = 0;
     int num_c_su_match = 0;
     int num_c_su_setid = 0;
     int read_names = 0;
     int some_names = 0;
     int nrow  = 0;

     while (fgets(textraw, maxtext, fpU) != NULL) { /* read a line */

       nrow++;

/* Stop this looping? Sometimes, it really will just loop through to last record. */
                   
       if((in_num_names>0 || read_names==-2) && num_to_sort_by>0 && lenid>0) break;

/* Remove all blanks and tabs because tparse is not designed to handle them.      */

       int tsize = 0;
       for(n=0; n<maxtext; n++) { /*   linux \n            windows \r  */ 
         if(textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') break;
         if(textraw[n] != ' ' && textraw[n] != '\t') {
           textbeg[tsize] = textraw[n];
           tsize++;
         }
       }
       
       if(lenid==0) { /* note the id itself is case-sensitive. */
         for(n=0; n<10; n++) textbeg[n] = tolower(textbeg[n]);
         if(strncmp(textbeg,"c_su_setid",10) == 0) {
           char *ids[99];
           int nids;
           textbeg[tsize] = '\0';
           tparse(textbeg, ',', ids, &nids) ; 
           Rid = ids[1];
           lenid = strlen(Rid);
         }   
       }   

       for(n=0; n<sizeof(textbeg); n++) textbeg[n] = tolower(textbeg[n]);
       if(strncmp(textbeg,"c_su_match",10) == 0) num_c_su_match++;
       if(strncmp(textbeg,"c_su_setid",10) == 0) num_c_su_setid++;
       if(strncmp(textbeg,"c_su_names",10) == 0) num_c_su_names++;

       if(num_to_sort_by<1) { 
         if(strncmp(textbeg,"c_su_match",9) == 0) {
           int num_found;
           textbeg[tsize] = '\0';
           tparse(textbeg, ',', match, &num_found) ; 
           num_to_sort_by = num_found;
           for (j=1; j<num_found; j++) { 
             if(strcmp(match[j],"null") == 0) {
               num_to_sort_by = j;
               break;
             }
             match[j-1] = match[j]; /* get rid of c_su_match returned in first element  */
           }
           num_to_sort_by--;
         }
       } 

       if(read_names==-1) {
         if(strncmp(textbeg,"c_su_more",9) == 0) {
           read_names = 1;
           names_more = 1;
         }
         else read_names = -2;
       }

       if(read_names>0) {
         textbeg[tsize] = '\0'; /* andre */
         tparse(textbeg, ',', names+num_names, &some_names) ; 
         num_names += some_names;
         read_names = -1;
       }

       if(in_num_names==0 && strncmp(textbeg,"c_su_names",10) == 0) read_names = 1;

     } /* end of while (fgets(textraw,..... */
  
     fseek(fpU, 0L, SEEK_SET); /* reposition file to beginning record */ 

     int ierr = 0;
     if(in_num_names<1) {
       if(num_c_su_names>1) {
         ierr = 1;
         warn("**** Error: Text file has more than one C_SU_NAMES parameter record.");
         warn("****        Remove duplicates or override with names= on command line.");
       }
       else if(num_c_su_names==0) {
         ierr = 1;
         warn("**** Error: Text file has no C_SU_NAMES parameter record.");
         warn("****        Add it, or specify names= on command line.");
       }
     }
   
     if(in_num_to_sort_by<1) {
       if(num_c_su_match>1) {
         ierr = 1;
         warn("**** Error: Text file has more than one C_SU_MATCH record.");
         warn("****        Remove duplicates or override with match= on command line.");
       }
       else if(num_c_su_match==0) {
         ierr = 1;
         warn("**** Error: Text file has no C_SU_MATCH record.");
         warn("****        Add it, or specify match= command line.");
       }
     }
   
     if(in_lenid==0) {
       if(num_c_su_setid>1) {
         ierr = 1;
         warn("**** Error: Text file has more than one C_SU_SETID record.");
         warn("****        Remove duplicates or override with setid= on command line.");
       }
       else if(num_c_su_setid==0) {
         ierr = 1;
         warn("**** Error: Text file has no C_SU_SETID record.");
         warn("****        Add it, or specify setid= on command line.");
       }
     }
  
     if(ierr>0) {
       err("**** Error: Text file has duplicate or missing C_SU_NAMES or C_SU_MATCH or C_SU_SETID records.");
     }

   } /* end of   if(num_names==0 || num_to_sort_by<1 || lenid==0) { */ 

/* Resolve setid options.  */

   if(lenid>3 && Rid[0]=='"' && Rid[lenid-1]=='"') {
     for(n=0; n<lenid-1; n++) Rid[n] = Rid[n+1];
     lenid -= 2;
   }
   else {
     for(n=0; n<lenid; n++) Rid[n] = toupper(Rid[n]);
   }

   if(lenid==3) {
     char text3[3]; 
     for(n=0; n<3; n++) text3[n] = tolower(Rid[n]);
     if(strncmp(text3,"any",3) == 0) lenid = 0;        
   }
   else if(lenid==4) {
     char text4[4]; 
     for(n=0; n<4; n++) text4[n] = tolower(Rid[n]);
     if(strncmp(text4,"none",4) == 0) {
       lenid = 0;         
     }
   }

/* --------------------------------------------- */
/* --------------------------------------------- */
/* --------------------------------------------- */

   if(strcmp(match[0],"sps2") == 0 || strcmp(match[0],"sps1") == 0) {
     if(strcmp(Rid,"X") == 0) {
       match[0] = ealloc1(4,1);
       strcpy(match[0],"fldr");

       match[1] = ealloc1(5,1);
       strcpy(match[1],"tracf");
   
       num_to_sort_by = 2;
     }
     else if(strcmp(Rid,"S") == 0) {
       match[0] = ealloc1(6,1);
       strcpy(match[0],"grnofr");

       match[1] = ealloc1(6,1);
       strcpy(match[1],"grnlof");
   
       num_to_sort_by = 2;
     }
     else if(strcmp(Rid,"R") == 0) {
       match[0] = ealloc1(6,1);
       strcpy(match[0],"grnors");

       match[1] = ealloc1(4,1);
       strcpy(match[1],"gaps");
   
       num_to_sort_by = 2;
     }
   }

/* Note sps2 and sps2all are same in this program since both null and numb are ignored. */

   if(strcmp(names[0],"sps2") == 0 || strcmp(names[0],"sps2all") == 0) {
     if(strcmp(Rid,"X") == 0) {
       names[0] = ealloc1(7,1);
       strcpy(names[0],"c_su_id");

       names[1] = ealloc1(9,1);
       strcpy(names[1],"2_null2_7");

       names[2] = ealloc1(11,1);
       strcpy(names[2],"8_match1_15");

       names[3] = ealloc1(11,1);
       strcpy(names[3],"16_null4_16");

       names[4] = ealloc1(11,1);
       strcpy(names[4],"17_null5_17");

       names[5] = ealloc1(12,1);
       strcpy(names[5],"18_grnofr_27");

       names[6] = ealloc1(12,1);
       strcpy(names[6],"28_grnlof_37");

       names[7] = ealloc1(11,1);
       strcpy(names[7],"38_null8_38");

       names[8] = ealloc1(16,1);
       strcpy(names[8],"39_matche1_cf_43");

       names[9] = ealloc1(16,1);
       strcpy(names[9],"44_matche1_ct_48");

       names[10] = ealloc1(16,1);
       strcpy(names[10],"49_matche1_ci_49");

       names[11] = ealloc1(12,1);
       strcpy(names[11],"50_grnors_59");

       names[12] = ealloc1(13,1);
       strcpy(names[12],"60_gaps_rf_69");

       names[13] = ealloc1(13,1);
       strcpy(names[13],"70_gaps_rt_79");

       names[14] = ealloc1(12,1);
       strcpy(names[14],"80_null15_80");

       num_names = 15;

     }
     else if(strcmp(Rid,"S") == 0) {
       names[0] = ealloc1(7,1);
       strcpy(names[0],"c_su_id");

       names[1] = ealloc1(11,1);
       strcpy(names[1],"2_grnofr_11");

       names[2] = ealloc1(12,1);
       strcpy(names[2],"12_grnlof_21");

       names[3] = ealloc1(11,1);
       strcpy(names[3],"22_null4_23");

       names[4] = ealloc1(11,1);
       strcpy(names[4],"24_null5_24");

       names[5] = ealloc1(11,1);
       strcpy(names[5],"25_null6_26");

       names[6] = ealloc1(11,1);
       strcpy(names[6],"27_sstat_30");

       names[7] = ealloc1(12,1);
       strcpy(names[7],"31_sdepth_34");

       names[8] = ealloc1(10,1);
       strcpy(names[8],"35_sdel_38");

       names[9] = ealloc1(9,1);
       strcpy(names[9],"39_sut_40");

       names[10] = ealloc1(11,1);
       strcpy(names[10],"41_swdep_46");

       names[11] = ealloc1(8,1);
       strcpy(names[11],"47_sx_55");

       names[12] = ealloc1(8,1);
       strcpy(names[12],"56_sy_65");

       names[13] = ealloc1(11,1);
       strcpy(names[13],"66_selev_71");

       names[14] = ealloc1(12,1);
       strcpy(names[14],"72_null15_74");

       names[15] = ealloc1(12,1);
       strcpy(names[15],"75_null16_76");

       names[16] = ealloc1(12,1);
       strcpy(names[16],"77_null17_78");

       names[17] = ealloc1(12,1);
       strcpy(names[17],"79_null18_80");

       num_names = 18;

     }
     else if(strcmp(Rid,"R") == 0) {
       names[0] = ealloc1(7,1);
       strcpy(names[0],"c_su_id");

       names[1] = ealloc1(11,1);
       strcpy(names[1],"2_grnors_11");

       names[2] = ealloc1(10,1);
       strcpy(names[2],"12_gaps_21");

       names[3] = ealloc1(11,1);
       strcpy(names[3],"22_null4_23");

       names[4] = ealloc1(11,1);
       strcpy(names[4],"24_null5_24");

       names[5] = ealloc1(11,1);
       strcpy(names[5],"25_null6_26");

       names[6] = ealloc1(11,1);
       strcpy(names[6],"27_gstat_30");

       names[7] = ealloc1(11,1);
       strcpy(names[7],"31_null8_34");

       names[8] = ealloc1(10,1);
       strcpy(names[8],"35_gdel_38");

       names[9] = ealloc1(9,1);
       strcpy(names[9],"39_gut_40");

       names[10] = ealloc1(11,1);
       strcpy(names[10],"41_gwdep_46");

       names[11] = ealloc1(8,1);
       strcpy(names[11],"47_gx_55");

       names[12] = ealloc1(8,1);
       strcpy(names[12],"56_gy_65");

       names[13] = ealloc1(11,1);
       strcpy(names[13],"66_gelev_71");

       names[14] = ealloc1(12,1);
       strcpy(names[14],"72_null15_74");

       names[15] = ealloc1(12,1);
       strcpy(names[15],"75_null16_76");

       names[16] = ealloc1(12,1);
       strcpy(names[16],"77_null17_78");

       names[17] = ealloc1(12,1);
       strcpy(names[17],"79_null18_80");

       num_names = 18;

     }
   } /* end of  if(strcmp(names[0],"sps2") == 0 || strcmp(names[0],"sps2all") == 0) { */

/* Note sps1 and sps1all are same in this program since both null and numb are ignored. */

   if(strcmp(names[0],"sps1") == 0 || strcmp(names[0],"sps1all") == 0) {
     if(strcmp(Rid,"X") == 0) {
       names[0] = ealloc1(7,1);
       strcpy(names[0],"c_su_id");

       names[1] = ealloc1(9,1);
       strcpy(names[1],"2_null2_7");

       names[2] = ealloc1(11,1);
       strcpy(names[2],"8_match1_11");

       names[3] = ealloc1(11,1);
       strcpy(names[3],"12_null4_12");

       names[4] = ealloc1(11,1);
       strcpy(names[4],"13_null5_13");

       names[5] = ealloc1(12,1);
       strcpy(names[5],"14_grnofr_29");

       names[6] = ealloc1(12,1);
       strcpy(names[6],"30_grnlof_37");

       names[7] = ealloc1(11,1);
       strcpy(names[7],"38_null8_38");

       names[8] = ealloc1(16,1);
       strcpy(names[8],"39_matche1_cf_42");

       names[9] = ealloc1(16,1);
       strcpy(names[9],"43_matche1_ct_46");

       names[10] = ealloc1(16,1);
       strcpy(names[10],"47_matche1_ci_47");

       names[11] = ealloc1(12,1);
       strcpy(names[11],"48_grnors_63");

       names[12] = ealloc1(13,1);
       strcpy(names[12],"64_gaps_rf_71");

       names[13] = ealloc1(13,1);
       strcpy(names[13],"72_gaps_rt_79");

       names[14] = ealloc1(12,1);
       strcpy(names[14],"80_null15_80");

       num_names = 15;

     }
     else if(strcmp(Rid,"S") == 0) {
       names[0] = ealloc1(7,1);
       strcpy(names[0],"c_su_id");

       names[1] = ealloc1(11,1);
       strcpy(names[1],"2_grnofr_17");

       names[2] = ealloc1(12,1);
       strcpy(names[2],"18_grnlof_25");

       names[3] = ealloc1(11,1);
       strcpy(names[3],"26_null4_26");

       names[4] = ealloc1(11,1);
       strcpy(names[4],"26_null5_26");

       names[5] = ealloc1(11,1);
       strcpy(names[5],"27_null6_28");

       names[6] = ealloc1(11,1);
       strcpy(names[6],"29_sstat_32");

       names[7] = ealloc1(12,1);
       strcpy(names[7],"33_sdepth_36");

       names[8] = ealloc1(10,1);
       strcpy(names[8],"37_sdel_40");

       names[9] = ealloc1(9,1);
       strcpy(names[9],"41_sut_42");

       names[10] = ealloc1(11,1);
       strcpy(names[10],"43_swdep_46");

       names[11] = ealloc1(8,1);
       strcpy(names[11],"47_sx_55");

       names[12] = ealloc1(8,1);
       strcpy(names[12],"56_sy_65");

       names[13] = ealloc1(11,1);
       strcpy(names[13],"66_selev_71");

       names[14] = ealloc1(12,1);
       strcpy(names[14],"72_null15_74");

       names[15] = ealloc1(12,1);
       strcpy(names[15],"75_null16_76");

       names[16] = ealloc1(12,1);
       strcpy(names[16],"77_null17_78");

       names[17] = ealloc1(12,1);
       strcpy(names[17],"79_null18_80");

       num_names = 18;

     }
     else if(strcmp(Rid,"R") == 0) {
       names[0] = ealloc1(7,1);
       strcpy(names[0],"c_su_id");

       names[1] = ealloc1(11,1);
       strcpy(names[1],"2_grnors_17");

       names[2] = ealloc1(10,1);
       strcpy(names[2],"18_gaps_25");

       names[3] = ealloc1(11,1);
       strcpy(names[3],"26_null4_26");

       names[4] = ealloc1(11,1);
       strcpy(names[4],"26_null5_26");

       names[5] = ealloc1(11,1);
       strcpy(names[5],"27_null6_28");

       names[6] = ealloc1(11,1);
       strcpy(names[6],"29_gstat_32");

       names[7] = ealloc1(11,1);
       strcpy(names[7],"33_null8_36");

       names[8] = ealloc1(10,1);
       strcpy(names[8],"37_gdel_40");

       names[9] = ealloc1(9,1);
       strcpy(names[9],"41_gut_42");

       names[10] = ealloc1(11,1);
       strcpy(names[10],"43_gwdep_46");

       names[11] = ealloc1(8,1);
       strcpy(names[11],"47_gx_55");

       names[12] = ealloc1(8,1);
       strcpy(names[12],"56_gy_65");

       names[13] = ealloc1(11,1);
       strcpy(names[13],"66_gelev_71");

       names[14] = ealloc1(12,1);
       strcpy(names[14],"72_null15_74");

       names[15] = ealloc1(12,1);
       strcpy(names[15],"75_null16_76");

       names[16] = ealloc1(12,1);
       strcpy(names[16],"77_null17_78");

       names[17] = ealloc1(12,1);
       strcpy(names[17],"79_null18_80");

       num_names = 18;

     }
   } /* end of  if(strcmp(names[0],"sps1") == 0 || strcmp(names[0],"sps1all") == 0) { */

   ilead = calloc(num_names,sizeof(int));
   if(ilead == NULL) err("**** Unable to allocate memory.");
   itrail = calloc(num_names,sizeof(int));
   if(itrail == NULL) err("**** Unable to allocate memory.");
   ncase = calloc(num_names,sizeof(int));
   if(ncase == NULL) err("**** Unable to allocate memory.");
   nspot = calloc(num_names,sizeof(int));
   if(nspot == NULL) err("**** Unable to allocate memory.");
   valmx = calloc(num_names,sizeof(double));
   if(valmx == NULL) err("**** Unable to allocate memory.");

   kcase = calloc(num_to_sort_by,sizeof(int));
   if(kcase == NULL) err("**** Unable to allocate memory.");
   ktol = calloc(num_to_sort_by,sizeof(int));
   if(ktol == NULL) err("**** Unable to allocate memory.");
   klocn = calloc(num_to_sort_by,sizeof(int));
   if(klocn == NULL) err("**** Unable to allocate memory.");
   dvals = calloc(num_to_sort_by,sizeof(double));
   if(dvals == NULL) err("**** Unable to allocate memory.");

/* --------------------------------------------- */
/* --------------------------------------------- */
/* --------------------------------------------- */

/* Get rid of  c_su_more from names= extension records. */

   int all_names = num_names;
   num_names = 0;
   for (n=0; n<all_names; n++) {
     if(strncmp(names[n],"c_su_more",9) != 0) {
       names[num_names] = names[n];
       num_names++;
     }
   }

/* Do not try to get c_su_id from record, it is S,R,X or whatever.      */

   if(strncmp(names[0],"c_su_id",7) == 0) names[0] = "null";

/* Parse the names to find character ranges and from,to,inc identifiers.           */
/* The names here could look like: gaps or gapis_rf or 55_gaps_58 or 55_gaps_rf_58 */
/* So, parse by underscore and figure out what is what.                            */

   int incomma = 0;
   int extra_parts = 0;
   for (n=0; n<num_names; n++) {

     cwp_String nparts[99]; /* really only 4 needed (currently) */
     int num_parts;

     strcpy(textbeg,names[n]);
     textbeg[strlen(names[n])] = '\0';
     tparse(textbeg, '_', nparts, &num_parts) ; 

     ilead[n]  = -1;
     itrail[n] = -1;
     namex[n] = "";
     if(num_parts==1) {
       names[n] = nparts[0];
     }
     else if(num_parts==2) {
       names[n] = nparts[0];
       namex[n] = nparts[1];
     }
     else if(num_parts==3) {
       ilead[n]  = atoi(nparts[0]);
       itrail[n] = atoi(nparts[2]);
       if(itrail[n]<ilead[n] || ilead[n]<1) { 
         err("**** Error: Your names= list has an entry with incorrect integers: %s",names[n]);
       }
       names[n]  = nparts[1];
     }
     else if(num_parts==4) {
       ilead[n]  = atoi(nparts[0]);
       itrail[n] = atoi(nparts[3]);
       if(itrail[n]<ilead[n] || ilead[n]<1) { 
         err("**** Error: Your names= list has an entry with incorrect integers: %s",names[n]);
       }
       names[n]  = nparts[1];
       namex[n]  = nparts[2];
     }
     else {
       err("**** Error: Your names= list has an entry that parses incorrectly: %s",names[n]);
     }

     if(strcmp(names[n],"null") != 0 && strcmp(names[n],"c_su_id") != 0 && strcmp(names[n],"numb") != 0) {
       if(num_parts<3) incomma = 1;
       if(num_parts==2 || num_parts==4) extra_parts++;
     }

   } /* end of   for (n=0; n<num_names; n++) {      */

   int lerr = 0;
   if(incomma==1 && irtype==0) {
     lerr = 1;
     warn("The rtype=fixed but at least one non-null name has no leading and trailing range.");
   }

   if(extra_parts!=0 && extra_parts!=5) {
     lerr = 1;
     warn("**** Error: Your names= list only has some of _cf _ct _ci _rf _rt. Need all 5 or none.");
   }

   if(incomma==1 && names_more==1) {
     lerr = 1;
     warn("**** Error: C_SU_MORE not permitted after C_SU_NAMES for comma-separated files.");
   }

/* substitute the actual match= key name for match1 and matche1 type specification */ 

   for (n=0; n<num_names; n++) {
     if(strncmp(names[n],"null",4) != 0 && strncmp(names[n],"numb",4) != 0) {

       if(strcmp(names[n],"match1") == 0 && num_to_sort_by>0) names[n] = match[0];
       if(strcmp(names[n],"match2") == 0 && num_to_sort_by>1) names[n] = match[1];
       if(strcmp(names[n],"match3") == 0 && num_to_sort_by>2) names[n] = match[2];
       if(strcmp(names[n],"match4") == 0 && num_to_sort_by>3) names[n] = match[3];
       if(strcmp(names[n],"match5") == 0 && num_to_sort_by>4) names[n] = match[4];
       if(strcmp(names[n],"match6") == 0 && num_to_sort_by>5) names[n] = match[5];
       if(strcmp(names[n],"match7") == 0 && num_to_sort_by>6) names[n] = match[6];
       if(strcmp(names[n],"match8") == 0 && num_to_sort_by>7) names[n] = match[7];
       if(strcmp(names[n],"match9") == 0 && num_to_sort_by>8) names[n] = match[8];
       if(strcmp(names[n],"matche1") == 0 && num_to_sort_by>0) names[n] = match[num_to_sort_by-1];
       if(strcmp(names[n],"matche2") == 0 && num_to_sort_by>1) names[n] = match[num_to_sort_by-2];
       if(strcmp(names[n],"matche3") == 0 && num_to_sort_by>2) names[n] = match[num_to_sort_by-3];
       if(strcmp(names[n],"matche4") == 0 && num_to_sort_by>3) names[n] = match[num_to_sort_by-4];
       if(strcmp(names[n],"matche5") == 0 && num_to_sort_by>4) names[n] = match[num_to_sort_by-5];
       if(strcmp(names[n],"matche6") == 0 && num_to_sort_by>5) names[n] = match[num_to_sort_by-6];
       if(strcmp(names[n],"matche7") == 0 && num_to_sort_by>6) names[n] = match[num_to_sort_by-7];
       if(strcmp(names[n],"matche8") == 0 && num_to_sort_by>7) names[n] = match[num_to_sort_by-8];
       if(strcmp(names[n],"matche9") == 0 && num_to_sort_by>8) names[n] = match[num_to_sort_by-9];
       if(strncmp(names[n],"match",5) == 0) {
         lerr = 1;
         warn("**** Error: Name  %s  could not be substituted from match= list.",names[n]);
       }

       for (m=n+1; m<num_names; m++) {
         if(strcmp(names[n],names[m]) == 0 && strcmp(namex[n],namex[m]) == 0) {  
           lerr = 1;
           warn("**** Error: Name  %s  exists at least twice in the names list.",names[n]);
         }
       }
     }
   }

   if(lerr>0) {
     err("**** Error: Related to names= or C_SU_NAMES record (details above).");
   }

/* Get key cases for switch.  */

   for(i=0; i<num_to_sort_by;i++) { 
     kcase[i] = GetCase(match[i]);
     if(kcase[i]<1) err("**** Error: a match name not recognized (or not allowed).");
   }

   int jscalel = 0;
   int jscalco = 0;
   int joffset = 0;
   int numcases = 0;
   for(i=0; i<num_names;i++) { 
     int iamcase = GetCase(names[i]);
     if(iamcase < 0) { /* name not recognized (and not null or numb either) */
       err("**** Error: Name  %s  in the names list is not recognized.",names[i]);
     }
     else if(iamcase > 0 && iamcase < 1000) { /* if not null or numb, will access/store this field. */
       ncase[numcases]  = iamcase;
       ilead[numcases]  = ilead[i];
       itrail[numcases] = itrail[i];
       names[numcases]  = names[i];
       namex[numcases]  = namex[i];
       if(irtype == 1) nspot[numcases] = i;
       else nspot[numcases] = numcases;

       char * ktype = hdtype(names[i]); 
       if(ktype[0] == 'i') valmx[numcases] = 2147483645.; /* giving leeway of 3 for odd compilors */
       else if(ktype[0] == 'h') valmx[numcases] = 32765.;
       else if(ktype[0] == 'u') valmx[numcases] = 65533.;
       else if(ktype[0] == 'f') valmx[numcases] = 3.4e37;  /* actually e38,  but giving leeway */    
       else                     valmx[numcases] = 1.7e307; /* actually e308, but giving leeway */

       if(strcmp(names[i],"gelev") == 0  || strcmp(names[i],"selev") == 0 ||
          strcmp(names[i],"sdepth") == 0 || strcmp(names[i],"gdel") == 0  ||
          strcmp(names[i],"sdel") == 0   || strcmp(names[i],"swdep") == 0 ||
          strcmp(names[i],"gwdep") == 0) {
         valmx[numcases] /= dscalel; /* if scaling by 10 on output, the input range is 10 smaller */
         jscalel = 1;
       }
       if(strcmp(names[i],"sx") == 0 || strcmp(names[i],"sy") == 0 ||
          strcmp(names[i],"gx") == 0 || strcmp(names[i],"gy") == 0) {
         valmx[numcases] /= dscalco;
         jscalco = 1;
       }
       if(strcmp(names[i],"offset") == 0) joffset = 1; 

       numcases++; 
     }
   }
   if(jscalel==0) {    
     if(iscaleldef==1) iscalel = 0;      /* defaulted scalel so do not set it */
     else err("**** Error: Cannot specify scalel= if not updating any of the 7 elevation related values.");
   }
   if(jscalco==0) {    
     if(iscalcodef==1) iscalco = 0;      /* defaulted scalco so do not set it */
     else err("**** Error: Cannot specify scalco= if not updating any of the 4 coordinate values.");
   }

/* if offset not in names list, recompute if some XY updating */
/* (but accept manual offset value, do not assume user wrong) */

   int ioffset = 0;
   if(joffset==0 && jscalco==1) ioffset = 1; 
  
/* Find the name with _cf _ct _ci and the name with _rf _rt appendices. */

   for(i=0; i<10;i++) mapx[i] = -1;

   int kerr = 0;

   if(extra_parts==5) {
     for(i=0; i<numcases;i++) { 
       if(strcmp(namex[i],"cf") == 0) {  
         if(mapx[0] != -1) {
           kerr = 1;
           warn("**** Error: Only one name with _cf appended is allowed.");
         }
         mapx[0] = i;
/* Do not set ncase[i] = 0 yet for this one. Want to include in sort. */
/* But set ncase[i] = 0 for _ct and _ci since we do not want to sort  */
/* by them or include them in update to output trace header.          */
       }
       else if(strcmp(namex[i],"ct") == 0) {  
         if(mapx[1] != -1) {
           kerr = 1;
           warn("**** Error: Only one name with _ct appended is allowed.");
         }
         mapx[1] = i;
         ncase[i] = 0; /* Set to NOT update the output header with this */
       }
       else if(strcmp(namex[i],"ci") == 0) {  
         if(mapx[2] != -1) {
           kerr = 1;
           warn("**** Error: Only one name with _ci appended is allowed.");
         }
         mapx[2] = i;
         ncase[i] = 0; /* Set to NOT update the output header with this */
       }
       else if(strcmp(namex[i],"rf") == 0) {  
         if(mapx[3] != -1) {
           kerr = 1;
           warn("**** Error: Only one name with _rf appended is allowed.");
         }
         mapx[3] = i;
/* Do not set ncase[i] = 0 for this. Will use to update header with result.*/
       }
       else if(strcmp(namex[i],"rt") == 0) {  
         if(mapx[4] != -1) {
           kerr = 1;
           warn("**** Error: Only one name with _rt appended is allowed.");
         }
         mapx[4] = i;
         ncase[i] = 0; /* Set to NOT update the output header with this */
       }
       else if(strcmp(names[i],"grnofr") == 0) { /* note names not namex */ 
         mapx[7] = i; /* if user does not use my standard name, gets no warning */
       }
       else if(strcmp(names[i],"grnlof") == 0) { /* note names not namex */
         mapx[8] = i; /* if user does not use my standard name, gets no warning */
       }

     }

     if(strcmp(names[mapx[0]],names[mapx[1]]) != 0 ||
        strcmp(names[mapx[0]],names[mapx[2]]) != 0) { 
       kerr = 1;
       warn("**** Error: _cf _ct _ci must be appended to the same name.");
     }
     if(strcmp(names[mapx[3]],names[mapx[4]]) != 0) { 
       kerr = 1;
       warn("**** Error: _rf _rt must be appended to the same name.");
     }

     if(kerr>0) {
       err("**** Error: Your _cf _ct _ci _rf _rt specification is incorrect.");
     }

     for(k=0; k<num_to_sort_by; k++) { /* will need the value from input header */ 
       if(strcmp(names[mapx[0]],match[k]) == 0) {
         mapx[5] = k;
         mapx[6] = kcase[k]; /* need this if this is a create run  */
       }
     }
     if(mapx[5] == -1) { /* should not be possible (as coded May 2021) */
       err("**** Error: Cannot find the name appended with _cf _ct and _ci in match= list.");
     }

   } /* end of  if(extra_parts==5)  */

/* Avoid precision issues by also storing sort fields as long long int. */
/* And then use them to sort/search.                                    */

   for(k=0; k<num_to_sort_by; k++) {
     int ifound = 0;
     for(n=0; n<numcases; n++) { /* non-null fields end-up in sequence at dfield[]     */
       if(ncase[n] == kcase[k]) {    /* fortunately, each name is unique case (as of now). */
         ktol[k] = n;                /* store the number where it ends-up in dfield[]      */
         ifound = 1; 
         break;
       }
     }
     if(ifound==0) err("**** Error: match= name not found in names=");
   }

/* For the create traces from nothing option, retain location of the match= values. */

   if(ncreate>0) {
     int jnum = 0;
     for(k=0; k<num_to_sort_by; k++) {
       for(n=0; n<numcases; n++) { 
         if(ncase[n] == kcase[k]) {
           klocn[jnum] = n; /* andre, simplify to klocn[k] = n; once you can test it */
           jnum++;
           break;
         }
       }
     }
   }

/* Stop the input match= name values from being updated to output. They have to   */
/* be stored from text file so they can be looked-up to find the correct record, */
/* so their value is already in trace header. This means the output switch will  */
/* cycle over null cases. This could be eliminated but it would require yet      */
/* another indexing array, so it is unlikely to be faster than the null cycling. */

   for(n=0; n<numcases; n++) { 
     for(k=0; k<num_to_sort_by; k++) {
       if(ncase[n] == kcase[k]) ncase[n] = 0;
     }
   }

   if(numR<1) {
     numR = countRec(fpR, textraw, maxtext, Rid, lenid, nicerecord);
     warn("Counted %d data records.",numR);
     if(numR == 0) err("**** No data records found. Wrong setid value? Wrong file?");
   }

   RecInfo = calloc(numR,sizeof(struct PointInfo)); 
   if(RecInfo == NULL) err("**** Unable to allocate memory for data values.");

/* The values from each record are going to be stored.              */ 
/* For quick "finding" we will sort as requested by user match=.    */
/* But we are only going to sort the pointers to the record values, */
/* not the record values themselves. The record values will stay    */
/* where they were stored during read-in.                           */ 
/*                                                                  */
/* We could just allocate a small chunk of memory for each record,  */
/* but by allocating contiguous memory and preserving the pointer   */ 
/* to that contiguous memory, we can still access them in their     */
/* input order if that becomes useful.                              */ 
/*                                                                  */ 
/* So, allocate contiguous block of doubles. Then pointer arithmatic*/
/* to divide that memory amoung the individual record pointers.     */

   double *dinput = calloc(numR*numcases,sizeof(double));
   if(dinput == NULL) {
     err("**** Unable to allocate memory for data fields per record.");
   }
   for(n=0; n<numR; n++) RecInfo[n].dfield = dinput + n*numcases;

   long long int *linput = calloc(numR*num_to_sort_by,sizeof(long long int));
   if(linput == NULL) {
     err("**** Unable to allocate memory for access fields per record.");
   }
   for(n=0; n<numR; n++) RecInfo[n].lfield = linput + n*num_to_sort_by;
 
/* Allocate for the record for use with bhigh function. */

   guy.dfield = calloc(numcases,sizeof(double));
   if(guy.dfield == NULL) {
     err("**** Unable to allocate memory for a record.");
   }

   guy.lfield = calloc(numcases,sizeof(long long int));
   if(guy.lfield == NULL) {
     err("**** Unable to allocate memory for a record.");
   }

/*  Get all record values into memory and perform checking.*/

   int count = 0;
   int ncount = 0;
   int lenerr = 0;
   int comerr = 0;
   int morerr = 0;
   int numerr = 0;
   int lrgerr = 0;
   int nblank = 0;
   int nextrow = 0;

/* these are for the unrepeat option */

   int nup = 0;
   long long int nlfield = -999999999999999ll;
   long long int incint = 100000000000000ll;
   long long int nrep = incint; /* always 1 here despite unrepeat value */

   while (fgets(textraw, maxtext, fpR) != NULL) { /*read a line*/
     ncount++;
     if(ncount>=nicerecord) {
       for(n=0; n<10; n++) textfront[n] = tolower(textraw[n]);
       if(strncmp(textfront,"c_su",4) == 0 || nextrow==1) {
         nextrow = 0;                       
         if(strncmp(textfront,"c_su_names",10) == 0 ||
            strncmp(textfront,"c_su_forms",10) == 0) nextrow = 1;
       }
       else {
         if(lenid<1 || strncmp(textraw,Rid,lenid) == 0) { /* Rid compare is case-sensitive */

           if(count<numR) {

             int lenraw = strlen(textraw);

             if(irtype == 0) {      /* input is fied format */

               int icomma = 0;
               for(ineed=0; ineed<numcases; ineed++) { 
                 if(itrail[ineed] >= lenraw) {
                   lenerr++;
                   if(lenerr<4) {
                     warn("Error at record %d   Record-too-short for requested fied ranges",ncount);
                     if(lenerr==3) {
                       warn("Have 3 Record-too-short warnings, no more will be printed.");
                     }
                   }
                   break;
                 }
                 strncpy(textbeg+icomma,textraw+ilead[ineed]-1,itrail[ineed]-ilead[ineed]+1);
                 icomma += itrail[ineed]-ilead[ineed]+1;
                 textbeg[icomma] = rdel;
                 icomma++;
               }
               strncpy(textraw,textbeg,icomma);
               textraw[icomma] = '\0';
             } /* else, it is already rdelimited, which getCSV handles (using nspot) */

             getCSV(textraw, textbeg, maxtext, rdel, 
                    RecInfo[count].dfield, nspot, numcases,
                    ncount, &comerr,&morerr,&numerr,&nblank);

             for(j=0; j<numcases; j++) { 
               if(RecInfo[count].dfield[j]<1.e308) { /* already counted as an error */
                 if(RecInfo[count].dfield[j]>valmx[j] || 0.-RecInfo[count].dfield[j]>valmx[j]) {
                   lrgerr = lrgerr + 1;
                   if(lrgerr<4) {
                     warn("Error at record %d number-too-large for SU name (%.2f)",
                     ncount,RecInfo[count].dfield[j]);
                     if(lrgerr==3) {
                       warn("Have 3 number-too-large-to-output warnings, no more will be printed.");
                     }
                   }
                 }
               }
             }

/* The channel range should be: _cf (from/lowest), _ct (to/highest), _ci (increment). */
/* Repair them, if they are not that way. This might not matter for the computation,  */
/* but the qsort includes _cf and bhigh finds the record with smallest _cf.           */
/* Note: This is not error-checking, that occurs after qsort.                         */

             if(mapx[3] > -1) { 
               double cfr = RecInfo[count].dfield[mapx[0]];
               double ctr = RecInfo[count].dfield[mapx[1]];
               double cir = RecInfo[count].dfield[mapx[2]];  
               double rfr = RecInfo[count].dfield[mapx[3]];
               double rtr = RecInfo[count].dfield[mapx[4]];
               if(cfr>ctr) {
                 RecInfo[count].dfield[mapx[0]] = ctr;
                 RecInfo[count].dfield[mapx[1]] = cfr;
                 RecInfo[count].dfield[mapx[3]] = rtr;
                 RecInfo[count].dfield[mapx[4]] = rfr;
               }
               if(fabs(RecInfo[count].dfield[mapx[0]] - RecInfo[count].dfield[mapx[1]]) < dtolh*2.) {
                 RecInfo[count].dfield[mapx[2]] = 1.;
               }
               else {
                 RecInfo[count].dfield[mapx[2]] = fabs(cir);
               }
             }

/* Cycle over the records, round-off and copy to lfield (used by qsort). */

             for(k=0; k<num_to_sort_by; k++) {
               RecInfo[count].lfield[k] = longt(RecInfo[count].dfield[ktol[k]],dtolh,dtol);
             }

/* Modify the first lfield for the unrepeat option. Example: if first match= */
/* is fldr and the first fldr has value of 5 then the longt just above here  */
/* with a dtol=100 has just stored 500 for that 5. The following code adds   */
/* large integers onto that, so 500 becomes 100000000000500. When there is   */
/* a reversal in fldr (such as 5,6,7,8,9,10,11,3,4,5,6 we increment large    */
/* integer. So the first 5 becomes 100000000000500 and the second 5 becomes  */
/* 200000000000500. Later, we do the same thing to the incomming trace fldr. */
/* The qsort and bhigh logic performs normally (which is exactly the idea).  */

             if(unrepeat > -2147483645) {
               if(nup>0) {
                 if(nlfield > RecInfo[count].lfield[0]) {
                   nup = -1;
                   nrep += incint;
                 }
               }
               else if(nup<0) {
                 if(nlfield < RecInfo[count].lfield[0]) {
                   nup = 1;
                   nrep += incint;
                 }
               }
               else { /* initialize nup on second record */
                 if(nlfield > -999999999999999ll) {
                   nup = 1;
                   if(nlfield > RecInfo[count].lfield[0]) nup = -1;
                 }
               }
               nlfield = RecInfo[count].lfield[0]; /* preserve unmodified (fldr) to compare next */
               RecInfo[count].lfield[0] += nrep;   /* modify this record (fldr) value for qsort  */
             }

           }
           count++;

         }  
       }
     }
   }

   if(nblank>0) {
     warn("Total all-blank fields: %d. Assumed zero for all.",nblank);
   }

   if(numerr>0) warn("Total errors for Field-unreadable as a number:        %d",numerr);
   if(morerr>0) warn("Total errors for Two-numbers in one field:            %d",morerr);
   if(lenerr>0) warn("Total errors for Record-too-short to get all values:  %d",lenerr);
   if(comerr>0) warn("Total errors for Not-enough-commas to get all values: %d",comerr);
   if(lrgerr>0) warn("Total errors for Number-too-large for SU name:        %d",lrgerr);
   if(lenerr>0 || comerr> 0 || morerr>0 || numerr>0 || lrgerr>0) {
     err("Errors detected while reading text file (see details above).");
   }

   if(unrepeat > -2147483645) {
     warn("For unrepeat option, the text-based incrementing integer ended at: %d ",nrep/incint);
   }
   warn("Have allocated memory to store values from %d records. Found %d records.",numR,count);
   if(count > numR) {
     err("Too many records read (%d) for your maxrecords= value (%d).",count,numR);
   }
   numR = count;

/* -----------------  */
/* -----------------  */

   qsort(RecInfo,numR,sizeof(struct PointInfo),compSort);

/* Error-check channel ranges: _cf (from/lowest), _ct (to/highest), _ci (increment). */
/* So we need to know which records are a set (typiclly the same fldr value).        */
/* But, it is whatever match= we sorted by except for the last key (the _cf number). */

   int lapover = 0;
   if(mapx[3] > -1) { 
     int l1verr = 0;
     int l2verr = 0;
     int l3verr = 0;
     int l4verr = 0;
     int l5verr = 0;
     int l6verr = 0;
     int l7verr = 0;
     int l8verr = 0;
     int l1same = 0;
     int l2same = 0;
     int l3same = 0;
     int l4same = 0;
     int l5same = 0;
     int l6same = 0;
     int l7same = 0;
     int l8same = 0;
     long long int ntop_grnofr = 0;
     long long int ntop_grnlof = 0;

     num_of_others = num_to_sort_by - 1; /* compOther uses num_of_others internally  */ 

/* You may have wondered what the channel increment was actually doing for a living. */
/* Well, this is it. Sometimes the channel numbers are deliberately overlapping.     */
/* For instance: from 1 to 11 increment 2 and from 2 to 12 increment 2.              */
/* This can put 1,2 on same receiver number and 3,4 on same receiver, and so on.     */
/* That was/is sometimes done for multi-component geophones/sensors, so you end-up   */
/* with 2 or more traces from the same geophone/sensor location.                     */
/* That possibility makes both checking and actual computation harder.               */

     int ntop = 0; 
     int npchan = 0;
     int npsegs = 0;
     double rpinc = 0.;

     for(n=1; n<=numR; n++) {

       if (n==numR || compOther(RecInfo+ntop,RecInfo+n) != 0) { 

         l1same = 0;
         l2same = 0;
         l3same = 0;
         l4same = 0;
         l5same = 0;
         l6same = 0;
         l7same = 0;
         l8same = 0;

         int nchan = 0;
         int nsegs = 0;

         for(m=ntop; m<n; m++) {

           long long int mcf = longt(RecInfo[m].dfield[mapx[0]],dtolh,dtol);
           long long int mct = longt(RecInfo[m].dfield[mapx[1]],dtolh,dtol);
           long long int mci = longt(RecInfo[m].dfield[mapx[2]],dtolh,dtol);  
           if(mci == 0) mci = 1;

           nchan += (mct - mcf) / mci + 1;
           nsegs++;

           double rinc = 1.;
           if(mct != mcf) rinc = (RecInfo[m].dfield[mapx[3]] - RecInfo[m].dfield[mapx[4]]) / ((mct - mcf) / mci);

           if(ntop>0 && l5same==0 && fabs(rinc-rpinc) > dtolh*2.) {
             l5same = 1;
             l5verr = l5verr + 1;
             if(l5verr<4) {
               int j = (int) (RecInfo[m].lfield[0]/dtol + 0.5); 
               warn("Warning: Receiver points per channel changed at %s= %d ",match[0],j);
               if(l5verr==3) {
                 warn("Have 3 Receiver-points-per-channel changed warnings, no more will be printed.");
               }
             }
           }
           rpinc = rinc;

           if(l1same==0 && mcf%mci != mct%mci) {
             l1same = 1;
             l1verr = l1verr + 1;
             if(l1verr<4) {
               int j = (int) (RecInfo[m].lfield[0]/dtol + 0.5);
               warn("Error: Layout ends do not conform to same increment at %s= %d ",match[0],j);
               if(l1verr==3) {
                 warn("Have 3 Layout-ends-do-not-conform-to-same-increment errors, no more will be printed.");
               }
             }
           }

           if(mapx[7]>-1) {
             if(m==ntop) {
               ntop_grnofr = longt(RecInfo[m].dfield[mapx[7]],dtolh,dtol);  
             }
             else {
               if(l7same==0 && ntop_grnofr != longt(RecInfo[m].dfield[mapx[7]],dtolh,dtol)) {
                 l7same = 1;
                 l7verr = l7verr + 1;
                 if(l7verr<4) {
                   int j = (int) (RecInfo[m].lfield[0]/dtol + 0.5);
                   warn("Error: The grnofr values are different for same shot at %s= %d ",match[0],j);
                   if(l7verr==3) {
                     warn("Have 3 grnofr-values-are-different-for-same-shot, no more will be printed.");
                   }
                 }
               }
             }
           } 

           if(mapx[8]>-1) {
             if(m==ntop) {
               ntop_grnlof = longt(RecInfo[m].dfield[mapx[8]],dtolh,dtol);  
             } 
             else {
               if(l8same==0 && ntop_grnlof != longt(RecInfo[m].dfield[mapx[8]],dtolh,dtol)) { 
                 l8same = 1;
                 l8verr = l8verr + 1;
                 if(l8verr<4) {
                   int j = (int) (RecInfo[m].lfield[0]/dtol + 0.5);
                   warn("Error: The grnlof values are different for same shot at %s= %d ",match[0],j);
                   if(l8verr==3) {
                     warn("Have 3 grnlof-values-are-different-for-same-shot, no more will be printed.");
                   }
                 }
               }
             }
           } 

           for(i=m+1; i<n-1; i++) {
             long long int icf = longt(RecInfo[i].dfield[mapx[0]],dtolh,dtol);
             long long int ici = longt(RecInfo[i].dfield[mapx[2]],dtolh,dtol);  
             if(icf > mct) break;     /* Current -from- greater than the other -to-*/ 
             if(lapover<1) {
               int j = (int) (RecInfo[i].lfield[0]/dtol + 0.5);
               warn("Warning: Overlapping channel range at %s= %d  Unusual, but not always an error.",match[0],j);
             }
             lapover++;

/* Different increments in overlapping ranges might work, but complex and more     */
/* likely just an error. */

             if(l2same==0 && ici != mci) {
               l2same = 1;
               l2verr = l2verr + 1;
               if(l2verr<4) {
                 int j = (int) (RecInfo[i].lfield[0]/dtol + 0.5);
                 warn("Error: Different increments in overlapping layout at %s= %d ",match[0],j);
                 if(l2verr==3) {
                   warn("Have 3 Different-increments-in-overlapping layout errors, no more will be printed.");
                 }
               }
             }
/* Going to hit the same numbers within the overlap?                               */
             if(l3same==0 && mcf%mci == icf%ici) {
               l3same = 1;
               l3verr = l3verr + 1;
               if(l3verr<4) {
                 int j = (int) (RecInfo[i].lfield[0]/dtol + 0.5);
                 warn("Error: Overlapping layout hits same channels at %s= %d ",match[0],j);
                 if(l3verr==3) {
                   warn("Have 3 Overlapping-layout-hits-same-channels errors, no more will be printed.");
                 }
               }
             }

           } /* end of  for(i=m+1; i<n-1; i++) { */

         } /* end of  for(m=ntop; m<n; m++) {  */


         if(ntop>0 && l4same==0 && nchan != npchan) {
           l4same = 1;
           l4verr = l4verr + 1;
           if(l4verr<4) {
             int j = (int) (RecInfo[ntop].lfield[0]/dtol + 0.5);
             warn("Warning: Number of channels in layout changed at %s= %d ",match[0],j);
             if(l4verr==3) {
               warn("Have 3 Number-of-channels-in-layout-changed warnings, no more will be printed.");
             }
           }
         }
         npchan = nchan;

         if(ntop>0 && l6same==0 && nsegs != npsegs) {
           l6same = 1;
           l6verr = l6verr + 1;
           if(l6verr<4) {
             int j = (int) (RecInfo[ntop].lfield[0]/dtol + 0.5); 
             if(nsegs==npsegs*2 || nsegs*2==npsegs) {
               warn("Warning: Number of segments in layout changed by ratio 2 at %s= %d (see unrepeat=).",match[0],j);
             }
             else {
               warn("Warning: Number of segments in layout changed at %s= %d ",match[0],j);
             }
             if(l6verr==3) {
               warn("Have 3 Number-of-segments-in-layout-changed warnings, no more will be printed.");
             }
           }
         }
         npsegs = nsegs;

         ntop = n;  

       } /* end of  if (compOther(RecInfo+ntop,RecInfo+n) != 0 || n==numR-1) */

     } /* end of  for(n=1; n<numR; n++) {  */


     if(lapover>0) warn("Total Overlapping-channel-ranges in layouts:        %d (very unusual)",lapover);
     if(l2verr>0) warn("Total Different-increments-in-overlapping layout:   %d (error-halt)",l2verr);
     if(l3verr>0) warn("Total Overlapping-layout-hits-same-channels:        %d (error-halt)",l3verr);
     if(l1verr>0) warn("Total Layout-ends-do-not-conform-to-same-increment: %d (error-halt)",l1verr);
     if(l4verr>0) warn("Total Number-of-channels-in-layout-changed:         %d (unusual)",l4verr);
     if(l5verr>0) warn("Total Receiver-points-per-channel changed:          %d (very unusual)",l5verr);
     if(l6verr>0) warn("Total Number-of-segments-in-layout-changed:         %d (unusual)",l6verr);
     if(l7verr>0) warn("Total grnofr-values-are-different-for-same-shot:    %d (error-halt)",l7verr);
     if(l8verr>0) warn("Total grnlof-values-are-different-for-same-shot:    %d (error-halt)",l8verr);

     if(l1verr>0 || l2verr>0 || l3verr>0 || l7verr>0 || l8verr>0) {
       err("Errors detected while reading text file (see details above).");
     }

   } /* end of  if(mapx[3] > -1) { */

   else {

     int l7verr = 0;
     for(n=1; n<numR; n++) { 
       if(compSort(RecInfo+n,RecInfo+n-1) == 0) { /* check for duplicate records */ 
         l7verr = l7verr + 1;
         if(l7verr<4) {
           if(num_to_sort_by==1) {
             int j = (int) (RecInfo[n].lfield[0]/dtol + 0.5); 
             warn("Error: Records have duplicate match values of %s=%d",match[0],j);
           }
           else if(num_to_sort_by==2) {
             int j = (int) (RecInfo[n].lfield[0]/dtol + 0.5); 
             int k = (int) (RecInfo[n].lfield[1]/dtol + 0.5); 
             warn("Error: Records have duplicate match values of %s=%d   %s=%d",match[0],j,match[1],k);
           }
           else if(num_to_sort_by==3) {
             int j = (int) (RecInfo[n].lfield[0]/dtol + 0.5); 
             int k = (int) (RecInfo[n].lfield[1]/dtol + 0.5); 
             int l = (int) (RecInfo[n].lfield[2]/dtol + 0.5); 
             warn("Error: Records have duplicate match values of %s=%d   %s=%d   %s=%d",
                  match[0],j,match[1],k,match[2],l);
           }
           else {
             int j = (int) (RecInfo[n].lfield[0]/dtol + 0.5); 
             int k = (int) (RecInfo[n].lfield[1]/dtol + 0.5); 
             int l = (int) (RecInfo[n].lfield[2]/dtol + 0.5); 
             int m = (int) (RecInfo[n].lfield[3]/dtol + 0.5); 
             warn("Error: Records have duplicate match values of %s=%d   %s=%d   %s=%d   %s=%d",
                  match[0],j,match[1],k,match[2],l,match[3],m);
           }
           if(l7verr==3) {
             warn("Have 3 Records-have-duplicate-match-values errors, no more will be printed.");
           }
         }
       } /* end of  if(compSort(RecInfo+n,RecInfo+n-1) == 0) { */
     }

     if(l7verr>0) err("Total  errors  for Records-have-duplicate-match-values: %d",l7verr);

     if(num_to_sort_by==1 || num_to_sort_by==2) { 
       int jn = 0;
       int jp = 0;
       int jd = 0;
       int jc = 1;
       int je = -1;
       int kn = 0;
       int kp = 0;
       int kd = 0;
       int kc = -1;

       for(n=0; n<numR; n++) { 

         if(num_to_sort_by==2) {
           jn = (int) (RecInfo[n].lfield[0]/dtol + 0.5); 
           kn = (int) (RecInfo[n].lfield[1]/dtol + 0.5); 
         }
         else kn = (int) (RecInfo[n].lfield[0]/dtol + 0.5); 

         if(n==0) {
           jp = jn;
           kp = kn;
         }

         if(jn != jp) {          /* new line                      */
           jc++;                 /* count of number of lines      */
           if(jd != jn-jp) je++; /* count line difference changes */
           jd = jn-jp;           /* reset line difference         */
         }
         else {                  /* not a new line                 */
           if(kn-kp != kd) {     /* check point difference         */
             kc++;               /* count point difference changes */
             kd = kn-kp;         /* reset point difference         */
           }
         }
         jp = jn;
         kp = kn;

       } /* end of  for(n=0; n<numR; n++) { */

       if(num_to_sort_by==2) {
         warn("Note: There are: %d sets of %s values (lines?).",jc,match[0]);
         if(je>1) warn("Warning: There are: %d irregular %s increments between the sets (missing lines?).",(je+1)/2,match[0]); 
         if(kc>1) warn("Warning: There are: %d irregular %s increments within the lines (missing points?).",(kc+1)/2,match[1]); 
       }
       else {
         if(kc>1) warn("Warning: There are: %d irregular %s increments within the line (missing points?).",(kc+1)/2,match[0]); 
       } 

     }

   }

   /* ==========================================  */

   /* Get info from first trace */ 

   if(ncreate<1 || ireplicate>0) { /* normal input run? or replicating one input only? */
     if(!gettr(&tr)) err("can't read first trace");
     for(n=1; n<ireplicate; n++) {
       if(!gettr(&tr)) err("can't read enough traces for spikes=n option.");
     }
   }
   else { /* going to create first trace from nothing and maybe put spikes on it */
     if(nspikes==0) {
       tr.ns = 1; /* samples */
       tr.dt = NINT(1000.0 * 4.);
       tr.data[0] = 0.;
     }
     else {
       float tmax = 0.;
       for(n=2; n<nspikes; n+=2) {
         if(spikes[n] > tmax) tmax = spikes[n];
       }
       int nsamps = NINT((tmax)/spikes[0]);
       tr.ns = nsamps; 
       tr.dt = NINT(1000.0 * spikes[0]);
       for(n=0; n<nsamps; n++) tr.data[n] = spikes[1];
       for(n=0; n<nspikes/2; n++) { /* n=0 includes the 0 sample in negative check */
         int m = NINT(spikes[n*2]/spikes[0]) - 1;
         if(m < 0) err("Error: Spike times less than first sample are not permitted.");
         tr.data[m] = spikes[n*2+1]; 
       }
     }

   }

   /* ==========================================  */
   /* ======  Main loop over traces  ===========  */
   /* ==========================================  */

   nup = 0;
   nrep = unrepeat * incint; /* yes, set to user value of unrepeat= here */
   nlfield = -999999999999999ll;

   int notfound = 0;
   int nmade = 0;
   double dmade = 0.;
   if(mapx[3] > -1) {
     dmade = RecInfo[0].dfield[mapx[0]] - RecInfo[0].dfield[mapx[2]];
   }

   int igo = 1; 
   do {

/* When creating traces from nothing, or replicating from a user      */
/* selected input trace, first update the created/selected header from*/
/* the match= values in the text input file and then proceed as if    */
/* that header had been input. So we put values in created/selected   */
/* header and then just extract them again? Yes. We do not bypass this*/
/* because we want things to proceed like the actual input situation. */
/* In particular, we do not want to use the double values directly,   */
/* we want to convert them into their header keys, then convert back. */

     if(ncreate>0) {

       if(nproct==ncreate) break;

       for(n=0;n<num_to_sort_by;n++) { 
         double dvalue = RecInfo[nmade].dfield[klocn[n]];
         tohead(&tr, kcase[n], dvalue);
       }

       if(mapx[3] < 0) { 
         nmade++;  
         if(nmade==numR) igo = 0;
       }
       else { /* special for SPS X-type files (and similar).*/
         dmade += RecInfo[nmade].dfield[mapx[2]]; 
         if(dmade>RecInfo[nmade].dfield[mapx[1]]) {
           nmade++;
           if(nmade==numR) break;
           for(n=0;n<num_to_sort_by;n++) { 
             double dvalue = RecInfo[nmade].dfield[klocn[n]];
             tohead(&tr, kcase[n], dvalue);
           }
           dmade = RecInfo[nmade].dfield[mapx[0]];
         }
         tohead(&tr, mapx[6], dmade);
       }

     } /* end of  if(ncreate>0) {  */

/* Get the match= values from the input header.          */

     for(n=0;n<num_to_sort_by;n++) { 
       dvals[n] = fromhead(tr, kcase[n]); 
       guy.lfield[n] = longt(dvals[n],dtolh,dtol);
     } 

/* -----------------  */

     if(unrepeat > -2147483645) {
       if(nup>0) {
         if(nlfield > guy.lfield[0]) {
           nup = -1;
           nrep += incint;
         }
       }
       else if(nup<0) {
         if(nlfield < guy.lfield[0]) {
           nup = 1;
           nrep += incint;
         }
       }
       else { /* initialize nup on second trace */
         if(nlfield > -999999999999999ll) {
           nup = 1;
           if(nlfield > guy.lfield[0]) nup = -1;
         }
       }
       nlfield = guy.lfield[0]; /* preserve unmodified (fldr) to compare next */
       guy.lfield[0] += nrep;   /* modify this traces (fldr) value for bhigh  */
     }

/* -----------------  */

/* Use the match= values to find where the text record values are for trace.      */
/* bhigh returns 1 more than exact match but also returns same idxR for exact+0.1 */ 
/* bhigh was intended to do this because X-record specification has _cf value as  */ 
/* the lowest end of that (channel) range and makes the following code simpler.   */ 

     idxR = bhigh(RecInfo, numR, &guy);

/* If zero retained next line, compSort will return != 0, so ifound will be <0 soon.*/ 
/* Similarly, if numR returned from bhigh, will be numR-1 and compOther checks it.  */ 
     if(idxR>0) idxR--; 

     int ifound = 1; 

     if(mapx[3] < 0) { /* not X-records */
       if (compSort(RecInfo+idxR,&guy) != 0) { 
         ifound = -1;
       }
     }
     else { /* is X-records or similar */

/* If there is no match to rest of key(s) values, then this trace is less than lowest _cf */
/* (from-channel) number amoung the set of records for this shot (same fldr, usually).   */
/* Or, of course, no such fldr number exists in text records.                             */
       if (compOther(RecInfo+idxR,&guy) != 0) ifound = -2;

/* If we do match, then is trace channel greater than _ct (to-channel) of this record - not */
/* just the highest record of the fldr set. Remember that channel numbers do not have to be */
/* contiguous, some can be left out of layout definition X-records and we need to make sure */
/* that this input trace is BETWEEN the from-to channel range.                              */
       else {
         if(lapover==0) {
           if(guy.lfield[num_to_sort_by-1] 
              > longt(RecInfo[idxR].dfield[mapx[1]],dtolh,dtol)) ifound = -3; 
         }

/* If overlap exists, grind thru the records with lower _cf (from-channel) numbers         */
/* and see which incrementing set of channels the input trace channel fits into.           */
         else { 
           ifound = -4; 
           for (m=idxR; m>=0; m--) {
             if (compOther(RecInfo+m,&guy) != 0) break;
             if(guy.lfield[num_to_sort_by-1] 
                > longt(RecInfo[m].dfield[mapx[1]],dtolh,dtol)) continue;
             long long int mcf = longt(RecInfo[m].dfield[mapx[0]],dtolh,dtol);
             long long int mci = longt(RecInfo[m].dfield[mapx[2]],dtolh,dtol);  
             if(guy.lfield[num_to_sort_by-1]%mci == mcf%mci) {
               ifound = 2;
               idxR = m;
               break;
             } 
           } /* end of  for (m=idxR; m>=0; m--) { */
         } 
       } 

     } /* end of else  is X-records or similar */

     if(ifound>0) { 

       if(iscalel!=0) tr.scalel = iscalel;
       if(iscalco!=0) {
         tr.scalco = iscalco;
         tr.counit = 1;
       }

/* Get the names= values and update the output header.   */

       for(n=0;n<numcases;n++) { 

         double dvalue = RecInfo[idxR].dfield[n];

         if(mapx[3] == n) { /* special computaion for SPS X-type files (and similar).*/
           double cfo = RecInfo[idxR].dfield[mapx[0]];
           double cto = RecInfo[idxR].dfield[mapx[1]];
/*         double cio = RecInfo[idxR].dfield[mapx[2]];  increment is unneeded (here) */
           double rfo = RecInfo[idxR].dfield[mapx[3]];
           double rto = RecInfo[idxR].dfield[mapx[4]];

           if(cto - cfo < dtolh*2.) {
             dvalue = rfo;
           }
           else {  
             dvalue = rfo + (dvals[mapx[5]] - cfo) * (rto - rfo) / (cto - cfo);
           }

         }

         if(ncase[n]>12 && ncase[n]<20) dvalue *= dscalel;
         if(ncase[n]>21 && ncase[n]<26) dvalue *= dscalco;

         if(iaction!=0) dvalue = iaction * dvalue + fromhead(tr, ncase[n]);

         tohead(&tr, ncase[n], dvalue);

        } /* end of   for(n=0,n<numcases,n++) {  */ 

        if(ioffset!=0) {
          tr.offset = sqrt((tr.sx-tr.gx)/dscalco*(tr.sx-tr.gx)/dscalco + (tr.sy-tr.gy)/dscalco*(tr.sy-tr.gy)/dscalco);
        }
        
      } /* end of  if(ifound>0) { */

      if(ifound<0) {
        notfound++;
        if(missing==-1) {
          err(" Error: Could not find a match= text record for trace and missing= option is error-halt.");
        }
      }

/* delete missing = 0  pass missing = 1  error-halt missing = -1   */
      if(ifound>0 || missing==1) {
        puttr(&tr);
        nproct++; /* count output traces */
      }

      if(ncreate<1) {
        if(!gettr(&tr)) igo = 0; 
      }

   } while (igo);

/* } while (gettr(&tr)); */

/* =========  End of Trace loop =============  */

   if(unrepeat > -2147483645) {
     warn("For unrepeat option, the trace-based incrementing integer ended at: %d ",nrep/incint);
   }

   warn("Number of traces output: %d ",nproct);
   if(notfound>0) {
     if(missing==0) {
       warn("Traces not output due to no match= in text file: %d (your missing= option)",notfound);
     }
     else {
       warn("Traces output anyway despite no match= in text file: %d (your missing= option)",notfound);
     }
   }

   return(CWP_Exit());
}

int countRec(FILE *ufile, char *textraw, int maxtext, char *Rid, int lenid, int nicerecord) {

/*  This function reads the *ufile looking for lines with */     
/*  1st field == Rid.   The file must be already opened   */
/*  and will be repositioned at start point at exit.      */

   int count = 0;
   int ncount = 0;
   int nextrow = 0;
   int n=0;
   char textfront[10];    

   while (fgets(textraw, maxtext, ufile) != NULL) { /* read a line */
     ncount++;
     if(ncount>=nicerecord) {
       for(n=0; n<10; n++) textfront[n] = tolower(textraw[n]);
       if(strncmp(textfront,"c_su",4) == 0 || nextrow==1) {
         nextrow = 0; 
         if(strncmp(textfront,"c_su_names",10) == 0 || 
            strncmp(textfront,"c_su_forms",10) == 0) nextrow = 1; 
       }
       else {
         if(lenid<1 || strncmp(textraw,Rid,lenid) == 0) count++;
       }
     }
   }

   fseek(ufile, 0L, SEEK_SET);    /* reposition input file */

   return count;

}

void getCSV(char *textraw, char *textbeg, int maxtext, char rdel, 
            double *dfield, int *nspot, int numcases,   
            int ncount, int *comerr,int *morerr,int *numerr,int *nblank) {

  int nbeg = -1;
  int nfield = 0;
  int ineed  = 0;
  int igot;
  int n=0;
  int m=0;
  double dval=0.0;

  for(n=0; n<maxtext; n++) {                         /* linux \n            windows \r */
    if(textraw[n] == rdel || textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') {
      if(nfield == nspot[ineed]) {
        dval = 1.1e308;
        int nb = -1;
        if(n-nbeg-1 > 0) {
          strncpy(textbeg,textraw+nbeg+1,n-nbeg-1);
          textbeg[n-nbeg-1] = '\0'; /* so sscanf knows where to stop */
          int ib = -1;
          for (m=0; m<n-nbeg-1; m++) {
            if(textbeg[m] != ' ') {
              nb = m;
              if(ib>-1) {
                *morerr = *morerr + 1;
                if(*morerr<4) {
                  warn("Error at record %d   two-numbers in field (%s)",ncount,textbeg);
                  if(*morerr==3) warn("Have 3 two-numbers in field warnings, no more will be printed.");
                }
                nb=-1;
                break;
              }
            }
            if(textbeg[m] == ' ' && nb>-1) ib = m;
          }  

          if(nb>-1) {
            igot = sscanf(textbeg,"%lf",&dval);  
            if(igot<1) {
              *numerr = *numerr + 1;
              if(*numerr<4) {
                warn("Error at record %d   field-unreadable as a number (%s)",ncount,textbeg);
                if(*numerr==3) warn("Have 3 field-unreadable warnings, no more will be printed.");
              }
              break;
            }
          }
        } /* end of  if(n-nbeg-1 > 0) { */
        if(nb<0) {
          *nblank = *nblank + 1;
          dval = 0.; 
        }
        dfield[ineed] = dval;
        ineed++;
        if(ineed>=numcases) break; 
      }
      if(textraw[n] == '\0' || textraw[n] == '\n' || textraw[n] == '\r') {
        if(ineed<numcases) {
          *comerr = *comerr + 1;
          if(*comerr<4) {
            warn("Error at record %d   Not-enough-commas in record to get all values",ncount);
            if(*comerr==3) warn("Have 3 Not-enough-comma warnings, no more will be printed.");
          }
        }
      }
      nbeg = n;
      nfield++;
    }
  } 
  return;
}

/* --------------------------- */

int GetCase(char* cbuf) {
   
       int ncase = -1;
   
       if(strncmp(cbuf,"null",4) == 0) ncase = 0;  /* any name starting with null */
       else if(strcmp(cbuf,"tracl") == 0) ncase = 1;
       else if(strcmp(cbuf,"tracr") == 0) ncase = 2;
       else if(strcmp(cbuf,"fldr" ) == 0) ncase = 3;
       else if(strcmp(cbuf,"tracf") == 0) ncase = 4;
       else if(strcmp(cbuf,"ep"   ) == 0) ncase = 5;
       else if(strcmp(cbuf,"cdp") == 0) ncase = 6;
       else if(strcmp(cbuf,"cdpt") == 0) ncase = 7;
       else if(strcmp(cbuf,"trid") == 0) ncase = 8;
       else if(strcmp(cbuf,"nvs") == 0) ncase = 9;
       else if(strcmp(cbuf,"nhs") == 0) ncase = 10;
       else if(strcmp(cbuf,"duse") == 0) ncase = 11;
       else if(strcmp(cbuf,"offset") == 0) ncase = 12;
       else if(strcmp(cbuf,"gelev") == 0) ncase = 13;
       else if(strcmp(cbuf,"selev") == 0) ncase = 14;
       else if(strcmp(cbuf,"sdepth") == 0) ncase = 15;
       else if(strcmp(cbuf,"gdel") == 0) ncase = 16;
       else if(strcmp(cbuf,"sdel") == 0) ncase = 17;
       else if(strcmp(cbuf,"swdep") == 0) ncase = 18;
       else if(strcmp(cbuf,"gwdep") == 0) ncase = 19;
       else if(strcmp(cbuf,"scalel") == 0) ncase = 20;
       else if(strcmp(cbuf,"scalco") == 0) ncase = 21;
       else if(strcmp(cbuf,"sx") == 0) ncase = 22;
       else if(strcmp(cbuf,"sy") == 0) ncase = 23;
       else if(strcmp(cbuf,"gx") == 0) ncase = 24;
       else if(strcmp(cbuf,"gy") == 0) ncase = 25;
       else if(strcmp(cbuf,"counit") == 0) ncase = 26;
       else if(strcmp(cbuf,"wevel") == 0) ncase = 27;
       else if(strcmp(cbuf,"swevel") == 0) ncase = 28;
       else if(strcmp(cbuf,"sut") == 0) ncase = 29;
       else if(strcmp(cbuf,"gut") == 0) ncase = 30;
       else if(strcmp(cbuf,"sstat") == 0) ncase = 31;
       else if(strcmp(cbuf,"gstat") == 0) ncase = 32;
       else if(strcmp(cbuf,"tstat") == 0) ncase = 33;
       else if(strcmp(cbuf,"laga") == 0) ncase = 34;
       else if(strcmp(cbuf,"lagb") == 0) ncase = 35;
       else if(strcmp(cbuf,"delrt") == 0) ncase = 36;
       else if(strcmp(cbuf,"muts") == 0) ncase = 37;
       else if(strcmp(cbuf,"mute") == 0) ncase = 38;
       else if(strcmp(cbuf,"ns") == 0) ncase = 39;
       else if(strcmp(cbuf,"dt") == 0) ncase = 40;
       else if(strcmp(cbuf,"gain") == 0) ncase = 41;
       else if(strcmp(cbuf,"igc") == 0) ncase = 42;
       else if(strcmp(cbuf,"igi") == 0) ncase = 43;
       else if(strcmp(cbuf,"corr") == 0) ncase = 44;
       else if(strcmp(cbuf,"sfs") == 0) ncase = 45;
       else if(strcmp(cbuf,"sfe") == 0) ncase = 46;
       else if(strcmp(cbuf,"slen") == 0) ncase = 47;
       else if(strcmp(cbuf,"styp") == 0) ncase = 48;
       else if(strcmp(cbuf,"stas") == 0) ncase = 49;
       else if(strcmp(cbuf,"stae") == 0) ncase = 50;
       else if(strcmp(cbuf,"tatyp") == 0) ncase = 51;
       else if(strcmp(cbuf,"afilf") == 0) ncase = 52;
       else if(strcmp(cbuf,"afils") == 0) ncase = 53;
       else if(strcmp(cbuf,"nofilf") == 0) ncase =54;
       else if(strcmp(cbuf,"nofils") == 0) ncase = 55;
       else if(strcmp(cbuf,"lcf") == 0) ncase = 56;
       else if(strcmp(cbuf,"hcf") == 0) ncase = 57;
       else if(strcmp(cbuf,"lcs") == 0) ncase = 58;
       else if(strcmp(cbuf,"hcs") == 0) ncase = 59;
       else if(strcmp(cbuf,"year") == 0) ncase = 60;
       else if(strcmp(cbuf,"day") == 0) ncase = 61;
       else if(strcmp(cbuf,"hour") == 0) ncase = 62;
       else if(strcmp(cbuf,"minute") == 0) ncase = 63;
       else if(strcmp(cbuf,"sec") == 0) ncase = 64;
       else if(strcmp(cbuf,"timbas") == 0) ncase = 65;
       else if(strcmp(cbuf,"trwf") == 0) ncase = 66;
       else if(strcmp(cbuf,"grnors") == 0) ncase = 67;
       else if(strcmp(cbuf,"grnofr") == 0) ncase = 68;
       else if(strcmp(cbuf,"grnlof") == 0) ncase = 69;
       else if(strcmp(cbuf,"gaps") == 0) ncase = 70;
       else if(strcmp(cbuf,"otrav") == 0) ncase = 71;
       else if(strcmp(cbuf,"d1") == 0) ncase = 72;
       else if(strcmp(cbuf,"f1") == 0) ncase = 73;
       else if(strcmp(cbuf,"d2") == 0) ncase = 74;
       else if(strcmp(cbuf,"f2") == 0) ncase = 75;
       else if(strcmp(cbuf,"ungpow") == 0) ncase = 76;
       else if(strcmp(cbuf,"unscale") == 0) ncase = 77;
       else if(strcmp(cbuf,"ntr") == 0) ncase = 78;
       else if(strcmp(cbuf,"mark") == 0) ncase = 79;
       else if(strncmp(cbuf,"numb",4) == 0) {
         ncase = 1000 + atoi(cbuf+4);
       }
  
   return ncase;

}

/* --------------------------- */
/* Convert from double to long long int with tolerance factor and multiplier. */

long long int longt (double dvalue, double dtolh, double dtol) {

  long long int ltint; 

  if(dvalue >= 0.0) ltint = (long long int) ((dvalue + dtolh) * dtol);
  else              ltint = (long long int) ((dvalue - dtolh) * dtol);

  return (ltint);

}

/* --------------------------- */
/* Specify compare function for qsort and my bhigh function.  */

int compSort (const void * q1, const void * q2) {
  
  int n=0;

  struct PointInfo* p1 = (struct PointInfo*) q1;
  struct PointInfo* p2 = (struct PointInfo*) q2;

  for(n=0; n<num_to_sort_by; n++) {  
    if(p1->lfield[n] < p2->lfield[n]) return (-1);
    if(p1->lfield[n] > p2->lfield[n]) return (1);
  }

  return (0);

}

/* --------------------------- */
/* Specify compare function for limited searching.            */

int compOther (const void * q1, const void * q2) {

  struct PointInfo* p1 = (struct PointInfo*) q1;
  struct PointInfo* p2 = (struct PointInfo*) q2;

  int n=0;

/* note num_of_others==0 returns with 0 (equal) which is */
/* nice for interpolation options since they work with   */
/* last match= and expect previous match= to be equal    */

  for(n=0; n<num_of_others; n++) {  
    if(p1->lfield[n] < p2->lfield[n]) return (-1);
    if(p1->lfield[n] > p2->lfield[n]) return (1);
  }

  return (0);

}

/* --------------------------- */
int bhigh(struct PointInfo *all, int last, struct PointInfo* guy) {

  int mid;
  int low = 0;
  int high = last;

/* This is just a standard binary search where return from   */ 
/* bhigh lands 1 above exact but stays there for exact+0.1   */ 
/* which is easier to deal with than the blow code because   */ 
/* of the X-record type specification where the _cf value is */ 
/* the lowest value in that (channel) range and qsort has    */ 
/* _cf include in the sort match=                            */ 

  while (low < high) {
    mid = low + (high - low) / 2;
    if (compSort(guy,all+mid) >= 0) low = mid +1; 
    else high = mid;                               
  }

/* blow lands at exact but then goes up 1 for exact+0.1     */
/*  if (compSort(guy,all+mid) <= 0) high = mid;             */ 
/*  else low = mid + 1;                                     */

  return low; 
}

/* --------------------------- */
/* expects a string with no blanks and no tabs \t   */
void tparse(char *tbuf, char d, char **fields, int *numfields) { 
  int nbeg = -1;
  *numfields = 0;
  int n=0;

  for(n=0; ; n++) {
    if(tbuf[n] == d || tbuf[n] == '\0') {
      if(n-nbeg-1 > 0) {
        fields[*numfields] = ealloc1(n-nbeg-1,1);
        strncpy(fields[*numfields],tbuf+nbeg+1,n-nbeg-1);
      }
      else {
        fields[*numfields] = ealloc1(4,1);
        fields[*numfields][0] = 'n';
        fields[*numfields][1] = 'u';
        fields[*numfields][2] = 'l';
        fields[*numfields][3] = 'l';
/*      strncpy(fields[*numfields],"null",4);  makes compilor unhappy */
      }
      nbeg = n;
      *numfields = *numfields + 1;
    }
    if(tbuf[n] == '\0') break;
  }
 
}

/* --------------------------- */
double fromhead(segy tr, int k) {

       double dval=0.0;

       switch (k) {
   
         case -1: 
/*       null, name not found? */
         break;
         case 0:  
/*       null   do not read from header */ 
         break;
         case 1:
           dval = tr.tracl;
         break;
         case 2:
           dval = tr.tracr;
         break;
         case 3:
           dval = tr.fldr;
         break;
         case 4:
           dval = tr.tracf;
         break;
         case 5:
           dval = tr.ep;
         break;
         case 6:
           dval = tr.cdp;
         break;
         case 7:
           dval = tr.cdpt;
         break;
         case 8:
           dval = tr.trid;
         break;
         case 9:
           dval = tr.nvs;
         break;
         case 10:
           dval = tr.nhs;
         break;
         case 11:
           dval = tr.duse;
         break;
         case 12:
           dval = tr.offset;
         break;
         case 13:
           dval = tr.gelev;
         break;
         case 14:
           dval = tr.selev;
         break;
         case 15:
           dval = tr.sdepth;
         break;
         case 16:
           dval = tr.gdel;
         break;
         case 17:
           dval = tr.sdel;
         break;
         case 18:
           dval = tr.swdep;
         break;
         case 19:
           dval = tr.gwdep;
         break;
         case 20:
           dval = tr.scalel;
         break;
         case 21:
           dval = tr.scalco;
         break;
         case 22:
           dval = tr.sx;
         break;
         case 23:
           dval = tr.sy;
         break;
         case 24:
           dval = tr.gx;
         break;
         case 25:
           dval = tr.gy;
         break;
         case 26:
           dval = tr.counit;
         break;
         case 27:
           dval = tr.wevel;
         break;
         case 28:
           dval = tr.swevel;
         break;
         case 29:
           dval = tr.sut;
         break;
         case 30:
           dval = tr.gut;
         break;
         case 31:
           dval = tr.sstat;
         break;
         case 32:
           dval = tr.gstat;
         break;
         case 33:
           dval = tr.tstat;
         break;
         case 34:
           dval = tr.laga;
         break;
         case 35:
           dval = tr.lagb;
         break;
         case 36:
           dval = tr.delrt;
         break;
         case 37:
           dval = tr.muts;
         break;
         case 38:
           dval = tr.mute;
         break;
         case 39:
           dval = tr.ns;
         break;
         case 40:
           dval = tr.dt;
         break;
         case 41:
           dval = tr.gain;
         break;
         case 42:
           dval = tr.igc;
         break;
         case 43:
           dval = tr.igi;
         break;
         case 44:
           dval = tr.corr;
         break;
         case 45:
           dval = tr.sfs;
         break;
         case 46:
           dval = tr.sfe;
         break;
         case 47:
           dval = tr.slen;
         break;
         case 48:
           dval = tr.styp;
         break;
         case 49:
           dval = tr.stas;
         break;
         case 50:
           dval = tr.stae;
         break;
         case 51:
           dval = tr.tatyp;
         break;
         case 52:
           dval = tr.afilf;
         break;
         case 53:
           dval = tr.afils;
         break;
         case 54:
           dval = tr.nofilf;
         break;
         case 55:
           dval = tr.nofils;
         break;
         case 56:
           dval = tr.lcf;
         break;
         case 57:
           dval = tr.hcf;
         break;
         case 58:
           dval = tr.lcs;
         break;
         case 59:
           dval = tr.hcs;
         break;
         case 60:
           dval = tr.year;
         break;
         case 61:
           dval = tr.day;
         break;
         case 62:
           dval = tr.hour;
         break;
         case 63:
           dval = tr.minute;
         break;
         case 64:
           dval = tr.sec;
         break;
         case 65:
           dval = tr.timbas;
         break;
         case 66:
           dval = tr.trwf;
         break;
         case 67:
           dval = tr.grnors;
         break;
         case 68:
           dval = tr.grnofr;
         break;
         case 69:
           dval = tr.grnlof;
         break;
         case 70:
           dval = tr.gaps;
         break;
         case 71:
           dval = tr.otrav;
         break;
         case 72:
           dval = tr.d1;
         break;
         case 73:
           dval = tr.f1;
         break;
         case 74:
           dval = tr.d2;
         break;
         case 75:
           dval = tr.f2;
         break;
         case 76:
           dval = tr.ungpow;
         break;
         case 77:
           dval = tr.unscale;
         break;
         case 78:
           dval = tr.ntr;
         break;
         case 79:
           dval = tr.mark;
         break;
         case 80:
           dval = tr.shortpad;
         break;
  
/*      default:                           */
/*         err("unknown type %s", type);   */
/*      break;                             */
  
        } /* end of   switch */ 
        
      return (dval);
}



/* --------------------------- */
void tohead(segy *tr, int k, double dvalue) {

       switch (k) {
  
         case -1: 
/*       null, name not found? */
         break;
         case 0:  
/*       null   do not write to header */ 
         break;
         case 1:
           tr->tracl = lrint(dvalue);  
         break;
         case 2:
           tr->tracr = lrint(dvalue);
         break;
         case 3:
           tr->fldr = lrint(dvalue);
         break;
         case 4:
           tr->tracf = lrint(dvalue);
         break;
         case 5:
           tr->ep = lrint(dvalue);
         break;
         case 6:
           tr->cdp = lrint(dvalue);
         break;
         case 7:
           tr->cdpt = lrint(dvalue);
         break;
         case 8:
           tr->trid = lrint(dvalue);
         break;
         case 9:
           tr->nvs = lrint(dvalue);
         break;
         case 10:
           tr->nhs = lrint(dvalue);
         break;
         case 11:
           tr->duse = lrint(dvalue);
         break;
         case 12:
           tr->offset = lrint(dvalue);
         break;
         case 13:
           tr->gelev = lrint(dvalue);
         break;
         case 14:
           tr->selev = lrint(dvalue);
         break;
         case 15:
           tr->sdepth = lrint(dvalue);
         break;
         case 16:
           tr->gdel = lrint(dvalue);
         break;
         case 17:
           tr->sdel = lrint(dvalue);
         break;
         case 18:
           tr->swdep = lrint(dvalue);
         break;
         case 19:
           tr->gwdep = lrint(dvalue);
         break;
         case 20:
           tr->scalel = lrint(dvalue);
         break;
         case 21:
           tr->scalco = lrint(dvalue);
         break;
         case 22:
           tr->sx = lrint(dvalue);
         break;
         case 23:
           tr->sy = lrint(dvalue);
         break;
         case 24:
           tr->gx = lrint(dvalue);
         break;
         case 25:
           tr->gy = lrint(dvalue);
         break;
         case 26:
           tr->counit = lrint(dvalue);
         break;
         case 27:
           tr->wevel = lrint(dvalue);
         break;
         case 28:
           tr->swevel = lrint(dvalue);
         break;
         case 29:
           tr->sut = lrint(dvalue);
         break;
         case 30:
           tr->gut = lrint(dvalue);
         break;
         case 31:
           tr->sstat = lrint(dvalue);
         break;
         case 32:
           tr->gstat = lrint(dvalue);
         break;
         case 33:
           tr->tstat = lrint(dvalue);
         break;
         case 34:
           tr->laga = lrint(dvalue);
         break;
         case 35:
           tr->lagb = lrint(dvalue);
         break;
         case 36:
           tr->delrt = lrint(dvalue);
         break;
         case 37:
           tr->muts = lrint(dvalue);
         break;
         case 38:
           tr->mute = lrint(dvalue);
         break;
         case 39:
           tr->ns = lrint(dvalue);
         break;
         case 40:
           tr->dt = lrint(dvalue);
         break;
         case 41:
           tr->gain = lrint(dvalue);
         break;
         case 42:
           tr->igc = lrint(dvalue);
         break;
         case 43:
           tr->igi = lrint(dvalue);
         break;
         case 44:
           tr->corr = lrint(dvalue);
         break;
         case 45:
           tr->sfs = lrint(dvalue);
         break;
         case 46:
           tr->sfe = lrint(dvalue);
         break;
         case 47:
           tr->slen = lrint(dvalue);
         break;
         case 48:
           tr->styp = lrint(dvalue);
         break;
         case 49:
           tr->stas = lrint(dvalue);
         break;
         case 50:
           tr->stae = lrint(dvalue);
         break;
         case 51:
           tr->tatyp = lrint(dvalue);
         break;
         case 52:
           tr->afilf = lrint(dvalue);
         break;
         case 53:
           tr->afils = lrint(dvalue);
         break;
         case 54:
           tr->nofilf = lrint(dvalue);
         break;
         case 55:
           tr->nofils = lrint(dvalue);
         break;
         case 56:
           tr->lcf = lrint(dvalue);
         break;
         case 57:
           tr->hcf = lrint(dvalue);
         break;
         case 58:
           tr->lcs = lrint(dvalue);
         break;
         case 59:
           tr->hcs = lrint(dvalue);
         break;
         case 60:
           tr->year = lrint(dvalue);
         break;
         case 61:
           tr->day = lrint(dvalue);
         break;
         case 62:
           tr->hour = lrint(dvalue);
         break;
         case 63:
           tr->minute = lrint(dvalue);
         break;
         case 64:
           tr->sec = lrint(dvalue);
         break;
         case 65:
           tr->timbas = lrint(dvalue);
         break;
         case 66:
           tr->trwf = lrint(dvalue);
         break;
         case 67:
           tr->grnors = lrint(dvalue);
         break;
         case 68:
           tr->grnofr = lrint(dvalue);
         break;
         case 69:
           tr->grnlof = lrint(dvalue);
         break;
         case 70:
           tr->gaps = lrint(dvalue);
         break;
         case 71:
           tr->otrav = lrint(dvalue);
         break;
         case 72:
           tr->d1 = dvalue;
         break;
         case 73:
           tr->f1 = dvalue;
         break;
         case 74:
           tr->d2 = dvalue;
         break;
         case 75:
           tr->f2 = dvalue;
         break;
         case 76:
           tr->ungpow = dvalue;
         break;
         case 77:
           tr->unscale = dvalue;
         break;
         case 78:
           tr->ntr = lrint(dvalue);
         break;
         case 79:
           tr->mark = lrint(dvalue);
         break;
         case 80:
           tr->shortpad = lrint(dvalue);
         break;

        } /* end of   switch                       */ 

}
